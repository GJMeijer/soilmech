#' Influence factor for stress under strip load
#'
#' @description
#' Calculates the influence factor for the total stress underneath
#' a strip load, based on the position in normalised coordinates
#'
#' @param xB x-coordinate, normalised by width of strip (array)
#' @param zB depth, normalised by width of strip (array)
#' @param direction direction of stress increase. `z` for vertical, `x` for
#'   horizontal away from strip, `y` for horizontal in direction of strip, and
#'   `xz` for shear stress
#' @param nu Posisson's ratio
#' @return array with influence factors
#' @export

influencefactor_uniformstrip <- function(xB, zB, direction = "z", nu = 0.3){
  delta <- atan2(xB - 0.5, zB)
  alphadelta <- atan2(xB + 0.5, zB)
  alpha <- alphadelta - delta
  if (direction == "z") {
    return((alpha + sin(alpha)*cos(alpha + 2*delta))/pi)
  } else if (direction == "x") {
    return((alpha - sin(alpha)*cos(alpha + 2*delta))/pi)
  } else if (direction == "y") {
    return((2*nu*alpha)/pi)
  } else if (direction == "xz") {
    return((sin(alpha)*sin(alpha + 2*delta))/pi)
  } else {
    return(rep(NA, length(xB)))
  }
}


#' Influence factor for stress under triangular load
#'
#' @description
#' Calculates the influence factor for the total stress underneath
#' a strip load, based on the position in normalised coordinates
#'
#' @param xB x-coordinate, normalised by width of strip (array). `xB = 0` is
#'   located where the triangle height is 0, and `xB = 1` where the triangular
#'   load is at it maximum (just before it drops to zero after `xB > 1`)
#' @param zB depth, normalised by width of strip (array)
#' @param direction direction of stress increase. `z` for vertical, `x` for
#'   horizontal away from strip, `y` for horizontal in direction of strip, and
#'   `xz` for shear stress
#' @param nu Posisson's ratio
#' @return array with influence factors
#' @export

influencefactor_triangularstrip <- function(xB, zB, direction = "z", nu = 0.3){
  delta <- atan2(xB - 1, zB)
  alphadelta <- atan2(xB, zB)
  alpha <- alphadelta - delta
  if (direction == "z") {
    return((xB*alpha - 0.5*sin(2*delta))/pi)
  } else if (direction == "x") {
    R1sq <- xB^2 + zB^2
    R2sq <- (xB - 1)^2 + zB^2
    logR1R2 <- log(R1sq/R2sq)
    logR1R2[!is.finite(logR1R2)] <- 0
    return((xB*alpha - zB*logR1R2 + 0.5*sin(2*delta))/pi)
  } else if (direction == "xz") {
    return((0.5 + 0.5*cos(2*delta) - zB*alpha)/pi)
  } else {
    return(rep(NA, length(xB)))
  }
}


#' Influence factor for stress under corner of rectangular uniform load
#'
#' @description
#' Calculates the influence factor for a point underneath the corner of a
#' rectangular uniform load
#'
#' @param xz x-coordinate, normalised by depth z
#' @param yz y-coordinate, normalised by depth z
#' @param direction direction of stress increase. `z` for vertical
#' @return array with influence factors
#' @export

influencefactor_cornerrectangle <- function(xz, yz, direction = "z"){
  if (direction == "z") {
    return(
      1/(2*pi) *
        (xz*yz*(xz^2 + yz^2+2)/((1 + xz^2)*(1 + yz^2)*sqrt(1 + xz^2 + yz^2)) +
           atan(xz*yz/sqrt(1 + xz^2 + yz^2)))
    )
  } else {
    return(rep(NA, length(xz)))
  }
}


#' ggplot influence factors for stress underneath uniform strip foundation
#'
#' @description
#' Creates a ggplot chart with contour lines and heat map (optional) for the
#' influence factors for stress underneath a uniform strip load.
#' Assumes elastic solutions are valid.
#'
#' @param direction the direction of the stress of interest. `z` for
#'   vertical, `x` for horizontal away from strip, `y` for horizontal
#'   in direction of strip, and `xz` for shear stress
#' @param nu Poisson's ratio, required for stress in y-direction
#' @param xlim 2-value array with limits for x-axis
#' @param ylim 2-value array with limits for depth axis
#' @param gridsize grid size for calculations
#' @param contour_width I-interval between contour lines
#' @param contour_label values of I for which to draw contour lines
#' @param showgrid if `TRUE`, heatmap of influence factors is shown in
#'   background of the plot
#' @param larrow length of load arrows (in normalised coordinates)
#' @param lsize thickness of load arrows
#' @param labelsize size of labels on contour lines
#' @param palette RColorBrewer palette for heat map
#' @param xB crosshairs annotation
#' @param zB crosshairs annotation
#' @param nround number of decimals in crosshairs I label
#' @importFrom magrittr `%>%`
#' @importFrom rlang .data
#' @return ggplot object
#' @examples
#' #plot empty chart
#' ggplot_stress_uniformstrip()
#'
#' #plot with annotation
#' ggplot_stress_uniformstrip(xB = 1, zB = 2)
#' @export

ggplot_stress_uniformstrip <- function(
  direction = "z",
  nu = 0.3,
  xlim = c(-2.5, 2.5),
  ylim = c(0, 7),
  gridsize = 0.05,
  contour_width = 0.1,
  contour_label = c(0, 0.1, 0.2, 0.4, 0.6, 0.8, 1.0),
  showgrid = TRUE,
  larrow = 0.3,
  lsize = 3.5,
  labelsize = 3,
  palette = "Blues",
  xB = NULL,
  zB = NULL,
  nround = 3
){
  #all points
  dl <- tidyr::expand_grid(
    xB = seq(xlim[1], xlim[2], gridsize),
    zB = seq(ylim[1], ylim[2], gridsize)
  )
  #calculate influence factors
  dl$I <- abs(influencefactor_uniformstrip(dl$xB, dl$zB, direction = direction, nu = nu))
  #contour line labels and positions of labels - try to place near bottom of contourline
  contour_label <- contour_label[contour_label <= max(dl$I)]
  dlab <- tibble::tibble(
    I = contour_label,
  ) %>% dplyr::mutate(
    purrr::map_dfr(
      contour_label,
      function(xx) {
        dl %>%
          dplyr::filter((I >= xx) & (.data$xB >= 0)) %>%
          dplyr::slice_max(.data$zB) %>%
          dplyr::slice_max(.data$I) %>%
          dplyr::summarize(zB = mean(.data$zB), xB = min(.data$xB))
      }
    )
  )
  dlab$xB[dlab$I == 0] <- 0.5*(xlim[2] + 0.5)
  dlab$zB[dlab$I == 0] <- 0
  dlab$label <- dlab$I
  dlab$label[dlab$I == 0.0] <- paste("I[", direction, "]==0.0")
  #plot
  plt <- ggplot2::ggplot() +
    theme_soilmech() +
    ggplot2::theme(legend.position = "none") +
    ggplot2::theme(
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank()
    )
  #add color to ground
  if (showgrid == TRUE){
    plt <- plt + ggplot2::geom_tile(
      data = dl,
      ggplot2::aes(x = .data$xB, y = .data$zB + 0.5*gridsize, fill = .data$I),
      show.legend = FALSE
    )
  }
  #manually plot a grid
  grid_size <- 0.5
  grid_x <- grid_size*round(seq(xlim[1], xlim[2], by = grid_size)/grid_size)
  grid_y <- grid_size*round(seq(ylim[1], ylim[2], by = grid_size)/grid_size)
  plt <- plt +
    ggplot2::geom_hline(yintercept = grid_y, size = 0.2, color = "grey90") +
    ggplot2::geom_vline(xintercept = grid_x, size = 0.2, color = "grey90") +
    ggplot2::geom_hline(yintercept = 0)
  #add contour lines
  plt <- plt +
    ggplot2::stat_contour(
      data = dl,
      ggplot2::aes(x = .data$xB, y = .data$zB, z = .data$I),
      binwidth = contour_width,
      color = "black",
      linetype = 1,
      show.legend = FALSE
    ) +
    ggplot2::xlab("x/B [-]") +
    ggplot2::ylab("z/B [-]") +
    ggplot2::scale_y_reverse() +
    ggplot2::coord_fixed(
      ratio = 1,
      xlim = xlim,
      ylim = c(ylim[2], ylim[1] - 2*larrow),
      expand = FALSE
    ) +
    ggplot2::scale_fill_distiller(palette = palette, direction = 1)
  #annotate load at surface by means of arrows
  darr <- tibble::tibble(
    x = seq(-0.5, 0.5, 0.2),
    xend = seq(-0.5, 0.5, 0.2),
    y = -larrow,
    yend = 0
  )
  plt <- plt +
    ggplot2::geom_segment(
      data = darr,
      ggplot2::aes(x = .data$x, xend = .data$xend, y = .data$y, yend = .data$yend),
      arrow = ggplot2::arrow(length = ggplot2::unit(0.03, "npc"))
    ) +
    ggplot2::annotate("segment", x = -0.5, y = -larrow, xend = 0.5, yend = -larrow)
  #add contour labels
  plt <- plt + ggplot2::geom_label(
    data = dlab,
    ggplot2::aes(x = .data$xB, y = .data$zB, label = .data$label),
    fill = "white",
    alpha = 0.7,
    size = labelsize,
    label.padding = ggplot2::unit(0.1, "lines"),
    label.r = ggplot2::unit(0.1, "lines"),
    parse = TRUE,
    hjust = 0.5,
    vjust = 0.5,
  )
  #add annotation
  if (!is.null(xB) & !is.null(zB)) {
    I <- influencefactor_uniformstrip(xB, zB, direction = direction, nu = nu)
    plt <- ggplot_addcrosshairs(
      plt,
      x = xB,
      y = zB,
      group = paste0("I[", direction, "]==", round(I, nround)),
      add_colour_scale = TRUE,
      label_parse = TRUE,
      xlim = xlim[1],
      ylim = ylim[2],
    )
  }
  #return
  plt
}


#' ggplot influence factors for stress underneath triangular strip foundation
#'
#' @description
#' Creates a ggplot chart with contour lines and heat map (optional) for the
#' influence factors for stress underneath a triangular strip load.
#' Assumes elastic solutions are valid
#'
#' @inheritParams ggplot_stress_uniformstrip
#' @importFrom magrittr `%>%`
#' @importFrom rlang .data
#' @return ggplot object
#' @examples
#' #plot empty chart
#' ggplot_stress_triangularstrip()
#'
#' #plot with annotation
#' ggplot_stress_triangularstrip(xB = 1, zB = 2)
#' @export

ggplot_stress_triangularstrip <- function(
  direction = "z",
  nu = 0.3,
  xlim = c(-1, 2),
  ylim = c(0, 3.5),
  gridsize = 0.05,
  contour_width = 0.1,
  contour_label = c(0, 0.1, 0.2, 0.4, 0.6, 0.8, 1.0),
  showgrid = TRUE,
  larrow = 0.3, #plot arrow length
  lsize = 3.5,
  labelsize = 3,
  palette = "Blues",
  xB = NULL,
  zB = NULL,
  nround = 3
){
  #all points
  dl <- tidyr::expand_grid(
    xB = seq(xlim[1], xlim[2], gridsize),
    zB = seq(ylim[1], ylim[2], gridsize)
  )
  #calculate influence factors
  dl$I <- abs(influencefactor_triangularstrip(dl$xB, dl$zB, direction = direction, nu = nu))
  #contour line labels and positions of labels - try to place near bottom of contourline
  contour_label <- contour_label[contour_label <= max(dl$I, na.rm = TRUE)]
  dlab <- tibble::tibble(
    I = contour_label,
  ) %>% dplyr::mutate(
    purrr::map_dfr(
      contour_label,
      ~(dl %>%
          dplyr::filter((I >= .x) & (xB >= 0)) %>%
          dplyr::slice_max(zB) %>%
          dplyr::slice_max(I) %>%
          dplyr::summarize(zB = mean(.data$zB), xB = min(.data$xB))
      )
    )
  )
  dlab$xB[dlab$I == 0] <- 0.5*(xlim[2] + 1)
  dlab$zB[dlab$I == 0] <- 0
  dlab$label <- dlab$I
  dlab$label[dlab$I == 0.0] <- paste0("I[", direction, "]==0.0")
  #plot
  plt <- ggplot2::ggplot() +
    theme_soilmech() +
    ggplot2::theme(legend.position='none') +
    ggplot2::theme(
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank()
    )
  #add heatmap
  if (showgrid == TRUE){
    plt <- plt + ggplot2::geom_tile(
      data = dl,
      ggplot2::aes(x = .data$xB, y = .data$zB + 0.5*gridsize, fill = .data$I)
    )
  }
  #manually plot a grid
  grid_size <- 0.5
  grid_x <- grid_size*round(seq(xlim[1], xlim[2], by = grid_size)/grid_size)
  grid_y <- grid_size*round(seq(ylim[1], ylim[2], by = grid_size)/grid_size)
  plt <- plt +
    ggplot2::geom_hline(yintercept = grid_y, size = 0.2, color = "grey90") +
    ggplot2::geom_vline(xintercept = grid_x, size = 0.2, color = "grey90") +
    ggplot2::geom_hline(yintercept = 0)
  #add countour lines and labels
  plt <- plt +
    ggplot2::stat_contour(
      data = dl,
      ggplot2::aes(x = .data$xB, y = .data$zB, z = .data$I),
      binwidth = 0.1,
      color = "black",
      linetype = 1
    ) +
    ggplot2::xlab("x/B [-]") +
    ggplot2::ylab("z/B [-]") +
    ggplot2::scale_y_reverse() +
    ggplot2::coord_fixed(
      ratio = 1,
      xlim = xlim,
      ylim = c(ylim[2], ylim[1] - 2*larrow),
      expand = FALSE
    ) +
    ggplot2::scale_fill_distiller(palette = palette, direction = 1)
  #add arrows to indicate load
  darr <- tibble::tibble(
    x = seq(0.2, 1, 0.2),
    xend = seq(0.2, 1, 0.2),
    y = -seq(0.2, 1, 0.2)*larrow,
    yend = 0
  )
  plt <- plt +
    ggplot2::geom_segment(
      data = darr,
      ggplot2::aes(x = .data$x, xend = .data$xend, y = .data$y, yend = .data$yend),
      arrow = ggplot2::arrow(length = ggplot2::unit(0.03, "npc"))
    ) +
    ggplot2::annotate("segment", x = 0.0, y = 0.0, xend = 1.0, yend = -larrow)
  #add contour labels
  plt <- plt + ggplot2::geom_label(
    data = dlab,
    ggplot2::aes(x = .data$xB, y = .data$zB, label = .data$label),
    fill = "white",
    alpha = 0.7,
    size = labelsize,
    label.padding = ggplot2::unit(0.1, "lines"),
    label.r = ggplot2::unit(0.1, "lines"),
    parse = TRUE
  )
  #add annotation
  if (!is.null(xB) & !is.null(zB)) {
    I <- influencefactor_triangularstrip(xB, zB, direction = direction, nu = nu)
    plt <- ggplot_addcrosshairs(
      plt,
      x = xB,
      y = zB,
      group = paste0("I[", direction, "]==", round(I, nround)),
      add_colour_scale = TRUE,
      label_parse = TRUE,
      xlim = xlim[1],
      ylim = ylim[2],
    )
  }
  #return
  plt
}


#' ggplot Fadum's chart
#'
#' @description
#' Plot Fadum's chart for the influence factors for a point underneath a
#' corner of a rectagular uniform load
#'
#' @param palette RColorBrewer color palette to use to more easily differentiate
#'   between lines.
#' @param Lz_inf L/z value associated with a strip, i.e L/z --> infinity. Choose
#'   very large
#' @param Bz,Lz arrays with values for which to draw crosshairs. If these are
#'   specified, all Fadum's lines are drawn as black, and the annotated crosshairs
#'   are plotted in color
#' @param nround number of digits to use in rounding in crosshair labels
#' @importFrom magrittr `%>%`
#' @importFrom rlang .data
#' @return ggplot object
#' @examples
#' #plot empty chart
#' ggplot_stress_fadum()
#'
#' #plot with annotation
#' ggplot_stress_fadum(Bz = 1, Lz = 1)
#' @export

ggplot_stress_fadum <- function(
  palette = NULL,
  Lz_inf = 1e6,
  Bz = NULL,
  Lz = NULL,
  nround = 3
){
  #colors
  if (is.null(palette)){
    colo <- rep("black", 4)
  } else {
    colo <- RColorBrewer::brewer.pal(4, palette)
  }
  #intervals
  d <- dplyr::bind_rows(
    tibble::tibble(yz = c(0.5, 1.0, Lz_inf), xz0 = 0.01, xz1 = 10, color = colo[1]),
    tibble::tibble(yz = c(0.1, 0.2, 0.3, 0.4), xz0 = 0.03, xz1 = 10, color = colo[2]),
    tibble::tibble(yz = c(0.6, 0.7, 0.8, 0.9, 1.2), xz0 = 0.1, xz1 = 10, color = colo[3]),
    tibble::tibble(yz = c(1.4, 1.6, 1.8, 2.0, 2.5, 3.0), xz0 = 1, xz1 = 10, color = colo[4])
  ) %>%
    dplyr::arrange(.data$yz)
  #for each trace, draw evenly spaced points on logscale
  dl <- d %>%
    tidyr::expand_grid(frac = seq(0, 1, l = 101)) %>%
    dplyr::mutate(
      xz = 10^(log10(.data$xz0) + .data$frac*(log10(.data$xz1) - log10(.data$xz0))),
      I = influencefactor_cornerrectangle(.data$xz, .data$yz, direction = "z")
    )
  #create a label for each trace
  marg <- 0.06 * c(rep(1, 12), 0.5, 1.0, 1.5, 0.5, 1.0, 1.5)
  d$xz_lab <- 10^(log10(min(d$xz0)) + (1.0 - marg)*(log10(max(d$xz1)) - log10(min(d$xz0))))
  d$yz_lab <- influencefactor_cornerrectangle(d$xz_lab, d$yz, direction = "z")
  d$text <- d$yz
  d$text[d$yz == Lz_inf] <- "infinity"
  d$text[d$yz == min(d$yz)] <- paste0("B/z==", min(d$yz))
  #plot
  if (!is.null(palette) & is.null(Bz) & is.null(Lz)) {
    plt <- ggplot2::ggplot(
      dl,
      ggplot2::aes(x = .data$xz, y = .data$I, color = as.factor(.data$yz))
    ) +
      ggplot2::scale_color_manual(values = as.character(d$color))
  } else {
    plt <- ggplot2::ggplot() +
      ggplot2::geom_line(
        data = dl,
        ggplot2::aes(x = .data$xz, y = .data$I, group = as.factor(.data$yz))
      )
  }
  plt <- plt +
    theme_soilmech() +
    ggplot2::theme(legend.position = "none") +
    ggplot2::geom_line(size = 0.3, show.legend = FALSE) +
    ggplot2::xlab("L/z [-]") +
    ggplot2::ylab(expression("Influence factor"~I[z]~"[-]")) +
    ggplot2::coord_cartesian(
      xlim = c(min(dl$xz), max(dl$xz)),
      ylim = c(0, 0.26),
      expand = FALSE
    ) +
    ggplot2::scale_y_continuous(
      breaks = seq(0, 0.26, by = 0.02),
      minor_breaks = seq(0, 0.26, by = 0.01)
    ) +
    ggplot2::scale_x_log10(
      breaks = c(0.01, 0.1, 1, 10),
      minor_breaks = get_log10_minorbreaks(c(0.01, 10))
    ) +
    ggplot2::geom_label(
      data = d,
      ggplot2::aes(x = .data$xz_lab, y = .data$yz_lab, label = .data$text),
      hjust = 0.5,
      vjust = 0.5,
      alpha = 1,
      size = 2,
      label.padding = ggplot2::unit(0.1, "lines"),
      parse = TRUE,
      show.legend = FALSE
    )
  #add crosshairs
  if (!is.null(Bz) & !is.null(Lz)){
    I <- influencefactor_cornerrectangle(Bz, Lz, direction = "z")
    lab <- paste("I[z]==", round(I, nround))
    plt <- ggplot_addcrosshairs(
      plt,
      Lz,
      I,
      group = lab,
      xlim = min(d$xz0),
      arrow_x = FALSE,
      label_parse = TRUE,
      add_colour_scale = TRUE
    )
  }
  #return
  plt
}


#' ggplot Giroud's chart
#'
#' @description
#' Plot Giroud's chart for the influence factors for a point underneath a
#' corner of a rectagular uniform load
#'
#' @param LB_inf L/B value associated with a strip, i.e L/B --> infinity. Choose
#'   very large
#' @param zB,LB arrays with values for which to draw crosshairs. If these are
#'   specified, all Giroud's lines are drawn as black, and the annotated crosshairs
#'   are plotted in color
#' @param nround number of digits to use in rounding in crosshair labels
#' @importFrom magrittr `%>%`
#' @importFrom rlang .data
#' @return ggplot object
#' @examples
#' #plot empty chart
#' ggplot_stress_giroud()
#'
#' #plot with annotation
#' ggplot_stress_giroud(zB = 1, LB = 1)
#' @export

ggplot_stress_giroud <- function(
  LB_inf = 1e6,
  zB = NULL,
  LB = NULL,
  nround = 3
){
  #intervals
  d <- tibble::tibble(
    zB0 = 0.1,
    zB1 = 10,
    LB = c(1, 2, 5, LB_inf)
  )
  #for each trace, draw evenly spaced points on logscale
  dl <- d %>%
    tidyr::expand_grid(frac = seq(0, 1, l = 101)) %>%
    dplyr::mutate(
      zB = 10^(log10(.data$zB0) + .data$frac*(log10(.data$zB1) - log10(.data$zB0))),
      I = ifelse(
        .data$zB == 0,
        0.25,
        influencefactor_cornerrectangle(1/.data$zB, .data$LB/.data$zB, direction = "z")
      )
    )
  #create a label for each trace
  d$I_lab <- c(0.05, 0.05, 0.05, 0.04)
  d$zB_lab <- c(2.8, 3.8, 5.2, 8.0)
  d$text <- c("square", "L/B==2", "L/B==5", "strip")
  #plot
  plt <- ggplot2::ggplot() +
    theme_soilmech() +
    ggplot2::theme(legend.position = "none") +
    ggplot2::geom_path(
      data = dl,
      ggplot2::aes(x = .data$I, y = log10(.data$zB), group = as.factor(.data$LB)),
      size = 0.3, show.legend = FALSE
    ) +
    ggplot2::xlab(expression("Influence factor"~I[z]~"[-]")) +
    ggplot2::ylab("z/B [-]") +
    ggplot2::coord_cartesian(
      xlim = c(0, 0.25),
      ylim = c(log10(10), log10(0.1)),
      expand = FALSE
    ) +
    ggplot2::scale_x_continuous(
      breaks = c(seq(0, 0.20, 0.04), 0.25),
      minor_breaks = c(seq(0, 0.24, 0.02), 0.25)
    ) +
    ggplot2::scale_y_reverse(
      breaks = log10(c(10, 8, 6, 5, 4, 3, 2, 1.6, 1.4, 1.2, 1.0, 0.8, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1)),
      labels = (c(10, 8, 6, 5, 4, 3, 2, 1.6, 1.4, 1.2, 1.0, 0.8, 0.6, 0.5, 0.4, 0.3, 0.2, 0.0)),
      minor_breaks = log10(get_log10_minorbreaks(c(0.01, 10))),
      position = "right"
    ) +
    ggplot2::geom_label(
      data = d,
      ggplot2::aes(x = .data$I_lab, y = log10(.data$zB_lab), label = .data$text),
      alpha = 0.7,
      hjust = 0.5,
      vjust = 0.5,
      size = 2,
      label.padding = ggplot2::unit(0.1, "lines"),
      parse = TRUE,
      show.legend = FALSE
    )
  #add crosshairs
  if (!is.null(zB) & !is.null(LB)){
    I <- ifelse(
      zB == 0,
      0.25,
      influencefactor_cornerrectangle(1/zB, LB/zB, direction = "z")
    )
    zB <- pmax(0.1, zB)
    lab <- paste0("I[z]==", round(I, nround))
    plt <- ggplot_addcrosshairs(
      plt,
      I,
      log10(zB),
      group = lab,
      xlim = 0.25,
      ylim = log10(max(d$zB1)),
      arrow_y = FALSE,
      add_colour_scale = TRUE,
      label_parse = TRUE,
      label_hjust = 1
    )
  }
  #return
  plt
}
