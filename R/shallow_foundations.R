#' Shallow foundation bearing capacity factors
#'
#' @description
#' Calculate shallow foundation bearing capacity factors, based on the
#' value of the friction angle, according to EC7
#'
#' @param phi angle of internal friction, in radians (array)
#' @return tibble with fields `phi` for angle of internal friction, and
#'   `Nc`, `Nq` and `Ngamma` for the three bearing capacity factors
#' @export

calculate_bearingcapacityfactors <- function(phi) {
  df <- tibble::tibble(phi = phi, Nc = 2 + pi)
  df$Nq <- (1 + sin(phi))/(1 - sin(phi))*exp(pi*tan(phi))
  df$Nc[phi > 0] <- (df$Nq[phi > 0] - 1)/tan(phi[phi > 0])
  df$Ngamma <- 2*(df$Nq - 1)*tan(phi)
  df
}


#' Plotly chart with foundation bearing capacity factors
#'
#' @description
#' Creates an interactive plotly chart with the shallow foundation bearing
#' capacity factors, varying with the angle of internal friction, according
#' to EC7
#'
#' @param phi_deg array with all values for the angle of internal friction,
#'   in degrees (array)
#' @param ylim 2-value array with y-axis limits
#' @param nround number of decimals to round angle and N values to
#' @param palette RCOlorBrewer color palette to use for line colors
#' @importFrom rlang .data
#' @return plotly object
#' @examples
#' plotly_bearingcapacityfactors()
#' @export

plotly_bearingcapacityfactors <- function(
  phi_deg = seq(0, 50, 1),
  ylim = c(1, 1000),
  nround = 1,
  palette = "Set1"
){
  #colors
  colo <- RColorBrewer::brewer.pal(3, palette)
  #generate bearing capacity factors
  df <- calculate_bearingcapacityfactors(phi_deg/180*pi)
  df$phi_deg <- phi_deg
  #generate plotly labels
  df <- dplyr::mutate(
    df,
    phi_label = paste0("\u03c6' = ", round(.data$phi_deg, nround), "\u00b0"),
    Nc_label = paste0("N<sub>c</sub> = ", round(.data$Nc, nround)),
    Nq_label = paste0("N<sub>q</sub> = ", round(.data$Nq, nround)),
    Ngamma_label = paste0("N<sub>\u03b3</sub> = ", round(.data$Ngamma, nround))
  )
  #generate empty plotly
  p <- plotly::plot_ly()
  #hidden trace for phi
  p <- plotly::add_trace(
    p,
    type = "scatter",
    mode = "lines",
    name = "phi",
    data = df,
    x = ~phi_deg,
    y = ylim[1],
    line = list(color = "black"),
    text = ~phi_label,
    hoverinfo = "text",
    showlegend = FALSE,
    opacity = 0
  )
  #add traces for N
  p <- plotly::add_trace(
    p,
    type = "scatter",
    mode = "lines",
    name = "N<sub>c</sub>",
    data = df,
    x = ~phi_deg,
    y = ~Nc,
    line = list(color = colo[1]),
    text = ~Nc_label,
    hoverinfo = "text"
  )
  p <- plotly::add_trace(
    p,
    type = "scatter",
    mode = "lines",
    name = "N<sub>q</sub>",
    data = df,
    x = ~phi_deg,
    y = ~Nq,
    line = list(color = colo[2]),
    text = ~Nq_label,
    hoverinfo = "text"
  )
  p <- plotly::add_trace(
    p,
    type = "scatter",
    mode = "lines",
    name = "N<sub>\u03b3</sub>",
    data = df,
    x = ~phi_deg,
    y = ~Ngamma,
    line = list(color = colo[3]),
    text = ~Ngamma_label,
    hoverinfo = "text"
  )
  #add layout
  plotly::layout(
    p,
    xaxis = list(
      range = c(min(phi_deg), max(phi_deg)),
      title = "\u03c6' [\u00b0]"
    ),
    yaxis = list(
      type = "log",
      range = log10(ylim),
      title = "N [-]"
    ),
    hovermode = "y_unified",
    showlegend = TRUE,
    template = "none"
  )
}


#' ggplot chart with foundation bearing capacity factors
#'
#' @description
#' Creates an ggplot chart with the shallow foundation bearing
#' capacity factors, varying with the angle of internal friction, according
#' to EC7
#'
#' @param phi_deg array with all values for the angle of internal friction,
#'   in degrees (array)
#' @param ylim 2-value array with y-axis limits
#' @param palette RCOlorBrewer color palette to use for line colors
#' @return ggplot object
#' @examples
#' ggplot_bearingcapacityfactors()
#' @export

ggplot_bearingcapacityfactors <- function(
  phi_deg = seq(0, 50, 1),
  ylim = c(1, 1000),
  palette = "Set1"
){
  #generate bearing capacity factors
  df <- calculate_bearingcapacityfactors(phi_deg/180*pi)
  df$phi_deg <- phi_deg
  #long data
  dl <- tidyr::pivot_longer(
    df,
    cols = c("Nq", "Nc", "Ngamma"),
    names_to = "parameter",
    values_to = "value"
  ) %>%
    dplyr::filter(.data$value > 0)
  #create and return ggplot
  ggplot2::ggplot(
    dl,
    ggplot2::aes(x = .data$phi_deg, y = .data$value, color = .data$parameter)
  ) +
    theme_soilmech() +
    ggplot2::geom_line() +
    ggplot2::scale_color_brewer(
      name = "",
      palette = palette,
      labels = c(expression(N[c]), expression(N[gamma]), expression(N[q]))
    ) +
    ggplot2::coord_cartesian(
      xlim = c(min(phi_deg), max(phi_deg)),
      ylim = ylim,
      expand = FALSE
    ) +
    ggplot2::scale_y_log10(
      breaks = c(1, 10, 100, 1000),
      minor_breaks = get_log10_minorbreaks(c(1, 1000))
    ) +
    ggplot2::annotation_logticks(side = "l") +
    ggplot2::xlab(expression("Angle of internal friction"~phi*minute~"["*degree*"]")) +
    ggplot2::ylab(expression("Bearing capacity factor"~N~"[-]"))
}


#' ggplot for block failure mechanism underneath shallow foundation
#'
#' @description
#' Creates a ggplot for shallow bearing capacity mechanism, assuming
#' blocks sliding across each other (upper bound mechanism)
#'
#' @param uB vertical displacement of foundation, relative to width B
#' @param phi_deg soil angle of internal friction, in degrees
#' @param nwedge number of wedges to use in the triangular section
#' @param symmetrical if `FALSE`, the bearing capacity failure mechanism is
#'   plotted on one size only. If `TRUE`, failure is symmetrical (both failing
#'   towards the right and the left simultaneously)
#' @param beta_passive_deg angle of the passive wedge, in degrees.
#'   Defaults to 45 + phi_deg/4
#' @param beta_active_deg angle of the active wedge. in degrees.
#'   Defaults to 45 - phi_deg/4
#' @param hB_foundation thickness of the plotted foundation
#' @param xlim user-defined x-axis limits
#' @param ylim user-defined y-axis limits
#' @param fill_soil fill color for soil
#' @param fill_soil2 fill color for moving soil
#' @param color_soil color for outline of soil
#' @param fill_foundation fill color for foundation
#' @param color_foundation color for outline of foundation
#' @param size_soil line thickness of soil polygon outlines
#' @param label_text text to plot in upper right hand corner. If set to
#'   `label_text = Nc_auto`, the Nc coefficient is calculated and
#'   displayed
#' @param label_parse determines whether the label should be parsed or not
#' @param label_size font size of label
#' @param label_nsignif = number of significant digits shown for Nc, in case
#'   `label_text = Nc_auto`
#' @param arrow_interval distance between load arrows (normalised by foundation
#'   width B)
#' @param arrow_load size of load arrow (normalised by foundation
#'   width B). If `arrow_load = 0`, no arrows are plotted
#' @param arrow_surcharge size of surcharge arrow (normalised by foundation
#'   width B). If `arrow_surcharge = 0`, no arrows are plotted
#' @param arrow_size thickness of arrows
#' @param arrow_headsize arrow tip size, in ggplot unit 'npc' (normalised plot
#'   size)
#' @param color_load color of load arrow and text
#' @param color_surcharge color of surcharge arrows and text
#' @importFrom magrittr `%>%`
#' @importFrom rlang .data
#' @return a ggplot object
#' @examples
#' #undrained - undisplaced blocks - symmetrical
#' ggplot_shallowfailuremechanism(uB = 0, phi_deg = 0, symmetrical = TRUE)
#'
#' #drained - displaced blocks - one-sided mechanism
#' ggplot_shallowfailuremechanism(uB = 0.4, phi_deg = 30, symmetrical = FALSE)
#' @export

ggplot_shallowfailuremechanism <- function(
  uB = 0.25,
  phi_deg = 30,
  nwedge = 10,
  symmetrical = FALSE,
  beta_passive_deg = NULL,
  beta_active_deg = NULL,
  hB_foundation = 0.15,
  xlim = c(NA, NA),
  ylim = c(-1, NA),
  fill_soil = "#d3bc5f",
  fill_soil2 = "#bca134",
  color_soil = "#65571d",
  fill_foundation = "black",
  color_foundation = "black",
  size_soil = 0.3,
  label_text = "",
  label_parse = FALSE,
  label_size = 3.5,
  label_nsignif = 4,
  arrow_interval = 0.2,
  arrow_load = 0.40,
  arrow_surcharge = 0.20,
  arrow_size = 0.35,
  arrow_headsize = 0.04,
  color_load = "black",
  color_surcharge = "#6a0dad"
){
  #friction angle
  phi <- phi_deg/180*pi
  #passive and active angles
  if (is.null(beta_passive_deg)){
    beta_passive <- pi/4 + phi/2
  } else {
    beta_passive <- beta_passive_deg/180*pi
  }
  if (is.null(beta_active_deg)){
    beta_active <- pi/4 - phi/2
  } else {
    beta_active <- beta_active_deg/180*pi
  }
  #length of sides of passsive and active wedge
  LAB <- 0.5/cos(beta_passive)
  LBC <- LAB*exp((pi - beta_passive - beta_active)*tan(phi))
  LC0 <- 2*LBC*cos(beta_active)
  #points on wedges
  dp <- tibble::tibble(
    #radii
    r = c(
      2*LAB*cos(beta_passive),
      LAB*exp(seq(0, (pi - beta_passive - beta_active), l = nwedge + 1)*tan(phi)),
      LC0
    ),
    #angles of separating lines
    alpha = c(
      pi,
      pi - beta_passive - seq(0, (pi - beta_passive - beta_active), l = nwedge + 1),
      0
    )
  ) %>% dplyr::mutate(
    #positions
    x0 = .data$r*cos(.data$alpha),
    y0 = .data$r*sin(.data$alpha),
    #angles of failure surfaces
    beta = c(atan2(diff(.data$y0), diff(.data$x0)), NA),
    #initiate displacements - values of passive wedge - assume uB = 1
    uBx = (1 / tan(beta_passive)) * (1 - symmetrical),
    uBy = 1
  )
  #get all displacements - loop
  for (i in 2:(nrow(dp) - 1)) {
    dp$uBx[i] <- (dp$uBy[i-1] - dp$uBx[i-1]*tan(dp$alpha[i])) / (tan(dp$beta[i]) - tan(dp$alpha[i]))
    dp$uBy[i] <- dp$uBx[i]*tan(dp$beta[i])
  }
  #create all polygons
  dpol <- tibble::tibble(
    x0 = 0.5 + c(
      utils::head(dp$x0, -1),
      utils::tail(dp$x0, -1),
      rep(0, nwedge + 2)
    ),
    y0 = c(
      utils::head(dp$y0, -1),
      utils::tail(dp$y0, -1),
      rep(0, nwedge + 2)
    ),
    wedge = rep(seq(nwedge + 2), 3)
  ) %>%
    dplyr::mutate(
      x = .data$x0 + uB * rep(utils::head(dp$uBx, -1), 3),
      y = .data$y0 + uB * rep(utils::head(dp$uBy, -1), 3)
    )
  #if symmetrical - add extra polygons on negative x-side
  if (symmetrical == TRUE) {
    dpol2 <- dpol %>%
      dplyr::filter(.data$wedge != 1) %>%
      dplyr::mutate(
        x0 = -.data$x0,
        x = -.data$x,
        wedge = .data$wedge + nwedge + 2
      )
    dpol <- dplyr::bind_rows(dpol, dpol2)
  }
  #axis limits
  xlim <- round_limits(
    1.05*c(dpol$x, dpol$x0),
    lower = xlim[1],
    upper = xlim[2]
  )
  ylim <- round_limits(
    1.10*c(dpol$y, dpol$y0, -hB_foundation - arrow_load + uB, -arrow_surcharge + uB*dp$uBy[nwedge + 2]),
    lower = ylim[1],
    upper = ylim[2]
  )
  #foundation polygon
  df <- tibble::tibble(
    x0 = c(-0.5, 0.5, 0.5, -0.5),
    y0 = c(-hB_foundation, -hB_foundation, 0, 0)
  ) %>%
    dplyr::mutate(
      x = .data$x0 + (uB/tan(beta_passive))*(1 - symmetrical),
      y = .data$y0 + uB
    )
  #soil polygon
  if (symmetrical == TRUE) {
    ds <- tibble::tibble(
      x = c(xlim[1], -rev(0.5 + dp$x0[2:nrow(dp)]), 0.5 + dp$x0[2:nrow(dp)], xlim[2], xlim[2], xlim[1]),
      y = c(0, rev(dp$y0[2:nrow(dp)]), dp$y0[2:nrow(dp)], 0, ylim[2], ylim[2])
    )
  } else {
    ds <- tibble::tibble(
      x = c(xlim[1], 0.5 + dp$x0, xlim[2], xlim[2], xlim[1]),
      y = c(0, dp$y0, 0, ylim[2], ylim[2])
    )
  }
  #Nc values
  if (label_text == "Nc_auto"){
    if (symmetrical == TRUE) {
      Nc <- dp %>% dplyr::mutate(
        L_beta = c(sqrt(diff(.data$x0)^2 + diff(.data$y0)^2), 0),
        L_alpha = c(sqrt(.data$x0[1:(nrow(dp)-1)]^2 + .data$y0[1:(nrow(dp)-1)]^2), 0),
        u_beta = c(0, sqrt(.data$uBx[2:nrow(dp)]^2 + .data$uBy[2:nrow(dp)]^2)),
        u_alpha = c(0, sqrt(diff(.data$uBx)^2 + diff(.data$uBy)^2))
      ) %>%
        dplyr::summarize(Nc = 2*sum(.data$L_beta*.data$u_beta + .data$L_alpha*.data$u_alpha, na.rm = TRUE)) %>%
        dplyr::pull(Nc)
    } else {
      Nc <- dp %>% dplyr::mutate(
        L_beta = c(sqrt(diff(.data$x0)^2 + diff(.data$y0)^2), 0),
        L_alpha = c(0, sqrt(.data$x0[2:(nrow(dp)-1)]^2 + .data$y0[2:(nrow(dp)-1)]^2), 0),
        u_beta = sqrt(.data$uBx^2 + .data$uBy^2),
        u_alpha = c(0, sqrt(diff(.data$uBx)^2 + diff(.data$uBy)^2))
      ) %>%
        dplyr::summarize(Nc = sum(.data$L_beta*.data$u_beta + .data$L_alpha*.data$u_alpha, na.rm = TRUE)) %>%
        dplyr::pull(Nc)
    }
    #Nc label
    label_text <- paste0("N[c]==", signif(Nc, label_nsignif))
    label_parse <- TRUE
  }
  #initiate plot
  plt <- ggplot2::ggplot() +
    theme_soilmech() +
    #plot non-moving soil
    ggplot2::geom_polygon(
      data = ds,
      ggplot2::aes(x = .data$x, y = .data$y),
      fill = fill_soil,
      color = color_soil,
      size = size_soil
    )
  #add load arrows
  if (arrow_load != 0) {
    dloadarrow <- tibble::tibble(
      x = c(
        -utils::tail(c(seq(0, 0.5, arrow_interval)), -1),
        c(seq(0, 0.5, arrow_interval))
      ) + (uB/tan(beta_passive))*(1 - symmetrical),
      y = -hB_foundation - arrow_load + uB,
      xend = .data$x,
      yend = -hB_foundation + uB
    )
    dloadline <- tibble::tibble(
      x = c(-0.5, 0.5) + (uB/tan(beta_passive))*(1 - symmetrical),
      y = -hB_foundation - arrow_load + uB
    )
    plt <- plt +
      ggplot2::geom_segment(
        data = dloadarrow,
        ggplot2::aes(x = .data$x, xend = .data$xend, y = .data$y, yend = .data$yend),
        color = color_load,
        size = arrow_size,
        arrow = ggplot2::arrow(length = ggplot2::unit(arrow_headsize, "npc"))
      ) +
      ggplot2::geom_path(
        data = dloadline,
        ggplot2::aes(x = .data$x, y = .data$y),
        color = color_load,
        size = arrow_size
      ) +
      ggplot2::annotate(
        "text",
        x = 0 + (uB/tan(beta_passive))*(1 - symmetrical),
        y = -hB_foundation - arrow_load + uB - 0.01*diff(ylim),
        label = ifelse(phi_deg == 0, "q[f]", "q*minute[f]"),
        parse = TRUE,
        hjust = 0.5,
        vjust = 0,
        size = label_size,
        color = color_load
      )
  }
  #add surcharge arrows
  if (arrow_surcharge != 0) {
    #surcharge that moves
    xsurchargearrow <- seq(0, 0.5 + LC0, arrow_interval)
    dsurchargearrow <- tibble::tibble(
      x = xsurchargearrow[xsurchargearrow > 0.5] + uB*dp$uBx[nwedge + 2],
      y = -arrow_surcharge + uB*dp$uBy[nwedge + 2],
      xend = xsurchargearrow[xsurchargearrow > 0.5] + uB*dp$uBx[nwedge + 2],
      yend = uB*dp$uBy[nwedge + 2]
    )
    dsurchargeline <- tibble::tibble(
      x = c(0.5, 0.5 + LC0) + uB*dp$uBx[nwedge + 2],
      y = c(0, 0) - arrow_surcharge + uB*dp$uBy[nwedge + 2]
    )
    #surcharge that remains stationary
    xsurchargearrow1 <- seq(0, xlim[2], arrow_interval)
    dsurchargearrow1 <- tibble::tibble(
      x = xsurchargearrow1[xsurchargearrow1 > (0.5 + LC0)],
      y = -arrow_surcharge,
      xend = xsurchargearrow1[xsurchargearrow1 > (0.5 + LC0)],
      yend = 0
    )
    dsurchargeline1 <- tibble::tibble(
      x = c(0.5 + LC0, xlim[2]),
      y = c(0, 0) - arrow_surcharge
    )
    xsurchargearrow2 <- seq(0, xlim[1], -arrow_interval)
    if (symmetrical == TRUE) {
      xsurchargearrow2 <- xsurchargearrow2[xsurchargearrow2 < (-0.5 - LC0)]
    } else {
      xsurchargearrow2 <- xsurchargearrow2[xsurchargearrow2 < -0.5]
    }
    dsurchargearrow2 <- tibble::tibble(
      x = xsurchargearrow2,
      y = -arrow_surcharge,
      xend = xsurchargearrow2,
      yend = 0
    )
    dsurchargeline2 <- tibble::tibble(
      x = c(-0.5 - LC0*as.double(symmetrical), xlim[1]),
      y = c(0, 0) - arrow_surcharge
    )
    #add to plots
    plt <- plt +
      ggplot2::geom_segment(
        data = dsurchargearrow,
        ggplot2::aes(x = .data$x, xend = .data$xend, y = .data$y, yend = .data$yend),
        color = color_surcharge,
        size = arrow_size,
        arrow = ggplot2::arrow(length = ggplot2::unit(arrow_headsize, "npc"))
      ) +
      ggplot2::geom_path(
        data = dsurchargeline,
        ggplot2::aes(x = .data$x, y = .data$y),
        color = color_surcharge,
        size = arrow_size
      ) +
      ggplot2::geom_segment(
        data = dsurchargearrow1,
        ggplot2::aes(x = .data$x, xend = .data$xend, y = .data$y, yend = .data$yend),
        color = color_surcharge,
        size = arrow_size,
        arrow = ggplot2::arrow(length = ggplot2::unit(arrow_headsize, "npc"))
      ) +
      ggplot2::geom_path(
        data = dsurchargeline1,
        ggplot2::aes(x = .data$x, y = .data$y),
        color = color_surcharge,
        size = arrow_size
      ) +
      ggplot2::geom_segment(
        data = dsurchargearrow2,
        ggplot2::aes(x = .data$x, xend = .data$xend, y = .data$y, yend = .data$yend),
        color = color_surcharge,
        size = arrow_size,
        arrow = ggplot2::arrow(length = ggplot2::unit(arrow_headsize, "npc"))
      ) +
      ggplot2::geom_path(
        data = dsurchargeline2,
        ggplot2::aes(x = .data$x, y = .data$y),
        color = color_surcharge,
        size = arrow_size
      ) +
      ggplot2::annotate(
        "text",
        x = 0.5 + 0.5*LC0 + uB*dp$uBx[nwedge + 2],
        y = -arrow_surcharge + uB*dp$uBy[nwedge + 2] - 0.01*diff(ylim),
        label = ifelse(phi_deg == 0, "q[0]", "q*minute[0]"),
        parse = TRUE,
        hjust = 0.5,
        vjust = 0,
        size = label_size,
        color = color_surcharge
      )
    if (symmetrical == TRUE) {
      plt <- plt +
        ggplot2::geom_segment(
          data = dsurchargearrow,
          ggplot2::aes(x = -.data$x, xend = -.data$xend, y = .data$y, yend = .data$yend),
          color = color_surcharge,
          size = arrow_size,
          arrow = ggplot2::arrow(length = ggplot2::unit(arrow_headsize, "npc"))
        ) +
        ggplot2::geom_path(
          data = dsurchargeline,
          ggplot2::aes(x = -.data$x, y = .data$y),
          color = color_surcharge,
          size = arrow_size
        )
    }
  }
  #plot moving wedges
  plt <- plt + ggplot2::geom_polygon(
    data = dpol,
    ggplot2::aes(x = .data$x, y = .data$y, group = as.factor(.data$wedge)),
    fill = fill_soil2,
    color = color_soil,
    size = size_soil
  ) +
  #plot foundation polygon
  ggplot2::geom_polygon(
    data = df,
    ggplot2::aes(x = .data$x, y = .data$y),
    fill = fill_foundation,
    color = color_foundation
  ) +
  #annotate label, if requested
  ggplot2::annotate(
    "text",
    x = xlim[1] + 0.95*(xlim[2] - xlim[1]),
    y = ylim[1] + 0.05*(ylim[2] - ylim[1]),
    label = label_text,
    hjust = 1,
    vjust = 1,
    parse = label_parse,
    size = label_size
  ) +
  #set axes
  ggplot2::scale_y_reverse() +
  ggplot2::coord_fixed(
    ratio = 1,
    xlim = xlim,
    ylim = rev(ylim),
    expand = FALSE
  ) +
  ggplot2::xlab("x/B [-]") +
  ggplot2::ylab("z/B [-]")
  #return
  return(plt)
}
