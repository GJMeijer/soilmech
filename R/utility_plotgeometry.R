#' Annotate water table marker to a ggplot object
#'
#' @description
#' Adds a standard water table marker (upside-down triangle over a number of
#' horizontal lines) to an existing ggplot
#'
#' @param plt ggplot object
#' @param xc,yc x and y-position where marker touches water table
#' @param theta rotation angle (in case label needs to be rotated, in rad)
#' @param scale scaling factor for marker. Default height is approx 1
#' @param linesize thickness of lines
#' @param fill_water fill colour for triangle
#' @param colour_water line colour for all lines
#' @return ggplot object with added water table marker
#' @examples
#' #generate simple plot
#' plt <- ggplot2::ggplot() +
#'   ggplot2::annotate("rect", xmin = -2, xmax = 2, ymin = -2, ymax = 2, fill = "yellow") +
#'   ggplot2::coord_fixed(ratio = 1)
#'
#' #add water table marker
#' ggplot_add_watermarker(plt, x = 1, y = 1, scale = 1, theta = 0.4)
#' @export

ggplot_add_watermarker <- function(
  plt,
  xc,
  yc,
  theta = 0,
  scale = 1,
  linesize = 0.5,
  fill_water = "#2a7fff",
  colour_water = "#000080"
){
  #polygon coordinates
  xpu <- scale*c(0, -0.5*tan(pi/6), 0.5*tan(pi/6))
  ypu <- scale*c(0, 0.5, 0.5)
  xp <- xc + xpu*cos(theta) - ypu*sin(theta)
  yp <- yc + xpu*sin(theta) + ypu*cos(theta)
  #line coordinates
  xlu1 <- 0.25*scale*c(-1, 1)
  xlu2 <- 0.15*scale*c(-1, 1)
  xlu3 <- 0.05*scale*c(-1, 1)
  ylu1 <- -0.15*scale*c(1, 1)
  ylu2 <- -0.30*scale*c(1, 1)
  ylu3 <- -0.45*scale*c(1, 1)
  xl1 <- xc + xlu1*cos(theta) - ylu1*sin(theta)
  xl2 <- xc + xlu2*cos(theta) - ylu2*sin(theta)
  xl3 <- xc + xlu3*cos(theta) - ylu3*sin(theta)
  yl1 <- yc + xlu1*sin(theta) + ylu1*cos(theta)
  yl2 <- yc + xlu2*sin(theta) + ylu2*cos(theta)
  yl3 <- yc + xlu3*sin(theta) + ylu3*cos(theta)
  #plot
  plt <- plt +
    ggplot2::annotate(
      "polygon",
      x = xp,
      y = yp,
      size = linesize,
      fill = fill_water,
      color = colour_water
    ) +
    ggplot2::annotate(
      "segment",
      x = xl1[1],
      xend = xl1[2],
      y = yl1[1],
      yend = yl1[2],
      size = linesize,
      color = colour_water
    ) +
    ggplot2::annotate(
      "segment",
      x = xl2[1],
      xend = xl2[2],
      y = yl2[1],
      yend = yl2[2],
      size = linesize,
      color = colour_water
    ) +
    ggplot2::annotate(
      "segment",
      x = xl3[1],
      xend = xl3[2],
      y = yl3[1],
      yend = yl3[2],
      size = linesize,
      color = colour_water
    )
  #return
  return(plt)
}


#' Annotate soil surface marker to a ggplot object
#'
#' @description
#' Adds a standard soil surface indicator icon to an existing ggplot object
#'
#' @param plt ggplot object
#' @param xc,yc x and y-position where marker touches the soil surface (middle)
#' @param theta rotation angle (in case label needs to be rotated, in rad)
#' @param scale scaling factor for marker. Default width is approx 1
#' @param linesize thickness of lines
#' @param colour_soil line colour for all lines
#' @param n number of lines to plot underneath the water table, as part of the symbol
#' @return ggplot object with added water table marker
#' @importFrom rlang .data
#' @examples
#' #generate simple plot
#' plt <- ggplot2::ggplot() +
#'   ggplot2::annotate("rect", xmin = -2, xmax = 2, ymin = -2, ymax = 2, fill = "yellow") +
#'   ggplot2::coord_fixed(ratio = 1)
#'
#' #add water table marker
#' ggplot_add_soilmarker(plt, x = 1, y = 1, scale = 1, theta = 0.5)
#' @export

ggplot_add_soilmarker <- function(
  plt,
  xc,
  yc,
  theta = 0,
  scale = 1,
  linesize = 0.5,
  colour_soil = "#65571d",
  n = 3
){
  #tibble with all line segments to be annotated
  df <- tibble::tibble(
    xu = scale*c(seq(-0.5, 0.5 - 0.5/n, l = 2*n), seq(-0.25, 0.25, l = (n + 1))),
    yu = scale*c(rep(0, 2*n), rep(-0.25, (n + 1))),
    xuend = scale*c(seq(-0.25, 0 - 0.25/n, l = n), seq(0.25, 0.5 - 0.25/n, l = n), seq(0, 0.25 - 0.25/n, l = n), 0.5),
    yuend = -scale*0.25*c(rep(seq(1, 0 + 1/n, l = n), 2), seq(0, 1 - 1/n, l = n), 0),
    x = xc + .data$xu*cos(theta) - .data$yu*sin(theta),
    y = yc + .data$xu*sin(theta) + .data$yu*cos(theta),
    xend = xc + .data$xuend*cos(theta) - .data$yuend*sin(theta),
    yend = yc + .data$xuend*sin(theta) + .data$yuend*cos(theta),
  )
  #add to existing ggplot
  plt <- plt + ggplot2::annotate(
    "segment",
    x = df$x,
    xend = df$xend,
    y = df$y,
    yend = df$yend,
    size = linesize,
    color = colour_soil
  )
  #return
  return(plt)
}


#' Add a series of polygons to a ggplot object
#'
#' @description
#' Function to make it faster to add a series of soil, water or structure
#' polygons to an existing ggplot object
#'
#' @param plt ggplot to be added to
#' @param df tibble with data. Should contain fields `x` and `y` for
#'   polygon point positions. May contain optional field `group` to
#'   differentiate between different polygons
#' @param type either `soil`, `water` or `structure`. This setting is
#'   used to determine the fill colour of the polygons
#' @param fill if specified, overrides the fill colour defined through the
#'   `settings` settings
#' @param settings dataframe with default settings for soil, water and
#'   structure elements
#' @return ggplot object
#' @export

ggplot_add_polygon <- function(
  plt,
  df,
  type = "soil",
  fill = NULL,
  settings = tibble::tibble(
    id = c("soil", "saturated", "water", "structure"),
    fill = c("#d3bc5f" ,"#aebab7", "#2a7fff", "#333333"),
    colour = c(NA, NA, NA, NA),
  )
) {
  #global settings: soil - water - structure
  sett <- dplyr::filter(settings, .data$id == type)
  if (!is.null(fill)) {sett$fill <- fill}
  #add grouping if required
  if (!("group" %in% colnames(df))) {
    df$group <- 1
  }
  #add to plot
  plt + ggplot2::geom_polygon(
    data = df,
    ggplot2::aes(x = .data$x, y = .data$y, group = as.factor(.data$group)),
    fill = sett$fill,
    color = sett$colour,
    show.legend = FALSE
  )
}


#' Add surfaces to a ggplot object
#'
#' @description
#' Function to make it faster to add a series of soil, water or structure
#' surfaces to an existing ggplot object
#'
#' @param plt ggplot to be added to
#' @param df tibble with data. Should contain fields `x` and `y` for
#'   polygon point positions. May contain optional field `group` to
#'   differentiate between different polygons
#' @param type either `soil`, `water` or `structure`. This setting is
#'   used to determine the fill colour of the polygons
#' @param colour if specied, overrides the colour settings defined through the
#'   `settings` settings
#' @param fill if specified, overrides the fill colour defined through the
#'   `settings` settings
#' @param linetype if specied, overrides the linetype settings defined through
#'   the `settings` settings
#' @param linesize thickness of surface lines
#' @param linesize_marker thickness of lines in surface marker elements
#' @param settings dataframe with default settings for soil, water and
#'   structure elements
#' @param linesize line size thickness
#' @importFrom rlang .data
#' @return ggplot object
#' @export

ggplot_add_surface <- function(
  plt,
  df,
  type = "soil",
  colour = NULL,
  fill = NULL,
  linetype = NULL,
  linesize = 0.5,
  linesize_marker = 0.3,
  settings = tibble::tibble(
    id = c("soil", "water", "structure"),
    fill = c("#d3bc5f" ,"#2a7fff", "#333333"),
    colour = c("#65571d", "#000080", "#000000"),
    linetype = c(1, 2, 1),
    pos_marker = c(1/3, 2/3, 1/2),
    scale_marker = c(0.075, 0.050, 0.050)
  )
) {
  #global settings: soil - water - structure
  sett <- dplyr::filter(settings, .data$id == type)
  if (!is.null(colour)) {sett$colour <- colour}
  if (!is.null(fill)) {sett$fill <- fill}
  if (!is.null(linetype)) {sett$linetype <- linetype}
  #add grouping if required
  if (!("group" %in% colnames(df))) {
    df$group <- 1
  }
  #add lines to plot
  plt <- plt + ggplot2::geom_path(
    data = df,
    ggplot2::aes(x = .data$x, y = .data$y, group = as.factor(.data$group)),
    color = sett$colour,
    linetype = sett$linetype,
    size = linesize,
    show.legend = FALSE
  )
  #current y-axis limits
  ylim <- ggplot2::ggplot_build(plt)$layout$panel_scales_y[[1]]$range$range
  #add soil surface markers
  if (sett$scale_marker > 0) {
    dfm <- df %>%
      dplyr::group_by(.data$group) %>%
      dplyr::slice_head(n = 2) %>%
      dplyr::summarise(
        xc = .data$x[1] + sett$pos_marker*(.data$x[2] - .data$x[1]),
        yc = .data$y[1] + sett$pos_marker*(.data$y[2] - .data$y[1]),
        theta = atan2(.data$y[2] - .data$y[1], .data$x[2] - .data$x[1]),
        scale = sett$scale_marker*diff(ylim)
      )
    if (type == "soil") {
      for (i in 1:nrow(dfm)) {
        plt <- ggplot_add_soilmarker(
          plt,
          xc = dfm$xc[i],
          yc = dfm$yc[i],
          theta = dfm$theta[i],
          scale = dfm$scale[i],
          linesize = linesize_marker,
          colour_soil = sett$colour
        )
      }
    } else if (type == "water") {
      for (i in 1:nrow(dfm)) {
        plt <- ggplot_add_watermarker(
          plt,
          xc = dfm$xc[i],
          yc = dfm$yc[i],
          theta = dfm$theta[i],
          scale = dfm$scale[i],
          linesize = linesize_marker,
          fill_water = sett$fill,
          colour_water = sett$colour
        )
      }
    }
  }
  #return
  plt
}


#' ggplot a soil + water + structure geometry
#'
#' @description
#' Wrapper function to plot a soil, water and structure geometry
#'
#' @param pol_water tibble with surface water polygons
#' @param pol_soil tibble with soil polygons
#' @param pol_saturated tibble with saturated soil polygons
#' @param pol_structure tibble with structure polygons
#' @param surf_soil tibble with soil surfaces
#' @param surf_water tibble with water tables
#' @param xlim x-axis limits (array with min and max)
#' @param ylim y-axis limits (array with min and max)
#' @param axes if `TRUE`, plot axis ticks and labels. If `FALSE`, no axis
#'   are plotted
#' @param xlab x-axis label
#' @param ylab y-axis label
#' @return ggplot object
#' @examples
#' #input
#' pol_soil = data.frame(x = c(0, 0, 2, 2), y = c(0, 1, 1.1, 0), group = 0)
#' pol_saturated = data.frame(x = c(0, 0, 1, 1), y = c(0, 0.5, 0.5, 0))
#' pol_water = data.frame(x = c(0, 0, 2, 2), y = c(1, 1.5, 1.2, 1))
#' pol_structure = data.frame(x = c(0.5, 1, 1), y = c(0.5, 0.5, 0.75))
#' surf_soil = data.frame(x = c(0, 2), y = c(1, 1.1))
#' surf_water = data.frame(x = c(0, 2), y = c(1.5, 1.2))
#'
#' #plot
#' ggplot_geometry(
#'   pol_water = pol_water,
#'   pol_soil = pol_soil,
#'   pol_saturated = pol_saturated,
#'   pol_structure = pol_structure,
#'   surf_soil = surf_soil,
#'   surf_water = surf_water
#' )
#' @export

ggplot_geometry <- function(
  pol_water = NULL,
  pol_soil = NULL,
  pol_saturated = NULL,
  pol_structure = NULL,
  surf_soil = NULL,
  surf_water = NULL,
  xlim = c(NA, NA),
  ylim = c(NA, NA),
  axes = TRUE,
  xlab = "x [m]",
  ylab = "z [m]"
){
  #get plot limits
  if (length(unique(c(pol_soil$x, pol_saturated$x, pol_water$x, pol_structure$x))) > 1) {
    xlim <- round_limits(
      c(pol_soil$x, pol_saturated$x, pol_water$x, pol_structure$x),
      lower = xlim[1],
      upper = xlim[2]
    )
  }
  if (length(unique(c(pol_soil$y, pol_saturated$y, pol_water$y, pol_structure$y))) > 1) {
    ylim <- round_limits(
      c(pol_soil$y, pol_saturated$y, pol_water$y, pol_structure$y),
      lower = ylim[1],
      upper = ylim[2]
    )
  }
  #initiate plot
  plt <- ggplot2::ggplot() +
    soilmech::theme_soilmech() +
    ggplot2::theme(
      panel.background = ggplot2::element_blank(),
      panel.border = ggplot2::element_blank(),
    ) +
    ggplot2::coord_fixed(ratio = 1, xlim = xlim, ylim = ylim, expand = FALSE) +
    ggplot2::xlab(xlab) +
    ggplot2::ylab(ylab)
  #remove axes if needed
  if (axes == FALSE) {
    plt <- plt + ggplot2::theme(
      axis.text.x = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank(),
      axis.title.x = ggplot2::element_blank(),
      axis.title.y = ggplot2::element_blank()
      #panel.grid.major = ggplot2::element_blank(),
      #panel.grid.minor = ggplot2::element_blank()
    )
  }
  #add water polygons
  if (!is.null(pol_water)) {
    plt <- ggplot_add_polygon(plt, pol_water, type = "water")
  }
  #add soil polygons
  if (!is.null(pol_soil)) {
    plt <- ggplot_add_polygon(plt, pol_soil, type = "soil")
  }
  #add saturated soil polygons
  if (!is.null(pol_saturated)) {
    plt <- ggplot_add_polygon(plt, pol_saturated, type = "saturated")
  }
  #add structure polygons
  if (!is.null(pol_structure)) {
    plt <- ggplot_add_polygon(plt, pol_structure, type = "structure")
  }
  #plot water tables
  if (!is.null(surf_water)) {
    plt <- ggplot_add_surface(plt, surf_water, type = "water")
  }
  #plot soil surfaces
  if (!is.null(surf_soil)) {
    plt <- ggplot_add_surface(plt, surf_soil, type = "soil")
  }
  #return
  plt
}


#' Add arrows to ggplot
#'
#' @description
#' Add some arrow objects and annotation to an existing ggplot
#'
#' @param plt existing ggplot object
#' @param x,y x,y positions of arrow tips
#' @param theta orientation of arrows (in rad)
#' @param scale arrow length, in fraction of y-range
#' @param length arrow length in x,y coordinate system. If specified,
#'   overrides `scale` setting
#' @param label labels for each arrow (string). Labels are plotted near the
#'   tail of each arrow
#' @param label_size label text size
#' @param label_parse if `TRUE`, parse the label input
#' @param label_offset small offset to create distance between arrow tails.
#'   Expressed in fraction of y-range
#' @param linewidth arrow line width
#' @param colour arrow colour
#' @param scale_arrow length of arrow head, as fraction of arrow length
#' @param scale_arrow_max maximum arrow head size, as fraction of y-range
#' @examples
#' #generate a plot
#' plt <- ggplot2::ggplot() +
#'   ggplot2::geom_rect(ggplot2::aes(xmin = 0, xmax = 2, ymin = 0, ymax = 2), fill = "yellow") +
#'   ggplot2::coord_fixed(ratio = 1)
#'
#' #arrow
#' x = c(1, 1, 1, 1)
#' y = c(1, 1, 1, 1)
#' theta = c(-pi/2, 0, pi/2, pi)
#' label = "test*theta"
#' ggplot_add_arrow(plt, x, y, theta = theta, label = label, label_parse = TRUE)
#' @export


ggplot_add_arrow <- function(
  plt,
  x,
  y,
  theta = -pi/2,
  scale = 0.1,
  length = NULL,
  label = NULL,
  label_size = 3.5,
  label_parse = FALSE,
  label_offset = 0.02,
  linewidth = 0.5,
  colour = "black",
  scale_arrow = 0.3,
  scale_arrow_max = 0.1
) {
  #get current plot y-axis limits
  ylim <- ggplot2::ggplot_build(plt)$layout$panel_scales_y[[1]]$range$range
  #if arrow length not defined, determine from scale
  if (is.null(length)) {
    if (is.double(scale)) {
      length <- scale*diff(ylim)
    } else {
      length <- 1
    }
  }
  scale <- length/diff(ylim)
  #create dataframe with begin and end positions
  df <- tibble::tibble(
    xend = x,
    yend = y,
    theta = theta %% (2*pi),
    x = xend - length*cos(.data$theta),
    y = yend - length*sin(.data$theta)
  )
  #add to ggplot
  plt <- plt + ggplot2::geom_segment(
    data = df,
    ggplot2::aes(x = .data$x, y = .data$y, xend = .data$xend, yend = .data$yend),
    color = colour,
    size = linewidth,
    arrow = ggplot2::arrow(length = ggplot2::unit(min(scale_arrow_max, scale*scale_arrow), "npc"))
  )
  #add labels
  if (!is.null(label)) {
    df <- df %>% dplyr::mutate(
      label = label,
      xlabel = x,
      ylabel = y,
      hjust = 0.5,
      vjust = 0.5
    )
    #bottom
    i2 <- (df$theta >= pi/4) & (df$theta < 3*pi/4)
    df$vjust[i2] <- 1
    df$ylabel[i2] < df$y[i2] - label_offset*diff(ylim)
    #right
    i1 <- (df$theta >= 3*pi/4) & (df$theta < 5*pi/4)
    df$hjust[i1] <- 0
    df$xlabel[i1] <- df$x[i1] + label_offset*diff(ylim)
    #top
    i4 <- (df$theta >= 5*pi/4) & (df$theta < 7*pi/4)
    df$vjust[i4] <- 0
    df$ylabel[i4] <- df$y[i4] + label_offset*diff(ylim)
    #right
    i3 <- as.logical(1 - i2 - i1 - i4)
    df$hjust[i3] <- 1
    df$xlabel[i3] <- df$x[i3] - label_offset*diff(ylim)
    #add text
    plt <- plt + ggplot2::geom_text(
      data = df,
      ggplot2::aes(
        x = .data$xlabel,
        y = .data$ylabel,
        label = .data$label,
        hjust = .data$hjust,
        vjust = .data$vjust
      ),
      color = colour,
      size = label_size,
      parse = label_parse
    )
  }
  #return
  return(plt)
}
