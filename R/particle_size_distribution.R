#' Plot empty particle size distribution chart
#'
#' @description
#' Function to plot an empty particle size distribution chart. Default
#' plots the range from clay all the way to boulders
#'
#' @param trace_diameter array with diameters for plotting traces
#' @param trace_passing array with fractions passing for plotting traces
#' @param trace_group array with grouping for plotting traces
#' @param crosshairs_diameter diameter values of crosshairs to add
#' @param crosshairs_passing fraction passing of crosshairs to add
#' @param crosshairs_group crosshairs grouping + value for labels
#' @param crosshairs_parse if `TRUE`, parse `crosshairs_group` string in label
#' @param palette RColorBrewer palette for plotting traces
#' @param legend_position position for traces legend. names of traces are
#'   taken from the field `group`
#' @param legend_title title of the legend for traces
#' @param add_points if `TRUE`, all points are plotted in addition to the
#'   lines of the traces
#' @param height_labelbar height of single label bar above plot, relative
#'   to the height of the plot
#' @param percentage if `TRUE`, plot percentages rather than fractions
#' @param box_thickness thickness of label boxes
#' @param textsize_primary size of labels, primary fraction
#' @param textsize_secondary size of labels, fine/medium/course
#' @param xlim user-defined x-axis limits
#' @importFrom magrittr `%>%`
#' @importFrom rlang .data
#' @return a ggplot object
#' @examples
#' #empty chart
#' ggplot_psdchart()
#'
#' #chart with two trace
#' ggplot_psdchart(
#'   trace_diameter = c(0.01, 1, 100, 0.1, 1, 100),
#'   trace_passing = c(0.1, 0.5, 0.9, 0.2, 0.6, 0.9),
#'   trace_group = c(rep("A", 3), rep("B", 3)),
#'   legend_position = "right"
#' )
#'
#' #chart with one trace and some annotations
#' ggplot_psdchart(
#'   trace_diameter = c(0.001, 0.01, 0.1, 1, 10, 100),
#'   trace_passing = c(0.05, 0.1, 0.2, 0.5, 0.6, 0.9),
#'   crosshairs_diameter = c(0.01, 1),
#'   crosshairs_passing = c(0.1, 0.5),
#'   crosshairs_group = c("D[10]==0.01~mm", "D[50]==1~mm"),
#'   crosshairs_parse = TRUE
#' )
#' @export

ggplot_psdchart <- function(
  trace_diameter = NULL,
  trace_passing = NULL,
  trace_group = NULL,
  crosshairs_diameter = NULL,
  crosshairs_passing = NULL,
  crosshairs_group = NULL,
  crosshairs_parse = NULL,
  palette = "Set1",
  legend_position = "none",
  legend_title = "",
  add_points = TRUE,
  height_labelbar = 0.075,
  percentage = TRUE,
  box_thickness = 0.25,
  textsize_primary = 2.5,
  textsize_secondary = 2,
  xlim = c(NA, NA)
){
  #Particle size distribution classes
  dc <- tibble::tibble(
    type = c("CLAY", "SILT", "SILT", "SILT", "SAND", "SAND", "SAND", "GRAVEL", "GRAVEL", "GRAVEL", "COBBLES", "BOULDERS"),
    subtype = c(NA, "Fine", "Medium", "Coarse", "Fine", "Medium", "Coarse", "Fine", "Medium", "Coarse", NA, NA),
    dmin = c(4e-05, 0.002, 0.006, 0.02, 0.06, 0.2, 0.6, 2, 6, 20, 60, 200),
    dmax = c(0.002, 0.006, 0.02, 0.06, 0.2, 0.6, 2, 6, 20, 60, 200, 2e+03)
  )
  #limit based on input - take subset of classes
  if (any(is.na(xlim))){
    xlim <- c(min(dc$dmin), max(dc$dmax))
  } else {
    dc <- dc[(dc$dmax > xlim[1]) & (dc$dmin < xlim[2]), ]
    dc$dmin <- pmax(dc$dmin, xlim[1])
    dc$dmax <- pmin(dc$dmax, xlim[2])
  }
  #add plot data to classes - major classes labels
  dc1 <- dc %>%
    dplyr::group_by(.data$type) %>%
    dplyr::summarise(
      dmin = min(.data$dmin),
      dmax = max(.data$dmax),
      ymin = ifelse(any(is.na(.data$subtype)), 1, 1 + height_labelbar)
    ) %>%
    dplyr::mutate(
      ymax = 1 + 2*height_labelbar,
      dlabel = 10^(0.5*(log10(.data$dmin) + log10(.data$dmax))),
      ylabel = 0.5*(.data$ymin + .data$ymax)
    )
  #add plot data to classes - minor classes labels
  dc2 <- dc %>%
    dplyr::filter(!is.na(.data$subtype)) %>%
    dplyr::group_by(.data$type, .data$subtype) %>%
    dplyr::summarise(
      dmin = min(.data$dmin),
      dmax = max(.data$dmax)
    ) %>%
    dplyr::mutate(
      ymin = 1,
      ymax = 1 + 1*height_labelbar,
      dlabel = 10^(0.5*(log10(.data$dmin) + log10(.data$dmax))),
      ylabel = 0.5*(.data$ymin + .data$ymax)
    )
  #Create multiplier if percentages are used
  if (percentage == TRUE){
    mult <- 100
  } else {
    mult <- 1
  }
  #plot empty plot
  plt <- ggplot2::ggplot() +
    theme_soilmech() +
    ggplot2::geom_rect(
      data = dc1,
      ggplot2::aes(xmin = .data$dmin, xmax = .data$dmax, ymin = mult*.data$ymin, ymax = mult*.data$ymax),
      fill = "white",
      color = "black",
      size = box_thickness
    ) +
    ggplot2::geom_rect(
      data = dc2,
      ggplot2::aes(xmin = .data$dmin, xmax = .data$dmax, ymin = mult*.data$ymin, ymax = mult*.data$ymax),
      fill = "white",
      color = "black",
      size = box_thickness
    ) +
    ggplot2::geom_text(
      data = dc1,
      ggplot2::aes(x = .data$dlabel, y = mult*.data$ylabel, label = .data$type),
      hjust = 0.5,
      vjust = 0.5,
      size = textsize_primary
    ) +
    ggplot2::geom_text(
      data = dc2,
      ggplot2::aes(x = .data$dlabel, y = mult*.data$ylabel, label = .data$subtype),
      hjust = 0.5,
      vjust = 0.5,
      size = textsize_secondary
    ) +
    ggplot2::scale_y_continuous(
      limits = c(0, mult*(1 + 2*height_labelbar)),
      breaks = mult*seq(0, 1, 0.2),
      expand = c(0, 0)
    ) +
    ggplot2::scale_x_log10(
      limits = xlim,
      breaks = scales::trans_breaks("log10", function(x) 10^x),
      labels = scales::trans_format("log10", scales::math_format(expr = ~10^.x)),
      minor_breaks = get_log10_minorbreaks(xlim),
      expand = c(0, 0)
    ) +
    ggplot2::annotation_logticks(side = "b")
  #add labels
  plt <- plt + ggplot2::xlab('Particle size [mm]')
  if (percentage == TRUE){
    plt <- plt + ggplot2::ylab('Percentage passing [%]')
  } else {
    plt <- plt + ggplot2::ylab('Fraction passing [-]')
  }
  #add traces - if specified in input
  if (!is.null(trace_diameter) & !is.null(trace_passing)) {
    dt <- tibble::tibble(
      diameter = trace_diameter,
      passing = mult*trace_passing,
    )
    if (is.null(trace_group)) {
      plt <- plt +
        ggplot2::geom_line(
          data = dt,
          ggplot2::aes(x = .data$diameter, y = .data$passing)
        )
      if (add_points == TRUE) {
        plt <- plt +
          ggplot2::geom_point(
            data = dt,
            ggplot2::aes(x = .data$diameter, y = .data$passing)
          )
      }
    } else {
      dt$group <- trace_group
      plt <- plt +
        ggplot2::geom_line(
          data = dt,
          ggplot2::aes(x = .data$diameter, y = .data$passing, color = as.factor(.data$group))
        ) +
        ggplot2::theme(legend.position = legend_position) +
        ggplot2::scale_color_brewer(name = legend_title, palette = palette)
      if (add_points == TRUE) {
        plt <- plt +
          ggplot2::geom_point(
            data = dt,
            ggplot2::aes(x = .data$diameter, y = .data$passing, color = as.factor(.data$group))
          )
      }
    }
  }
  #add crosshairs - if specified in input
  if (!is.null(crosshairs_diameter) & !is.null(crosshairs_passing)) {
    dc <- tibble::tibble(
      diameter = crosshairs_diameter,
      passing = mult*crosshairs_passing,
      group = crosshairs_group
    )
    plt <- ggplot_addcrosshairs(
      plt,
      x = dc$diameter,
      y = dc$passing,
      group = dc$group,
      add_colour_scale = TRUE,
      xlim = xlim[1],
      arrow_y = FALSE,
      label_parse = crosshairs_parse
    )
  }
  #return
  plt
}
