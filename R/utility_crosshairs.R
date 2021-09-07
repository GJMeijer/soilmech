#' Add datapoints/crosshairs to an existing plot
#'
#' @description
#' Function takes an existing ggplot and add points and/or crosshairs to the
#' plot to show the position of points
#'
#' @param plt a ggplot object
#' @param x x-coordinates for all points
#' @param y y-coordinates for all points
#' @param group names of group, for all points - this is plotted in the label
#' @param legend_title title for legend
#' @param legend_position legend position. If `NULL`, no legend is drawn
#' @param crosshairs if `TRUE`, plot crosshairs for each data point
#' @param add_colour_scale if `TRUE`, add new color scale for different crosshairs
#' @param xlim start point of x crosshairs
#' @param ylim start point of y crosshairs
#' @param arrow_x,arrow_y if `TRUE`, point x or y-arrows towards the point.
#'   If `FALSE`, points towards the axis
#' @param arrow_length length of arrow heads
#' @param palette RColorBrewer color palette to use
#' @param label_size text size in label
#' @param label_hjust horizontal justification of labels
#' @param label_vjust vertical justification of labels
#' @param label_nudge_x x-offset of labels
#' @param label_nudge_y y-offset of labels
#' @param label_parse if `TRUE`, parse label text
#' @importFrom rlang .data
#' @export

ggplot_addcrosshairs <- function(
  plt,
  x,
  y,
  group = NULL,
  legend_title = NULL,
  legend_position = "none",
  crosshairs = TRUE,
  add_colour_scale = FALSE,
  xlim = 0,
  ylim = 0,
  arrow_x = TRUE,
  arrow_y = TRUE,
  arrow_length = 0.1,
  palette = "Set1",
  label_size = 4,
  label_hjust = 0,
  label_vjust = 0,
  label_nudge_x = 0,
  label_nudge_y = 0,
  label_parse = FALSE
){
  #create tibble with all points
  dp <- tibble::tibble(x = x, y = y)
  if (!is.null(group)) {
    dp$group <- group
  } else {
    dp$group <- as.character(seq(length(x)))
  }
  #create line segments
  if (crosshairs == TRUE) {
    if (arrow_x == TRUE){
      xx <- c(rep(xlim, nrow(dp)), dp$x, dp$x, dp$x)
    } else {
      xx <- c(dp$x, rep(xlim, nrow(dp)), dp$x, dp$x)
    }
    if (arrow_y == TRUE){
      yy <- c(dp$y, dp$y, rep(ylim, nrow(dp)), dp$y)
    } else {
      yy <- c(dp$y, dp$y, dp$y, rep(ylim, nrow(dp)))
    }
    dl <- tibble::tibble(
      x = xx,
      y = yy,
      group = rep(dp$group, 4),
      group_segment = rep(c(1, 2), each = 2*nrow(dp)),
      group_plot = paste0(.data$group, "-", .data$group_segment)
    )
    #add to plot
    if (add_colour_scale == TRUE) {
      plt <- plt +
        ggplot2::geom_path(
          data = dl,
          ggplot2::aes(x = .data$x, y = .data$y, color = as.factor(.data$group), group = .data$group_plot),
          arrow = ggplot2::arrow(length = ggplot2::unit(arrow_length, "inches")),
          show.legend = FALSE
        )
    } else {
      plt <- plt +
        ggplot2::geom_path(
          data = dl,
          ggplot2::aes(x = .data$x, y = .data$y, group = .data$group_plot),
          arrow = ggplot2::arrow(length = ggplot2::unit(arrow_length, "inches")),
          show.legend = FALSE
        )
    }
  }
  #add labels
  if (!is.null(group)) {
    if (add_colour_scale == TRUE) {
      plt <- plt +
        ggplot2::geom_label(
          data = dp,
          ggplot2::aes(x = .data$x, y = .data$y, label = .data$group, color = as.factor(.data$group)),
          hjust = label_hjust,
          vjust = label_vjust,
          nudge_x = label_nudge_x,
          nudge_y = label_nudge_y,
          parse = label_parse,
          label.padding = ggplot2::unit(0.1, "lines"),
          show.legend = FALSE
        )
    } else {
      plt <- plt +
        ggplot2::geom_label(
          data = dp,
          ggplot2::aes(x = .data$x, y = .data$y, label = .data$group, group = as.factor(.data$group)),
          hjust = label_hjust,
          vjust = label_vjust,
          nudge_x = label_nudge_x,
          nudge_y = label_nudge_y,
          parse = label_parse,
          label.padding = ggplot2::unit(0.1, "lines"),
          show.legend = FALSE
        )
    }
  }
  #add points
  if (add_colour_scale == TRUE) {
    plt <- plt +
      ggplot2::geom_point(
        data = dp,
        ggplot2::aes(x = .data$x, y = .data$y, color = as.factor(.data$group))
      ) +
      ggplot2::scale_color_brewer(name = legend_title, palette = palette) +
      ggplot2::scale_linetype_discrete(name = legend_title) +
      ggplot2::scale_shape_discrete(name = legend_title)
  } else {
    plt <- plt +
      ggplot2::geom_point(
        data = dp,
        ggplot2::aes(x = .data$x, y = .data$y, group = as.factor(.data$group))
      )
  }
  #add legend
  plt <- plt + ggplot2::theme(legend.position = legend_position)
  #return
  return(plt)
}
