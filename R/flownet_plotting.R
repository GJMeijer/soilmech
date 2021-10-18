#' Decide on number of equipotential intervals
#'
#' @description
#' This function takes a finite difference solution for both head `h` and flow
#' potential function `psi` to get the best number of equipotential intervals
#' `Nd` based on a given number of flow channels `Nf`. This function assumes
#' the permeabilities in all domains is the same
#'
#' @param df tibble with information about domain
#' @param dp tibble with head and flow potential results for each node
#' @param Nf number of requested flow channels
#' @return the optimal number of equipotential intervals `Nd`
#' @importFrom magrittr `%>%`
#' @importFrom rlang .data
#' @examples
#' df <- flownet_geometry_rectangular()
#' dp <- flownet_solve_rectangular(df)
#' amount_equipotentialintervals(df, dp, Nf = 5)
#' @export

amount_equipotentialintervals <- function(df, dp, Nf = 5) {
  #total head difference
  dh <- max(dp$h) - min(dp$h)
  #total flow per second
  Q <- max(dp$psi) - min(dp$psi)
  #k_avg
  k_avg <- df$dom %>%
    dplyr::summarize(k = mean(sqrt(.data$kx*.data$ky))) %>%
    dplyr::pull(.data$k)
  #Nd
  Nd <- round(k_avg*dh*Nf/Q)
  #return
  return(Nd)
}


#' Check if series of gridded domains is rectangular
#'
#' @description
#' Check if a flownet solution is given on a rectangular grid, i.e. nodes
#' are organised in horizontal rows and vertical columns. The distances
#' between rows and columns may vary.
#'
#' For each domain in the grid, the function checks if the number of nodes
#' is equal to the product of unique x and y-positions. If this is true for
#' all domains, the grid is regular and the function returns `TRUE`.
#'
#' @param dp tibble with flow net solution. Should contain fields `x` and `y`
#'   for positions of the nodes
#' @return `TRUE` is the grid is rectangular, `FALSE` is not
#' @importFrom magrittr `%>%`
#' @importFrom rlang .data
#' @examples
#' #rectangular grid
#' df <- flownet_geometry_rectangular()
#' dp <- flownet_solve_rectangular(df)
#' is_grid_rectangular(dp)
#'
#' #quadrilateral grid
#' df <- flownet_geometry_quadrilateral()
#' dp <- flownet_solve_quadrilateral(df)
#' is_grid_rectangular(dp)
#' @export

is_grid_rectangular <- function(dp) {
  return(
    dp %>%
      dplyr::group_by(.data$domain) %>%
      dplyr::summarize(
        n = dplyr::n(),
        nx = length(unique(.data$x)),
        ny = length(unique(.data$y)),
      ) %>%
      dplyr::mutate(rect = (.data$n == .data$nx*.data$ny)) %>%
      dplyr::pull(.data$rect) %>%
      all()
  )
}

#' Add flownet to an existing ggplot
#'
#' @description
#' Function to add a flow net to an existing ggplot. If the plot does not
#' exist yet, a new one if created
#'
#' @param df tibble with domains/problem description
#' @param dp flow net solution tibble
#' @param di flow net solution rectangular interpolation tibble. This tibble
#'   is used for plotting contour lines (flow lines and equipotential lines).
#'   If the grid is rectangular, this is simply equal to `dp`. If not, an
#'   interpolation function is used. If not defined, it is calculated from
#'   `dp`
#' @param Nf number of requested flow paths
#' @param linewidth line thickness for all lines
#' @param colour_flowline colour of flow lines
#' @param colour_equipotentialline colour of equipotential lines
#' @param fill fill colour for (saturated) soil
#' @param plt a ggplot object. If not defined, a new plot is generated
#' @param xlab,ylab x and y-axis label - if `plt` not defined
#' @param xlim,ylim user-defined axis limits - if `plt` not defined
#' @param axes if `FALSE` no axes are plotted - if `plt` not defined
#' @return a ggplot object
#' @importFrom magrittr `%>%`
#' @importFrom rlang .data
#' @examples
#' #rectangular grid
#' df <- flownet_geometry_rectangular()
#' dp <- flownet_solve_rectangular(df)
#' ggplot_add_flownet(df, dp, axes = FALSE)
#' ggplot_add_flownet(df, dp, axes = TRUE, fill = NA)
#'
#' #quadrilateral grid
#' df <- flownet_geometry_quadrilateral()
#' dp <- flownet_solve_quadrilateral(df)
#' di <- flownet_interpolate_quadrilateral(df, dp)
#' ggplot_add_flownet(df, dp, di = di, axes = FALSE)
#' @export

ggplot_add_flownet <- function(
  df,
  dp,
  di = NULL,
  Nf = 5,
  linewidth = 0.5,
  colour_flowline = "black",
  colour_equipotentialline = "black",
  fill = "#aebab7",
  plt = NULL,
  xlab = "x [m]",
  ylab = "z [m]",
  xlim = c(NA, NA),
  ylim = c(NA, NA),
  axes = TRUE
){
  #check if grid rectangular
  grid_rectangular <- is_grid_rectangular(dp)
  #interpolate rectangular grid if required - and get soil polygons
  if (is.null(di)) {
    if (grid_rectangular == TRUE) {
      di <- dp
    } else {
      di <- flownet_interpolate_quadrilateral(df, dp)
    }
  }
  #if plot does not exist, generate a new one
  if (is.null(plt)) {
    plt <- ggplot_geometry(xlim = xlim, ylim = ylim, axes = axes, xlab = xlab, ylab = ylab)
  }
  #calculate number of equamoutipotential intervals
  Nd <- amount_equipotentialintervals(df, dp, Nf = Nf)
  Nd_levels <- seq(min(dp$h), max(dp$h), l = (Nd + 1))[2:Nd]
  Nf_levels <- seq(min(dp$psi), max(dp$psi), l = (Nf + 1))
  #plot soil polygons
  if (is.character(fill)) {
    if (grid_rectangular == TRUE) {
      dpol <- df$dom %>%
        dplyr::mutate(
          x4 = .data$x1,
          x3 = .data$x1,
          x2 = .data$x0,
          x1 = .data$x0,
          y4 = .data$y0,
          y3 = .data$y1,
          y2 = .data$y1,
          y1 = .data$y0
        ) %>%
        tidyr::pivot_longer(
          cols = c("x1", "x2", "x3", "x4", "y1", "y2", "y3", "y4"),
          names_to = c(".value", "set"),
          names_pattern = "(.+)(.+)"
        )
    } else {
      dpol <- tidyr::pivot_longer(
        df$dom,
        cols = c("x1", "x2", "x3", "x4", "y1", "y2", "y3", "y4"),
        names_to = c(".value", "set"),
        names_pattern = "(.+)(.+)"
      )
    }
    plt <- plt + ggplot2::geom_polygon(
      data = dpol,
      ggplot2::aes(x = .data$x, y = .data$y, group = as.factor(.data$domain)),
      fill = fill,
      color = NA
    )
  }
  #add flow and equipotential lines using contour lines
  if (is.character(colour_equipotentialline)) {
    plt <- plt + ggplot2::geom_contour(
      data = di,
      ggplot2::aes(x = .data$x, y = .data$y, z = .data$h),
      breaks = Nd_levels,
      color = colour_equipotentialline,
      size = linewidth
    )
  }
  if (is.character(colour_flowline)) {
    plt <- plt + ggplot2::geom_contour(
      data = di,
      ggplot2::aes(x = .data$x, y = .data$y, z = .data$psi),
      breaks = Nf_levels,
      color = colour_flowline,
      size = linewidth
    )
  }
  #return
  return(plt)
}


#' Create polygons from gridded points
#'
#' @description
#' Get the node indices of the four corner points of a series of grids
#'
#' @param df list with flownet problem definition
#' @param dp list with flownet results, per node. This dataframe should
#'   contain a field with the name `val` that is used to calculate a
#'   value for each polygon
#' @return tibble with polygon data. Unique polygons are differentiated
#'   using the unique identifier in the field `id`. Interpolated values
#'   are stored in field `val`, and the positions of the corner points
#'   in `x` and `y`
#' @importFrom magrittr `%>%`
#' @importFrom rlang .data
#' @examples
#' #calculate
#' df <- flownet_geometry_rectangular(grid_size = 2)
#' dp <- flownet_solve_rectangular(df)
#' dpol <- grid2polygon(df, dp)
#'
#' #plot
#' ggplot2::ggplot() +
#'   ggplot2::geom_polygon(
#'     data = dpol,
#'     ggplot2::aes(x = x, y = y, fill = h, group = id),
#'     color = NA
#'   )
#' @export

grid2polygon <- function(df, dp){
  #get indices of polygons
  dpol <- df$dom %>%
    dplyr::mutate(i0_real = c(0, utils::head(cumsum(.data$nx*.data$ny), -1))) %>%
    dplyr::select(.data$nx, .data$ny, .data$i0_real) %>%
    purrr::pmap_dfr(
      function(nx, ny, i0_real) {
        tibble::tibble(
          i1 = i0_real + rep(seq(1, nx - 1), ny - 1) + nx*rep(seq(0, ny - 2), each = nx - 1),
          i2 = i0_real + rep(seq(1, nx - 1), ny - 1) + nx*rep(seq(1, ny - 1), each = nx - 1),
          i3 = i0_real + rep(seq(2, nx), ny - 1) + nx*rep(seq(1, ny - 1), each = nx - 1),
          i4 = i0_real + rep(seq(2, nx), ny - 1) + nx*rep(seq(0, ny - 2), each = nx - 1)
        )
      }
    ) %>%
    #get properties and positions
    dplyr::mutate(
      id = seq(dplyr::n()),
      x1 = dp$x[.data$i1],
      x2 = dp$x[.data$i2],
      x3 = dp$x[.data$i3],
      x4 = dp$x[.data$i4],
      y1 = dp$y[.data$i1],
      y2 = dp$y[.data$i2],
      y3 = dp$y[.data$i3],
      y4 = dp$y[.data$i4]
    )
  #average properties per field
  if ("h" %in% colnames(dp)) {
    dpol$h <- (dp$h[dpol$i1] + dp$h[dpol$i2] + dp$h[dpol$i3] + dp$h[dpol$i4])/4
  }
  if ("qx" %in% colnames(dp)) {
    dpol$qx <- (dp$qx[dpol$i1] + dp$qx[dpol$i2] + dp$qx[dpol$i3] + dp$qx[dpol$i4])/4
  }
  if ("qy" %in% colnames(dp)) {
    dpol$qy <- (dp$qy[dpol$i1] + dp$qy[dpol$i2] + dp$qy[dpol$i3] + dp$qy[dpol$i4])/4
  }
  if ("psi" %in% colnames(dp)) {
    dpol$psi <- (dp$psi[dpol$i1] + dp$psi[dpol$i2] + dp$psi[dpol$i3] + dp$psi[dpol$i4])/4
  }
  if ("val" %in% colnames(dp)) {
    dpol$val <- (dp$val[dpol$i1] + dp$val[dpol$i2] + dp$val[dpol$i3] + dp$val[dpol$i4])/4
  }
  #convert to long form and return
  return(
    dpol %>%
      tidyr::pivot_longer(
        cols = c("x1", "x2", "x3", "x4", "y1", "y2", "y3", "y4"),
        names_to = c(".value", "set"),
        names_pattern = "(.+)(.+)"
      )
  )
}


#' ggplot 2D map of hydraulic heads and pressures
#'
#' @description
#' Plots a 2D map of the pressure or head distributions in a flow net problem.
#' The plot shows both a heat map and contour lines. Can be added to an existing
#' ggplot if required
#'
#' @inheritParams ggplot_add_flownet
#' @param type parameter to plot. `type = "h"` for hydraulic heads,
#'   `type = "hz"` for elevation head, `type = `hb` for pressure head
#'   or `type = "u"` for pore water pressures
#' @param binwidth contour line distance
#' @param label_size contour line text label size
#' @param label_n number of times contour label appears on each line
#' @param gamma_w unit weight of water, in kN/m3 (required for calculating
#'   pore water pressures)
#' @param palette RColorBrewer color map for shading
#' @param palette_direction RColorBrewer color map direction (`1` or `-1`)
#' @param legend_position position of the legend for colours. use
#'   `legend_position = "none"` to disable the legend
#' @param linewidth contour line thickness
#' @param colour_contour colour of contour lines
#' @return a ggplot object
#' @examples
#' #rectangular flow net
#' df <- flownet_geometry_rectangular()
#' dp <- flownet_solve_rectangular(df)
#' #plot heads
#' ggplot_add_hydraulichead(
#'   df,
#'   dp,
#'   type = "h",
#'   binwidth = 0.25
#' )
#' #plot elevation head
#' ggplot_add_hydraulichead(
#'   df,
#'   dp = dp,
#'   type = "hz",
#'   binwidth = 2.5
#' )
#' #plot pressure head
#' ggplot_add_hydraulichead(
#'   df,
#'   dp = dp,
#'   type = "hb",
#'   binwidth = 2.5
#' )
#' #plot pore water pressure
#' ggplot_add_hydraulichead(
#'   df,
#'   dp = dp,
#'   type = "u",
#'   binwidth = 25
#' )
#'
#' #quadrilateral flow net
#' df <- flownet_geometry_quadrilateral()
#' dp <- flownet_solve_quadrilateral(df)
#' #plot heads
#' ggplot_add_hydraulichead(
#'   df,
#'   dp,
#'   type = "h",
#'   binwidth = 0.5
#' )
#' @export

ggplot_add_hydraulichead <- function(
  df,
  dp,
  di = NULL,
  type = "h",  #h, hz, hb, u
  linewidth = 0.5,
  colour_contour = "black",
  binwidth = 0.25,
  gamma_w = 10,
  label_size = 3,
  label_n = 1,
  palette = "Blues",  #YlGnBu
  palette_direction = 1,
  legend_position = "right",
  plt = NULL,
  xlab = "x [m]",
  ylab = "z [m]",
  xlim = c(NA, NA),
  ylim = c(NA, NA),
  axes = TRUE
){
  #check if grid rectangular
  grid_rectangular <- is_grid_rectangular(dp)
  #if plot does not exist, generate a new one
  if (is.null(plt)) {
    plt <- ggplot_geometry(xlim = xlim, ylim = ylim, axes = axes, xlab = xlab, ylab = ylab)
  }
  #interpolate rectangular grid if required - and get soil polygons
  if (is.null(di)) {
    if (grid_rectangular == TRUE) {
      di <- dp
    } else {
      di <- flownet_interpolate_quadrilateral(df, dp)
    }
  }
  #select hydraulic head profiles
  if (type == "h") {
    dp$val <- dp$h
    di$val <- di$h
    label <- expression(h~"[m]")
  } else if (type == "hz") {
    dp$val <- dp$y
    di$val <- di$y
    label <- expression(h[z]~"[m]")
  } else if (type == "hb") {
    dp$val <- dp$h - dp$y
    di$val <- di$h - di$y
    label <- expression(h[b]~"[m]")
  } else if (type == "u") {
    dp$val <- (dp$h - dp$y)*gamma_w
    di$val <- (di$h - di$y)*gamma_w
    label <- expression(u~"[kPa]")
  }
  #plot colourmap
  if (is.character(palette) == TRUE) {
    #generate polygons with averages of the four corner points
    dpol <- grid2polygon(df, dp)
    #add to plot
    plt <- plt + ggplot2::geom_polygon(
      data = dpol,
      ggplot2::aes(x = .data$x, y = .data$y, group = .data$id, fill = .data$val)
    ) +
      ggplot2::scale_fill_distiller(
        name = label,
        palette = palette,
        direction = palette_direction
      ) +
      ggplot2::theme(legend.position = legend_position)
  }
  #plot contourlines
  if (is.double(binwidth) == TRUE) {
    plt <- plt + ggplot2::geom_contour(
      data = di,
      ggplot2::aes(x = .data$x, y = .data$y, z = .data$val),
      binwidth = binwidth,
      color = colour_contour,
      size = linewidth
    ) +
      metR::geom_text_contour(
        data = di,
        binwidth = binwidth,
        ggplot2::aes(x = .data$x, y = .data$y, z = .data$val),
        label.placement = metR::label_placement_n(label_n),
        size = label_size
      )
  }
  #return
  return(plt)
}
