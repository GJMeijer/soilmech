#' Double-drained consolidation at point in time and space
#'
#' @description
#' Calculates the degree of consolidation Uv in a double-drained
#' soil layer at a certain (normalised) times: Tv = c_v*t/d^2 and
#' at certain (normalised) position in the layer (zd = z/d), where
#' 'd' is half the thickness of the open layer
#'
#' Initial excess pore pressure is assumed to be uniform across the
#' entire soil layer
#'
#' Solution uses the fourier series method to solve the consolidation
#' equation, detailed in Craig's Soil Mechanics, among others
#'
#' @param zd normalised depth (array)
#' @param Tv normalised time (array)
#' @param nfourier number of fourier steps to use
#' @importFrom magrittr `%>%`
#' @importFrom rlang .data
#' @return array with degrees of consolidation (between 0 and 1)
#' @export

calculate_consolidation_openlayer <- function(
  zd,
  Tv,
  nfourier = 50
){
  #solution of consolidation equation using Fourier expansion
  tibble::tibble(
    zd = zd,
    Tv = Tv
  ) %>%
    tidyr::expand_grid(m = seq(nfourier)) %>%
    dplyr::mutate(
      M = pi/2*(2*.data$m + 1),
      ue = 2/.data$M*sin(.data$M*.data$zd)*exp(-.data$M^2*.data$Tv)) %>%
    dplyr::group_by(.data$zd, .data$Tv) %>%
    dplyr::summarize(Uv = 1 - sum(.data$ue)) %>%
    dplyr::pull(.data$Uv)
}


#' ggplot consolidation in double-drained layer in time and space
#'
#' @description
#' Plots a ggplot with development of consolidation in a double-drained
#' layer, in space and time.
#'
#' @param zd normalised depth for crosshair annotations
#' @param Tv normalised time for crosshair annotation
#' @param zd_calc normalised depth used to plot lines
#' @param Tv_calc normalised time steps for which to plot lines
#' @param label_zd depth position at which to plot time labels on traces
#' @param label_size size of the text in time labels
#' @param label_alpha opacity of the time labels
#' @param nround number of decimals to use in rounding values in crosshair
#'   annotations
#' @importFrom magrittr `%>%`
#' @importFrom rlang .data
#' @return a ggplot object
#' @examples
#' #empty chart
#' ggplot_consolidation_doubledrained()
#'
#' #annotated chart
#' ggplot_consolidation_doubledrained(zd = 0.6, Tv = 0.3)
#' @export

ggplot_consolidation_doubledrained <- function(
  zd = NULL,
  Tv = NULL,
  zd_calc = seq(0, 2, l = 101),
  Tv_calc = c(seq(0.05, 0.20, 0.05), seq(0.3, 0.9, 0.1)),
  label_zd = 0.7,
  label_size = 3,
  label_alpha = 1,
  nround = 2
){
  #all combinations
  df <- tidyr::expand_grid(
    zd = zd_calc,
    Tv = Tv_calc
  ) %>%
    dplyr::mutate(
      Uv = calculate_consolidation_openlayer(.data$zd, .data$Tv)
    )
  #add zero time line
  df <- dplyr::bind_rows(
    df,
    tibble::tibble(zd = c(0, 2), Tv = 0, Uv = 0.001)
  )
  #generate labels
  zd_closest <- zd_calc[which.min(abs(zd_calc - label_zd))]
  dl <- df %>%
    dplyr::filter(.data$zd == zd_closest) %>%
    dplyr::mutate(label = .data$Tv)
  #generate plot
  plt <- ggplot2::ggplot() +
    theme_soilmech() +
    ggplot2::geom_path(
      data = df,
      ggplot2::aes(x = .data$Uv, y = .data$zd, group = as.factor(.data$Tv))
    ) +
    ggplot2::geom_label(
      data = dl,
      ggplot2::aes(x = .data$Uv, y = .data$zd, group = as.factor(.data$Tv), label = .data$label),
      hjust = 0.5,
      vjust = 0.5,
      label.padding = ggplot2::unit(0.1, "lines"),
      alpha = label_alpha,
      size = label_size
    ) +
    ggplot2::coord_cartesian(
      xlim = c(0, 1),
      ylim = c(2, 0),
      expand = FALSE
    ) +
    ggplot2::geom_segment(
      ggplot2::aes(
        x = 0.05,
        xend = 0,
        y = 0.2,
        yend = 0.25,
        group = as.factor(0)
      ),
      arrow = ggplot2::arrow(length = ggplot2::unit(0.4, "lines"))
    ) +
    ggplot2::geom_label(
      ggplot2::aes(
        x = 0.05,
        y = 0.2,
        label = "T[v]==0",
        group = as.factor(0)
      ),
      parse = TRUE,
      hjust = 0,
      vjust = 0.5,
      label.padding = ggplot2::unit(0.1, "lines"),
      alpha = label_alpha,
      size = label_size
    ) +
    ggplot2::scale_y_reverse(
      breaks = seq(0, 2, 0.5),
      minor_breaks = seq(0, 2, 0.1)
    ) +
    ggplot2::scale_x_continuous(
      breaks = seq(0, 1, 0.2),
      minor_breaks = seq(0, 1, 0.05)
    ) +
    ggplot2::xlab(expression("Degree of consolidation"~U[v]~"[-]")) +
    ggplot2::ylab(expression(z/d~"[-]")) +
    ggplot2::theme(legend.position = "none")
  #annotate crosshairs if requested
  if (!is.null(zd) & !is.null(Tv)){
    Uv <- calculate_consolidation_openlayer(zd, Tv)
    plt <- ggplot_addcrosshairs(
      plt,
      Uv,
      zd,
      group = paste0("U[v]==", round(Uv, nround)),
      label_parse = TRUE,
      add_colour_scale = TRUE,
      ylim = 2,
      arrow_y = FALSE
    )
  }
  #return
  return(plt)
}


#' Calculate average consolidation in layer as function of time
#'
#' @description
#' calculate average consolidation in a layer as function of time.
#' Calculations approximate the solution using a Taylor series,
#' similar function `calculate_consolidation_openlayer()`.
#'
#' Three cases are distinguished:
#' * `trace = 1` drainage in open layers with thickness 2d (uniform or
#'     triangular), or in half-closed layer with uniform initial pressure
#'     and thickness d
#' * `trace = 2` drainage in half-closed layer with triangular distribution
#'     of initial excess pore water pressure. Larger pressure near closed
#'     end
#' * `trace = 3` drainage in half-closed layer with triangular distribution
#'     of initial excess pore water pressure. Larger pressure near open
#'     end
#'
#' @param Tv normalised time for crosshair annotation
#' @param trace trace selection (1, 2 or 3) for crosshair annotation
#' @param nfourier number of fourier steps to use
#' @importFrom magrittr `%>%`
#' @importFrom rlang .data
#' @return an array with average degree of consolidation Uv
#' @export

calculate_consolidation_averagetime <- function(
  Tv,
  trace = 1,
  nfourier = 50
){
  #generate fourier expansion points
  df <- tidyr::expand_grid(
    Tv = Tv,
    n = seq(1, nfourier)
  )
  #traces
  if (trace == 1) {
    df$Dn <- 2/(pi*df$n)*(1 - cos(pi*df$n))
  } else if (trace == 3) {
    df$Dn <- 2 * ((2/(pi*df$n) - 4/(pi^2*df$n^2)*sin(pi*df$n/2)) - 2/(pi^2*df$n^2)*(2*sin(pi*df$n/2) - 2*sin(pi*df$n) + pi*df$n*cos(pi*df$n)))
  } else if (trace == 2) {
    df$Dn <- 2 * ((4/(pi^2*df$n^2)*sin(pi*df$n/2) - 2/(pi*df$n)*cos(pi*df$n/2)) + 2/(pi^2*df$n^2)*(2*sin(pi*df$n/2) - 2*sin(pi*df$n) + pi*df$n*cos(pi*df$n/2)))
  }
  #calculate degree of consolidation
  df %>%
    dplyr::mutate(u = .data$Dn*exp(-pi^2*.data$n^2*.data$Tv/4)*(1 - cos(pi*.data$n))/(pi*.data$n)) %>%
    dplyr::group_by(.data$Tv) %>%
    dplyr::summarise(Uv = 1 - sum(.data$u)) %>%
    dplyr::pull(.data$Uv)
}


#' ggplot for average consolidation in layer as function of time
#'
#' @description
#' plots a ggplot for the average consolidation in a layer as function of
#' time.
#'
#' Three traces are plotted
#' * `trace = 1` drainage in open layers with thickness 2d (uniform or
#'     triangular), or in half-closed layer with uniform initial pressure
#'     and thickness d
#' * `trace = 2` drainage in half-closed layer with triangular distribution
#'     of initial excess pore water pressure. Larger pressure near closed
#'     end
#' * `trace = 3` drainage in half-closed layer with triangular distribution
#'     of initial excess pore water pressure. Larger pressure near open
#'     end
#'
#' @param Tv normalised time for crosshair annotation
#' @param trace trace selection (1, 2 or 3) for crosshair annotation
#' @param label_Tv time at which to plot time labels on traces
#' @param label_size size of the text in time labels
#' @param label_alpha opacity of the time labels
#' @param nround number of decimals to use in rounding values in crosshair
#'   annotations
#' @importFrom magrittr `%>%`
#' @return a ggplot object
#' #empty chart
#' ggplot_consolidation_averagetime()
#'
#' #annotated chart
#' ggplot_consolidation_averagetime(Tv = 0.5, trace = 1)
#' @export

ggplot_consolidation_averagetime <- function(
  Tv = NULL,
  trace = NULL,
  label_Tv = 0.03,
  label_size = 3,
  label_alpha = 1,
  nround = 2
){
  #data for labels
  dl <- tibble::tibble(
    Tv = label_Tv,
    Uv = c(
      calculate_consolidation_averagetime(label_Tv, trace = 1),
      calculate_consolidation_averagetime(label_Tv, trace = 2),
      calculate_consolidation_averagetime(label_Tv, trace = 3)
    ),
    trace = seq(3)
  )
  #all data for traces
  df <- tibble::tibble(
    Tv = lseq(0.001, 5, l = 101),
    Uv1 = calculate_consolidation_averagetime(Tv, trace = 1),
    Uv2 = calculate_consolidation_averagetime(Tv, trace = 2),
    Uv3 = calculate_consolidation_averagetime(Tv, trace = 3)
  ) %>%
    tidyr::pivot_longer(
      cols = c("Uv1", "Uv2", "Uv3"),
      names_to = "trace",
      values_to = "Uv"
    )
  #ggplot
  plt <- ggplot2::ggplot() +
    theme_soilmech() +
    ggplot2::geom_line(
      data = df,
      ggplot2::aes(x = .data$Tv, y = .data$Uv, group = as.factor(.data$trace))
    ) +
    ggplot2::geom_label(
      data = dl,
      ggplot2::aes(x = .data$Tv, y = .data$Uv, label = .data$trace, group = as.factor(.data$trace)),
      hjust = 0.5,
      vjust = 0.5,
      label.padding = ggplot2::unit(0.1, "lines"),
      alpha = label_alpha,
      size = label_size
    ) +
    ggplot2::scale_y_reverse(
      lim = c(1, 0),
      expand = c(0, 0),
      breaks = seq(0, 1, 0.2),
      minor_breaks = seq(0, 1, 0.1)
    ) +
    ggplot2::scale_x_log10(
      lim = c(0.001, 5),
      expand = c(0, 0),
      breaks = c(0.001, 0.01, 0.1, 1),
      minor_breaks = get_log10_minorbreaks(c(0.001, 5))
    ) +
    ggplot2::annotation_logticks(side = "b") +
    ggplot2::xlab(expression("Dimensionless time factor"~T[v])) +
    ggplot2::ylab(expression("Average degree of consolidation"~bar(U)[v]))
  #add crosshairs
  if (!is.null(Tv) & !is.null(trace)) {
    Uv <- calculate_consolidation_averagetime(Tv, trace)
    plt <- ggplot_addcrosshairs(
      plt,
      Tv,
      Uv,
      group = paste0("bar(U)[v]==", round(Uv, nround)),
      label_parse = TRUE,
      arrow_x = FALSE,
      legend_position = "none",
      add_colour_scale = TRUE,
      xlim = 0.001,
      ylim = 1
    )
  }
  #return
  return(plt)
}


#' ggplot profile types for average layer consolidation times
#'
#' @description
#' Creates a ggplot object with the various profiles of initial excess
#' pore water pressure (uniform, triangular) in open and half-closed
#' soil layers. The average consolidation with time can be plotted with
#' the function `ggplot_consolidation_averagetime()`
#'
#' @param xmarg relative horizontal offset of pressure profiles
#' @param ymarg relative thickness of soil layers on bottom and top
#' @param ymarg2 relative distance between top and pressure arrow
#' @param fill_pressure fill color of pore water pressure profiles
#' @param color_layer text color of soil layer type annotation
#'   (impermeable, permeable)
#' @param fill_layer fill color of soil layer annotations
#'   (impermeable, permeable)
#' @param label_size text font size for all labels
#' @param case_names names for cases: A, B and C
#' @param drain_names names for drainage types: impermeable and permeable
#' @param layer_names layer types: half-closed layer and open layer
#' @param arrow_size length of arrow heads, in unit 'lines'
#' @importFrom magrittr `%>%`
#' @return a ggplot object
#' @examples
#' ggplot_consolidation_averageprofiles()
#' @export

ggplot_consolidation_averageprofiles <- function(
  xmarg = 0.2,
  ymarg = 0.15,
  ymarg2 = 0.30,
  fill_pressure = "#2a7fff",
  color_layer = c("white", "black"),
  fill_layer = c("grey20", "#d3bc5f"),
  label_size = 3.5,
  case_names = c("Case A", "Case B", "Case C"),
  drain_names = c("Permeable", "Impermeable"),
  layer_names = c("Open layer", "Half-closed layer"),
  arrow_size = 0.4
){
  #top: permeable
  dl1 <- tibble::tibble(
    x = c(0, 1, 1, 0),
    y = c(1, 1, 1 - ymarg, 1 - ymarg),
    case1 = layer_names[1],
    type = drain_names[1],
    side = "top"
  )
  dl2 <- tibble::tibble(
    x = c(0, 1, 1, 0),
    y = c(1, 1, 1 - ymarg, 1 - ymarg),
    case1 = layer_names[2],
    type = drain_names[1],
    side = "top"
  )
  #bottom: permeable
  dl3 <- tibble::tibble(
    x = c(0, 1, 1, 0),
    y = c(0, 0, ymarg, ymarg),
    case1 = layer_names[1],
    type = drain_names[1],
    side = "bottom"
  )
  #bottom: closed
  dl4 <- tibble::tibble(
    x = c(0, 1, 1, 0),
    y = c(0, 0, ymarg, ymarg),
    case1 = layer_names[2],
    type = drain_names[2],
    side = "bottom"
  )
  #layers - merge
  dl <- dplyr::bind_rows(dl1, dl2, dl3, dl4)
  #profiles
  dp <- tibble::tibble(
    x = rep(c(
        c(xmarg, 1 - xmarg, 1 - xmarg, xmarg),
        c(xmarg, 1 - xmarg, xmarg, xmarg),
        c(xmarg, xmarg, 1 - xmarg, xmarg)
      ), 2),
    y = rep(c(1 - ymarg, 1 - ymarg, ymarg, ymarg), 6),
    case1 = c(
      rep(layer_names[1], 4*3),
      rep(layer_names[2], 4*3)
    ),
    case2 = rep(rep(case_names, each = 4), 2)
  )
  #labels
  db <- tibble::tibble(
    x = rep(c(0.5, xmarg + (1 - 2*xmarg)/3, xmarg + (1 - 2*xmarg)/3), 2),
    y = rep(c(0.5, ymarg + 2*(1 - 2*ymarg)/3, ymarg + (1 - 2*ymarg)/3), 2),
    label = c("1", "1", "1", "1", "3", "2"),
    case1 = rep(layer_names, each = 3),
    case2 = rep(case_names, 2)
  )
  #layer names
  dl2 <- dl %>%
    dplyr::group_by(.data$case1, .data$type, .data$side) %>%
    dplyr::summarize(
      x = mean(.data$x),
      y = mean(.data$y),
      label = .data$type[1]
    )
  #annotation: pressure arrows
  darr <- tibble::tibble(
    x = xmarg,
    xend = 1 - xmarg,
    y = 1 - ymarg2,
    yend = 1 - ymarg2,
    case1 = layer_names[2],
    case2 = case_names[1]
  )
  #thickness arrows
  darr2 <- tibble::tibble(
    x = rep(0.5*xmarg, 3),
    xend = rep(0.5*xmarg, 3),
    y = c(ymarg, ymarg, 0.5),
    yend = c(1 - ymarg, 0.5, 1 - ymarg),
    case1 = c(layer_names[2], rep(layer_names[1], 2)),
    case2 = case_names[1]
  )
  #text labels
  dtext <- tibble::tibble(
    x = c(rep(0.75*xmarg, 3), 0.5),
    y = c(0.5, 0.5*(ymarg + 0.5), 1 - 0.5*(ymarg + 0.5), 1 - 0.5*(ymarg + ymarg2)),
    label = c("d", "d", "d", "u[e]"),
    case1 = c(layer_names[2], rep(layer_names[1], 2), layer_names[2]),
    case2 = case_names[1]
  )
  #plot
  ggplot2::ggplot() +
    ggplot2::theme_void() +
    #pressure
    ggplot2::geom_polygon(
      data = dp,
      ggplot2::aes(x = .data$x, y = .data$y),
      fill = fill_pressure
    ) +
    #layers top and bottom
    ggplot2::geom_polygon(
      data = dl,
      ggplot2::aes(x = .data$x, y = .data$y, fill = as.factor(.data$type))
    ) +
    ggplot2::geom_text(
      data = dl2,
      ggplot2::aes(x = .data$x, y = .data$y, label = .data$label, color = as.factor(.data$type)),
      hjust = 0.5,
      vjust = 0.5,
      size = label_size
    ) +
    #type - in middle of pressure
    ggplot2::geom_label(
      data = db,
      ggplot2::aes(x = .data$x, y = .data$y, label = .data$label),
      hjust = 0.5,
      vjust = 0.5,
      size = label_size,
      label.padding = ggplot2::unit(0.15, "lines")
    ) +
    #annotations - pressure arrow
    ggplot2::geom_segment(
      data = darr,
      ggplot2::aes(x = .data$x, xend = .data$xend, y = .data$y, yend = .data$yend),
      arrow = ggplot2::arrow(ends = "both", length = ggplot2::unit(arrow_size, "lines"))
    ) +
    #annotation - thickness arrows
    ggplot2::geom_segment(
      data = darr2,
      ggplot2::aes(x = .data$x, xend = .data$xend, y = .data$y, yend = .data$yend),
      arrow = ggplot2::arrow(ends = "both", length = ggplot2::unit(arrow_size, "lines"))
    ) +
    #annotation - text
    ggplot2::geom_text(
      data = dtext,
      ggplot2::aes(x = .data$x, y = .data$y, label = .data$label),
      hjust = 0.5,
      vjust = 0.5,
      size = label_size,
      parse = TRUE
    ) +
    #facets and axes
    ggplot2::facet_grid(.data$case1 ~ .data$case2) +
    ggplot2::theme(legend.position = "none") +
    ggplot2::scale_fill_manual(values = fill_layer) +
    ggplot2::scale_color_manual(values = color_layer)
}
