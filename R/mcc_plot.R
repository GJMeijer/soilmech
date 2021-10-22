#' MCC yield surface coordinates (positive side only)
#'
#' @description
#' Create p and q positions for a MCC yield surface
#'
#' @param pc preconsolidation pressure
#' @param M MCC M-parameters
#' @param t array with angles at which to plot point
#' @return tibble with `p` and `q`
#' @export

mcc_yieldsurface_pq <- function(pc, M, t = seq(0, pi, l = 101)) {
  return(tibble::tibble(
    p = 0.5*pc*(1 + cos(t)),
    q = 0.5*M*sin(t)*pc
  ))
}


#' ggplot MCC axial strain versus deviatoric stress
#'
#' @description
#' Plot axial strain versus deviatoric stress
#'
#' @param df tibble with all MCC predictions
#' @param palette RColorBrewer palette for annotation colours
#' @param e1_max maximum strain to plot
#' @return a ggplot object
#' @examples
#' #drained
#' df <- mcc_solve_triax_drained()
#' ggplot_mcc_eps1_q(df)
#' ggplot_mcc_eps1_q(df, e1_max = 0.04)
#'
#' #undrained
#' df <- mcc_solve_triax_undrained()
#' ggplot_mcc_eps1_q(df)
#' @export

ggplot_mcc_eps1_q <- function(df, palette = "Set1", e1_max = NULL) {
  #generate some colors
  colo <- RColorBrewer::brewer.pal(3, palette)
  #subset data - only below certain strain
  if (is.null(e1_max)) {
    dff <- df
  } else {
    dff <- dplyr::filter(df, .data$e1 <= e1_max)
  }
  #plot
  ggplot2::ggplot() +
    theme_soilmech() +
    #add trace
    ggplot2::geom_path(
      data = dff,
      ggplot2::aes(x = .data$e1, y = .data$q)
    ) +
    #add starting point
    ggplot2::annotate(
      "point",
      x = dff$e1[1],
      y = dff$q[1],
      shape = 15
    ) +
    #add yield point
    ggplot2::geom_point(
      data = dff %>%
        dplyr::filter(.data$plastic == TRUE) %>%
        dplyr::slice_min(.data$e1),
      ggplot2::aes(x = .data$e1, y = .data$q),
      color = colo[2],
      shape = 17
    ) +
    #add final point
    ggplot2::geom_point(
      data = dff %>%
        dplyr::slice_max(.data$e1),
      ggplot2::aes(x = .data$e1, y = .data$q),
      color = colo[1]
    ) +
    #axes
    ggplot2::xlab(expression("Axial strain"~epsilon[1]~"[-]")) +
    ggplot2::ylab(expression("Deviatoric stress"~q~"[kPa]")) +
    ggplot2::coord_cartesian(
      xlim = round_limits(df$e1, lower = 0),
      ylim = round_limits(df$q, lower = 0),
      expand = FALSE
    )
}


#' ggplot MCC axial strain versus volumetric strain
#'
#' @description
#' Plot axial strain versus volumetric strain
#'
#' @param df tibble with all MCC predictions
#' @param palette RColorBrewer palette for annotation colours
#' @param e1_max maximum strain to plot
#' @return a ggplot object
#' @examples
#' #drained
#' df <- mcc_solve_triax_drained()
#' ggplot_mcc_eps1_epsv(df)
#' ggplot_mcc_eps1_epsv(df, e1_max = 0.04)
#'
#' #undrained
#' df <- mcc_solve_triax_undrained()
#' ggplot_mcc_eps1_epsv(df)
#' @export

ggplot_mcc_eps1_epsv <- function(df, palette = "Set1", e1_max = NULL) {
  #generate some colors
  colo <- RColorBrewer::brewer.pal(3, palette)
  #subset data - only below certain strain
  if (is.null(e1_max)) {
    dff <- df
  } else {
    dff <- dplyr::filter(df, .data$e1 <= e1_max)
  }
  #plot
  ggplot2::ggplot() +
    theme_soilmech() +
    #add trace
    ggplot2::geom_path(
      data = dff,
      ggplot2::aes(x = .data$e1, y = -.data$ev)
    ) +
    #add starting point
    ggplot2::annotate(
      "point",
      x = dff$e1[1],
      y = -dff$ev[1],
      shape = 15
    ) +
    #add yield point
    ggplot2::geom_point(
      data = dff %>%
        dplyr::filter(.data$plastic == TRUE) %>%
        dplyr::slice_min(.data$e1),
      ggplot2::aes(x = .data$e1, y = -.data$ev),
      color = colo[2],
      shape = 17
    ) +
    #add final point
    ggplot2::geom_point(
      data = dff %>%
        dplyr::slice_max(.data$e1),
      ggplot2::aes(x = .data$e1, y = -.data$ev),
      color = colo[1]
    ) +
    #axes
    ggplot2::xlab(expression("Axial strain"~epsilon[1]~"[-]")) +
    ggplot2::ylab(expression("Volumetric strain"~epsilon[v]~"[-]")) +
    ggplot2::coord_cartesian(
      xlim = round_limits(df$e1, lower = 0),
      ylim = round_limits(-df$ev),
      expand = FALSE
    )
}


#' ggplot MCC p-q stress paths
#'
#' @description
#' Plot (effective) p-q stress paths
#'
#' @param df tibble with all MCC predictions
#' @param M MCC M-parameter
#' @param pc initial preconsolidation pressure
#' @param palette RColorBrewer palette for annotation colours
#' @param e1_max maximum strain to plot
#' @return a ggplot object
#' @examples
#' M <- 1.35
#' pc <- 100
#'
#' #drained
#' df <- mcc_solve_triax_drained(M = M, pc = pc)
#' ggplot_mcc_p_q(df, M = M, pc = pc)
#' ggplot_mcc_p_q(df, e1_max = 0.01)
#'
#' #undrained
#' df <- mcc_solve_triax_undrained(M = M, pc = pc, p0 = 80, nu = 0.3)
#' ggplot_mcc_p_q(df, M = M, pc = pc)
#' @export

ggplot_mcc_p_q <- function(
  df,
  M = 1.35,
  pc = 100,
  palette = "Set1",
  e1_max = NULL
){
  #generate some colors
  colo <- RColorBrewer::brewer.pal(4, palette)
  #subset data - only below certain strain
  if (is.null(e1_max)) {
    dff <- df
  } else {
    dff <- dplyr::filter(df, .data$e1 <= e1_max)
  }
  #current point
  dffc <- dplyr::slice_max(dff, .data$e1)
  #yield surfaces
  dy1 <- mcc_yieldsurface_pq(pc, M)
  dy2 <- mcc_yieldsurface_pq(dffc$pc, M)
  #plot
  ggplot2::ggplot() +
    theme_soilmech() +
    #plot current yield surface
    ggplot2::geom_polygon(
      data = dy2,
      ggplot2::aes(x = .data$p, y = .data$q),
      color = colo[3],
      fill = colo[3],
      alpha = 0.10,
      linetype = 1
    ) +
    #plot initial yield surface
    ggplot2::geom_path(
      data = dy1,
      ggplot2::aes(x = .data$p, y = .data$q),
      color = colo[3],
      linetype = 2
    ) +
    #plot M-line
    ggplot2::annotate(
      "segment",
      x = 0, xend = c(max(dy1$p, dy2$p)),
      y = 0, yend = M*c(max(dy1$p, dy2$p)),
      color = colo[4],
      linetype = 1
    ) +
    #add trace
    ggplot2::geom_path(
      data = dff,
      ggplot2::aes(x = .data$p, y = .data$q)
    ) +
    #add starting point
    ggplot2::annotate(
      "point",
      x = dff$p[1],
      y = dff$q[1],
      shape = 15
    ) +
    #add yield point
    ggplot2::geom_point(
      data = dff %>%
        dplyr::filter(.data$plastic == TRUE) %>%
        dplyr::slice_min(.data$e1),
      ggplot2::aes(x = .data$p, y = .data$q),
      color = colo[2],
      shape = 17
    ) +
    #add final point
    ggplot2::geom_point(
      data = dffc,
      ggplot2::aes(x = .data$p, y = .data$q),
      color = colo[1]
    ) +
    #axes
    ggplot2::xlab(expression("Isotropic effective stress"~p*minute~"[kPa]")) +
    ggplot2::ylab(expression("Deviatoric stress"~q~"[kPa]")) +
    ggplot2::coord_cartesian(
      xlim = round_limits(c(dy1$p, dy2$p), lower = 0),
      ylim = round_limits(c(dy1$q, dy2$q), lower = 0),
      expand = FALSE
    )
}


#' ggplot MCC ln(p)-v stress paths
#'
#' @description
#' Plot (effective) ln(p)-v stress paths
#'
#' @param df tibble with all MCC predictions
#' @param Gamma MCC CSL intercept specific volume
#' @param lambda,kappa MCC compression parameters
#' @param palette RColorBrewer palette for annotation colours
#' @param e1_max maximum strain to plot
#' @return a ggplot object
#' @examples
#' Gamma <- 2
#' lambda <- 0.115
#' kappa <- 0.015
#' df <- mcc_solve_triax_drained(Gamma = Gamma, lambda = lambda, kappa = kappa)
#' ggplot_mcc_lnp_v(df, Gamma = Gamma, lambda = lambda, kappa = kappa)
#' ggplot_mcc_lnp_v(df, Gamma = Gamma, lambda = lambda, kappa = kappa, e1_max = 0.01)
#' @export

ggplot_mcc_lnp_v <- function(
  df,
  Gamma = 2,
  lambda = 0.115,
  kappa = 0.015,
  palette = "Set1",
  e1_max = NULL
){
  #generate some colors
  colo <- RColorBrewer::brewer.pal(4, palette)
  #subset data - only below certain strain
  if (is.null(e1_max)) {
    dff <- df
  } else {
    dff <- dplyr::filter(df, .data$e1 <= e1_max)
  }
  #current point
  dffc <- dplyr::slice_max(dff, .data$e1)
  #limits
  xlim <- round_limits(c(df$p, df$p), log = TRUE)
  ylim <- round_limits(c(df$v, df$v))
  #CSL and ICL
  dcsl <- tibble::tibble(
    p = xlim,
    v = Gamma - lambda*log(xlim)
  )
  dicl <- tibble::tibble(
    p = xlim,
    v = Gamma - lambda*log(xlim) + log(2)*(lambda - kappa)
  )
  #plot
  ggplot2::ggplot() +
    theme_soilmech() +
    #add critical state line
    ggplot2::geom_path(
      data = dcsl,
      ggplot2::aes(x = .data$p, y = .data$v),
      color = colo[3],
      linetype = 1
    ) +
    #add isotropic consolidation line
    ggplot2::geom_path(
      data = dicl,
      ggplot2::aes(x = .data$p, y = .data$v),
      color = colo[3],
      linetype = 2
    ) +
    #add trace
    ggplot2::geom_path(
      data = dff,
      ggplot2::aes(x = .data$p, y = .data$v)
    ) +
    #add starting point
    ggplot2::annotate(
      "point",
      x = dff$p[1],
      y = -dff$v[1],
      shape = 15
    ) +
    #add yield point
    ggplot2::geom_point(
      data = dff %>%
        dplyr::filter(.data$plastic == TRUE) %>%
        dplyr::slice_min(.data$e1),
      ggplot2::aes(x = .data$p, y = .data$v),
      color = colo[2],
      shape = 17
    ) +
    #add final point
    ggplot2::geom_point(
      data = dffc,
      ggplot2::aes(x = .data$p, y = .data$v),
      color = colo[1]
    ) +
    #axes
    ggplot2::xlab(expression("Isotropic effective stress"~p*minute~"[kPa]")) +
    ggplot2::ylab(expression("Specific volume"~v~"[-]")) +
    ggplot2::coord_cartesian(
      xlim = xlim,
      ylim = ylim,
      expand = FALSE
    ) +
    ggplot2::scale_x_log10()
}
