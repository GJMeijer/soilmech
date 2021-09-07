#' ggplot Mohr-Coulomb failure criterion based on Mohr circles
#'
#' @description
#' Plot an drained or undrained mohr-coulomb fit based on measured
#' principal stresses.
#'
#' Mohr-Coulomb fits are obtained by minimising the sum of square distances
#' between failure surface and Mohr circles
#'
#' @param sig1 major principal stresses (array)
#' @param sig3 minor principal stresses (array)
#' @param drained if `TRUE`, stresses are assumed to be effective, and an
#'   cohesional/frictional fit is returned. If `FALSE`, stresses are assumed
#'   to be total stresses, and a undrained shear strength fit is returned
#' @param cohesionless if `TRUE`, assume no cohesion when analysing a drained
#'   problem
#' @param palette RColorBrewer palette for plotting circles
#' @param nc number of points to use to plot a circle
#' @param sig_angle_line relative x-length of horizontal line for angle plot
#' @param sig_angle_arrow relative x-position of phi arrow
#' @param sig_angle_label relative x-position of phi label
#' @param sig_cohesion_arrow relative x-position of cohesion arrow
#' @param sig_cohesion_label relative x-position of cohesion label
#' @importFrom magrittr `%>%`
#' @importFrom rlang .data
#' @return ggplot object
#' @examples
#' #drained example - cohesional/frictional
#' sig1 = c(100, 130, 160)
#' sig3 = c(50, 70, 85)
#' ggplot_triaxial_mohrcoulomb(sig1, sig3, drained = TRUE)
#'
#' #drained example - cohesionless
#' sig1 = c(100, 130, 160)
#' sig3 = c(50, 70, 85)
#' ggplot_triaxial_mohrcoulomb(sig1, sig3, drained = TRUE, cohesionless = TRUE)
#'
#' #undrained example
#' sig1 = c(100, 130, 160)
#' sig3 = c(50, 75, 102)
#' ggplot_triaxial_mohrcoulomb(sig1, sig3, drained = FALSE)
#' @export

ggplot_triaxial_mohrcoulomb <- function(
  sig1 = c(75, 120, 160),
  sig3 = c(20, 40, 55),
  drained = TRUE,
  cohesionless = FALSE,
  palette = "Set1",
  nc = 181,
  sig_angle_line = 0.20,
  sig_angle_arrow = 0.18,
  sig_angle_label = 0.19,
  sig_cohesion_arrow = 0.02,
  sig_cohesion_label = 0.03
){
  #generate (half) Mohr-circles
  dc <- tidyr::expand_grid(
    tibble::tibble(sig1 = sig1, sig3 = sig3, id = seq(length(sig1))),
    theta = seq(0, pi, l = nc)
  ) %>%
    dplyr::mutate(
      sig = 0.5*(.data$sig1 + .data$sig3) + 0.5*(.data$sig1 - .data$sig3)*cos(.data$theta),
      q = 0.5*(.data$sig1 - .data$sig3)*sin(.data$theta)
    )
  #fit a trend line
  if (drained == TRUE) {
    if (cohesionless == TRUE) {
      #drained behaviour - cohesionless
      sol_fit <- stats::optimise(
        f = function(par, p, q){sum((q - p*sin(par))^2)},
        interval = c(0, 50/180*pi),
        p = (sig1 + sig3)/2,
        q = (sig1 - sig3)/2
      )
      coh <- 0
      phi <- sol_fit$minimum
    } else {
      #drained behaviour - cohesional/frictional
      sol_fit <- stats::optim(
        par = c(0.5*mean(sig1 - sig3), 0.4),
        fn = function(par, p, q){sum((q - ((p*par[2] + par[1])/(sqrt(1 + par[2]^2))))^2)},
        p = (sig1 + sig3)/2,
        q = (sig1 - sig3)/2
      )
      coh <- sol_fit$par[1]
      phi <- atan(sol_fit$par[2])
    }
    coh_label <- paste0("c*minute==", signif(coh, 3), "~kPa")
    phi_label <- paste0("phi*minute==", signif(phi*180/pi, 3), "*degree")
  } else {
    #undrained behaviour
    coh <- 0.5*mean(sig1 - sig3)
    phi <- 0
    coh_label <- paste0("c[u]==", signif(coh, 3), "~kPa")
    phi_label <- paste0("phi==0*degree")
  }
  #generate a fit line
  dfit <- tibble::tibble(
    sig = c(0, max(sig1)),
    q = c(coh, coh + tan(phi)*max(sig1))
  )
  #plot
  plt <- ggplot2::ggplot() +
    theme_soilmech() +
    ggplot2::geom_path(
      data = dc,
      ggplot2::aes(x = .data$sig, y = .data$q, color = as.factor(.data$id))
    ) +
    ggplot2::geom_line(
      data = dfit,
      ggplot2::aes(x = .data$sig, y = .data$q),
      color = "black",
      linetype = 2
    ) +
    ggplot2::coord_fixed(
      ratio = 1,
      xlim = round_limits(1.1*max(sig1), lower = 0),
      ylim = round_limits(1.3 * 0.5*max(sig1 - sig3), lower = 0),
      expand = FALSE
    ) +
    ggplot2::scale_color_brewer(palette = palette) +
    ggplot2::theme(legend.position = "none") +
    ggplot2::ylab("q [kPa]")
  #annotations for cohesion
  if (!dplyr::near(coh, 0)) {
    plt <- plt +
      ggplot2::annotate(
        "segment",
        x = sig_cohesion_arrow*max(sig1),
        xend = sig_cohesion_arrow*max(sig1),
        y = 0,
        yend = coh,
        arrow = ggplot2::arrow(length = ggplot2::unit(0.25, "lines"), ends = "both")
      ) +
      ggplot2::annotate(
        "text",
        x = (sig_cohesion_label)*max(sig1),
        y = 0.5*coh,
        label = coh_label,
        parse = TRUE,
        hjust = 0,
        vjust = 0.5
      )
  }
  #annotations for friction angle
  if (phi > 0) {
    plt <- plt +
      ggplot2::annotate(
        "segment",
        x = 0,
        xend = sig_angle_line*max(sig1),
        y = coh,
        yend = coh
      ) +
      ggplot2::geom_path(
        ggplot2::aes(
          x = sig_angle_arrow*max(sig1)*cos(seq(0, phi, l = 25)),
          y = coh + sig_angle_arrow*max(sig1)*sin(seq(0, phi, l = 25))
        ),
        arrow = ggplot2::arrow(length = ggplot2::unit(0.25, "lines"))
      ) +
    ggplot2::annotate(
      "text",
      x = sig_angle_label*max(sig1)*cos(0.5*phi),
      y = coh + sig_angle_label*max(sig1)*sin(0.5*phi),
      label = phi_label,
      hjust = 0,
      vjust = 0.5,
      parse = TRUE
    )
  } else {
    plt <- plt +
      ggplot2::annotate(
        "text",
        x = sig_angle_label*max(sig1),
        y = 1.05*coh,
        label = phi_label,
        hjust = 0,
        vjust = 0,
        parse = TRUE
      )
  }
  #add x-axis
  if (drained == TRUE) {
    plt <- plt + ggplot2::xlab(expression(sigma*"' [kPa]"))
  } else {
    plt <- plt + ggplot2::xlab(expression(sigma~"[kPa]"))
  }
  #return
  return(plt)
}
