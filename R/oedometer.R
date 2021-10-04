#' ggplot for Casagrande oedometer preconsolidation pressure
#'
#' @description
#' Function to create a ggplot showing how Casagrande's method for determining
#' the preconsolidation pressure from oedometer data works.
#'
#' The point of maximum curvature is determined by fitting a bilinear curve
#' with a hyperbolic transition area.
#'
#' @md
#' @param sigma_v array with oedometer stresses
#' @param e array with void ratios
#' @param xlim,ylim 2-parameters arrays for x and y axes limits. If not
#'   defined, these limits are determined automatically
#' @param palette RColorBrewer pallette for plot colors
#' @param label_line array with six character strings, which will be parsed
#'
#'   - string for tangent line at maximum curvature (1)
#'   - string for horizontal line (2)
#'   - string for tangent line virgin behaviour (3)
#'   - string for bisection line between lines (2) and (3)
#'   - string for vertical line down from intersection lines (3) and (4)
#'   - string for preconsolidation pressure label
#'
#' @param stages plot fitting stages (array with numbers)
#'
#'   - `1`: adds smooth trace through measurement points
#'   - `2`: adds trangent at max curvature
#'   - `3`: adds horizontal line at max curvature
#'   - `4`: adds tangent at virgin compression
#'   - `5`: adds bisection line between 3 and 4
#'   - `6`: adds vertical at intersection of 4 and 5
#'
#' @examples
#' #all stages
#' ggplot_casagrande_preconsolidation()
#'
#' #stage by stage
#' ggplot_casagrande_preconsolidation(stages = c())
#' ggplot_casagrande_preconsolidation(stages = seq(1))
#' ggplot_casagrande_preconsolidation(stages = seq(2))
#' ggplot_casagrande_preconsolidation(stages = seq(3))
#' ggplot_casagrande_preconsolidation(stages = seq(4))
#' ggplot_casagrande_preconsolidation(stages = seq(5))
#' ggplot_casagrande_preconsolidation(stages = seq(6))
#' @export

ggplot_casagrande_preconsolidation <- function(
  sigma_v = c(5, 10, 20, 50, 100, 200, 500, 1000, 2000),
  e = c(0.94, 0.93, 0.92, 0.90, 0.88, 0.85, 0.80, 0.75, 0.71),
  xlim = c(1, NA),
  ylim = c(NA, NA),
  palette = "Set1",
  label_line = c("A", "B", "C", "D", "E", 'sigma*minute[v*","*0]'),
  stages = seq(6)
) {
  #colors
  colo <- RColorBrewer::brewer.pal(5, palette)
  #fit to get point of maximum curvature
  ft <- stats::nls(
    e ~  e0 - Cs*log10(sigma_v) - 0.5*(b*log(cosh(log10(sigma_v/sigma_vc)/b)) + log10(sigma_v))*(Cc - Cs),
    sigma_v = sigma_v,
    algorithm = "port",
    start = list(e0 = min(e), Cs = 0.015, Cc = 0.15, sigma_vc = 10^mean(log10(sigma_v)), b = 0.2),
    lower = list(e0 = -Inf, Cs = 0, Cc = 0, sigma_vc = min(sigma_v), b = -Inf),
    upper = list(e0 = Inf, Cs = Inf, Cc = Inf, sigma_vc = max(sigma_v), b = Inf)
  )
  #assign coefficients
  e0 <- as.double(coef(ft)[1])
  Cs <- as.double(coef(ft)[2])
  Cc <- as.double(coef(ft)[3])
  sigma_vc <- as.double(coef(ft)[4])
  b <- as.double(coef(ft)[5])
  #generate fitting line
  dfit <- tibble::tibble(sigma_v = 10^seq(log10(min(sigma_v)), log10(max(sigma_v)), l = 101))
  dfit$e <- e0 - Cs*log10(dfit$sigma_v) - 0.5*(b*log(cosh(log10(dfit$sigma_v/sigma_vc)/b)) + log10(dfit$sigma_v))*(Cc - Cs)
  #void ratio and gradient at max curvature
  e_c <- e0 - Cs*log10(sigma_vc) - 0.5*(log10(sigma_vc))*(Cc - Cs)
  Cp <- Cs + 0.5*(Cc - Cs)
  #void ratio and gradient at maximum stress
  e_m <- utils::tail(dfit$e, 1)
  Cm <- Cs + 0.5*(Cc - Cs)*(1 + tanh(log10(max(sigma_v)/sigma_vc)/b))
  #half gradient
  Ch <- tan(0.5*atan2(Cp, 1))
  #intersection of half-angle line with Cc line
  sigma_v0 <- 10^((e_m - e_c + Cm*log10(max(sigma_v)) - Ch*log10(sigma_vc))/(Cm - Ch))
  e_v <- e_c - Ch*log10(sigma_v0/sigma_vc)
  #plot limits
  xlim <- round_limits(1.1*sigma_v, lower = xlim[1], upper = xlim[2])
  ylim <- round_limits(c(0.95*e, 1.05*e), lower = ylim[1], upper = ylim[2])
  #position of labels
  sigma_vl <- 10^(log10(min(sigma_v)) + (0.95*log10(max(sigma_v)/min(sigma_v))))
  #plot
  plt <- ggplot2::ggplot() +
    #add individual points - measured
    ggplot2::geom_point(
      ggplot2::aes(sigma_v, y = e)
    ) +
    #axis etc
    theme_soilmech() +
    ggplot2::coord_cartesian(xlim = xlim, ylim = ylim, expand = FALSE) +
    ggplot2::scale_x_log10(
      breaks = 10^seq(0, 6),
      labels = scales::trans_format("log10", scales::math_format(expr = ~10^.x)),
      minor_breaks = get_log10_minorbreaks(xlim),
    ) +
    ggplot2::annotation_logticks(side = "b") +
    ggplot2::xlab(expression("Vertical effective stress"~sigma*minute[v]~"[kPa]")) +
    ggplot2::ylab(expression("Void ratio"~e~"[-]"))
  #plot fitting line
  if (1 %in% stages) {
    plt <- plt +
      ggplot2::geom_line(
        data = dfit,
        ggplot2::aes(x = .data$sigma_v, y = .data$e),
        color = "black",
        linetype = 2
      )
  }
  #plot max gradient point and line
  if (2 %in% stages) {
    plt <- plt +
      ggplot2::annotate(
        "point",
        x = sigma_vc,
        y = e_c,
        color = colo[1]
      ) +
      ggplot2::annotate(
        "segment",
        x = xlim[1],
        xend = xlim[2],
        y = e_c - Cp*log10(xlim[1]/sigma_vc),
        yend = e_c - Cp*log10(xlim[2]/sigma_vc),
        color = colo[1]
      ) +
      ggplot2::annotate(
        "label",
        x = sigma_vl,
        y = e_c - Cp*log10(sigma_vl/sigma_vc),
        label = label_line[1],
        hjust = 0.5,
        vjust = 0.5,
        color = colo[1]
      )
  }
  #plot horizontal line at gradient
  if (3 %in% stages) {
    plt <- plt +
      ggplot2::annotate(
        "segment",
        x = sigma_vc,
        xend = xlim[2],
        y = e_c,
        yend = e_c,
        color = colo[2]
      ) +
      ggplot2::annotate(
        "label",
        x = sigma_vl,
        y = e_c,
        label = label_line[2],
        hjust = 0.5,
        vjust = 0.5,
        color = colo[2]
      )
  }
  #plot virgin line
  if (4 %in% stages) {
    plt <- plt +
      ggplot2::annotate(
        "segment",
        x = xlim[1],
        xend = xlim[2],
        y = e_m + Cm*log10(max(sigma_v)/xlim[1]),
        yend = e_m + Cm*log10(max(sigma_v)/xlim[2]),
        color = colo[3]
      ) +
      ggplot2::annotate(
        "label",
        x = sigma_vl,
        y = e_m - Cm*log10(sigma_vl/max(sigma_v)),
        label = label_line[3],
        hjust = 0.5,
        vjust = 0.5,
        color = colo[3]
      )
  }
  #plot half-line
  if (5 %in% stages) {
    plt <- plt +
      ggplot2::annotate(
        "segment",
        x = sigma_vc,
        xend = xlim[2],
        y = e_c,
        yend = e_c - tan(0.5*atan2(Cp, 1))*log10(xlim[2]/sigma_vc),
        color = colo[4]
      ) +
      ggplot2::annotate(
        "label",
        x = sigma_vl,
        y = e_c - Ch*log10(sigma_vl/sigma_vc),
        label = label_line[4],
        hjust = 0.5,
        vjust = 0.5,
        color = colo[4]
      )
  }
  #plot vertical and preconsolidation pressure
  if (6 %in% stages) {
    plt <- plt +
      ggplot2::annotate(
        "point",
        x = sigma_v0,
        y = e_v,
        color = colo[5]
      ) +
      ggplot2::annotate(
        "segment",
        x = sigma_v0,
        xend = sigma_v0,
        y = e_v,
        yend = min(ylim),
        arrow = ggplot2::arrow(length = ggplot2::unit(0.1, "inches")),
        color = colo[5]
      ) +
      ggplot2::annotate(
        "text",
        x = 10^(1.01*log10(sigma_v0)),
        y = min(ylim) + 0.25*(e_v - min(ylim)),
        label = label_line[6],
        hjust = 0,
        vjust = 0.5,
        parse = TRUE,
        color = colo[5]
      ) +
      ggplot2::annotate(
        "label",
        x = sigma_v0,
        y = 0.5*(e_v + ylim[1]),
        label = label_line[5],
        hjust = 0.5,
        vjust = 0.5,
        color = colo[5]
      )
  }
  #return
  return(plt)
}
