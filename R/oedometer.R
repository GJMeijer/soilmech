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
#'   - string for bisection line between lines (1) and (2)
#'   - string for tangent line virgin behaviour (3)
#'   - string for vertical line down from intersection lines (3) and (4)
#'   - string for preconsolidation pressure label
#'
#' @param stages plot fitting stages (array with numbers)
#'
#'   - `1`: adds smooth trace through measurement points
#'   - `2`: adds trangent at max curvature
#'   - `3`: adds horizontal line at max curvature
#'   - `4`: adds bisection line between 2 and 3
#'   - `5`: adds tangent at virgin compression
#'   - `6`: adds vertical at intersection of 4 and 5
#'
#' @param include_value if `TRUE`, the value of the preconsolidation pressure
#'   is added to the plotted label
#' @param nsignif the number of significant digits to use in the (rounded)
#'   value of the show preconsolidation pressure (if `include_value = TRUE`)
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
  xlim = c(4, NA),
  ylim = c(NA, NA),
  palette = "Set1",
  label_line = c("1", "2", "3", "4", "5", 'sigma*minute[v*","*c]'),
  stages = seq(6),
  include_value = TRUE,
  nsignif = 2
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
  e0 <- as.double(stats::coef(ft)[1])
  Cs <- as.double(stats::coef(ft)[2])
  Cc <- as.double(stats::coef(ft)[3])
  sigma_vc <- as.double(stats::coef(ft)[4])
  b <- as.double(stats::coef(ft)[5])
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
  #plot bisection line
  if (4 %in% stages) {
    plt <- plt +
      ggplot2::annotate(
        "segment",
        x = sigma_vc,
        xend = xlim[2],
        y = e_c,
        yend = e_c - tan(0.5*atan2(Cp, 1))*log10(xlim[2]/sigma_vc),
        color = colo[3]
      ) +
      ggplot2::annotate(
        "label",
        x = sigma_vl,
        y = e_c - Ch*log10(sigma_vl/sigma_vc),
        label = label_line[3],
        hjust = 0.5,
        vjust = 0.5,
        color = colo[3]
      )
  }
  #plot virgin line
  if (5 %in% stages) {
    plt <- plt +
      ggplot2::annotate(
        "segment",
        x = xlim[1],
        xend = xlim[2],
        y = e_m + Cm*log10(max(sigma_v)/xlim[1]),
        yend = e_m + Cm*log10(max(sigma_v)/xlim[2]),
        color = colo[4]
      ) +
      ggplot2::annotate(
        "label",
        x = sigma_vl,
        y = e_m - Cm*log10(sigma_vl/max(sigma_v)),
        label = label_line[4],
        hjust = 0.5,
        vjust = 0.5,
        color = colo[4]
      )
  }
  #plot vertical and preconsolidation pressure
  if (6 %in% stages) {
    #label for pressure
    if (include_value == TRUE) {
      label_line[6] <- paste0(
        label_line[6], "%~~%", signif(sigma_v0, nsignif), "~kPa"
      )
    }
    #plot
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


#' Create oedometer time-sample height data for single load step
#'
#' @description
#' Create some time - sample height oedometer data, given known soil
#' properties, for a single load step
#'
#' @param t array with time steps. If not defined, a wide range is assumed
#' @param h0 initial sample height [m]
#' @param dsigmav increase in vertical stress [Pa]
#' @param mv coefficient of compressibility [m^2/N]
#' @param cv coefficient of consolidation [m^2/s]
#' @param Calpha creep strain per log10-cycle of time. Creep is only added
#'   after primary consolidation is (almost) finished. A small smoothing
#'   is applied to smoothen the sudden jump in gradient of settlement
#' @param Tvalpha value of normalised strain at which creep stats to be
#'   included
#' @return tibble with field for time (`t`, in seconds) and current sample
#'   height (`h`, in meters)
#' @examples
#' #get data
#' df <- oedometer_create_data()
#'
#' #plot data - log scale for time
#' ggplot2::ggplot() +
#'   ggplot2::geom_line(
#'     data = df,
#'     ggplot2::aes(x = t, y = h)
#'   ) +
#'   ggplot2::scale_x_log10()
#' @export

oedometer_create_data <- function(
  t = NULL,
  h0 = 20 * 1e-3,
  dsigmav = 20e3 * 1e6,
  mv = 5e-6,
  cv = 1e-8,
  Calpha = 0.01,
  Tvalpha = 10^0.2
){
  #if time not specified, generate generic trace based on standard open layer
  #consolidation
  if (is.null(t)) {
    Tv <- lseq(0.005, 10, l = 101)
    t <- Tv*(0.5*h0)^2/cv
  } else {
    Tv <- cv*t/((0.5*h0)^2)
  }
  #degree of consolidation - open layer - simplified curve
  # Tv = pi/4*Uv^2               when Uv < 0.6
  # Tv = -0.933*log10(1-Uv)-0.085  when Uv > 0.6
  Uv <- sqrt(4*Tv/pi)
  i <- (Uv > 0.6)
  Uv[i] <- 1 - 10^(-(Tv[i] + 0.085)/0.933)
  #sample height
  h <- h0 - h0*mv*dsigmav*Uv
  #add some creep
  h <- h - Calpha*h0*0.5*(log10(Tv) + sqrt(log10(Tv)^2 + log10(Tvalpha)^2))
  #return
  return(tibble::tibble(t = t, h = h))
}


#' ggplot cv determination log method
#'
#' @description
#' ggplot showing how to obtain t50 from oedometer data using the log-time
#' method
#'
#' @inheritParams oedometer_create_data
#' @param t time points for individual measurement points
#' @param t1 time point for getting initial height
#' @param xlim,ylim manual x and y-axis limits
#' @param palette RColorBrewer palette to use for colours of annotations
#' @param line linestyles for lines in the three steps
#' @param shape marker shape for the three steps
#' @return ggplot object
#' @examples
#' #default example
#' ggplot_cv_log()
#' @export

ggplot_cv_log <- function(
  t = c(seq(10, 60, 10), c(1, 2, 4, 8, 15, 30)*60, c(1, 2, 4, 8, 24)*60*60),
  t1 = 20,
  h0 = 0.0195 * 1e3,
  dsigmav = 20e3 * 1e-6,
  mv = 5e-6 * 1e6,
  cv = 1e-8 * 1e6,
  Calpha = 0.005,
  Tvalpha = 10^0.1,
  xlim = c(5, NA),
  ylim = c(NA, NA),
  palette = "Set1",
  line = c(2, 4, 5),
  shape = c(15, 15, 15)
){
  ## Get data
  #get trace data
  df <- oedometer_create_data(
    t = lseq(min(t), max(t), l = 101),
    h0 = h0, dsigmav = dsigmav, mv = mv, cv = cv, Calpha = Calpha, Tvalpha = Tvalpha
  )
  #get point data
  dp <- oedometer_create_data(
    t = t,
    h0 = h0, dsigmav = dsigmav, mv = mv, cv = cv, Calpha = Calpha, Tvalpha = Tvalpha
  )
  #get data at t1 and t2
  dt <- oedometer_create_data(
    t = c(t1, 4*t1),
    h0 = h0, dsigmav = dsigmav, mv = mv, cv = cv, Calpha = Calpha, Tvalpha = Tvalpha
  )
  #plot limits
  xlim <- round_limits(df$t, lower = xlim[1], upper = xlim[2], log = TRUE)
  ylim <- round_limits(df$h, lower = ylim[1], upper = ylim[2])

  ## Get intersection points etc.
  #get initial height
  h0 <- 2*dt$h[1] - dt$h[2]
  #get max gradient line
  tgrad <- (pi/4*0.6^2)*(0.5*h0)^2/cv
  hgrad <- h0 - 0.6*h0*dsigmav*mv
  htgrad <- -h0*mv*dsigmav*0.5*log(10)*sqrt(16*cv*tgrad/(pi*h0^2)) #dh/dlog10(t)
  hgradline <- ylim[1] + c(0, 0.9*diff(ylim))
  tgradline <- 10^(log10(tgrad) + (hgradline - hgrad)/htgrad)
  #get creep gradient line
  tcreep <- Tvalpha*(0.5*h0)^2/cv
  hcreep <- h0 - h0*dsigmav*mv
  htcreep <- -h0*Calpha
  tcreepline <- 10^(log10(xlim[1]) + diff(log10(xlim))*c(0.5, 1))
  hcreepline <- hcreep + htcreep*log10(tcreepline/tcreep)
  #intersection
  t100 <- 10^((hcreep - hgrad + htgrad*log10(tgrad) - htcreep*log10(tcreep))/(htgrad - htcreep))
  h100 <- hcreep + htcreep*(log10(t100/tcreep))
  #half deformation - and time
  h50 <- (h0 + h100)/2
  t50 <- pi/4*0.5^2*(0.5*h0)^2/cv

  ## Plot
  #colors
  colo <- RColorBrewer::brewer.pal(3, palette)
  #initiate plot
  plt <- ggplot2::ggplot() +
    theme_soilmech() +
    #plot trace and points
    ggplot2::geom_path(
      data = df,
      ggplot2::aes(x = .data$t, y = .data$h)
    ) +
    ggplot2::geom_point(
      data = dp,
      ggplot2::aes(x = .data$t, y = .data$h)
    ) +
    #plot axes
    ggplot2::scale_y_continuous(
      limits = ylim
    ) +
    ggplot2::scale_x_log10(
      limits = xlim,
      breaks = scales::trans_breaks("log10", function(x) 10^x),
      labels = scales::trans_format("log10", scales::math_format(expr = ~10^.x)),
      minor_breaks = get_log10_minorbreaks(xlim),
      expand = c(0, 0)
    ) +
    ggplot2::xlab(expression("Time"~t~"[s]")) +
    ggplot2::ylab(expression("Sample height"~h~"[mm]"))
  #add lines
  plt <- plt +
    #first step - find h0
    ggplot2::geom_vline(
      xintercept = dt$t,
      color = colo[1], linetype = line[1]
    ) +
    ggplot2::geom_hline(
      yintercept = dt$h,
      color = colo[1], linetype = line[1]
    ) +
    ggplot2::geom_hline(
      yintercept = 2*dt$h[1] - dt$h[2],
      color = colo[1], linetype = line[1]
    ) +
    #second step - find h100
    ggplot2::annotate(
      "segment",
      x = tgradline[1], xend = tgradline[2],
      y = hgradline[1], yend = hgradline[2],
      color = colo[2], linetype = line[2]
    ) +
    ggplot2::annotate(
      "segment",
      x = tcreepline[1], xend = tcreepline[2],
      y = hcreepline[1], yend = hcreepline[2],
      color = colo[2], linetype = line[2]
    ) +
    ggplot2::geom_hline(
      yintercept = h100,
      color = colo[2], linetype = line[2]
    ) +
    #third step - find h50 and t50
    ggplot2::geom_hline(
      yintercept = h50,
      color = colo[3], linetype = line[3]
    ) +
    ggplot2::geom_vline(
      xintercept = t50,
      color = colo[3], linetype = line[3]
    )
  #add points
  plt <- plt +
    #first step - find h0
    ggplot2::annotate(
      "point",
      x = dt$t, y = dt$h,
      color = colo[1], shape = shape[1]
    ) +
    #second step - find h100
    ggplot2::annotate(
      "point",
      x = t100, y = h100,
      color = colo[2], shape = shape[2]
    ) +
    #third step - find h50 and t50
    ggplot2::annotate(
      "point",
      x = t50, y = h50,
      color = colo[3], shape = shape[3]
    )
  #add labels
  plt <- plt +
    #first step - find h0
    ggplot2::annotate(
      "label",
      x = dt$t, y = ylim[1], label = c("t[1]", "4*t[1]"),
      parse = TRUE, color = colo[1], hjust = 0.5, vjust = 0
    ) +
    ggplot2::annotate(
      "label",
      x = xlim[2], y = h0, label = "h[0]",
      parse = TRUE, color = colo[1], hjust = 1, vjust = 0.5) +
    #second step - find h100
    ggplot2::annotate(
      "label",
      x = xlim[2], y = h100, label = "h[100]",
      parse = TRUE, color = colo[2], hjust = 1, vjust = 0.5
    ) +
    #third step - find h50 and t50
    ggplot2::annotate(
      "label",
      x = xlim[2], y = h50, label = "h[50]",
      parse = TRUE, color = colo[3], hjust = 1, vjust = 0.5
    ) +
    ggplot2::annotate(
      "label",
      x = t50, y = ylim[1], label = "t[50]",
      parse = TRUE, color = colo[3], hjust = 0.5, vjust = 0
    )
  #return
  return(plt)
}


#' ggplot cv determination sqrt method
#'
#' @description
#' ggplot showing how to obtain t50 from oedometer data using the sqrt
#' method
#'
#' @inheritParams ggplot_cv_log
#' @return ggplot object
#' @examples
#' #default example
#' ggplot_cv_sqrt()
#' @export

ggplot_cv_sqrt <- function(
  t = c(seq(10, 60, 10), c(1, 2, 4, 8, 15, 30)*60, c(1, 2, 4, 8, 24)*60*60),
  h0 = 0.0195 * 1e3,
  dsigmav = 20e3 * 1e-6,
  mv = 5e-6 * 1e6,
  cv = 1e-8 * 1e6,
  Calpha = 0.005,
  Tvalpha = 10^0.1,
  xlim = c(NA, 150),
  ylim = c(NA, NA),
  palette = "Set1",
  line = c(2, 4, 5),
  shape = c(15, 15, 15)
){
  ## Get data
  #get trace data
  df <- oedometer_create_data(
    t = lseq(min(t), max(t), l = 101),
    h0 = h0, dsigmav = dsigmav, mv = mv, cv = cv, Calpha = Calpha, Tvalpha = Tvalpha
  )
  #get point data
  dp <- oedometer_create_data(
    t = t,
    h0 = h0, dsigmav = dsigmav, mv = mv, cv = cv, Calpha = Calpha, Tvalpha = Tvalpha
  )
  #plot limits
  xlim <- round_limits(sqrt(df$t), lower = xlim[1], upper = xlim[2])
  ylim <- round_limits(df$h, lower = ylim[1], upper = ylim[2])

  ## Plot elements
  #gradient line
  hgrad <- c(h0, ylim[1])
  htgrad <- -4*mv*dsigmav*sqrt(cv/pi)
  tgrad <- ((hgrad - h0)/htgrad)^2
  #offsetted gradient line (15%)
  hgrad2 <- hgrad
  tgrad2 <- 1.15^2*tgrad
  #intersection
  sol <- stats::uniroot(
    function(h) {4*cv/h0^2*((h0 - h)*1.15/htgrad)^2 - (-0.933*log10(1 - (h0 - h)/(mv*dsigmav*h0)) - 0.085)},
    h0 - c(1, 0)*mv*dsigmav*h0
  )
  hinters <- sol$root
  tinters <- ((h0 - hinters)*1.15/htgrad)^2

  ## Plot
  #colors
  colo <- RColorBrewer::brewer.pal(3, palette)
  #initiate plot
  plt <- ggplot2::ggplot() +
    theme_soilmech() +
    #plot trace and points
    ggplot2::geom_path(
      data = df,
      ggplot2::aes(x = sqrt(.data$t), y = .data$h)
    ) +
    ggplot2::geom_point(
      data = dp,
      ggplot2::aes(x = sqrt(.data$t), y = .data$h)
    ) +
    #plot axes
    ggplot2::scale_x_continuous(lim = xlim, expand = c(0, 0)) +
    ggplot2::scale_y_continuous(lim = ylim) +
    ggplot2::xlab(expression(sqrt(Time)~sqrt(t)~"["*s^{0.5}*"]")) +
    ggplot2::ylab(expression("Sample height "*h*" [mm]"))
  #add gradient lines
  plt <- plt +
    ggplot2::annotate(
      "segment",
      x = sqrt(tgrad[1]), xend = sqrt(tgrad[2]),
      y = hgrad[1], yend = hgrad[2],
      color = colo[1],
      linetype = line[1]
    ) +
    ggplot2::annotate(
      "segment",
      x = sqrt(tgrad2[1]), xend = sqrt(tgrad2[2]),
      y = hgrad2[1], yend = hgrad2[2],
      color = colo[2],
      linetype = line[2]
    ) +
    ggplot2::annotate(
      "point",
      x = tgrad[1], y = hgrad[1],
      color = colo[1],
      shape = shape[1]
    ) +
    ggplot2::geom_vline(
      xintercept = sqrt(tinters),
      color = colo[3],
      linetype = line[3]
    ) +
    ggplot2::annotate(
      "point",
      x = sqrt(tinters), y = hinters,
      color = colo[3],
      shape = shape[3]
    ) +
    ggplot2::annotate(
      "label",
      x = sqrt(tinters), y = ylim[2],
      label = "t[90]",
      parse = TRUE,
      hjust = 0.5, vjust = 1,
      color = colo[3]
    )
  #return
  return(plt)
}
