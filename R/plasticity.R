#' Classify soil based on plasticity chart
#'
#' @description
#' Get shorthand for soil major type and plasticity level. Output is given
#' in shorthand, i.e. MH for a high plasticity silt etc. Classification is
#' based on the plasticity chart
#'
#' @param wL liquid limits, in \%
#' @param Ip plasticity index, in \%
#' @return array with shorthand classifications
#' @examples
#' classify_plasticity(wL = 60, Ip = 20)
#' @export

classify_plasticity <- function(wL, Ip) {
  #get soil type, based on A-line
  soil_type <- rep("M", length(wL))
  soil_type[Ip >= pmax(6.3, 0.73*(wL - 20))] <- "C"
  #get plasticity level
  plasticity_level <- rep("L", length(wL))
  plasticity_level[(wL >= 35) & (wL < 50)] <- "I"
  plasticity_level[(wL >= 50) & (wL < 70)] <- "H"
  plasticity_level[(wL >= 70) & (wL < 90)] <- "V"
  plasticity_level[(wL >= 90)] <- "E"
  #return
  paste0(soil_type, plasticity_level)
}


#' Plot plasticity chart
#'
#' @description
#' Function to plot a plasticity chart containing both the A and B lines
#'
#' @param wL liquid limits for crosshairs (in \%)
#' @param Ip plasticity index for crosshairs (in \%)
#' @param add_legend if `TRUE`, return a plot with a legend
#' @param textsize_class font size of text in behaviour types
#' @param textsize_legend font size of text in legend
#' @return a ggplot object
#' @examples
#' #plot empty chart
#' ggplot_plasticitychart()
#'
#' #plot chart with annotation
#' ggplot_plasticitychart(wL = 60, Ip = 20)
#' @export

ggplot_plasticitychart <- function(
  wL = NULL,
  Ip = NULL,
  add_legend = TRUE,
  textsize_class = 3,
  textsize_legend = 2.5
) {
  #Plasticity chart limits [%]
  wLmin <- 0
  wLmax <- 120
  Ipmin <- 0
  Ipmax <- 80
  #Plasticity classification
  dlim <- data.frame(
    Symbol = c("L", "I", "H", "V", "E"),
    Text = c("Low plasticity", "Medium plasticity", "High plasticity", "Very high plasticity", "Extremely high plasticity")
  )
  #limits for regions
  dlim$wLmin <- c(wLmin, 35, 50, 70, 90)
  dlim$wLmax <- c(35, 50, 70, 90, wLmax)
  dlim$wLavg <- 0.5*(dlim$wLmin + dlim$wLmax)
  #A-B-line
  dAB <- data.frame(wL = c(20, 20 + 6.3/0.73, wLmax, 15, wLmax))
  dAB$Type <- c("A", "A", "A", "B", "B")
  dAB$Ip <- pmax(6.3, 0.73*(dAB$wL - 20))
  dAB$Ip[dAB$Type == "B"] <- pmax(6.3, 0.9*(dAB$wL[dAB$Type == "B"] - 8))
  #Symbols
  dsym <- data.frame(PSymbol = rep(dlim$Symbol, 2))
  dsym$SSymbol <- c(rep("M", 5), rep("C", 5))
  dsym$Symbol <- paste0(dsym$SSymbol, dsym$PSymbol)
  dsym$wL <- pmax(30, rep(dlim$wLavg, 2))
  dsym$IpA <- pmax(6.3, 0.73*(dsym$wL - 20))
  dsym$IpB <-  pmax(6.3, 0.9*(dsym$wL - 8))
  dsym$Ip <- 0.5*dsym$IpA
  dsym$Ip[dsym$SSymbol == "C"] <- 0.5*(dsym$IpA[dsym$SSymbol == "C"] + pmin(Ipmax, dsym$IpB[dsym$SSymbol == "C"]))
  #legend label
  lab <- "Predominant behaviour\nmaterial <0.425 mm\n- C = Clay\n- M = Silt\n\nPlasticity:\n- L = Low\n- I = Intermediate\n- H = High\n- V = Very High\n- E = Extremely high"
  #Create plot
  plt <- ggplot2::ggplot() +
    theme_soilmech() +
    ggplot2::theme(legend.position = "none") +
    ggplot2::geom_line(
      data = dAB,
      ggplot2::aes(x = .data$wL, y = .data$Ip, linetype = .data$Type),
      show.legend = FALSE
    ) +
    ggplot2::geom_vline(
      data = dlim,
      ggplot2::aes(xintercept = wLmin),
      size = 0.5,
      linetype = 3,
      show.legend = FALSE
    ) +
    ggplot2::scale_x_continuous(breaks = seq(wLmin, wLmax, 20)) +
    ggplot2::scale_y_continuous(breaks = seq(Ipmin, Ipmax, 20)) +
    ggplot2::coord_cartesian(
      xlim = c(wLmin, wLmax),
      ylim = c(Ipmin, Ipmax),
      expand = FALSE
    ) +
    ggplot2::xlab(expression(Liquid ~limit~w[L]~"[%]")) +
    ggplot2::ylab(expression(Plasticity~index~I[p]~"[%]")) +
    ggplot2::geom_text(
      data = dsym,
      ggplot2::aes(x = .data$wL,y = .data$Ip, label = .data$Symbol),
      show.legend = FALSE
    ) +
    ggplot2::annotate(
      "label",
      x = 105,
      y = 0.73*(105 - 20),
      label = "'A'-line",
      size = textsize_class
    ) +
    ggplot2::annotate(
      "label",
      x = 15,
      y = 0.90*(15 - 8),
      label = "'B'-line",
      hjust = 1,
      size = textsize_class
    )
  if (add_legend == TRUE){
    plt <- plt + ggplot2::annotate(
      "label",
      x = wLmin + 0.01*(wLmax - wLmin),
      y = 0.99*Ipmax,
      label = lab,
      size = textsize_legend,
      hjust = 0,
      vjust = 1
    )
  }
  #add crosshairs is requested
  if (!is.null(wL) & !is.null(Ip)) {
    plt <- ggplot_addcrosshairs(
      plt,
      wL,
      Ip,
      group = classify_plasticity(wL, Ip),
      add_colour_scale = TRUE
    )
  }
  #return
  plt
}


#' ggplotly plasticity fall cone chart: water content - penetration depth
#'
#' @description
#' Plots a chart with measured water content versus fall cone penetration
#' depth in a fall cone test. The best linear fit is plotted and the
#' liquid limit (rounded) is indicated
#'
#' @param w measured water contents, in \%
#' @param u measured penetration depths, in mm
#' @param nround number of decimals in rounded liquid limit
#' @param palette RColorBrewer color palette to use for colors
#' @return ggplot object
#' @examples
#' w = c(9, 13, 16, 20)
#' u = c(12, 17, 21, 25)
#' ggplot_fallcone(w, u, nround = 0)
#' @export

ggplot_fallcone <- function(
  w = c(9, 13, 16, 20),
  u = c(12, 17, 21, 25),
  nround = 0,
  palette = "Set1"
){
  #generate colors
  colo <- RColorBrewer::brewer.pal(3, palette)
  #get liquid limit through linear fit
  fit <- stats::lm(u ~ w)
  wL <- (20 - stats::coef(fit)[1])/stats::coef(fit)[2]
  #fit line
  dft <- tibble::tibble(
    w = c(min(w), max(w)),
    u = stats::coef(fit)[1] + stats::coef(fit)[2]*c(min(w), max(w))
  )
  #plot limits
  xlim <- round_limits(c(0.95*min(c(w, wL)), 1.05*max(c(w, wL))))
  ylim <- round_limits(c(0.95*min(c(u, 20)), 1.05*max(c(u, 20))))
  #plot
  plt <- ggplot2::ggplot() +
    theme_soilmech() +
    ggplot2::geom_line(
      data = dft,
      ggplot2::aes(x = .data$w, y = .data$u),
      color = colo[3],
      linetype = 2
    ) +
    ggplot2::geom_point(
      ggplot2::aes(x = w, y = u),
      color = colo[2]
    ) +
    ggplot2::coord_cartesian(
      xlim = xlim,
      ylim = ylim,
      expand = FALSE
    ) +
    ggplot2::xlab("Water content w [%]") +
    ggplot2::ylab("Penetration depth u [mm]")
  #add crosshairs and liquid limit label
  ggplot_addcrosshairs(
    plt,
    wL,
    20,
    add_colour_scale = TRUE,
    group = paste0("W[L]==", round(wL, nround), "*'%'"),
    arrow_y = FALSE,
    ylim = ylim[1],
    legend_position = "none",
    label_parse = TRUE,
    label_nudge_x = 0.2,
    label_nudge_y = 0.5*(ylim[1] - 20),
    label_vjust = 0.5
  )
}
