#' Calculate Ic from known Fr and Qtn, in zone 2 - 7
#'
#' @description
#' Calculate Ic from known Fr and Qtn, in zone 2 - 7 in Robertson's
#' chart
#'
#' @param Fr normalised friction ratio, in \%
#' @param Qtn normalised tip resistance
#' @return Robertson's Ic index
#' @export

FrQtn2Ic <- function(Fr, Qtn) {
  sqrt((3.47 - log10(Qtn))^2 + (1.22 + log10(Fr))^2)
}


#' Calculate Fr from known Qtn and Ic, in zone 2 - 7
#'
#' @description
#' Calculate Fr from known Qtn and Ic, in zone 2 - 7 in Robertson's
#' chart
#'
#' @param Ic Robertson's Ic index
#' @param Qtn normalised tip resistance
#' @return normalised friction ratio, in \%
#' @export

IcQtn2Fr <- function(Ic, Qtn) {
  10^(sqrt(Ic^2 - (3.47 - log10(Qtn))^2) - 1.22)
}


#' Calculate Qtn from known Ic and Fr, in zone 2 - 7
#'
#' @description
#' Calculate Ic from known Ic and Fr, in zone 2 - 7 in Robertson's
#' chart
#'
#' @param Ic Robertson's Ic index
#' @param Fr normalised friction ratio, in \%
#' @return normalised tip resistance Qtn
#' @export

IcFr2Qtn <- function(Ic, Fr) {
  10^(3.47 - sqrt(Ic^2 - (1.22 + log10(Fr))^2))
}


#' Fr-Qtn relationship for boundary zone 1
#'
#' @description
#' Calculate Qtn based on known Fr for line seperating zone 1 from the
#' other zones
#'
#' @param Fr Friction ratio, in percentage
#' @return Qtn values
#' @export

QtnA <- function(Fr) {
  12*exp(-1.4*Fr)
}


#' Fr-Qtn relationship for boundary zone 8 and 9
#'
#' @description
#' Calculate Qtn based on known Fr for line seperating zone 8 and 9 from the
#' other zones
#'
#' @param Fr Friction ratio, in percentage
#' @return Qtn values
#' @export

QtnG <- function(Fr) {
  1/(0.005*(Fr - 1) - 0.0003*(Fr - 1)^2 - 0.002)
}


#' Qtn-Fr relationship for boundary between zone 8 and 9
#'
#' @description
#' Calculate Fr based on known Qtn for line seperating zones 8 and 9
#'
#' @return Fr values
#' @export

FrH <- function() {
  4.5
}


#' Create polygons for each soil zone in Robertson's chart
#'
#' @description
#' Generate polygons defined by (Frn, Qtn) points for each soil zone in
#' Robertson's chart
#'
#' @param nFr number of points on curved polygon edges. Fr interval is split
#'   linearly
#' @param Fr_lim Fr lower and upper limit of plot
#' @param Qtn_lim Qtn lower and upper limit of plot
#' @return tibble with Fr (field `Fr`), Qtn (field `Qtn`), zone numbers
#'   (field `zone`) and soil descriptions (field `soil`)
#' @export

robertson_polygons <- function(
  nFr = 21,
  Fr_lim = c(0.1, 10),
  Qtn_lim = c(1, 1000)
){
  #Fr coordinates of intersections between zones - solve functions
  Fr012 <- stats::uniroot(function(Fr){QtnA(Fr) - Qtn_lim[1]}, c(1, 2))$root
  Fr123 <- stats::uniroot(function(Fr){QtnA(Fr) - IcFr2Qtn(3.60, Fr)}, c(1.4, 1.6))$root
  Fr134 <- stats::uniroot(function(Fr){QtnA(Fr) - IcFr2Qtn(2.95, Fr)}, c(0.4, 0.8))$root
  Fr145 <- stats::uniroot(function(Fr){QtnA(Fr) - IcFr2Qtn(2.60, Fr)}, c(0.2, 0.3))$root
  Fr459 <- stats::uniroot(function(Fr){IcFr2Qtn(2.60, Fr) - QtnG(Fr)}, c(5.0, 7.0))$root
  Fr568 <- stats::uniroot(function(Fr){IcFr2Qtn(2.05, Fr) - QtnG(Fr)}, c(2.0, 3.0))$root
  Fr067 <- stats::uniroot(function(Fr){IcFr2Qtn(1.31, Fr) - Qtn_lim[2]}, c(0.9, 1.1))$root
  Fr068 <- stats::uniroot(function(Fr){QtnG(Fr) - Qtn_lim[2]}, c(1.5, 2.0))$root
  Fr089 <- FrH()
  #zones - always start on left, than lower points, clockwise
  pol1 <- tibble::tibble(
    Fr = c(
      Fr_lim[1],
      lseq(Fr_lim[1], Fr012, l = nFr)
    ),
    Qtn = c(
      Qtn_lim[1],
      QtnA(lseq(Fr_lim[1], Fr012, l = nFr))
    ),
    zone = 1,
    soil = "Sensitive clays and silts"
  )
  pol2 <- tibble::tibble(
    Fr = c(
      lseq(Fr123, Fr_lim[2], l = nFr),
      Fr_lim[2],
      lseq(Fr012, Fr123, l = nFr)
    ),
    Qtn = c(
      IcFr2Qtn(3.60, lseq(Fr123, Fr_lim[2], l = nFr)),
      Qtn_lim[1],
      QtnA(lseq(Fr012, Fr123, l = nFr))
    ),
    zone = 2,
    soil = "Organic soils"
  )
  pol3 <- data.frame(
    Fr = c(
      lseq(Fr134, Fr_lim[2], l = nFr),
      lseq(Fr_lim[2], Fr123, l = nFr),
      lseq(Fr123, Fr134, l = nFr)
    ),
    Qtn = c(
      IcFr2Qtn(2.95, lseq(Fr134, Fr_lim[2], l = nFr)),
      IcFr2Qtn(3.60, lseq(Fr_lim[2], Fr123, l = nFr)),
      QtnA(lseq(Fr123, Fr134, l = nFr))
    ),
    zone = 3,
    soil = "Clays"
  )
  pol4 <- tibble::tibble(
    Fr = c(
      lseq(Fr145,Fr459,l=nFr),
      lseq(Fr459,Fr_lim[2],l=nFr),
      lseq(Fr_lim[2],Fr134,l=nFr),
      lseq(Fr134,Fr145,l=nFr)
    ),
    Qtn = c(
      IcFr2Qtn(2.60, lseq(Fr145, Fr459, l = nFr)),
      QtnG(lseq(Fr459, Fr_lim[2], l = nFr)),
      IcFr2Qtn(2.95, lseq(Fr_lim[2], Fr134, l = nFr)),
      QtnA(lseq(Fr134, Fr145, l = nFr))
    ),
    zone = 4,
    soil = "Silty mixtures"
  )
  pol5 <- tibble::tibble(
    Fr = c(
      lseq(Fr_lim[1], Fr568, l = nFr),
      lseq(Fr568, Fr459, l = nFr),
      lseq(Fr459, Fr145, l = nFr),
      lseq(Fr145, Fr_lim[1], l = nFr)
    ),
    Qtn = c(
      IcFr2Qtn(2.05, lseq(Fr_lim[1], Fr568, l = nFr)),
      QtnG(lseq(Fr568, Fr459, l = nFr)),
      IcFr2Qtn(2.60, lseq(Fr459, Fr145, l = nFr)),
      QtnA(lseq(Fr145, Fr_lim[1], l = nFr))
    ),
    zone = 5,
    soil = "Sandy mixtures"
  )
  pol6 <- tibble::tibble(
    Fr = c(
      lseq(Fr_lim[1], Fr067, l = nFr),
      lseq(Fr068, Fr568, l = nFr),
      lseq(Fr568, Fr_lim[1], l = nFr)
    ),
    Qtn = c(
      IcFr2Qtn(1.31, lseq(Fr_lim[1], Fr067, l = nFr)),
      QtnG(lseq(Fr068, Fr568, l = nFr)),
      IcFr2Qtn(2.05, lseq(Fr568, Fr_lim[1], l = nFr))
    ),
    zone = 6,
    soil = "Sands"
  )
  pol7 <- tibble::tibble(
    Fr = c(
      Fr_lim[1],
      lseq(Fr067, Fr_lim[1], l = nFr)
    ),
    Qtn = c(
      Qtn_lim[2],
      IcFr2Qtn(1.31, lseq(Fr067, Fr_lim[1], l = nFr))
    ),
    zone = 7,
    soil = "Gravelly sands"
  )
  pol8 <- tibble::tibble(
    Fr = c(
      Fr089,
      lseq(Fr089, Fr068, l = nFr)
    ),
    Qtn = c(
      Qtn_lim[2],
      QtnG(lseq(Fr089, Fr068, l = nFr))
    ),
    zone = 8,
    soil = "Very stiff OC sands to clayey sands"
  )
  pol9 <- tibble::tibble(
    Fr = c(
      Fr089,
      Fr_lim[2],
      lseq(Fr_lim[2], Fr089, l = nFr)
    ),
    Qtn = c(
      Qtn_lim[2],
      Qtn_lim[2],
      QtnG(lseq(Fr_lim[2], Fr089, l = nFr))
    ),
    zone = 9,
    soil = "Very stiff OC clays to silts"
  )
  #Bind together, and return
  dplyr::bind_rows(pol1, pol2, pol3, pol4, pol5, pol6, pol7, pol8, pol9)
}


#' Assign Robertson's soil type based in Fr and Qtn values
#'
#' @description
#' Assign Robertson's soil type based in Fr and Qtn values. Function returns
#' the number of the soil type, i.e. 1 - 9
#'
#' @param Fr Friction ratio, in \%
#' @param Qtn Normalised tip resistance
#' @return array with soil type numbers
#' @export

robertson_assign_zone <- function(Fr, Qtn){
  #initiate output
  zone <- rep(NA, length(Fr))
  #check if below lines (or to left, in case of H)
  bA <- (Qtn <= QtnA(Fr))
  bB <- (Qtn <= IcFr2Qtn(3.60, Fr))
  bC <- (Qtn <= IcFr2Qtn(2.95, Fr))
  bD <- (Qtn <= IcFr2Qtn(2.60, Fr))
  bE <- ifelse((Fr >= 10^(2.05 - 1.22)), TRUE, (Qtn <= IcFr2Qtn(2.05, Fr)))
  bF <- ifelse((Fr >= 10^(1.31 - 1.22)), TRUE, (Qtn <= IcFr2Qtn(1.31, Fr)))
  bG <- ifelse(Fr <= (28/3 - sqrt(2260/36)), TRUE, (Qtn <= QtnG(Fr)))
  lH <- (Fr <= FrH())
  #assign zones
  zone[(bA == TRUE)] <- 1
  zone[(bA == FALSE) & (bG == TRUE) & (bB == TRUE)] <- 2
  zone[(bA == FALSE) & (bG == TRUE) & (bB == FALSE) & (bC == TRUE)] <- 3
  zone[(bA == FALSE) & (bG == TRUE) & (bC == FALSE) & (bD == TRUE)] <- 4
  zone[(bA == FALSE) & (bG == TRUE) & (bD == FALSE) & (bE == TRUE)] <- 5
  zone[(bA == FALSE) & (bG == TRUE) & (bE == FALSE) & (bF == TRUE)] <- 6
  zone[(bA == FALSE) & (bG == TRUE) & (bF == FALSE)] <- 7
  zone[(bG == FALSE) & (lH == TRUE)] <- 8
  zone[(bG == FALSE) & (lH == FALSE)] <- 9
  #return
  return(zone)
}


#' ggplot Robertson's chart
#'
#' @description
#' Plot an empty Robertson's chart
#'
#' @inheritParams robertson_polygons
#' @param Fr normalised friction values for crosshairs, in \%
#' @param Qtn normalised tip resitance values for crosshairs
#' @param group labels for crosshairs
#' @param palette colormap to use for coloring different zones
#' @param alpha transparancy of soil polygon fills
#' @param legend_position position of legend
#' @importFrom magrittr `%>%`
#' @importFrom rlang .data
#' @return ggplot object
#' @examples
#' #empty chart
#' ggplot_robertsonschart()
#'
#' #annotated chart
#' ggplot_robertsonschart(Fr = c(0.5, 2), Qtn = c(200, 50), group = c("A", "B"))
#' @export

ggplot_robertsonschart <- function(
  Fr = NULL,
  Qtn = NULL,
  group = NULL,
  nFr = 25,
  Fr_lim = c(0.1, 10),
  Qtn_lim = c(1, 1000),
  palette = "Set1",
  alpha = 0.25,
  legend_position = "right"
){
  #create polygons
  rob_polygon <- robertson_polygons(nFr = nFr, Fr_lim = Fr_lim, Qtn_lim = Qtn_lim)
  #create labels to plot in zones
  rob_labels <- rob_polygon %>%
    dplyr::group_by(.data$zone, .data$soil) %>%
    dplyr::summarize(10^polygon_centroid(log10(.data$Fr), log10(.data$Qtn))) %>%
    dplyr::rename(Fr = .data$x, Qtn = .data$y)
  #ggplot
  plt <- ggplot2::ggplot() +
    theme_soilmech() +
    ggplot2::theme(
      legend.position = legend_position
    ) +
    ggplot2::geom_polygon(
      data = rob_polygon,
      ggplot2::aes(x = .data$Fr, y = .data$Qtn, fill = paste0(.data$zone, ": ", .data$soil)),
      color = "black",
      alpha = alpha
    ) +
    ggplot2::geom_text(
      data = rob_labels,
      ggplot2::aes(x = .data$Fr, y = .data$Qtn, label = .data$zone),
      color = "black"
    ) +
    ggplot2::coord_cartesian(xlim = Fr_lim, ylim = Qtn_lim, expand = FALSE) +
    ggplot2::scale_x_log10(
      minor_breaks = c(seq(0.1, 1, 0.1), seq(2, 10, 1)),
      breaks = c(0.1, 1, 10),
      labels = scales::trans_format("log10", scales::math_format(expr = ~10^.x))
    ) +
    ggplot2::scale_y_log10(
      minor_breaks = c(seq(1, 10, 1), seq(20, 100, 10), seq(200, 1000, 100)),
      breaks = c(1, 10, 100, 1000),
      labels = scales::trans_format("log10", scales::math_format(expr = ~10^.x))
    ) +
    ggplot2::annotation_logticks(sides = "trbl") +
    ggplot2::xlab(expression("Normalised friction"~F[r]~"[%]")) +
    ggplot2::ylab(expression("Normalised tip resistance"~Q[tn]~"[-]")) +
    ggplot2::scale_fill_brewer(name = "Soil type:", palette = palette)
  #add annotations
  if (!is.null(Fr) & !is.null(Qtn)) {
    zone <- robertson_assign_zone(Fr, Qtn)
    if (is.null(group)) {
      label <- as.character(zone)
    } else {
      label <- paste0(group, ": ", zone)
    }
    plt <- ggplot_addcrosshairs(
      plt,
      x = Fr,
      y = Qtn,
      group = label,
      add_colour_scale = TRUE,
      xlim = Fr_lim[1],
      ylim = Qtn_lim[1],
      label_parse = FALSE
    )
  }
  #return
  plt
}
