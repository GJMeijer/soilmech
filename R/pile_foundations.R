#' plotly Berezantsev's chart for pile tip bearing capacity factor
#'
#' @description
#' Create a plotly chart for pile tip bearing capacity factor Nq as function
#' of angle of internal friction for a sandy soil, according to Berezantsev.
#'
#' This function uses data from a digitised plot, so may not be 100% accurate
#'
#' @param palette RColorBrewer pallete for line colors
#' @return plotly object
#' @examples
#' plotly_piletip_berezantsev()
#' @export

plotly_piletip_berezantsev <- function(
  palette = "Set1"
){
  #file path
  file_path <- system.file(
    "extdata",
    "Nq_berezantsev_datapoints_digitised.csv",
    package = "soilmech"
  )
  #load data from file
  df_raw <- readr::read_csv(
    file_path,
    col_types = readr::cols()
  )
  #interpolate at every degree
  df <- tibble::tibble(
    phi = seq(25, 45),
    Nq = stats::approx(df_raw$phi, df_raw$Nq, xout = .data$phi)$y,
    hover_label = paste0("\u03c6' = ", .data$phi, "<br>N<sub>q</sub> = ", round(.data$Nq, 1))
  )
  #colors
  colo <- RColorBrewer::brewer.pal(3, palette)
  #create plotly object
  plt <- plotly::plot_ly(
    data = df,
    x = ~phi,
    y = ~Nq,
    type = "scatter",
    mode = "lines",
    line = list(color = colo[1]),
    text = ~hover_label,
    hoverinfo = "text"
  )
  #add layout
  plotly::layout(
    plt,
    xaxis = list(
      range = c(min(df$phi), max(df$phi)),
      title = "\u03c6' [\u00b0]"
    ),
    yaxis = list(
      type = "log",
      range = log10(c(10, 1000)),
      title = "N<sub>q</sub> [-]"
    )
  )
}


#' ggplot Berezantsev's chart for pile tip bearing capacity factor
#'
#' @description
#' Create a ggplot chart for pile tip bearing capacity factor Nq as function
#' of angle of internal friction for a sandy soil, according to Berezantsev.
#'
#' This function uses data from a digitised plot, so may not be 100% accurate
#'
#' @param palette RColorBrewer pallete for line colors
#' @param phi list of angles of internal friction (in degrees) for which to
#'   plot crosshairs and show Nq values
#' @param nround number of decimals to use in displayed Nq values in
#'   crosshairs
#' @return ggplot object
#' @examples
#' #empty chart
#' ggplot_piletip_berezantsev()
#'
#' #annotated chart
#' ggplot_piletip_berezantsev(phi = 35)
#' @export

ggplot_piletip_berezantsev <- function(
  palette = "Set1",
  phi = NULL,
  nround = 1
){
  #file path
  file_path <- system.file(
    "extdata",
    "Nq_berezantsev_datapoints_digitised.csv",
    package = "soilmech"
  )
  #load data from file
  df_raw <- readr::read_csv(
    file_path,
    col_types = readr::cols()
  )
  #interpolate at every degree
  df <- tibble::tibble(
    phi = seq(25, 45),
    Nq = stats::approx(df_raw$phi, df_raw$Nq, xout = phi)$y,
    hover_label = paste0("\u03c6' = ", phi, "<br>N<sub>q</sub> = ", round(.data$Nq, 1))
  )
  #axis limits
  xlim <- c(min(df$phi), max(df$phi))
  ylim <- c(10, 1000)
  #create empty ggplot
  plt <- ggplot2::ggplot() +
    theme_soilmech() +
    ggplot2::geom_line(
      data = df,
      ggplot2::aes(x = .data$phi, y = .data$Nq)
    ) +
    ggplot2::scale_x_continuous(
      lim = xlim,
      expand = c(0, 0)
    ) +
    ggplot2::scale_y_log10(
      lim = ylim,
      breaks = c(10, 100, 1000),
      minor_breaks = get_log10_minorbreaks(ylim),
      expand = c(0, 0)) +
    ggplot2::annotation_logticks(sides = "l") +
    ggplot2::xlab(expression(phi*"'"~"["*degree*"]")) +
    ggplot2::ylab(expression(N[q]~"[-]"))
  #annotate with crosshairs
  if (!is.null(phi)){
    da <- tibble::tibble(
      phi = phi,
      Nq = stats::approx(df_raw$phi, df_raw$Nq, xout = phi)$y,
      group = paste0("N[q]==", round(.data$Nq, nround))
    )
    plt <- ggplot_addcrosshairs(
      plt,
      x = da$phi,
      y = da$Nq,
      group = da$group,
      add_colour_scale = TRUE,
      xlim = 25,
      ylim = 10,
      arrow_x = FALSE,
      label_parse = TRUE,
      legend_position = "none"
    )
  }
  #return
  return(plt)
}


#' plotly Poulos's chart for friction along pile in sand
#'
#' @description
#' Create a plotly chart for Ks*tan(delta) factors for resistance along pile
#' shafts in sand, according to Poulos
#'
#' This function uses data from a digitised plot, so may not be 100% accurate
#'
#' @param palette RColorBrewer pallete for line colors
#' @param nround number of decimals to round values in hoverlabels to
#' @importFrom rlang .data
#' @return plotly object
#' @examples
#' plotly_pileshaft_poulos()
#' @export

plotly_pileshaft_poulos <- function(
  palette = "Set1",
  nround = 2
){
  #file path
  file_path <- system.file(
    "extdata",
    "Kstandelta_poulos_datapoints_digitised.csv",
    package = "soilmech"
  )
  #load data from file
  df_raw <- readr::read_csv(
    file_path,
    col_types = readr::cols()
  )
  #interpolate at every degree
  df <- df_raw %>%
    dplyr::group_by(.data$Type) %>%
    dplyr::summarise(
      Kstand = stats::approx(.data$phi, .data$Kstand, xout = seq(32.5, 40, 0.25))$y,
      phi = seq(32.5, 40, 0.25)
    ) %>%
    tidyr::drop_na() %>%
    dplyr::mutate(
      phi2 = ifelse(.data$Type == "Jacked", .data$phi, NA),
      Kstand_label = paste0(.data$Type, ": K<sub>s</sub>tan\u3b4 = ", round(.data$Kstand, nround)),
      phi_label = paste0("\u03c6' = ", round(.data$phi2, nround), "\u00b0"),
    )
  #plotly
  plt <- plotly::plot_ly(
    data = df,
    x = ~phi,
    y = ~Kstand,
    type = "scatter",
    mode = "lines",
    color = ~Type,
    colors = RColorBrewer::brewer.pal(3, palette),
    text = ~Kstand_label,
    hoverinfo = "text"
  )
  #add empty trace for phi
  plt <- plotly::add_trace(
    plt,
    type = "scatter",
    mode = "lines",
    name = "phi",
    data = df,
    x = ~phi2,
    y = 0,
    line = list(color = "black"),
    text = ~phi_label,
    hoverinfo = "text",
    showlegend = FALSE,
    opacity = 0
  )
  #add layout
  plotly::layout(
    plt,
    xaxis = list(
      range = c(32, 40),
      title = "\u03c6' [\u00b0]"
    ),
    yaxis = list(
      range = c(0, 1.5),
      title = "K<sub>s</sub>tan\u3b4 [-]"
    ),
    hovermode = "y_unified",
    showlegend = TRUE,
    template = "none"
  )
}


#' ggplot Poulos's chart for friction along pile in sand
#'
#' @description
#' Create a ggplot chart for Ks*tan(delta) factors for resistance along pile
#' shafts in sand, according to Poulos
#'
#' This function uses data from a digitised plot, so may not be 100% accurate
#'
#' @param palette RColorBrewer pallete for line colors
#' @param phi value of angle of internal friction for crosshairs
#' @param type pile type for crosshairs. Value should be`Bored`, `Driven` or
#'   `Jacked`
#' @param nround number of digits to use for rounding Ks*tan(delta)
#' @return plotly object
#' @examples
#' #empty chart
#' ggplot_pileshaft_poulos()
#'
#' #annotated chart
#' ggplot_pileshaft_poulos(phi = 36, type = "Driven")
#' @export

ggplot_pileshaft_poulos <- function(
  palette = "Set1",
  phi = NULL,
  type = NULL,
  nround = 2
){
  #file path
  file_path <- system.file(
    "extdata",
    "Kstandelta_poulos_datapoints_digitised.csv",
    package = "soilmech"
  )
  #load data from file
  df_raw <- readr::read_csv(
    file_path,
    col_types = readr::cols()
  )
  #create ggplot
  plt <- ggplot2::ggplot() +
    theme_soilmech() +
    ggplot2::geom_line(
      data = df_raw,
      ggplot2::aes(x = .data$phi, y = .data$Kstand, color = .data$Type)
    ) +
    ggplot2::coord_cartesian(
      xlim = c(32, 40),
      ylim = c(0, 1.5),
      expand = FALSE
    ) +
    ggplot2::xlab(expression(phi*minute~"["*degree*"]")) +
    ggplot2::ylab(expression(K[s]*tan*delta~"[-]")) +
    ggplot2::scale_color_brewer(name = "", palette = palette)
  #annotate with crosshairs
  if (!is.null(phi)) {
    da <- tibble::tibble(
      phi = phi,
      Type = type,
    ) %>%
      dplyr::mutate(
        Kstand = purrr::map2_dbl(
          .data$phi, .data$Type,
          ~approx(
            df_raw$phi[df_raw$Type == .y],
            df_raw$Kstand[df_raw$Type == .y],
            xout = .x
          )$y
        ),
        group = paste0("K[s]*tan*delta==", round(.data$Kstand, nround))
      )
    plt <- ggplot_addcrosshairs(
      plt,
      x = da$phi,
      y = da$Kstand,
      group = da$group,
      add_colour_scale = FALSE,
      xlim = 32,
      arrow_x = FALSE,
      label_parse = TRUE
    )
  }
  #return
  return(plt)
}


#' ggplot Poulos's chart for friction along pile in sand
#'
#' @description
#' Create a ggplot chart for the stiffness under the pile tip as used in the
#' Hyperbolic method by Fleming (1992). Data digitised from Craig's Soil
#' Mechanics
#'
#' This function uses data from a digitised plot, so may not be 100% accurate
#'
#' @return ggplot object
#' @examples
#' #plot chart
#' ggplot_hyperbolicmethod_stiffness()
#' @export

ggplot_hyperbolicmethod_stiffness <- function(){
  #file path
  file_path <- system.file(
    "extdata",
    "Hyperbolic_method_soil_stiffness_digitised.csv",
    package = "soilmech"
  )
  #load data from file
  df_raw <- readr::read_csv(
    file_path,
    col_types = readr::cols()
  )
  #change factor levels - to get correct order in plot
  soiltype_unique <- unique(df_raw$soiltype)
  df_raw$soiltype <- as.factor(df_raw$soiltype)
  df_raw$soiltype <- factor(df_raw$soiltype, levels = soiltype_unique)
  consistency_unique <- unique(df_raw$consistency)
  df_raw$consistency <- as.factor(df_raw$consistency)
  df_raw$consistency <- factor(df_raw$consistency, levels = consistency_unique)
  #create ggplot   - plt
  plt <- ggplot2::ggplot() +
    theme_soilmech() +
    ggplot2::geom_errorbar(
      data = df_raw,
      ggplot2::aes(x = .data$consistency, ymin = .data$`Eb_min [MPa]`, ymax = .data$`Eb_max [MPa]`),
      width = 0.5
    ) +
    ggplot2::facet_wrap(~.data$soiltype, scales = "free_x") +
    ggplot2::xlab(NULL) +
    ggplot2::ylab(expression(E[b]~"[MPa]")) +
    ggplot2::scale_y_continuous(lim = c(0, 200), expand = c(0, 0)) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1))
  #return
  return(plt)
}


#' Calculate pile displacements using the hyperbolic method
#'
#' @description
#' Calculates the settlement of a rigid pile using Fleming's (1992)
#' hyperbolic method for pile displacements.
#'
#' @param Q vector of loads applied to pile
#' @param Qs pile shaft capacity
#' @param Qb pile base capacity
#' @param Eb soil Young's modulus underneath the pile tip
#' @param Ds pile shaft diameter
#' @param Db pile base diameter. If not defined, this is assumed to be equal
#'   to the shaft diameter
#' @param Ms empirical factor, assumed as Ms = 0.001 by default
#' @examples
#' # Input parameters
#' Qs <- 1000  # [kN]
#' Qb <- 1500  # [kN]
#' Q <- seq(0, Qs + Qb, l = 101)  # [kN]
#' Eb <- 20*1e3  # [kPa]
#' Ds <- 0.5  # [m]
#'
#' # Calculate settlement
#' s <- hyperbolicmethod_settlement(Q, Qb, Qs, Eb, Ds)  # [m]
#'
#' # plot
#' plot(Q, s)
#' @export

hyperbolicmethod_settlement <- function(Q, Qb, Qs, Eb, Ds, Db = NULL, Ms = 0.001) {
  # pile not underreamed
  if (is.null(Db)) {
    Db <- Ds
  }
  # calculate parameters
  alp <- Qs
  bet <- Db*Qb*Eb
  del <- 0.6*Qb
  lam <- Ms*Ds
  eta <- Db*Eb
  a <- eta*(Q - alp) - bet
  b <- Q*(del + lam*eta) - alp*del - bet*lam
  c <- lam*del*Q
  # solve quadratic equation
  D <- sqrt(b^2 - 4*a*c)
  s <- ifelse(D/a > 0, (-b + D)/(2*a), (-b - D)/(2*a))
  # return settlement
  return(s)
}


#' ggplot exampe pile settlement using the hyperbolic method
#'
#' @description
#' Returns an example pile load-settlement graph for a rigid pile using
#' Fleming's (1992) hyperbolic method for pile displacements.
#'
#' @param Q Example load applied to pile
#' @param Qs pile shaft capacity
#' @param Qb pile base capacity
#' @param Eb soil Young's modulus underneath the pile tip
#' @param Ds pile shaft diameter
#' @param Db pile base diameter. If not defined, this is assumed to be equal
#'   to the shaft diameter
#' @param Ms empirical factor, assumed as Ms = 0.001 by default
#' @param ylim vector with upper and lower plot settlement limits
#' @param palette RColorBrewer color palette to use for annotations
#' @examples
#' # Input parameters
#' Qs <- 1000  # [kN]
#' Qb <- 1500  # [kN]
#' Q <- 2000  # [kN]
#' Eb <- 50e3  # [kPa]
#' Ds <- 0.5  # [m]
#'
#' # plot
#' ggplot_hyperbolicmethod(Qb = Qb, Qs = Qs, Q = Q, Eb = Eb, Ds = Ds)
#' @export

ggplot_hyperbolicmethod <- function(
  Qb = 1500,
  Qs = 1000,
  Q = 2000,
  Eb = 50e3,
  Ds = 0.50,
  Db = NULL,
  Ms = 0.001,
  ylim = c(0, 0.5),
  palette = "Set1"
) {
  # no loads applied - create vector
  Qp <- seq(0, Qb + Qs, l = 251)[1:250]
  # calculate settlements
  sp <- hyperbolicmethod_settlement(Qp, Qb, Qs, Eb, Ds, Db = Db, Ms = Ms)
  # colors
  colo <- RColorBrewer::brewer.pal(3, palette)
  # plot
  plt <- ggplot2::ggplot() +
    soilmech::theme_soilmech() +
    ggplot2::geom_path(ggplot2::aes(x = Qp, y = sp)) +
    ggplot2::xlab("Pile load Q [kN]") +
    ggplot2::ylab("Displacement s [m]") +
    ggplot2::scale_x_continuous(lim = c(0, Qb + Qs), expand = ggplot2::expansion(mult = c(0, 0.025))) +
    ggplot2::coord_cartesian(ylim = rev(ylim)) +
    ggplot2::scale_y_reverse(expand = c(0, 0))
  # add max load
  plt <- plt +
    ggplot2::geom_vline(xintercept = Qb + Qs, color = colo[2], linetype = 2) +
    ggplot2::annotate(
      "text",
      x = 0.99*(Qb + Qs), y = 0.5*(ylim[1] + ylim[2]),
      label = "Q[u]", parse = TRUE,
      hjust = 1, vjust = 0.5, color = colo[2]
    )
  # add example load
  if (!is.null(Q)) {
    # settlement at example load
    s <- hyperbolicmethod_settlement(Q, Qb, Qs, Eb, Ds, Db = Db, Ms = Ms)
    # add annotation for applied load
    plt <- ggplot_addcrosshairs(
        plt,
        x = Q, y = s, ylim = ylim[2],
        arrow_x = FALSE, add_colour_scale = TRUE, palette = palette
      ) +
      ggplot2::annotate(
        "text",
        x = Q - 0.01*(Qb + Qs), y = 0.5*(s + ylim[2]),
        label = "Q", parse = TRUE,
        hjust = 1, vjust = 0.5, color = colo[1]
      ) +
      ggplot2::annotate(
        "text",
        x = 0.5*Q, y = s + 0.01*(ylim[2] - ylim[1]),
        label = "s", parse = TRUE,
        hjust = 0.5, vjust = 1, color = colo[1]
      )
  }
  # return plot
  return(plt)
}
