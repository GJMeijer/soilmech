#' Calculate profiles of horizontal effective stresses
#'
#' @description
#' Function takes two stress profiles, one corresponding with the largest
#' historical vertical effective stress and one with the current stresses,
#' and generates profiles with horizontal stresses based on soil properties
#'
#' @param z_soil depth of tops of current soil layer,
#' @param z_watertable current water table level,
#' @param q current surcharge,
#' @param gamma_b current bulk unit weights
#' @param z_soil0 historic depths of tops of soil layers
#' @param z_watertable0 historic water table
#' @param q0 historic surcharge
#' @param gamma_b0 historic bulk unit weights
#' @param z_max maximum depth to consider
#' @param z_interval depth interval at which to calculate stresses
#' @param gamma_w unit weight of water
#' @param phi_deg angle of internal friction for current soil layers
#'   (in degrees)
#' @param description type of current soil layers, either `sand` or `clay`.
#'   used to determine K0
#' @importFrom magrittr `%>%`
#' @importFrom rlang .data
#' @return a tibble with many values, e.g. `sigma'_v` and `sigma'_v0` for
#'   current and historic vertical effective stresses, `OCR` for
#'   overconsolidation ratio, `K0NC` and `K0` for normally consolidated and
#'   current coefficient of lateral earth pressure, and `sigma'_h` and
#'   `sigma'_hNC` for current and normally consolidated horizontal effective
#'   stress
#' @examples
#' horizontal_stress_profile()
#' @export

horizontal_stress_profile <- function(
  z_soil = c(0, 4, 6),
  z_watertable = 4,
  q = 0,
  gamma_b = c(18, 20, 22),
  z_soil0 = NULL,
  z_watertable0 = NULL,
  q0 = NULL,
  gamma_b0 = NULL,
  z_max = 10,
  z_interval = 0.25,
  gamma_w = 10,
  phi_deg = c(30, 20, 25),
  description = "sand"
){
  #assume soil profile same as current, if no input provided
  if (is.null(z_soil0)) {
    z_soil0 <- z_soil
  }
  if (is.null(z_watertable0)) {
    z_watertable0 <- z_watertable
  }
  if (is.null(q0)) {
    q0 <- q
  }
  if (is.null(gamma_b0)) {
    gamma_b0 <- gamma_b
  }
  #properties for soil
  dp <- tibble::tibble(
    z0 = pmin(z_max, z_soil),
    z1 = pmin(z_max, c(z_soil[2:length(z_soil)], z_max)),
    phi_deg = phi_deg,
    description = description
  ) %>%
    dplyr::filter(.data$z0 < z_max) %>%
    dplyr::mutate(
      layer = seq(dplyr::n()),
      K0NC = ifelse(
        stringr::str_detect(
          stringr::str_to_lower(.data$description),
          "clay$"
        ),
        0.95 - sin(.data$phi_deg/180*pi),
        1.00 - sin(.data$phi_deg/180*pi)
      )
    )
  #calculate vertical effective stresses
  dsig0 <- vertical_stress_profile(
    z_soil,
    z_watertable = z_watertable0,
    z_max = z_max,
    gamma_b = gamma_b0,
    gamma_w = gamma_w,
    q = q0
  )
  dsig <- vertical_stress_profile(
    z_soil,
    z_watertable = z_watertable,
    z_max = z_max,
    gamma_b = gamma_b,
    gamma_w = gamma_w,
    q = q,
  )
  #get all depth levels + K0NC
  ds <- tibble::tibble(
    z = unique(c(seq(z_soil[1], z_max, z_interval), z_max)),
    layer = floor(stats::approx(dp$z0, dp$layer, xout = .data$z, yright = max(dp$layer))$y)
  ) %>% dplyr::bind_rows(
    tibble::tibble(
      z = c(utils::tail(z_soil[z_soil < z_max], -1), z_max),
      layer = seq(nrow(dp))
    )
  ) %>%
    unique() %>%
    dplyr::arrange(.data$layer, .data$z) %>%
    dplyr::left_join(
      dp %>% dplyr::select(.data$layer, .data$K0NC, .data$phi_deg, .data$description),
      by = "layer"
    ) %>%
    dplyr::mutate(
      `sigma'_v0` = stats::approx(dsig0$z, dsig0$`sigma'_v`, xout = .data$z)$y,
      `sigma'_v` = stats::approx(dsig$z, dsig$`sigma'_v`, xout = .data$z)$y,
      OCR = pmax(1, .data$`sigma'_v0`/.data$`sigma'_v`),
      K0 = .data$K0NC*.data$OCR^sin(.data$phi_deg/180*pi),
      `sigma'_h` = ifelse(dplyr::near(.data$`sigma'_v`, 0), 0, .data$K0*.data$`sigma'_v`),
      `sigma'_hNC` = .data$K0NC*.data$`sigma'_v`
    )
  #return
  return(ds)
}


#' Plotly historic and current vertical stress traces
#'
#' @description
#' Plots two traces of vertical effective stress in a single plotly plot
#'
#' @param z array width depths
#' @param sigma_v0 array with previously largests vertical effective stress
#'   at depths `z`
#' @param sigma_v array with current vertical effective stress at depths `z`
#' @param ylab depth axis label
#' @param xlim manually override stress axis limits
#' @param ylim manually override depth axis limits
#' @param mode plotly plotting mode, e.g. `lines` or `lines+markers`
#' @param line_colour array with two colors
#' @param line_type array with two plotly line styles
#' @param line_name names for the two traces
#' @param line_width line width
#' @param nround number of digits for rounding hoverlabels
#' @importFrom magrittr `%>%`
#' @return a plotly object
#' @examples
#' ds <- horizontal_stress_profile(q0 = 50)
#' plotly_stressprofile_horizontal_sigmav(
#'   ds$z,
#'   ds$`sigma'_v0`,
#'   ds$`sigma'_v`
#' )
#' @export

plotly_stressprofile_horizontal_sigmav <- function(
  z,
  sigma_v0,
  sigma_v,
  ylab = "Depth [m]",
  xlim = c(0, NA),
  ylim = c(NA, NA),
  mode = "lines",
  line_colour = c("#ff0000", "#ff0000"),
  line_type = c("dash", "solid"),
  line_name = c("Historic", "Current"),
  line_width = 1.5,
  nround = 2
){
  plotly::plot_ly() %>%
    plotly::add_trace(
      type = "scatter",
      mode = "lines",
      x = xlim[1],
      y = z,
      line = list(
        color = "white"
      ),
      text = ~paste0("depth = ", round(z, nround), " m"),
      hoverinfo = "text",
      showlegend = FALSE,
      opacity = 0
    ) %>%
    plotly::add_trace(
      name = line_name[1],
      type = "scatter",
      mode = mode,
      x = sigma_v0,
      y = z,
      line = list(
        color = line_colour[1],
        dash = line_type[1],
        width = line_width
      ),
      text = paste0(
        line_name[1],
        "<br>\u03c3'<sub>v,0</sub> = ",
        round(sigma_v0, nround),
        " kPa"
      ),
      hoverinfo = "text",
      showlegend = FALSE
    ) %>%
    plotly::add_trace(
      name = line_name[2],
      type = "scatter",
      mode = mode,
      x = sigma_v,
      y = z,
      line = list(
        color = line_colour[2],
        dash = line_type[2],
        width = line_width
      ),
      text = paste0(
        line_name[2],
        "<br>\u03c3'<sub>v</sub> = ",
        round(sigma_v, nround),
        " kPa"
      ),
      hoverinfo = "text",
      showlegend = FALSE
    ) %>%
    plotly::add_annotations(
      x = sigma_v0[round(0.5*length(z))],
      y = z[round(0.5*length(z))],
      text = line_name[1],
      xref = "x",
      yref = "y",
      showarrow = TRUE,
      arrowhead = 7,
      ax = 20,
      ay = -40
    ) %>%
    plotly::add_annotations(
      x = sigma_v[round(0.5*length(z))],
      y = z[round(0.5*length(z))],
      text = line_name[2],
      xref = "x",
      yref = "y",
      showarrow = TRUE,
      arrowhead = 7,
      ax = -20,
      ay = 40
    ) %>%
    plotly::layout(
      xaxis = list(
        range = round_limits(c(sigma_v, sigma_v0), lower = xlim[1], upper = xlim[2]),
        title = "\u03c3'<sub>v</sub> [kPa]",
        side = "top",
        visible = TRUE
      ),
      yaxis = list(
        range = rev(round_limits(z, lower = ylim[1], upper = ylim[2])),
        title = ylab
      ),
      hovermode = "y"
    )
}


#' Plotly OCR trace with depth
#'
#' @description
#' Plots trace of OCR with depth, and add line for OCR=1 (normally
#' consolidated soil)
#'
#' @param z array width depths
#' @param OCR array with overconsolidation ratio's for each depth `z`
#' @param ylab depth axis label
#' @param xlim manually override stress axis limits
#' @param ylim manually override depth axis limits
#' @param mode plotly plotting mode, e.g. `lines` or `lines+markers`
#' @param line_colour array with two colors
#' @param line_type array with two plotly line styles
#' @param line_name names for the two traces
#' @param line_width line width
#' @param nround number of digits for rounding hoverlabels
#' @importFrom magrittr `%>%`
#' @return a plotly object
#' @examples
#' ds <- horizontal_stress_profile(q0 = 50)
#' plotly_stressprofile_horizontal_ocr(
#'   ds$z,
#'   ds$OCR
#' )
#' @export

plotly_stressprofile_horizontal_ocr <- function(
  z,
  OCR,
  ylab = "Depth [m]",
  xlim = c(0, NA),
  ylim = c(NA, NA),
  mode = "lines",
  line_colour = c("#FF0000", "#707070"),
  line_type = c("solid", "dash"),
  line_name = c("Current", "Normally consolidated"),
  line_width = 1.5,
  nround = 2
){
  plotly::plot_ly() %>%
    plotly::add_trace(
      type = "scatter",
      mode = "lines",
      x = xlim[1],
      y = z,
      line = list(
        color = "white"
      ),
      text = ~paste0("depth = ", round(z, nround), " m"),
      hoverinfo = "text",
      showlegend = FALSE,
      opacity = 0
    ) %>%
    plotly::add_trace(
      name = line_name[1],
      type = "scatter",
      mode = mode,
      x = OCR,
      y = z,
      line = list(
        color = line_colour[1],
        dash = line_type[1],
        width = line_width
      ),
      text = paste0(
        line_name[1],
        "<br>OCR = ",
        round(OCR, nround)
      ),
      hoverinfo = "text",
      showlegend = FALSE
    ) %>%
    plotly::add_trace(
      name = line_name[2],
      type = "scatter",
      mode = mode,
      x = 1,
      y = z,
      line = list(
        color = line_colour[2],
        dash = line_type[2],
        width = line_width
      ),
      text = paste0(
        line_name[2],
        "<br>OCR = 1"
      ),
      hoverinfo = "text",
      showlegend = FALSE
    ) %>%
    plotly::add_annotations(
      x = OCR[round(0.5*length(z))],
      y = z[round(0.5*length(z))],
      text = line_name[1],
      xref = "x",
      yref = "y",
      showarrow = TRUE,
      arrowhead = 7,
      ax = 40,
      ay = 20
    ) %>%
    #plotly::add_annotations(
    #  x = 1,
    #  y = z[round(0.5*length(z))],
    #  text = line_name[2],
    #  xref = "x",
    #  yref = "y",
    #  showarrow = TRUE,
    #  arrowhead = 7,
    #  ax = 40,
    #  ay = 40
    #) %>%
    plotly::layout(
      xaxis = list(
        range = round_limits(pmin(10, c(1, OCR)), lower = xlim[1], upper = xlim[2]),
        title = "OCR [-]",
        side = "top",
        visible = TRUE
      ),
      yaxis = list(
        range = rev(round_limits(z, lower = ylim[1], upper = ylim[2])),
        title = ylab
      ),
      hovermode = "y"
    )
}


#' Plotly K0 traces with depth
#'
#' @description
#' Plots traces of K0 (current coefficient of lateral earth pressure) and K0NC
#' (coefficient of lateral earth pressure in the normally consolidated case)
#' with depth in a single plotly plot
#'
#' @param z array width depths
#' @param K0 current coefficient of lateral pressure for each depth `z`
#' @param K0NC normally consolidated coefficient of lateral pressure for each
#'   depth `z`
#' @param ylab depth axis label
#' @param xlim manually override stress axis limits
#' @param ylim manually override depth axis limits
#' @param mode plotly plotting mode, e.g. `lines` or `lines+markers`
#' @param line_colour array with two colors
#' @param line_type array with two plotly line styles
#' @param line_name names for the two traces
#' @param line_width line width
#' @param nround number of digits for rounding hoverlabels
#' @importFrom magrittr `%>%`
#' @return a plotly object
#' @examples
#' ds <- horizontal_stress_profile(q0 = 50)
#' plotly_stressprofile_horizontal_k0(
#'   ds$z,
#'   ds$K0,
#'   ds$K0NC
#' )
#' @export

plotly_stressprofile_horizontal_k0 <- function(
  z,
  K0,
  K0NC,
  ylab = "Depth [m]",
  xlim = c(0, NA),
  ylim = c(NA, NA),
  mode = "lines",
  line_colour = c("#FF0000", "#707070"),
  line_type = c("solid", "dash"),
  line_name = c("Current", "Normally consolidated"),
  line_width = 1.5,
  nround = 3
){
  plotly::plot_ly() %>%
    plotly::add_trace(
      type = "scatter",
      mode = "lines",
      x = xlim[1],
      y = z,
      line = list(
        color = "white"
      ),
      text = ~paste0("depth = ", round(z, nround), " m"),
      hoverinfo = "text",
      showlegend = FALSE,
      opacity = 0
    ) %>%
    plotly::add_trace(
      name = line_name[1],
      type = "scatter",
      mode = mode,
      x = K0,
      y = z,
      line = list(
        color = line_colour[1],
        dash = line_type[1],
        width = line_width
      ),
      text = paste0(
        line_name[1],
        "<br>K<sub>0</sub> = ",
        round(K0, nround)
      ),
      hoverinfo = "text",
      showlegend = FALSE
    ) %>%
    plotly::add_trace(
      name = line_name[2],
      type = "scatter",
      mode = mode,
      x = K0NC,
      y = z,
      line = list(
        color = line_colour[2],
        dash = line_type[2],
        width = line_width
      ),
      text = paste0(
        line_name[2],
        "<br>K<sub>0,NC</sub> = ",
        round(K0NC, nround)
      ),
      hoverinfo = "text",
      showlegend = FALSE
    ) %>%
    plotly::add_annotations(
      x = K0[round(0.5*length(z))],
      y = z[round(0.5*length(z))],
      text = line_name[1],
      xref = "x",
      yref = "y",
      showarrow = TRUE,
      arrowhead = 7,
      ax = 40,
      ay = 20
    ) %>%
    #plotly::add_annotations(
    #  x = 1,
    #  y = z[round(0.5*length(z))],
    #  text = line_name[2],
    #  xref = "x",
    #  yref = "y",
    #  showarrow = TRUE,
    #  arrowhead = 7,
    #  ax = 40,
    #  ay = 40
    #) %>%
  plotly::layout(
    xaxis = list(
      range = round_limits(c(K0, K0NC), lower = xlim[1], upper = xlim[2]),
      title = "K<sub>0</sub> [-]",
      side = "top",
      visible = TRUE
    ),
    yaxis = list(
      range = rev(round_limits(z, lower = ylim[1], upper = ylim[2])),
      title = ylab
    ),
    hovermode = "y"
  )
}


#' Plotly horizontal stress profile
#'
#' @description
#' Plots two traces of horizontal effective stress in a single plotly plot.
#' One trace is for current state (may be overconsolidated), and one for the
#' normally consolidated case.
#'
#' @param z array width depths
#' @param sigma_h current horizontal effective stress for each depth `z`
#' @param sigma_hNC normally consolidated horizontal effective stress for
#'   each depth `z`
#' @param ylab depth axis label
#' @param xlim manually override stress axis limits
#' @param ylim manually override depth axis limits
#' @param mode plotly plotting mode, e.g. `lines` or `lines+markers`
#' @param line_colour array with two colors
#' @param line_type array with two plotly line styles
#' @param line_name names for the two traces
#' @param line_width line width
#' @param nround number of digits for rounding hoverlabels
#' @importFrom magrittr `%>%`
#' @return a plotly object
#' @examples
#' ds <- horizontal_stress_profile(q0 = 50)
#' plotly_stressprofile_horizontal_sigmah(
#'   ds$z,
#'   ds$`sigma'_h`,
#'   ds$`sigma'_hNC`
#' )
#' @export

plotly_stressprofile_horizontal_sigmah <- function(
  z,
  sigma_h,
  sigma_hNC,
  ylab = "Depth [m]",
  xlim = c(0, NA),
  ylim = c(NA, NA),
  mode = "lines",
  line_colour = c("#FF0000", "#707070"),
  line_type = c("solid", "dash"),
  line_name = c("Current", "Normally consolidated"),
  line_width = 1.5,
  nround = 2
){
  plotly::plot_ly() %>%
    plotly::add_trace(
      type = "scatter",
      mode = "lines",
      x = xlim[1],
      y = z,
      line = list(
        color = "white"
      ),
      text = ~paste0("depth = ", round(z, nround), " m"),
      hoverinfo = "text",
      showlegend = FALSE,
      opacity = 0
    ) %>%
    plotly::add_trace(
      name = line_name[1],
      type = "scatter",
      mode = mode,
      x = sigma_h,
      y = z,
      line = list(
        color = line_colour[1],
        dash = line_type[1],
        width = line_width
      ),
      text = paste0(
        line_name[1],
        "<br>\u03c3'<sub>h</sub> = ",
        round(sigma_h, nround),
        " kPa"
      ),
      hoverinfo = "text",
      showlegend = FALSE
    ) %>%
    plotly::add_trace(
      name = line_name[2],
      type = "scatter",
      mode = mode,
      x = sigma_hNC,
      y = z,
      line = list(
        color = line_colour[2],
        dash = line_type[2],
        width = line_width
      ),
      text = paste0(
        line_name[2],
        "<br>\u03c3'<sub>h</sub> = ",
        round(sigma_hNC, nround),
        " kPa"
      ),
      hoverinfo = "text",
      showlegend = FALSE
    ) %>%
    plotly::add_annotations(
      x = sigma_h[round(0.5*length(z))],
      y = z[round(0.5*length(z))],
      text = line_name[1],
      xref = "x",
      yref = "y",
      showarrow = TRUE,
      arrowhead = 7,
      ax = 40,
      ay = -20
    ) %>%
    #plotly::add_annotations(
    #  x = 1,
    #  y = z[round(0.5*length(z))],
    #  text = line_name[2],
    #  xref = "x",
    #  yref = "y",
    #  showarrow = TRUE,
    #  arrowhead = 7,
    #  ax = 40,
    #  ay = 40
    #) %>%
  plotly::layout(
    xaxis = list(
      range = round_limits(c(sigma_h, sigma_hNC), lower = xlim[1], upper = xlim[2]),
      title = "\u03c3'<sub>h</sub> [kPa]",
      side = "top",
      visible = TRUE
    ),
    yaxis = list(
      range = rev(round_limits(z, lower = ylim[1], upper = ylim[2])),
      title = ylab
    ),
    hovermode = "y"
  )
}


#' Plotly calculation process of current horizontal effective stresses
#'
#' @description
#' Function takes data for two soil profiles (current and historical when
#' effective stress was largests) and draws soil profiles, profiles with
#' vertical effective stress, OCR, K0 and horizontal effective stress with
#' depth
#'
#' @inheritParams horizontal_stress_profile
#' @param ylab depth axis label
#' @param mode plotly trace mode, e.g. `lines`
#' @param nrow number of rows for subplots
#' @param width width of plot
#' @param height height of plot
#' @return a plotly object
#' @examples
#' #removal of surcharge
#' plotly_soilstressprofile_horizontal(q = 0, q0 = 100)
#'
#' #raising of the water table
#' plotly_soilstressprofile_horizontal(q = 0, q0 = 0, z_watertable = 0, z_watertable0 = 6)
#' @export

plotly_soilstressprofile_horizontal <- function(
  z_soil = c(0, 4),
  z_watertable = 4,
  q = 0,
  gamma_b = c(18, 20),
  z_soil0 = NULL,
  z_watertable0 = NULL,
  q0 = 100,
  gamma_b0 = NULL,
  z_max = 10,
  z_interval = 0.1,
  gamma_w = 10,
  phi_deg = c(35, 15),
  description = c("sand", "clay"),
  ylab = "Depth [m]",
  mode = "lines",
  nrow = 2,
  width = 800,
  height = 600
){
  #assume soil profile same as current, if no input provided
  if (is.null(z_soil0)) {
    z_soil0 <- z_soil
  }
  if (is.null(z_watertable0)) {
    z_watertable0 <- z_watertable
  }
  if (is.null(q0)) {
    q0 <- q
  }
  if (is.null(gamma_b0)) {
    gamma_b0 <- gamma_b
  }
  #calculate stresses
  ds <- horizontal_stress_profile(
    z_soil = z_soil,
    z_watertable = z_watertable,
    q = q,
    gamma_b = gamma_b,
    z_soil0 = z_soil0,
    z_watertable0 = z_watertable0,
    q0 = q0,
    gamma_b0 = gamma_b0,
    z_max = z_max,
    z_interval = z_interval,
    gamma_w = gamma_w,
    phi_deg = phi_deg,
    description = description
  )
  #generate plots
  plt_prof0 <- plotly_soilprofile(
    z_soil0,
    z_watertable = z_watertable0,
    z_max = z_max,
    q = q0,
    gamma_b = gamma_b0
  )
  plt_prof <- plotly_soilprofile(
    z_soil,
    z_watertable = z_watertable,
    z_max = z_max,
    q = q,
    gamma_b = gamma_b,
    fields_hover = c("description", "gamma_b", "phi_deg", "thickness"),
    fields_label = c("description", "gamma_b", "phi_deg"),
    description = description,
    phi_deg = phi_deg
  )
  plt_sigv <- plotly_stressprofile_horizontal_sigmav(
    ds$z,
    ds$`sigma'_v0`,
    ds$`sigma'_v`,
    mode = mode,
    ylab = ylab
  )
  plt_OCR <- plotly_stressprofile_horizontal_ocr(
    ds$z,
    ds$OCR,
    mode = mode,
    ylab = ylab
  )
  plt_K0 <- plotly_stressprofile_horizontal_k0(
    ds$z,
    ds$K0,
    ds$K0NC,
    mode = mode,
    ylab = ylab
  )
  plt_sigh <- plotly_stressprofile_horizontal_sigmah(
    ds$z,
    ds$`sigma'_h`,
    ds$`sigma'_hNC`,
    mode = mode,
    ylab = ylab
  )
  #merge plots together
  plt <- plotly::subplot(
    plt_prof0, plt_prof, plt_sigv, plt_OCR, plt_K0, plt_sigh,
    shareY = TRUE,
    titleX = TRUE,
    margin = 0.08,
    nrows = nrow
  )
  #add titles for profiles
  annotations = list(
    list(
      x = 1/6,
      y = 1.0,
      text = "Historic",
      xref = "paper",
      yref = "paper",
      xanchor = "center",
      yanchor = "bottom",
      showarrow = FALSE
    ),
    list(
      x = 1/2,
      y = 1,
      text = "Current",
      xref = "paper",
      yref = "paper",
      xanchor = "center",
      yanchor = "bottom",
      showarrow = FALSE
    )
  )
  plt <- plotly::layout(
    plt,
    annotations = annotations
  )
  #set plot sizes
  plt$width <- width
  plt$height <- height
  #return
  return(plt)
}
