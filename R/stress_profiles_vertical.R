#' Calculate verticals stresses using layers
#'
#' @description
#' Function takes a soil profile and returns the vertical stresses at all
#' layer interfaces
#'
#' @param z_soil depth of top of each layer
#' @param z_watertable depth of the water table
#' @param z_max maximum depth of interest
#' @param z_interval additional interval at which to calculate stresses
#' @param gamma_b unit weight of soil, in kN/m3. Should be a scalar or a
#'   vector with the same length as `z_soil`
#' @param gamma_w unit weight of water: 10 kN/m3
#' @param q surcharge at soil surface
#' @param ... additional fields to be added to the returned data
#' @importFrom magrittr `%>%`
#' @importFrom rlang .data
#' @return a tibble with depth `z`, vertical total stress `sigma_v`, vertical
#'   effective stress `sigma'_v` and pore pressure `u`. All fields in `...`
#'   are added too
#' @examples
#' #only calculate stresses at layer interfaces and water table
#' vertical_stress_profile(
#'   c(0, 4, 6),
#'   z_watertable = 2,
#'   gamma_b = c(18, 20, 22),
#'   q = 10
#' )
#'
#' #add stresses at regular intervals
#' vertical_stress_profile(
#'   c(0, 4, 6),
#'   z_watertable = 2,
#'   gamma_b = c(18, 20, 22),
#'   q = 10,
#'   z_interval = 1
#' )
#' @export

vertical_stress_profile <- function(
  z_soil = 0,
  z_watertable = 0,
  z_max = 10,
  z_interval = NULL,
  gamma_b = 19,
  gamma_w = 10,
  q = 0,
  ...
){
  #no water table
  if (is.null(z_watertable) | is.na(z_watertable)) {
    z_watertable <- z_max + 10
  }
  #create tibble - total stresses
  df <- tibble::tibble(
    z = c(z_soil, z_max),
    sigma_v = q + c(0, cumsum(gamma_b*diff(c(z_soil, z_max)))) +
      max(0, z_soil[1] - z_watertable[1])*gamma_w
  )
  #add intervals
  if (!is.null(z_interval)) {
    z_add <- seq(z_soil[1] + z_interval, z_max, z_interval)
    z_add <- z_add[!(z_add %in% df$z)]
    df <- dplyr::bind_rows(
      df,
      tibble::tibble(
        z = z_add,
        sigma_v = stats::approx(df$z, df$sigma_v, xout = z_add)$y
      )
    ) %>%
      dplyr::arrange(.data$z)
  }
  #add additional fields if inputted
  if (length(list(...)) > 0) {
    df <- dplyr::bind_cols(df, tibble::as_tibble(list(...)))
  }
  #add total stress at water table
  if (!(z_watertable %in% df$z)) {
    df <- dplyr::bind_rows(
      df,
      tibble::tibble(
        z = z_watertable,
        sigma_v = stats::approx(df$z, df$sigma_v, xout = z_watertable, yleft = 0, yright = 0)$y
      )
    )
  }
  #add extra point for surcharge if needed
  if ((q != 0) & (z_watertable < z_soil[1])) {
    df <- dplyr::bind_rows(
      df,
      tibble::tibble(
        z = z_soil[1],
        sigma_v = max(0, z_soil[1] - z_watertable[1])*gamma_w
      )
    )
  }
  #filter, arrange and add pore pressures and effective stress
  df <- df %>%
    dplyr::filter(.data$z <= z_max) %>%
    dplyr::arrange(.data$z, .data$sigma_v) %>%
    dplyr::mutate(
      u = pmax(0, .data$z - z_watertable)*gamma_w,
      `sigma'_v` = .data$sigma_v - .data$u
    )
  #return
  return(df)
}


#' Generate bookdown table for stress calculatinos
#'
#' @description
#' Generates a tibble that can be converted using knitr to a nicely looking
#' bookdown table
#'
#' @param ds tibble with stress values. Output from function
#'   `vertical_stress_profile()`
#' @param output if `values`, only return stress values. If `calculations`,
#'   it also return how to calculate the stresses
#' @param traces array with traces to be shown. May contain `sigma_v` for
#'   the total vertical stress, `u` for the pore water pressure and `sigma'_v`
#'   for the vertical effective stress
#' @param nround number of decimals to round stress values to
#' @return a tibble
#' @examples
#' ds <- vertical_stress_profile(z_watertable = -2, q = 10)
#' tabulate_vertical_stress_profile(ds)
#' @export

tabulate_vertical_stress_profile <- function(
  ds,
  output = "values",
  traces = c("sigma_v", "u", "sigma'_v"),
  nround = 2
){
  #get columns, and round values
  dt <- dplyr::select(ds, dplyr::all_of(c("z", traces)))
  if (output == "values") {
    #labels
    dt$sigma_v_label <- round(dt$sigma_v, nround)
    dt$u_label <- round(dt$u, nround)
    dt$`sigma'_v_label` <- round(dt$`sigma'_v`, nround)
  } else if (output == "calculations") {
    #initiate labels
    dt$sigma_v_label <- NA
    dt$u_label <- NA
    dt$`sigma'_v_label` <- NA
    #calculations - total stress
    if (dplyr::near(dt$sigma_v[1], 0)) {
      dt$sigma_v_label[1] <- 0
    } else {
      dt$sigma_v_label[1] <- paste(
        "$\\sigma_v = q = ",
        round(dt$sigma_v[1], nround),
        "$"
      )
    }
    for (i in 2:nrow(dt)) {
      if (dplyr::near(dt$z[i - 1], dt$z[i])) {
        dt$sigma_v_label[i] <- paste(
          "$",
          round(dt$sigma_v[i - 1], nround),
          "+",
          round(dt$sigma_v[i] - dt$sigma_v[i - 1], nround),
          "=",
          round(dt$sigma_v[i], nround),
          "$"
        )
      } else {
        dt$sigma_v_label[i] <- paste(
          "$",
          round(dt$sigma_v[i - 1], nround),
          "+",
          round(diff(dt$z[(i - 1):i]), nround),
          "\\cdot",
          round(diff(dt$sigma_v[(i - 1):i])/diff(dt$z[(i - 1):i]), nround),
          "=",
          round(dt$sigma_v[i], nround),
          "$"
        )
      }
    }
    #get unit weight of water
    if (max(dt$u) == 0) {
      gamma_w <- 10
    } else {
      gamma_w <- (max(ds$u) - min(ds$u)) / (ds$z[which.max(ds$u)] - max(ds$z[ds$z == min(ds$z)]))
    }
    #pore water pressures
    dt$u_label <- ifelse(
      dplyr::near(dt$u, 0),
      0,
      paste(
        "$",
        round(dt$u / gamma_w, nround),
        "\\cdot",
        round(gamma_w, nround),
        "=",
        round(dt$u, nround),
        "$"
      )
    )
    #effective stresses
    dt$`sigma'_v_label` <- paste(
      "$",
      round(dt$sigma_v, nround),
      "-",
      round(dt$u, nround),
      "=",
      round(dt$`sigma'_v`, nround),
      "$"
    )
  }
  #select rename columns
  dt <- dt %>%
    dplyr::select(dplyr::all_of(c("z", paste0(traces, "_label")))) %>%
    dplyr::rename(
      "z [m]" = .data$z,
      "$\\sigma_v$ [kPa]" = .data$sigma_v_label,
      "$u$ [kPa]" = .data$u_label,
      "$\\sigma'_v$ [kPa]" = .data$`sigma'_v_label`
    )
  #return
  return(dt)
}


#' Plotly vertical stress profiles
#'
#' @description
#' Plots an interactive graph using plotly to show how soil stresses change
#' with depth
#'
#' @param df dataframe with all depths, stresses etc, outputted from the
#'   function `vertical_stress_profile()`
#' @param traces array with traces to be plotted. May contain `sigma_v` for
#'   the total vertical stress, `u` for the pore water pressure and `sigma'_v`
#'   for the vertical effective stress
#' @param line_colour array with colors for total stress, pore water pressure
#'   and effective stress traces
#' @param line_type array with linetypes for total stress, pore water pressure
#'   and effective stress traces
#' @param line_width line thickness for trace
#' @param mode plotly argument `mode` for each trace. Use `lines+markers` to
#'   plot both lines and markers, or `lines` to only plot lines etc
#' @param labelmode setting controlling what hoverlabels to plot. If
#'   `stressonly`, only the values of stresses are displaced. If
#'   `calculations`, it is also shown how to calculate the stress based on the
#'   stress value in other points
#' @param hovermode value of `hovermode` setting to use in plotly
#' @param ylab y-axis label for depth
#' @param xlim manual override of x-axis limits (stresses)
#' @param ylim manual override of y-axis limits (depths)
#' @param nround number of digits to use in rounding stress values in plotly
#'   hoverlabels
#' @param xlim stress axis limits (manual override)
#' @param ylim depth axis limits (manual override)
#' @param hovermode hovermode to use in plotly
#' @examples
#' df <- vertical_stress_profile(
#'   z_soil = c(0, 4),
#'   gamma_b = c(17, 20),
#'   q = 10,
#'   z_watertable = 2
#' )
#' plotly_stressprofile_vertical(df)
#' @export

plotly_stressprofile_vertical <- function(
  df,
  traces = c("sigma_v", "u", "sigma'_v"),
  line_colour = c("#000000", "#0000ff", "#ff0000"),
  line_type = c("solid", "dash", "dot"),
  line_width = 1.5,
  mode = "lines+markers",
  labelmode = "calculations",
  hovermode = "y",
  ylab = "Depth [m]",
  xlim = c(0, NA),
  ylim = c(NA, NA),
  nround = 2
){
  #initiate plot
  plt <- plotly::plot_ly()
  #hover labels
  if (labelmode == "stressonly") {
    #total stress
    df$sigma_v_label = paste0("\u03c3<sub>v</sub> = ", round(df$sigma_v, nround)," kPa")
    #pore water pressure
    df$u_label = paste0("u = ", round(df$u, nround)," kPa")
    #effective stress
    df$`sigma'_v_label` = paste0("\u03c3'<sub>v</sub> = ", round(df$`sigma'_v`, nround)," kPa")
  } else if (labelmode == "calculations") {
    #total stress
    df$sigma_v_label <- c(
      ifelse(
        (df$sigma_v[1] == 0),
        "\u03c3<sub>v</sub> = 0 kPa",
        paste0("\u03c3<sub>v</sub> = q = ",  round(df$sigma_v[1], nround), " kPa")
      ),
      ifelse(
        (diff(df$z) == 0),
        paste0(
          "\u03c3<sub>v</sub> = ",
          round(df$sigma_v[1:(nrow(df) - 1)], nround),
          " + ",
          round(df$sigma_v[2:nrow(df)] - df$sigma_v[1:(nrow(df) - 1)], nround),
          " = ",
          round(df$sigma_v[2:nrow(df)], nround),
          " kPa"
        ),
        paste0(
          "\u03c3<sub>v</sub> = ",
          round(df$sigma_v[1:(nrow(df) - 1)], nround),
          " + ",
          diff(df$z),
          "\u00b7",
          round(diff(df$sigma_v) / diff(df$z), nround),
          " = ",
          round(df$sigma_v[2:nrow(df)], nround),
          " kPa"
        )
      )
    )
    #pore water pressure
    z_watertable <- max(df$z[dplyr::near(df$u, 0)])
    gamma_w <- max(diff(df$u) / diff(df$z), na.rm = TRUE)
    df$u_label <- ifelse(
      dplyr::near(df$u, 0),
      "u = 0 kPa",
      paste0(
        "u = ",
        df$z - z_watertable,
        "\u00b7",
        round(gamma_w, nround),
        " = ",
        round(df$u, nround),
        " kPa"
      )
    )
    #effective stress
    df$`sigma'_v_label` <- paste0(
      "\u03c3'<sub>v</sub> = \u03c3<sub>v</sub> - u = ",
      round(df$`sigma'_v`, nround)
    )
  }
  #depth label
  df$z_label <- paste0("Depth = ", round(df$z, nround), " m")
  #axis limits
  xlim <- round_limits(max(df[, traces]), lower = xlim[1], upper = xlim[2])
  ylim <- round_limits(df$z, lower = ylim[1], upper = ylim[2])
  #plot invisible trace with depths
  plt <- plotly::add_trace(
    plt,
    type = "scatter",
    mode = "lines",
    name = "phi",
    data = df,
    x = xlim[1],
    y = ~z,
    line = list(color = "white"),
    text = ~z_label,
    hoverinfo = "text",
    showlegend = FALSE,
    opacity = 0
  )
  #plot total stress
  if (("sigma_v" %in% colnames(df)) & ("sigma_v" %in% traces)) {
    plt <- plotly::add_trace(
      plt,
      name = "Vertical total stress",
      type = "scatter",
      mode = mode,
      data = df,
      x = ~`sigma_v`,
      y = ~z,
      hoverinfo = "text",
      text = ~sigma_v_label,
      marker = list(
        color = line_colour[1]
      ),
      line = list(
        width = line_width,
        color = line_colour[1],
        dash = line_type[1]
      )
    )
  }
  #plot pore water pressure
  if (("u" %in% colnames(df)) & ("u" %in% traces)) {
    #plot
    plt <- plotly::add_trace(
      plt,
      name = "Pore water pressure",
      type = "scatter",
      mode = mode,
      data = df,
      x = ~u,
      y = ~z,
      hoverinfo = "text",
      text = ~u_label,
      marker = list(
        color = line_colour[2]
      ),
      line = list(
        width = line_width,
        color = line_colour[2],
        dash = line_type[2]
      )
    )
  }
  #plot effective stress
  if (("sigma'_v" %in% colnames(df)) & ("sigma'_v" %in% traces)) {
    #plot
    plt <- plotly::add_trace(
      plt,
      name = "Vertical effective stress",
      type = "scatter",
      mode = mode,
      data = df,
      x = ~`sigma'_v`,
      y = ~z,
      hoverinfo = "text",
      text = ~`sigma'_v_label`,
      marker = list(
        color = line_colour[3]
      ),
      line = list(
        width = line_width,
        color = line_colour[3],
        dash = line_type[3]
      )
    )
  }
  #axis settings
  if (length(traces) == 1) {
    if (traces == "sigma_v") {
      xlab <- "Vertical total stress [kPa]"
    } else if (traces == "u") {
      xlab <- "Pore water pressure [kPa]"
    } else {
      xlab <- "Vertical effective stress [kPa]"
    }
  } else {
    xlab <- "Stress [kPa]"
  }
  #add layout
  plt <- plotly::layout(
    plt,
    xaxis = list(
      range = xlim,
      title = xlab,
      side = "top",
      visible = TRUE
    ),
    yaxis = list(
      range = rev(ylim),
      title = ylab
    ),
    legend = list(
      x = 1,
      y = 1,
      xanchor = "right",
      yanchor = "top"
    ),
    hovermode = hovermode,
    template = "none"
  )
  #return
  return(plt)
}


#' plotly soil profiles
#'
#' @description
#' Creates a plotly object with the soil profile and soil properties
#' @inheritParams vertical_stress_profile
#' @param description soil description: will be shown in hoverlabel
#' @param phi_deg angle of internal friction: will be shown in text and
#'   hoverlabels
#' @param ylim depth axis limits
#' @param ylab depth axis label
#' @param fill_soildry fill colour for dry soil
#' @param fill_soilwet fill colour for wet soil
#' @param fill_water fill colour for ponding water
#' @param colour_soil line color for soil polygons
#' @param colour_water line color for water table
#' @param colour_surcharge color of surcharge arrows and label
#' @param line_water line type for water table
#' @param line_width line tickness
#' @param nround number of decimals in plotly hoverlabels
#' @param fields_hover field names to show in hoverlabels
#' @param fields_label field names to plot in permanent text labels
#' @param arrow_size relative length of surcharge arrows
#' @param arrow_n number of surcharge arrows to use
#' @param title plot title
#' @param ... additional fields used for plotting labels
#' @examples
#' #ponding water
#' plotly_soilprofile(
#'   z_soil = 0,
#'   z_max = 5,
#'   z_watertable = -5,
#'   gamma_b = 20
#' )
#'
#' #water table within soil
#' plotly_soilprofile(
#'   z_soil = c(0, 3, 7),
#'   gamma_b = c(15, 16, 20),
#'   z_watertable = 2,
#'   description = c("sand", "clay", "gravel"),
#'   phi_deg = c(0, 10, 30),
#'   q = 10,
#'   title = "Test"
#' )
#' @export

plotly_soilprofile <- function(
  z_soil = c(0, 3, 6),
  z_watertable = 3,
  z_max = 10,
  gamma_b = 19,
  gamma_w = 10,
  q = 0,
  description = NULL,
  phi_deg = NULL,
  ylim = c(NA, NA),
  ylab = "Depth [m]",
  fill_soildry = "#d3bc5f",
  fill_soilwet = "#aebab7",
  fill_water = "#2a7fff",
  colour_soil = "#65571d",
  colour_water = "#2a7fff",
  colour_surcharge = "#6a0dad",
  line_water = "dash",
  line_width = 2.0,
  nround = 1,
  fields_hover = c("description", "gamma_b", "phi_deg", "thickness"),
  fields_label = c("gamma_b", "phi_deg"),
  arrow_size = 0.1,
  arrow_n = 6,
  title = NULL,
  ...
){
  #axis limits
  xlim <- c(0, 1)
  ylim <- round_limits(c(z_soil, z_max, min(0, min(z_max, z_watertable))), lower = ylim[1], upper = ylim[2])
  #polygons for all layers
  ds <- tibble::tibble(
    z0 = z_soil,
    z1 = pmin(z_max, c(z_soil[2:length(z_soil)], z_max)),
    gamma_b = gamma_b
  ) %>%
    tidyr::drop_na() %>%
    dplyr::filter(.data$z1 > .data$z0)
  #labels
  if (!is.null(description)) {
    ds$description_label <- description
  } else {
    ds$description_label <- NULL
  }
  if (!is.null(phi_deg)) {
    ds$phi_deg_label <- paste0("\u03c6 = ", round(phi_deg, nround), "\u00b0")
  } else {
    ds$phi_deg_label <- NULL
  }
  ds$description <- description
  ds$gamma_b_label <- paste0("\u03B3<sub>b</sub> = ", ds$gamma_b, " kN/m\u00b3")
  ds$thickness_label <- paste0("thickness = ", ds$z1 - ds$z0, " m")
  #create permanent labels and hoverlabels for each layer
  fields_hover <- paste0(fields_hover, "_label")
  fields_label <- paste0(fields_label, "_label")
  ds <- tidyr::unite(
    ds,
    "hover_label",
    dplyr::all_of(fields_hover[fields_hover %in% colnames(ds)]),
    sep = "<br>",
    remove = FALSE,
    na.rm = TRUE
  ) %>%
    tidyr::unite(
      "label",
      dplyr::all_of(fields_label[fields_label %in% colnames(ds)]),
      sep = "<br>",
      remove = FALSE,
      na.rm = TRUE
    )
  #filter
  ds <- dplyr::filter(ds, .data$z0 < z_max)
  #initiate plot
  plt <- plotly::plot_ly()
  #water polygon - for ponding water on top
  if (z_watertable < z_soil[1]) {
    plt <- plotly::add_trace(
      plt,
      type = "scatter",
      mode = "lines",
      x = c(xlim, rev(xlim)),
      y = c(z_watertable, z_watertable, z_soil[1], z_soil[1]),
      fill = "toself",
      fillcolor = fill_water,
      line = list("width" = 0.0),
      hoverinfo = "none",
      showlegend = FALSE
    )
  }
  #dry soil polygon
  if (z_watertable > z_soil[1]) {
    plt <- plotly::add_trace(
      plt,
      type = "scatter",
      mode = "lines",
      x = c(xlim, rev(xlim)),
      y = c(rep(z_soil[1], 2), rep(min(z_max, z_watertable), 2)),
      fill = "toself",
      fillcolor = fill_soildry,
      line = list("width" = 0.0),
      hoverinfo = "none",
      showlegend = FALSE
    )
  }
  #wet soil polygon
  if (z_watertable < z_max) {
    plt <- plotly::add_trace(
      plt,
      type = "scatter",
      mode = "lines",
      x = c(xlim, rev(xlim)),
      y = c(rep(max(z_soil[1], z_watertable), 2), rep(z_max, 2)),
      fill = "toself",
      fillcolor = fill_soilwet,
      line = list("width" = 0.0),
      hoverinfo = "none",
      showlegend = FALSE
    )
  }
  #plot soil outline polygons
  for (i in 1:nrow(ds)) {
    plt <- plotly::add_trace(
      plt,
      type = "scatter",
      mode = "lines",
      x = c(xlim, rev(xlim)),
      y = c(rep(ds$z0[i], 2), rep(ds$z1[i], 2)),
      line = list(
        width = line_width,
        color = colour_soil
      ),
      fill = "none",
      hoveron = "fills",
      text = ds$hover_label[i],
      hoverinfo = "text",
      showlegend = FALSE
    )
  }
  #add water table
  if (z_watertable < z_max) {
    plt <- plotly::add_trace(
      plt,
      type = "scatter",
      mode = "lines",
      x = c(xlim),
      y = c(z_watertable, z_watertable),
      line = list(
        width = line_width,
        color = colour_water,
        dash = line_water
      ),
      hoverinfo = "none",
      showlegend = FALSE
    )
  }
  #add text in soil boxes
  plt <- plt %>% plotly::add_trace(
    type = "scatter",
    mode = "text",
    x = 0.5*(xlim[2] - xlim[1]),
    y = 0.5*(ds$z0 + ds$z1),
    text = ds$label,
    textposition = "middle center",
    showlegend = FALSE,
    hoverinfo = "skip"
  )
  if (z_watertable < z_soil[1]) {
    plt <- plt %>% plotly::add_trace(
      type = "scatter",
      mode = "text",
      x = 0.5*(xlim[2] - xlim[1]),
      y = 0.5*(z_watertable + z_soil[1]),
      text = paste0("\u03B3<sub>w</sub> = ", gamma_w, " kN/m\u00b3"),
      textposition = "middle center",
      showlegend = FALSE,
      hoverinfo = "skip"
    )
  }
  #surcharge arrows
  if (!dplyr::near(q, 0)) {
    darr <- tibble::tibble(
      x0 = seq(xlim[1], xlim[2] - 1/arrow_n, l = arrow_n) + 0.5/arrow_n,
      x1 = seq(xlim[1], xlim[2] - 1/arrow_n, l = arrow_n) + 0.5/arrow_n,
      y0 = z_soil[1],
      y1 = z_soil[1] - (max(ds$z1) - min(ds$z0))*arrow_size
    )
    #add to plot
    plt <- plotly::add_annotations(
      plt,
      x = darr$x0,
      y = darr$y0,
      xref = "x", yref = "y",
      axref = "x", ayref = "y",
      text = "",
      showarrow = TRUE,
      arrowcolor = colour_surcharge,
      ax = darr$x1,
      ay = darr$y1,
      showlegend = FALSE
    )
    #add surcharge text
    plt <- plotly::add_trace(
      plt,
      type = "scatter",
      mode = "text",
      x = 0.5*(xlim[1] + xlim[2]),
      y = darr$y1[1],
      text = paste0("q = ", round(q, nround), " kPa"),
      textfont = list(
        color = colour_surcharge
      ),
      textposition = "top middle",
      hoverinfo = "skip",
      showlegend = FALSE
    )
    #change axis limits
    ylim[1] <- min(ylim[1], darr$y1[1] - (darr$y0[1] - darr$y1[1]))
  }
  #add layout
  plt <- plt %>% plotly::layout(
    title = title,
    xaxis = list(
      range = xlim,
      title = "position",
      side = "top",
      visible = FALSE
    ),
    yaxis = list(
      range = rev(ylim),
      title = ylab
    ),
    showlegend = FALSE,
    template = "none"
  )
  #return
  return(plt)
}


#' Plotly soil profile and vertical stress profile side-by-side
#'
#' @description
#' Creates a plotly plot with the soil profile on the left-hand side,
#' and the stress profiles on the right-hand side
#'
#' @inheritParams vertical_stress_profile
#' @inheritParams plotly_stressprofile_vertical
#' @param width width of the output plot
#' @param height height of the output plot
#' @examples
#' #dry soil - total stress only
#' plotly_soilstressprofile_vertical(
#'   z_soil = 0,
#'   z_watertable = 20,
#'   z_max = 6,
#'   z_interval = 1,
#'   traces = c("sigma_v")
#' )
#'
#' #partly saturated soil
#' plotly_soilstressprofile_vertical(
#'   z_soil = c(0, 4, 8),
#'   z_watertable = 2,
#'   z_max = 10,
#'   q = 20,
#'   gamma_b = c(18, 20, 15),
#'   traces = c("sigma_v", "u", "sigma'_v"),
#'   labelmode = "calculations"
#' )
#' @return a plotly object
#' @export

#plot together
plotly_soilstressprofile_vertical <- function(
  z_soil,
  z_watertable = 0,
  z_max = 10,
  z_interval = NULL,
  q = 0,
  gamma_b = 18,
  traces = c("sigma_v", "u", "sigma'_v"),
  labelmode = "stressonly",
  width = 600,
  height = 300
){
  #calculate stresses
  ds <- vertical_stress_profile(
    z_soil,
    z_watertable = z_watertable,
    z_max = z_max,
    z_interval = z_interval,
    q = q,
    gamma_b = gamma_b
  )
  #plot for soil profile
  plt1 <- plotly_soilprofile(
    z_soil,
    z_watertable = z_watertable,
    q = q,
    z_max = z_max,
    gamma_b = gamma_b
  )
  #plot stress profile
  plt2 <- plotly_stressprofile_vertical(
    ds,
    traces = traces
  )
  #combine plots into single plotly plot
  plt <- plotly::subplot(
    plt1, plt2,
    shareY = TRUE,
    titleX = TRUE,
    widths = c(0.3, 0.7),
    margin = 0.05
  )
  #set plot sizes
  plt$width <- width
  plt$height <- height
  #return
  return(plt)
}
