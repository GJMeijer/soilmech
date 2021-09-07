#' Calculate vertical soil stresses with depth and generate plotly labels
#'
#' @description
#' Function generates a dataframe with soil stresses with depth that can be
#' used for subsequent plotting (plotly) or displaying in a bookdown document.
#'
#' @param z array with depths (top of layers)
#' @param q surcharge at the soil surface. If array with length `z` is
#'   inputted, these surcharges are applied at the corresponding depth
#'   of `z`
#' @param saturated if `TRUE`, layers are fully saturated, and if `FALSE`,
#'   layers are assumed fully dry. If an array with length `z` is inputted,
#'   saturations can be set independently for each layer
#' @param gammab Bulk unit weight of the layers. If an array with length
#'   `z` is inputted, unit weghts can be set independently for each layer
#' @param zmax maximum depth to analyse
#' @param zmarg small depth inverval, used for suddenly jumping stress
#' @param zinterval if specified, add points at regular interval. parameter
#'   specifies depth interval distance
#' @param labelmode type of plotly hover labelling. can be `"stressonly`",
#'   showing only the stresses at that depth, or `"calculations"`, showing
#'   how the stress is calculated using values at depth level above current
#'   depth
#' @param gammaw unit weight of water
#' @param output type of ouput. If `output = "plot"`, the dataframe required
#'   for a plot is returned. If `output = "bookdown"`, a list is returned
#'   containing a dataframe and a vector with row headers

calculate_stressprofile <- function(
  z,
  q = 0,
  saturated = FALSE,
  gammab = 20,
  zmax = 10,
  zmarg = 1e-6,
  zinterval = NULL,
  labelmode = "stressonly",
  gammaw = 10,
  output = "plot"
){
  #increase zmax if required
  zmax <- max(zmax, z)
  #generate tibble
  d <- tibble::tibble(
    z = z,
    q = c(q, rep(0, length(z) - length(q))),
    saturated = saturated,
    gammab = gammab
  )
  #add points on regular interval if required
  if (!is.null(zinterval)) {
    #create points with default values
    dadd <- tibble::tibble(
      z = seq(seq(min(d$z), max(d$z), zinterval)),
      q = 0,
      saturated = TRUE,
      gammab = 20
    )
    #remove points already defined in <d>
    dadd <- dadd[!dadd$z %in% d$z, ]
    #change saturation and unit weight when required
    for (i in 1:dim(dadd)[1]){
      ind <- min(which(d$z >= dadd$z[i])[1] - 1, dim(d)[1], na.rm = TRUE)
      dadd$saturated[i] <- d$saturated[ind]
      dadd$gammab[i] <- d$gammab[ind]
    }
    #merge and sort
    d <- rbind(d, dadd)
    d <- d[order(d$z), ]
  }
  #get dataframe with stresses - excluding points just before surcharge
  ds1 <- tibble::tibble(
    layer = 1 + seq(0, nrow(d)),
    z = c(d$z, zmax),
    sigma = cumsum(c(0, diff(c(d$z, zmax))*d$gammab)) + cumsum(c(d$q, 0)),
    u = cumsum(c(0, d$saturated*diff(c(d$z, zmax))*gammaw))
  )
  #increase with respect to last point
  ds1$dgammab <- c(0, d$gammab)
  ds1$dq <- c(d$q, 0)
  ds1$dz <- c(0, diff(c(d$z, zmax)))
  ds1$sigma_last <- c(0, utils::head(ds1$sigma, -1))
  ds1$u_last <- c(0, utils::head(ds1$u, -1))
  #get dataframe with stresses - just before surcharge
  ds0 <- ds1[c(d$q != 0, FALSE), ]
  ds0$z <- ds0$z - zmarg
  ds0$sigma <- ds0$sigma - d$q[d$q!=0]
  #increase with espect to last point
  if (dim(ds0)[1] > 0){
    ds0$dq <- 0
  }
  #merge together and sort
  ds <- rbind(ds1, ds0[ds0$z >= min(ds1$z), ])
  ds <- ds[order(ds$z), ]
  #effective stress
  ds$sigmad <- with(ds, sigma - u)
  #hover labels
  ds$sigma_label <- NA
  ds$u_label <- NA
  ds$sigmad_label <- NA
  if (labelmode == "stressonly"){
    ds$sigma_label <- with(ds, paste("\u03c3<sub>v</sub> = ", ds$sigma, " kPa"))
    ds$u_label <- with(ds, paste("u = ", ds$sigma, " kPa"))
    ds$sigmad_label <- with(ds, paste("\u03c3'<sub>v</sub> = ", ds$sigmad, " kPa"))
  } else if (labelmode == "calculations"){
    #labels - first point
    if (ds$dq[1]==0) {
      ds$sigma_label[1] <- "\u03c3<sub>v</sub> = 0 kPa"
    } else {
      ds$sigma_label[1] <- paste0("\u03c3<sub>v</sub> = q = ", ds$dq[1], " kPa")
    }
    ds$u_label[1] <- "u = 0 kPa"
    ds$sigmad_label[1] <- paste0("\u03c3'<sub>v</sub> = \u03c3<sub>v</sub> - u = ", ds$sigma[1], " - ", ds$u[1], " = ", ds$sigmad[1], " kPa")
    # other labels - total
    ind0 <- which((ds$dq == 0) & (ds$layer > 1))
    ds$sigma_label[ind0] <- with(ds[ind0, ], paste0("\u03c3<sub>v</sub> = ", sigma_last , " + ", dz, "\u00b7", dgammab, " = ", sigma, " kPa"))
    ind1 <- which((ds$dq !=0) & (ds$layer > 1))
    ds$sigma_label[ind1] <- with(ds[ind1, ], paste0("\u03c3<sub>v</sub> = ", ds$sigma[ind1 - 1], " + ", dq, " = ", sigma, " kPa"))
    # labels - pore water pressure
    ds$u_label[ds$u == 0] <- "u = 0 kPa"
    ds$u_label[ds$u > 0] <- with(ds[ds$u>0, ], paste0("u = ", u_last, " + ", dz, "\u00b710 = ", u, " kPa"))
    # labels - effective stress
    ds$sigmad_label <- with(ds, paste0("\u03c3'<sub>v</sub> = \u03c3<sub>v</sub> - u = ", sigmad, " kPa"))
  }
  #return
  if (output == "bookdown"){
    return(list(
      ds[, c("z", "sigma", "u", "sigmad")],
      c("Depth [m]", "$\\sigma_v$ [kPa]", "$u$ [kPa]", "$\\sigma'_v$ [kPa]")
    ))
  } else {
    return(ds)
  }
}


#' Plotly vertical stress profiles with depth
#'
#' @description
#' Function generates a plotly object showing how vertical stresses in the
#' soil change with depth
#'
#' @param ds dataframe with all depths, stresses etc, outputted from the
#'   function `calculate_stressprofile()`
#' @param traces array with traces to be plotted. May contain `sigma` for
#'   the total vertical stress, `u` for the pore water pressure and `sigmad`
#'   for the vertical effective stress
#' @param linewidth line thickness for traces
#' @param cols array with colors for total stress, pore water pressure
#'   and effective stress traces
#' @param linetype array with linetypes for total stress, pore water pressure
#'   and effective stress traces
#' @param xlim stress axis limits (manual override)
#' @param ylim depth axis limits (manual override)
#' @param hovermode hovermode to use in plotly
#' @importFrom magrittr `%>%`
#' @return plotly object
#' @export

plotly_stressprofile <- function(
  ds,
  traces = c("sigma", "u", "sigmad"),
  linewidth = 1.5,
  cols = c("#000000", "#0000ff", "#ff0000"),
  linetype = c("solid", "dash", "dot"),
  xlim = NA,
  ylim = NA,
  hovermode = "y"
){
  #axis limits
  if (is.na(xlim[1])){
    xlim = round_limits(max(ds$sigma), lower = 0)
  } else {
    xlim <- sort(xlim, decreasing = TRUE)
  }
  if (is.na(ylim[1])){
    ylim = c(max(ds$z), min(ds$z))
  } else {
    ylim <- sort(ylim, decreasing = TRUE)
  }
  #initiate plot
  p <- plotly::plot_ly(ds)
  #add traces
  if ("sigma" %in% traces){
    p <- p %>% plotly::add_trace(
      x = ~sigma,
      y = ~z,
      name = "Total stress",
      type = "scatter",
      mode = "lines+markers",
      hoverinfo = "text",
      text = ~sigma_label,
      marker = list(
        color=cols[1]
      ),
      line = list(
        width = linewidth,
        color = cols[1],
        dash = linetype[1]
      )
    )
  }
  if ("u" %in% traces){
    p <- p %>% plotly::add_trace(
      x = ~u,
      y = ~z,
      name = "Pore water pressure",
      type = "scatter",
      mode = "lines+markers",
      hoverinfo = "text",
      text = ~u_label,
      marker = list(
        color = cols[2]
      ),
      line = list(
        width = linewidth,
        color = cols[2],
        dash = linetype[2]
      )
    )
  }
  if ("sigmad" %in% traces){
    p <- p %>% plotly::add_trace(
      x = ~sigmad,
      y = ~z,
      name = "Effective stress",
      type = "scatter",
      mode = "lines+markers",
      hoverinfo = "text",
      text = ~sigmad_label,
      marker = list(
        color = cols[3]
      ),
      line = list(
        width = linewidth,
        color = cols[3],
        dash = linetype[3]
      )
    )
  }
  #add layout
  p <- p %>% plotly::layout(
    xaxis = list(
      range = xlim,
      title = "Stress [kPa]",
      side = "top",
      visible = TRUE
    ),
    yaxis = list(
      range = ylim,
      title = "Depth [m]"
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
  #changes axis titles if only one plot
  if (length(traces) == 1){
    if (traces == "sigma"){
      p <- plotly::layout(p, xaxis = list(title = "Total stress [kPa]"))
    } else if (traces == "u"){
      p <- plotly::layout(p, xaxis = list(title = "Pore water pressure [kPa]"))
    } else if (traces == "sigmad"){
      p <- plotly::layout(p, xaxis = list(title = "Effective stress [kPa]"))
    }
  }
  #return
  return(p)
}


#' Plotly soil profile with depth
#'
#' @description
#' Function generates a plotly object showing the soil layering and properties
#' with depth
#'
#' @param z array with depths (top of layers)
#' @param q surcharge at the soil surface. If array with length `z` is
#'   inputted, these surcharges are applied at the corresponding depth
#'   of `z`
#' @param saturated if `TRUE`, layers are fully saturated, and if `FALSE`,
#'   layers are assumed fully dry. If an array with length `z` is inputted,
#'   saturations can be set independently for each layer
#' @param gammab Bulk unit weight of the layers. If an array with length
#'   `z` is inputted, unit weghts can be set independently for each layer
#' @param zmax maximum depth to analyse
#' @param gammaw unit weight for water
#' @param col_soildry face color for dry soil
#' @param col_soilsat face color for saturated soil
#' @param col_soilline line color for soil layer boundaries
#' @param col_water face color for ponding water
#' @param col_waterline line color for water table
#' @param linetype_soil line type for soil layer boundaries
#' @param linetype_water line type for water table
#' @param linewidth line thickness for soil layer boundaries and water table
#' @param xlim x-axis limits
#' @param ylim y-axis limits
#' @param arrow_size length of surcharge arrows
#' @param arrow_n number of surcharge arrows over width
#' @importFrom magrittr `%>%`
#' @return plotly object
#' @export

plotly_soilprofile <- function(
  z,
  q = 0,
  saturated = FALSE,
  gammab = 20,
  zmax = 10,
  gammaw = 10,
  col_soildry = "#d3bc5f",
  col_soilsat = "#aebab7",
  col_soilline = "#65571d",
  col_water = "#2a7fff",
  col_waterline = "#2a7fff",
  linetype_soil = "solid",
  linetype_water = "dash",
  linewidth = 2.0,
  xlim = NA,
  ylim = NA,
  arrow_size = 0.075,
  arrow_n = 6
){
  #increase zmax if required
  zmax <- max(zmax, z)
  #generate tibble with layering
  d <- tibble::tibble(
    z = z,
    q = c(q, rep(0, length(z) - length(q))),
    saturated = saturated,
    gammab = gammab
  )
  #axis limits
  if (is.na(ylim[1])){
    ylim <- c(max(d$z, zmax), min(d$z))
  }
  if (is.na(xlim[1])){
    xlim <- c(0, 1)
  }
  #initiate plot
  p <- plotly::plot_ly()
  #add layers + annotation
  for (i in 1:nrow(d)) {
    #get coordinates
    if (i == nrow(d)) {
      z0 <- d$z[i]
      z1 <- zmax
    } else {
      z0 <- d$z[i]
      z1 <- d$z[i + 1]
    }
    if (z1 > z0){
      #hover labels for soil masses
      label1 <- paste0("\u03B3<sub>b</sub> = ", d$gammab[i], " kN/m\u00b3")
      label2 <- paste0("h = ", c(d$z[-1], zmax) - d$z, " m")
      label3 <- rep("Fully saturated soil", dim(d)[1])
      label3[d$saturated == FALSE] <- "Dry soil"
      label3[(d$saturated == TRUE) & (d$gammab == gammaw)] <- "Water"
      label <- paste(label3, label1, label2, sep = '<br>')
      #get color
      if (d$saturated[i] == FALSE){
        soil_color <- col_soildry
      } else if ((d$saturated[i] == TRUE) & (d$gammab[i] == gammaw)) {
        soil_color <- col_water
      } else {
        soil_color <- col_soilsat
      }
      #plot rectangle
      p <- p %>% plotly::add_trace(
        type = "scatter",
        mode = "lines",
        x = c(xlim, rev(xlim)),
        y = c(z0, z0, z1, z1),
        fill = "toself",
        fillcolor = soil_color,
        line = list("width" = 0.0),
        hoveron = "fills",
        text = label[i],
        hoverinfo = "text",
        showlegend = FALSE
      )
      #add text
      p <- p %>% plotly::add_trace(
        type = "scatter",
        x = 0.5*(xlim[2] - xlim[1]),
        y = 0.5*(z0 + z1),
        mode = "text",
        text = paste0("\u03B3<sub>b</sub> = ", d$gammab[i], " kN/m\u00b3"),
        textposition = "middle center",
        showlegend = FALSE,
        hoverinfo = "skip"
      )
      #add interface
      if (c(d$saturated[1], diff(d$saturated))[i] == 1){
        p <- p %>% plotly::add_trace(
          type = "scatter",
          mode = "lines",
          x = xlim,
          y = c(z0, z0),
          line = list(
            width = linewidth,
            color = col_waterline,
            dash = linetype_water
          ),
          hoverinfo = "skip",
          showlegend = FALSE
        )
      } else {
        p <- p %>% plotly::add_trace(
          type = "scatter",
          mode = "lines",
          x = xlim,
          y = c(z0, z0),
          line = list(
            width = linewidth,
            color = col_soilline,
            dash = linetype_soil
          ),
          hoverinfo = "skip",
          showlegend = FALSE
        )
      }
    }
    #add overburden arrows
    if (d$q[i] != 0){
      #arrow data
      da <- data.frame(
        x0 = seq(xlim[1], xlim[2] - 1/arrow_n, l = arrow_n) + 0.5/arrow_n,
        x1 = seq(xlim[1], xlim[2] - 1/arrow_n, l = arrow_n) + 0.5/arrow_n,
        y0 = z0,
        y1 = z0 - (max(d$z) - min(d$z))*arrow_size
      )
      #add to plot
      p <- p %>% plotly::add_annotations(
        x = da$x0,
        y = da$y0,
        xref = "x", yref = "y",
        axref = "x", ayref = "y",
        text = "",
        showarrow = TRUE,
        ax = da$x1,
        ay = da$y1
      )
      #add overburden text
      p <- p %>% plotly::add_trace(
        type= "scatter",
        mode = "text",
        x = xlim[1],
        y = da$y1[1],
        text = paste0("q = ", d$q[i], " kPa"),
        textposition = "top right",
        hoverinfo = "skip",
        showlegend = FALSE
      )
      #change axis limits
      ylim[2] <- min(ylim[2], z0 - 2*(max(d$z) - min(d$z))*arrow_size)
    }
  }
  #add layout
  p <- p %>% plotly::layout(
    xaxis = list(
      range = xlim,
      title = "position",
      side = "top",
      visible = FALSE
    ),
    yaxis = list(
      range = ylim,
      title = "Depth [m]"
    ),
    showlegend = FALSE,
    template = "none"
  )
  #return
  return(p)
}


#' Plotly soil and stress profiles side-by-side
#'
#' @description
#' Generate a plotly object showing the soil profile (on the left) and the
#' stress profiles (on the right) for a particular soil layering
#'
#' @inheritParams plotly_soilprofile
#' @inheritParams plotly_stressprofile
#' @param zinterval equal interval for extra points
#' @param labelmode governs what to plot in hoverlabels, e.g `"stressonly"``
#' @param width plot width, in pixels?
#' @param height plot height, in pixels?
#' @return a plotly object
#' @export

#plot together
plotly_soilstressprofile <- function(
  z,
  q = 0,
  saturated = FALSE,
  gammab = 20,
  zmax = 10,
  traces = c("sigma", "u", "sigmad"),
  zinterval = NULL,
  labelmode = "stressonly",
  width = 600,
  height = 300
){
  #calculate stresses
  ds <- calculate_stressprofile(
    z, q = q, saturated = saturated, gammab = gammab,
    zmax = zmax, zinterval = zinterval, labelmode = labelmode
  )
  #generate individual plots
  p1 <- plotly_soilprofile(z, q = q, saturated = saturated, gammab = gammab, zmax = zmax)
  p2 <- plotly_stressprofile(ds, traces = traces)
  #combine plots into single plotly plot
  p <- plotly::subplot(
    p1, p2,
    shareY = TRUE, titleX = TRUE,
    widths = c(0.3, 0.7),
    margin = 0.05
  )
  #set plot sizes
  p$width <- width
  p$height <- height
  #return
  return(p)
}


#' ggplot vertical stress profiles with depth
#'
#' @description
#' Function generates a ggplot object showing how vertical stresses in the
#' soil change with depth
#'
#' @param ds dataframe with all depths, stresses etc, outputted from the
#'   function `calculate_stressprofile()`
#' @param traces array with traces to be plotted. May contain `sigma` for
#'   the total vertical stress, `u` for the pore water pressure and `sigmad`
#'   for the vertical effective stress
#' @param linewidth line thickness for traces
#' @param cols array with colors for total stress, pore water pressure
#'   and effective stress traces
#' @param linetype array with linetypes for total stress, pore water pressure
#'   and effective stress traces
#' @param xlim stress axis limits (manual override)
#' @param ylim depth axis limits (manual override)
#' @param layers if `TRUE`, plot layer interfaces as horizontal lines
#' @return ggplot object
#' @export

#ggplot profile
ggplot_stressprofile <- function(
  ds,
  traces = c("sigma", "u", "sigmad"),
  linewidth = 0.5,
  cols = c("#000000", "#0000ff", "#ff0000"),
  linetype = c(1, 2, 3),
  xlim = NA,
  ylim = NA,
  layers = TRUE
){
  #initiate plot
  p <- ggplot2::ggplot()
  #convert to long data
  dsl <- tidyr::pivot_longer(
    ds,
    cols = dplyr::all_of(.data$traces),
    names_to = "stresstype",
    values_to = "value"
  )
  #add data to plot
  p <- p + ggplot2::geom_path(
    data = dsl,
    ggplot2::aes(y = .data$z, x = .data$value, color = .data$stresstype, linetype = .data$stresstype)
  )
  if (layers == TRUE){
    p <- p + ggplot2::geom_hline(
      data = ds,
      ggplot2::aes(yintercept = .data$z),
      color = "grey50", linetype = 4
    )
  }
  #select axis labels
  if (length(traces) == 1){
    if (traces == "sigma"){
      xlab <- "Total stress [kPa]"
    } else if (traces == "u") {
      xlab <- "Pore water pressure [kPa]"
    } else if (traces == "sigmad") {
      xlab <- "Effective stress [kPa]"
    } else {
      xlab <- "Stress [kPa]"
    }
    legend <- FALSE
  } else {
    xlab <- "Stress [kPa]"
    legend <- TRUE
  }
  #axis extend
  if (is.na(xlim)) {
    xlim <- round_limits(max(dsl$value), lower = 0)
  }
  if (is.na(ylim)) {
    ylim <- c(max(ds$z), min(ds$z))
  }
  #change axes
  p <- p +
    ggplot2::theme_bw() +
    ggplot2::xlab(xlab) +
    ggplot2::ylab("Depth [m]") +
    ggplot2::scale_x_continuous(lim = xlim, expand = c(0, 0), position = "top") +
    ggplot2::scale_y_reverse(lim = ylim, expand = c(0, 0))
  #selected traces
  sel <- sort(match(traces, c("sigma", "u", "sigmad")))
  #change colors
  p <- p +
    ggplot2::scale_linetype_manual(name = "", values = linetype[sel], labels = expression(sigma[v], u, sigma*"'"[v])[sel]) +
    ggplot2::scale_color_manual(name = "", values = cols[sel], labels = expression(sigma[v], u, sigma*"'"[v])[sel])
  #legend
  if (legend == TRUE){
    p <- p + ggplot2::theme(
      legend.position = c(0.95, 0.95),
      legend.justification = c(1, 1),
      legend.title = ggplot2::element_blank()
    )
  } else {
    p <- p + ggplot2::theme(legend.position = "none")
  }
  #return
  return(p)
}
