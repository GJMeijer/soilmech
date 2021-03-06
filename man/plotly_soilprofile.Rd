% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/stress_profiles_vertical.R
\name{plotly_soilprofile}
\alias{plotly_soilprofile}
\title{plotly soil profiles}
\usage{
plotly_soilprofile(
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
  line_width = 2,
  nround = 1,
  fields_hover = c("description", "gamma_b", "phi_deg", "thickness"),
  fields_label = c("gamma_b", "phi_deg"),
  arrow_size = 0.1,
  arrow_n = 6,
  title = NULL,
  ...
)
}
\arguments{
\item{z_soil}{depth of top of each layer}

\item{z_watertable}{depth of the water table}

\item{z_max}{maximum depth of interest}

\item{gamma_b}{unit weight of soil, in kN/m3. Should be a scalar or a
vector with the same length as `z_soil`}

\item{gamma_w}{unit weight of water: 10 kN/m3}

\item{q}{surcharge at soil surface}

\item{description}{soil description: will be shown in hoverlabel}

\item{phi_deg}{angle of internal friction: will be shown in text and
hoverlabels}

\item{ylim}{depth axis limits}

\item{ylab}{depth axis label}

\item{fill_soildry}{fill colour for dry soil}

\item{fill_soilwet}{fill colour for wet soil}

\item{fill_water}{fill colour for ponding water}

\item{colour_soil}{line color for soil polygons}

\item{colour_water}{line color for water table}

\item{colour_surcharge}{color of surcharge arrows and label}

\item{line_water}{line type for water table}

\item{line_width}{line tickness}

\item{nround}{number of decimals in plotly hoverlabels}

\item{fields_hover}{field names to show in hoverlabels}

\item{fields_label}{field names to plot in permanent text labels}

\item{arrow_size}{relative length of surcharge arrows}

\item{arrow_n}{number of surcharge arrows to use}

\item{title}{plot title}

\item{...}{additional fields used for plotting labels}
}
\description{
Creates a plotly object with the soil profile and soil properties
}
\examples{
#ponding water
plotly_soilprofile(
  z_soil = 0,
  z_max = 5,
  z_watertable = -5,
  gamma_b = 20
)

#water table within soil
plotly_soilprofile(
  z_soil = c(0, 3, 7),
  gamma_b = c(15, 16, 20),
  z_watertable = 2,
  description = c("sand", "clay", "gravel"),
  phi_deg = c(0, 10, 30),
  q = 10,
  title = "Test"
)
}
