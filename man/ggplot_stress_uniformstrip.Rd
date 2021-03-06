% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/elastic_solutions.R
\name{ggplot_stress_uniformstrip}
\alias{ggplot_stress_uniformstrip}
\title{ggplot influence factors for stress underneath uniform strip foundation}
\usage{
ggplot_stress_uniformstrip(
  direction = "z",
  nu = 0.3,
  xlim = c(-2.5, 2.5),
  ylim = c(0, 7),
  gridsize = 0.05,
  contour_width = 0.1,
  contour_label = c(0, 0.1, 0.2, 0.4, 0.6, 0.8, 1),
  showgrid = TRUE,
  larrow = 0.3,
  lsize = 3.5,
  labelsize = 3,
  palette = "Blues",
  xB = NULL,
  zB = NULL,
  nround = 3
)
}
\arguments{
\item{direction}{the direction of the stress of interest. `z` for
vertical, `x` for horizontal away from strip, `y` for horizontal
in direction of strip, and `xz` for shear stress}

\item{nu}{Poisson's ratio, required for stress in y-direction}

\item{xlim}{2-value array with limits for x-axis}

\item{ylim}{2-value array with limits for depth axis}

\item{gridsize}{grid size for calculations}

\item{contour_width}{I-interval between contour lines}

\item{contour_label}{values of I for which to draw contour lines}

\item{showgrid}{if `TRUE`, heatmap of influence factors is shown in
background of the plot}

\item{larrow}{length of load arrows (in normalised coordinates)}

\item{lsize}{thickness of load arrows}

\item{labelsize}{size of labels on contour lines}

\item{palette}{RColorBrewer palette for heat map}

\item{xB}{crosshairs annotation}

\item{zB}{crosshairs annotation}

\item{nround}{number of decimals in crosshairs I label}
}
\value{
ggplot object
}
\description{
Creates a ggplot chart with contour lines and heat map (optional) for the
influence factors for stress underneath a uniform strip load.
Assumes elastic solutions are valid.
}
\examples{
#plot empty chart
ggplot_stress_uniformstrip()

#plot with annotation
ggplot_stress_uniformstrip(xB = 1, zB = 2)
}
