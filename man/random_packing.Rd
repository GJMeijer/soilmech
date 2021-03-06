% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/grain_packing.R
\name{random_packing}
\alias{random_packing}
\title{Create a random packing of particles}
\usage{
random_packing(
  porosity = 0.3,
  dmax = 0.4,
  dmin = 0.02,
  width = 1,
  height = 1,
  clip = TRUE,
  nc = 91,
  npos = 100^2,
  npar = 1000,
  fill_soil = "#d3bc5f",
  color_soil = "#65571d",
  size_soil = 0.5,
  size_border = 1
)
}
\arguments{
\item{porosity}{target porosity (fraction)}

\item{dmax}{maximum particle diameter}

\item{dmin}{minimum particle diameter}

\item{width}{soil element width}

\item{height}{soil element height}

\item{clip}{if `TRUE`, particles partially sticking out of sample are
cropped to the soil element}

\item{nc}{number of points on each particle circle}

\item{npos}{max number of grid positions to choose from}

\item{npar}{max number of particles}

\item{fill_soil}{color of fill of soil cube}

\item{color_soil}{color of outline of soil cube}

\item{size_soil}{thickness of particle border}

\item{size_border}{thickness of soil element border}
}
\value{
a ggplot object
}
\description{
Creates a random plot of particles
}
\examples{
#loose sample
random_packing(porosity = 0.35)

#dense sample
random_packing(porosity = 0.15)
}
