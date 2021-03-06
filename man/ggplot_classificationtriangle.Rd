% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/classification_triangle.R
\name{ggplot_classificationtriangle}
\alias{ggplot_classificationtriangle}
\title{ggplot classification triangle}
\usage{
ggplot_classificationtriangle(
  gravel = 0,
  sand = 0,
  silt = 0,
  clay = 0,
  group = NULL
)
}
\arguments{
\item{gravel}{mass of gravel for crosshair annotation}

\item{sand}{mass of sand for crosshair annotation}

\item{silt}{mass of silt for crosshair annotation}

\item{clay}{mass of clay for crosshair annotation}

\item{group}{grouping for crosshair annotation (label)}
}
\value{
a ggplot object
}
\description{
ggplot both parts of the BS classification triangle together into a single
ggplot object
}
\examples{
#empty chart
ggplot_classificationtriangle()

#example with annotations
ggplot_classificationtriangle(30,20,15,10)
}
