% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plasticity.R
\name{classify_plasticity}
\alias{classify_plasticity}
\title{Classify soil based on plasticity chart}
\usage{
classify_plasticity(wL, Ip)
}
\arguments{
\item{wL}{liquid limits, in \%}

\item{Ip}{plasticity index, in \%}
}
\value{
array with shorthand classifications
}
\description{
Get shorthand for soil major type and plasticity level. Output is given
in shorthand, i.e. MH for a high plasticity silt etc. Classification is
based on the plasticity chart
}
\examples{
classify_plasticity(wL = 60, Ip = 20)
}
