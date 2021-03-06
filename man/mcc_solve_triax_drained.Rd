% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mcc_calculations.R
\name{mcc_solve_triax_drained}
\alias{mcc_solve_triax_drained}
\title{Predict full drained triaxial compression test - Modified Cam-Clay}
\usage{
mcc_solve_triax_drained(
  p0 = 10,
  pc = 100,
  M = 1.35,
  Gamma = 2,
  lambda = 0.115,
  kappa = 0.015,
  nu = 0.3,
  e1 = seq(0, 0.2, l = 101)
)
}
\arguments{
\item{p0}{initial isotropic effective stress}

\item{pc}{initial preconsolidation pressure}

\item{M}{M-parameter}

\item{Gamma}{CSL specific volume intercept}

\item{kappa, lambda}{MCC compression parameters}

\item{nu}{Poisson's ratio}

\item{e1}{values of axial strain at which to calculate response}
}
\value{
a tibble with stress and strain properties at each level of
  axial strain `e1`
}
\description{
Predict the full MCC behaviour for a drained triaxial compression
test
}
\examples{
mcc_solve_triax_drained()
}
