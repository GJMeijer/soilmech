% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/flownet_quadrilateral.R,
%   R/flownet_quadrilateral_old.R
\name{flownet_solve_quadrilateral}
\alias{flownet_solve_quadrilateral}
\title{Solve a flow net problem with quadrilateral domains}
\usage{
flownet_solve_quadrilateral(df)

flownet_solve_quadrilateral(df)
}
\arguments{
\item{df}{list with geometry and boundary conditions of the problem. See
function `flownet_geometry_quadrilateral()` for more information}
}
\value{
a tibble with domains (`domain`) and positions of all nodes
  (`a`, `b`, `x`, `y`) and solutions for head (`h`), flow rates (`qx`,
  `qy`) and flow potential (`psi`)

a tibble with domains (`domain`) and positions of all nodes
  (`a`, `b`, `x`, `y`) and solutions for head (`h`), flow rates (`qx`,
  `qy`) and flow potential (`psi`)
}
\description{
Take a quadrilateral flow net problem and calculate results for the
hydraulic head, flow rates and flow potential, using finite differences

Take a quadrilateral flow net problem and calculate results for the
hydraulic head, flow rates and flow potential, using finite differences
}
\examples{
df <- flownet_geometry_quadrilateral(grid_size = 0.1)
flownet_solve_quadrilateral(df)
df <- flownet_geometry_quadrilateral(grid_size = 0.1)
flownet_solve_quadrilateral(df)
}
