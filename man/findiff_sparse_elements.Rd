% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/flownet_finitedifferences.R
\name{findiff_sparse_elements}
\alias{findiff_sparse_elements}
\title{Get sparse matrix elements for 1st or 2nd order differentiation}
\usage{
findiff_sparse_elements(nx, ny, direction, i0 = 0, multiplier = 1, ...)
}
\arguments{
\item{nx, ny}{number of real nodes in x and y-directions. There need to be
at least 3 (real or ghost) nodes in the direction of differentiation,
or 4 in case 2nd order differentiation is requested with
\code{real_only == TRUE}}

\item{direction}{direction of differentiation
\itemize{
\item \code{direction = "x"}: 1st order in x-direction (d/dx)
\item \code{direction = "y"}: 1st order in y-direction (d/dy)
\item \code{direction = "xx"}: 2nd order in x-direction (d^2/dx^2)
\item \code{direction = "xy"}: 2nd order in y-direction (d^2/dy^2)
\item \code{direction = "xy"}: 1st order in x and y-direction (d^2/dxdy)
}}

\item{i0}{optional index offset for node offset}

\item{multiplier}{optional multipliers for all values in the finite
difference matrix to differentiate real nodes}

\item{...}{potential extra arguments}
}
\value{
a tibble with node indices for row (\code{row}), column (\code{col}) and
value (\code{val}) columns of non-zero entries in the matrix
}
\description{
Function generates rows, column and values for sparse matrix to get a
first order approximation of values for all real nodes in a regularly
spaced, rectangular grid in either the x or y direction.

The function assumes a regular, square grid where nodes are numbered
row-first. The grid contains ghost nodes, and has a size 1*1.

The grid contains 'ghost' nodes on the side of all edges. For a 4*3 grid,
the numbering is ('g' indicates a ghost node):\tabular{ccccccc}{
    \tab ix=0 \tab ix=1 \tab ix=2 \tab ix=3 \tab ix=4 \tab ix=5 \cr
   iy=4 \tab  \tab 23 (g) \tab 24 (g) \tab 25 (g) \tab 26 (g) \tab  \cr
   iy=3 \tab 17 (g) \tab 18 \tab 19 \tab 20 \tab 21 \tab 22 (g) \cr
   iy=2 \tab 11 (g) \tab 12 \tab 13 \tab 14 \tab 15 \tab 16 (g) \cr
   iy=1 \tab 5 (g) \tab 6 \tab 7 \tab 8 \tab 9 \tab 10 (g) \cr
   iy=0 \tab  \tab 1 (g) \tab 2 (g) \tab 3 (g) \tab 4 (g) \tab  \cr
}


The finite difference approximation is of the second order (central
differences). When only real nodes are used and the node lies on an edge,
a second order approximation is used (forwards or backwards, depending
on the position)
}
\examples{
nx <- 4
ny <- 3
findiff_sparse_elements(nx, ny, "x")
findiff_sparse_elements(nx, ny, "yy")
findiff_sparse_elements(nx, ny, "xy")
}
