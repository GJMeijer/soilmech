% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/consolidation_time.R
\name{ggplot_consolidation_averagetime}
\alias{ggplot_consolidation_averagetime}
\title{ggplot for average consolidation in layer as function of time}
\usage{
ggplot_consolidation_averagetime(
  Tv = NULL,
  trace = NULL,
  label_Tv = 0.03,
  label_size = 3,
  label_alpha = 1,
  nround = 2
)
}
\arguments{
\item{Tv}{normalised time for crosshair annotation}

\item{trace}{trace selection (1, 2 or 3) for crosshair annotation}

\item{label_Tv}{time at which to plot time labels on traces}

\item{label_size}{size of the text in time labels}

\item{label_alpha}{opacity of the time labels}

\item{nround}{number of decimals to use in rounding values in crosshair
annotations}
}
\value{
a ggplot object
#empty chart
ggplot_consolidation_averagetime()

#annotated chart
ggplot_consolidation_averagetime(Tv = 0.5, trace = 1)
}
\description{
plots a ggplot for the average consolidation in a layer as function of
time.

Three traces are plotted
* `trace = 1` drainage in open layers with thickness 2d (uniform or
    triangular), or in half-closed layer with uniform initial pressure
    and thickness d
* `trace = 2` drainage in half-closed layer with triangular distribution
    of initial excess pore water pressure. Larger pressure near closed
    end
* `trace = 3` drainage in half-closed layer with triangular distribution
    of initial excess pore water pressure. Larger pressure near open
    end
}
