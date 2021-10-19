#' Custom ggplot theme
#'
#' @description
#' ggplot theme for all ggplots
#'
#' @return a ggplot theme object
#' @export

theme_soilmech <- function() {
  ggplot2::theme_bw(
    base_size = 12
    #base_family = "Arial"
  ) +
    ggplot2::theme(
      text = ggplot2::element_text(
        color = "black"
      ),
      axis.text.y = ggplot2::element_text(
        color = "black"
      ),
      axis.text.x = ggplot2::element_text(
        color = "black"
      )
    )
}


#' Get all minor log tick positions in log10 axis
#'
#' @description
#' Function returns a list of minor axis ticks for a certain range of values
#'
#' @param x array with minimum and maximum axis limits
#' @return an array with minor tick positions on a log-scale
#' @export

get_log10_minorbreaks <- function(x){
  #minimum and maximum order of magnitude
  minx <- floor(min(log10(x), na.rm = TRUE)) - 1
  maxx <- ceiling(max(log10(x), na.rm = TRUE)) + 1
  #number of orders of magnitude
  n_major <- maxx - minx + 1
  #get major breaks
  major_breaks <- seq(minx, maxx, by = 1)
  #get minor breaks
  minor_breaks <- 10^(
    rep(log10(seq(1, 9, by = 1)), times = n_major) +
    rep(major_breaks, each = 9)
  )
  #return only those in range
  minor_breaks[(minor_breaks >= min(x, na.rm = TRUE)) & (minor_breaks <= max(x, na.rm = TRUE))]
}


#' Round values in an array to a 'nice' nearby numbers in plots
#'
#' @description
#' Function to round an array of values to a 'nice' nearby numbers rounded
#' numbers. This can be useful for setting plot limits. Function returns
#' a bounding range
#'
#' @param vals values to be bounded in array
#' @param lower user-defined lower limit, which overrides limit found
#' @param upper user-defined upper limit, which overrides limit found
#' @param log if `log == FALSE`, a linear scale is assumed. If `log == TRUE`,
#'   a log10 scale is assumed
#' @param vals_round an array with potential rounding values, after values
#'   have been scaled by their order of magnitude. Optional values have to be
#'   on the domain `1 <= vals_round < 10`
#' @return two-value array with lower and upper limit
#' @examples
#' x <- c(51, 321)
#' round_limits(x)
#' round_limits(x, lower = 0)
#' round_limits(x, log = TRUE)
#' @export

round_limits <- function(
  vals,
  lower = NA,
  upper = NA,
  log = FALSE,
  vals_round = c(1, 2, 5, 10)
){
  #if single point, create very small range
  if (length(vals) == 1) {
    vals <- c(0.99, 1.01)*vals
  }
  #generate axis limits
  if (log == FALSE){
    #linear scale - range
    range <- c(min(vals, na.rm = TRUE), max(vals, na.rm = TRUE))
    range_width <- range[2] - range[1]
    #order of magnitude
    oom <- floor(log10(range_width))
    #round to nearest...
    round_to <- utils::head(vals_round[vals_round >= (range_width/(10^oom))], 1) * 10^(oom - 1)
    #round
    lims <- c(
      round_to*floor(range[1]/round_to),
      round_to*ceiling(range[2]/round_to)
    )
  } else {
    #log scale - orders of magnitude
    oom1 <- floor(log10(min(vals)))
    oom2 <- ceiling(log10(max(vals)))
    #values
    lims <- c(10^oom1, 10^oom2)
  }
  #override with user-specified limits
  if (!is.na(lower) & !is.null(lower)) {
    lims[1] <- lower
  }
  if (!is.na(upper) & !is.null(upper)) {
    lims[2] <- upper
  }
  #return
  return(lims)
}
