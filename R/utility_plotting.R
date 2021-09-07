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
    oom2 <- floor(log10(max(vals)))
    #values
    lims <- c(
      utils::tail(vals_round[vals_round <= (min(vals)/(10^oom1))], 1) * 10^oom1,
      utils::head(vals_round[vals_round >= (max(vals)/(10^oom2))], 1) * 10^oom2
    )
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


#############################
### DEPRECIATED FUNCTIONS ###
#############################

#' Round values in an array to a 'nice' nearby number
#'
#' @description
#' Function to round an array of values to a 'nice' nearby number rounded
#' numbers. This can be useful for setting plot limits. Values can be
#' rounded upwards or downwards depending in the user choice
#'
#' @param x array with values to be rounded
#' @param xl array with accepted values on a range [1,10].
#' @param upper if `upper == TRUE`, values are rounded upwards.
#'   If `upper = FALSE`, values are rounded downwards
#' @return array with rounded values. Same size as input `x`
#' @examples
#' x <- 100 * runif(5)
#' xf <- depreciated_round_limits(x, upper = TRUE)
#' print(x)
#' print(xf)
#' @export

depreciated_round_limits <- function(
  x,
  xl = c(1, 1.2, 1.5, 2, 2.5, 3, 4, 5, 6, 8, 10),
  upper = TRUE
){
  if (upper == TRUE){
    return(purrr::map_dbl(
      x,
      function(xi) ifelse((xi ==0), 0, utils::head(xl[xl>=(xi/(10^floor(log10(xi))))],1)*(10^floor(log10(xi))))
    ))
  } else {
    return(purrr::map_dbl(
      x,
      function(xi) ifelse((xi == 0), 0, utils::tail(xl[xl<=(xi/(10^floor(log10(xi))))],1)*(10^floor(log10(xi))))
    ))
  }
}


#' Annotate water table marker to a ggplot object
#'
#' @param plt ggplot object
#' @param x,y x and y-position where marker touches water table
#' @param scale scaling factor for marker. Default height is approx 1
#' @param line_size thickness of lines
#' @param fill_water fill colour for triangle
#' @param colour_water line colour for all lines
#' @return ggplot object with added water table marker
#' @examples
#' #generate simple plot
#' plt <- ggplot2::ggplot() +
#'   ggplot2::annotate("rect", xmin = -2, xmax = 2, ymin = -2, ymax = 2, fill = "yellow") +
#'   ggplot2::coord_fixed(ratio = 1)
#'
#' #add water table marker
#' ggplot_add_watermarker(plt, x = 1, y = 1, scale = 1)
#' @export

ggplot_add_watermarker <- function(
  plt,
  x,
  y,
  scale = 1,
  line_size = 0.5,
  fill_water = "#2a7fff",
  colour_water = "#000080"
){
  plt +
    ggplot2::annotate(
      "polygon",
      x = x + scale*c(0, -0.5*tan(pi/6), 0.5*tan(pi/6)),
      y = y + scale*c(0, 0.5, 0.5),
      size = line_size,
      fill = fill_water,
      color = colour_water
    ) +
    ggplot2::annotate(
      "segment",
      x = x - scale*0.25,
      xend = x + scale*0.25,
      y = y - scale*0.15,
      yend = y - scale*0.15,
      size = line_size,
      color = colour_water
    ) +
    ggplot2::annotate(
      "segment",
      x = x - scale*0.15,
      xend = x + scale*0.15,
      y = y - scale*0.30,
      yend = y - scale*0.30,
      size = line_size,
      color = colour_water
    ) +
    ggplot2::annotate(
      "segment",
      x = x - scale*0.05,
      xend = x + scale*0.05,
      y = y - scale*0.45,
      yend = y - scale*0.45,
      size = line_size,
      color = colour_water
    )
}


#' Annotate soil surface marker to a ggplot object
#'
#' @param plt ggplot object
#' @param x,y x and y-position where marker touches the soil surface (middle)
#' @param scale scaling factor for marker. Default width is approx 1
#' @param line_size thickness of lines
#' @param colour_soil line colour for all lines
#' @return ggplot object with added water table marker
#' @examples
#' #generate simple plot
#' plt <- ggplot2::ggplot() +
#'   ggplot2::annotate("rect", xmin = -2, xmax = 2, ymin = -2, ymax = 2, fill = "yellow") +
#'   ggplot2::coord_fixed(ratio = 1)
#'
#' #add water table marker
#' ggplot_add_soilmarker(plt, x = 1, y = 1, scale = 1)
#' @export

ggplot_add_soilmarker <- function(
  plt,
  x,
  y,
  scale = 1,
  line_size = 0.5,
  colour_soil = "#65571d",
  n = 3
){
  df <- tibble::tibble(
    x = scale*c(seq(-0.5, 0.5 - 0.5/n, l = 2*n), seq(-0.25, 0.25, l = (n + 1))),
    y = scale*c(rep(0, 2*n), rep(-0.25*tan(pi/3), 4)),
    xend = scale*c(seq(-0.25, 0 - 0.25/n, l = n), seq(0.25, 0.5 - 0.25/n, l = n), seq(0, 0.25 - 0.25/n, l = n), 0.5),
    yend = -scale*0.25*tan(pi/3)*c(rep(seq(1, 0 + 1/n, l = n), 2), seq(0, 1 - 1/n, l = n), 0)
  )
  plt + ggplot2::annotate(
    "segment",
    x = df$x,
    xend = df$xend,
    y = df$y,
    yend = df$yend,
    size = line_size,
    color = colour_soil
  )
}
