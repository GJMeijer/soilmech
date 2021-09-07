#' Plot simple shear element with contractive soil
#'
#' @description
#' Function returns a plot with a simple shear-deformed element filled
#' with soil particles, starting off in a loose state and densifying
#' with deformation
#'
#' @param u relative shear displacement, between 0 (no deformation, loose)
#'   and 1 (max deformation, dense)
#' @param width width of the sample relative to its height
#' @param n number of grains to use over the sample height
#' @param nc number of vertices for plotting circles for grains
#' @param crop if `TRUE`, only parts of particles that are present within
#'   the element of soil are plotted. If `FALSE`, also the parts of the
#'   particles that are sticking out are plotted
#' @param marg plot margins, relative to the height of the soil sample
#' @param color_water color for soil particles
#' @param color_soil color for water
#' @param size_water line thickness for water/soil element
#' @param size_soil line thickness for soil particles
#' @importFrom rlang .data
#' @return a ggplot object
#' @examples
#' ggplot_simpleshear_contraction(u = 0, n = 3)
#' ggplot_simpleshear_contraction(u = 1, n = 3)
#' @export

ggplot_simpleshear_contraction <- function(
  u = 1,
  width = 1,
  n = 4,
  nc = 91,
  crop = TRUE,
  marg = 0.02,
  color_water = "#2a7fff",
  color_soil = "#d3bc5f",
  size_water = 1,
  size_soil = 0.5
){
  #radius of grain
  r <- 0.5/n
  #max displacement angle - in radians
  t <- min(1, max(u, 0))*30/180*pi
  #number of grains in x and y direction
  nx <- ceiling(width*n)
  ny <- ceiling(n)
  #positions of undeformated grains
  d <- tidyr::expand_grid(
    x0 = 0.5*width + (seq(0, nx)- 0.5*nx) * 2*r,
    y0 = 0.5 + (seq(0, ny)- 0.5*ny) * 2*r
  )
  d$group <- seq(nrow(d))
  #deformed position
  d$x1 <- d$x0 + d$y0*sin(t)
  d$y1 <- d$y0*cos(t)
  #generate polygons for grains
  dc <- tidyr::expand_grid(d, t = seq(-pi, pi, l = nc))
  dc$x <- dc$x1 + r*cos(dc$t)
  dc$y <- dc$y1 + r*sin(dc$t)
  #crop
  if (crop == TRUE){
    #dc$y[dc$y < 0] <- 0
    #dc$y[dc$y > cos(t)] <- cos(t)
    #dc$x[dc$x < dc$y*tan(t)] <- dc$y[dc$x < dc$y*tan(t)]*tan(t)
    #dc$x[dc$x > (width + dc$y*tan(t))] <- width + dc$y[dc$x > (width + dc$y*tan(t))]*tan(t)
  }
  #create polygons for outline soil
  do <- tibble::tibble(
    x = c(0, sin(t), width + sin(t), width),
    y = c(0, cos(t), cos(t), 0)
  )
  #create polygons for outline water
  do2 <- tibble::tibble(
    x = c(0, tan(t), width + tan(t), width),
    y = c(0, 1, 1, 0)
  )
  #plot limits
  if (crop == TRUE){
    xlim <- c(0, width + tan(30*pi/180)) + c(-1, 1)*marg
    ylim <- c(0, 1) + c(-1, 1)*marg
  } else {
    xlim <- c(min(dc$x0), max(dc$x0)) + c(-1, 1)*r + c(0, tan(30*pi/180)) + c(-1, 1)*marg
    ylim <- c(min(dc$y0), max(dc$y0)) + c(-1, 1)*r + c(-1, 1)*marg
  }
  #plot
  ggplot2::ggplot() +
    ggplot2::theme_void() +
    ggplot2::geom_polygon(
      data = do2,
      ggplot2::aes(x = .data$x, y = .data$y),
      color = "black",
      fill = color_water,
      linetype = 2,
      size = size_water
    ) +
    ggplot2::geom_polygon(
      data = dc,
      ggplot2::aes(x = .data$x, y = .data$y, group = .data$group),
      color = "black",
      fill = color_soil,
      size = size_soil
    ) +
    ggplot2::geom_polygon(
      data = do,
      ggplot2::aes(x = .data$x, y = .data$y),
      fill = NA,
      colour = "black",
      size = size_water
    ) +
    ggplot2::coord_fixed(xlim = xlim, ylim = ylim, expand = FALSE)
}


#' Plot simple shear element with dilating soil
#'
#' @description
#' Function returns a plot with a simple shear-deformed element filled
#' with soil particles, starting off in a dense state and dilating
#' with deformation
#'
#' @param u relative shear displacement, between 0 (no deformation, dense)
#'   and 1 (max deformation, loose)
#' @param width width of the sample relative to its height
#' @param n number of grains to use over the sample height
#' @param nc number of vertices for plotting circles for grains
#' @param crop if `TRUE`, only parts of particles that are present within
#'   the element of soil are plotted. If `FALSE`, also the parts of the
#'   particles that are sticking out are plotted
#' @param marg plot margins, relative to the height of the soil sample
#' @param color_water color for soil particles
#' @param color_soil color for water
#' @param size_water line thickness for water/soil element
#' @param size_soil line thickness for soil particles
#' @importFrom rlang .data
#' @return a ggplot object
#' @examples
#' ggplot_simpleshear_dilation(u = 0, n = 3)
#' ggplot_simpleshear_dilation(u = 1, n = 3)
#' @export

ggplot_simpleshear_dilation <- function(
  u = 1,
  width = 1,
  n = 4,
  nc = 91,
  crop = TRUE,
  marg = 0.02,
  color_water = "#2a7fff",
  color_soil = "#d3bc5f",
  size_water = 1,
  size_soil = 0.5
){
  #initial angle
  t0 <- 30/180*pi
  #max displacement angle - in radians
  t <- min(1, max(u, 0))*t0
  #radius of grain
  r <- 0.5/n/cos(t0)
  #number of grains in x and y direction
  ny <- 1 + 2*ceiling((n-1)/2)
  nx <- 1 + 2*width*(ceiling((n-1)/2*cos(t0)))
  #positions of undeformated grains - always force a grain to be present in middle
  d <- tidyr::expand_grid(
    row = seq(-(1 + floor(ny/2)), (1 + floor(ny/2))),
    col = seq(-(1 + floor(nx/2)), floor(nx/2))
  )
  d$x0 <- 2*r*d$col + (d$row%%2)*sin(t0)*2*r
  d$y0 <- 2*r*d$row * cos(t0)
  d$group <- seq(nrow(d))
  #remove grains that are too far away
  d <- d[(d$x0 >= (-0.5*(width+2*r))) & (d$x0 <= (0.5*(width+2*r))), ]
  #deformed position
  d$y1 <- d$y0*cos(t0-t)/cos(t0)
  d$x1 <- d$x0 + d$y0*(sin(t0) - sin(t0-t))/cos(t0)
  #generate polygons for grains
  dc <- tidyr::expand_grid(d, t = seq(-pi, pi, l = nc))
  dc$x <- dc$x1 + r*cos(dc$t)
  dc$y <- dc$y1 + r*sin(dc$t)
  #crop
  if (crop == TRUE){
    ymax <- 0.5*cos(t0-t)/cos(t0)
    xmin <- -0.5*width + dc$y*(sin(t0) - sin(t0-t))/cos(t0-t)
    xmax <- 0.5*width + dc$y*(sin(t0) - sin(t0-t))/cos(t0-t)
    dc$y[dc$y <= (-ymax)] <- -ymax
    dc$y[dc$y >= ymax] <- ymax
    dc$x[dc$x <= xmin] <- xmin[dc$x <= xmin]
    dc$x[dc$x >= xmax] <- xmax[dc$x >= xmax]
  }
  #create polygons for outline soil
  do <- tibble::tibble(
    x = c(-0.5, -0.5, 0.5, 0.5)*width + c(-0.5, 0.5, 0.5, -0.5)*(sin(t0) - sin(t0-t))/cos(t0),
    y = c(-0.5, 0.5, 0.5, -0.5) * cos(t0-t)/cos(t0)
  )
  #create polygon for outline soil (no volume change)
  do2 <- tibble::tibble(
    x = c(-0.5, -0.5, 0.5, 0.5)*width + c(-0.5, 0.5, 0.5, -0.5)*(sin(t0) - sin(t0-t))/cos(t0-t),
    y = c(-0.5, 0.5, 0.5, -0.5)
  )
  #create polygons for outline water
  porosity0 <- 1 - pi/(4*cos(30*pi/180))
  porosity1 <- cos(t0-t)/cos(t0) - (1 - porosity0)
  ymaxwater1 <- porosity0/porosity1
  dow <- tibble::tibble(
    x = c(-0.5, -0.5, 0.5, 0.5)*width + ymaxwater1*c(-0.5, 0.5, 0.5, -0.5)*(sin(t0) - sin(t0-t))/cos(t0),
    y = c(-0.5, 0.5, 0.5, -0.5) * ymaxwater1*cos(t0-t)/cos(t0)
  )
  #plot limits
  if (crop == TRUE){
    xlim <- c(0, 1)*width + c(0, 1)*tan(t0) + c(-1, 1)*marg
    ylim <- c(0, 1)/cos(t0) + c(-1, 1)*marg
  } else {
    xlim <- c(min(d$x0) - r, max(d$x0) + r) + 0.5*width + c(0, tan(t0)) + c(-1, 1)*marg
    ylim <- c(min(d$y0)/cos(t0) - r, max(d$y0)/cos(t0) + r) + 0.5/cos(t0) + c(-1, 1)*marg
  }
  #add offset to deformed positions
  dc$x <- dc$x - min(do$x)
  dc$y <- dc$y - min(do$y)
  do2$x <- do2$x - min(do2$x)
  do2$y <- do2$y - min(do2$y)
  dow$y <- dow$y - min(dow$y)
  dow$x <- dow$x - min(dow$x)
  do$x <- do$x - min(do$x)
  do$y <- do$y - min(do$y)
  #plot
  ggplot2::ggplot() +
    ggplot2::theme_void() +
    ggplot2::geom_polygon(
      data = dow,
      ggplot2::aes(x = .data$x, y = .data$y),
      color = NA,
      fill = color_water
    ) +
    ggplot2::geom_polygon(
      data = dc,
      ggplot2::aes(x = .data$x, y = .data$y, group = .data$group),
      color = "black",
      fill = color_soil,
      size = size_soil
    ) +
    ggplot2::geom_polygon(
      data = do,
      ggplot2::aes(x = .data$x, y = .data$y),
      fill = NA,
      colour = "black",
      size = size_water
    ) +
    ggplot2::geom_polygon(
      data = do2,
      ggplot2::aes(x = .data$x, y = .data$y),
      color = "black",
      fill = NA,
      linetype = 2,
      size = size_water
    ) +
    ggplot2::coord_fixed(xlim = xlim, ylim = ylim, expand = FALSE)
}
