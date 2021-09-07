#' Create a random packing of particles
#'
#' @description
#' Creates a random plot of particles
#'
#' @param porosity target porosity (fraction)
#' @param dmax maximum particle diameter
#' @param dmin minimum particle diameter
#' @param width soil element width
#' @param height soil element height
#' @param clip if `TRUE`, particles partially sticking out of sample are
#'   cropped to the soil element
#' @param nc number of points on each particle circle
#' @param npos max number of grid positions to choose from
#' @param npar max number of particles
#' @param fill_soil color of fill of soil cube
#' @param color_soil color of outline of soil cube
#' @param size_soil thickness of particle border
#' @param size_border thickness of soil element border
#' @importFrom magrittr `%>%`
#' @importFrom rlang .data
#' @return a ggplot object
#' @examples
#' #loose sample
#' random_packing(porosity = 0.35)
#'
#' #dense sample
#' random_packing(porosity = 0.15)
#' @export

random_packing <- function(
  porosity = 0.3,
  dmax = 0.40,
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
){
  #draw many random positions
  dpos <- tidyr::expand_grid(
    x = seq(0 - 0.499*dmax, width + 0.499*dmax, l = round(sqrt(npos))),
    y = seq(0 - 0.499*dmax, height + 0.499*dmax, l = round(sqrt(npos))),
  ) %>%
    dplyr::mutate(
      valid = TRUE,
      id = seq(round(sqrt(npos))^2)
    )
  #initiate particles
  dpar <- tibble::tibble(
    id = seq(npar),
    x = NA,
    y = NA,
    d = NA,
    A = NA,
  )
  #loop
  for (i in 1:npar){
    #draw a random position
    ind <- sample(dpos$id[dpos$valid == TRUE], 1)
    #assign position
    dpar$x[i] <- dpos$x[ind]
    dpar$y[i] <- dpos$y[ind]
    #get diameter
    if (i == 1) {
      dpar$d[i] <- dmax
    } else {
      #distance between new centrepoint and all existing ones
      distance <- sqrt((dpar$x[1:(i - 1)] - dpar$x[i])^2 + (dpar$y[1:(i - 1)] - dpar$y[i])^2)
      #smallest clearances
      clearance <- min(distance - 0.5*dpar$d[1:(i - 1)])
      #diameter
      dpar$d[i] <- min(dmax, 2*clearance)
    }
    #assign area
    dpar$A[i] <- circle_area_positive(
      min(dpar$x[i], width - dpar$x[i]),
      min(dpar$y[i], height - dpar$y[i]),
      dpar$d[i]
    )
    #total area
    A_total <- sum(dpar$A[1:i])
    #update list of positions
    dpos$valid[((dpos$x - dpar$x[i])^2 + (dpos$y - dpar$y[i])^2) < (0.5*dpar$d[i] + 0.5*dmin)^2] <- FALSE
    #check if porosity reached
    if ((A_total >= (height*width*(1 - porosity))) | all(dpos$valid == FALSE)) {
      break
    }
  }
  #generate circles
  if (clip == TRUE) {
    dc <- circles_polygon(
      dpar$x,
      dpar$y,
      0.5*dpar$d,
      xp = c(0, 0, width, width),
      yp = c(0, height, height, 0)
    )
  } else {
    dc <- circles_polygon(
      dpar$x,
      dpar$y,
      0.5*dpar$d
    )
  }
  #plot
  ggplot2::ggplot() +
    ggplot2::theme_void() +
    ggplot2::geom_polygon(
      data = dc,
      ggplot2::aes(x = .data$x, y = .data$y, group = .data$id),
      color = color_soil,
      fill = fill_soil,
      size = size_soil
    ) +
    ggplot2::geom_polygon(
      ggplot2::aes(
        x = c(0, 0, width, width),
        y = c(0, height, height, 0)
      ),
      color = "black",
      fill = NA,
      size = size_border
    ) +
    ggplot2::coord_fixed(ratio = 1)
}

