#' Get position of centroid of polygon
#'
#' @description
#' Calculate the centroid position of a polygon based on the coordinates
#' of the polygon
#'
#' @param x,y arrays with x and y coordinates of the polygon
#' @return list with x and y positions
#' @export

polygon_centroid <- function(x, y){
  #second point for each vertex
  x2 <- c(utils::tail(x,-1), x[1])
  y2 <- c(utils::tail(y,-1), y[1])
  #area and centroid of triangle vertex + origin
  A <- 0.5*(x2*y - x*y2)
  xc <- (x + x2) / 3
  yc <- (y + y2) / 3
  #centroid of polygon
  xpc <- sum(A*xc)/sum(A)
  ypc <- sum(A*yc)/sum(A)
  #return
  return(data.frame(x = xpc, y = ypc))
}


#' Get cross-sectional area of polygon
#'
#' @description
#' Calculate the cross-sectional area of a polygon based on the coordinates
#' of the polygon
#'
#' @param x,y arrays with x and y coordinates of the polygon
#' @return list with x and y positions
#' @export

polygon_area <- function(x, y){
  #second point for each vertex
  x2 <- c(utils::tail(x,-1), x[1])
  y2 <- c(utils::tail(y,-1), y[1])
  #area
  abs(sum(0.5*(x2*y - x*y2)))
}


#' Get area of circle on positive side of a cartesian axis system
#'
#' @description
#' Calculate the area of a circle that lies on the positive side of a
#' two-dimensional cartesian coordinate system
#'
#' @param x,y x, y position of the circle midpoint (scalar)
#' @param d diameter of the circle
#' @return area (scalar)
#' @examples
#' #complete on positive sides
#' circle_area_positive(1, 1, 2)
#'
#' #crossing x-axis only
#' circle_area_positive(0, 1, 2)
#' circle_area_positive(-0.5, 1, 2)
#'
#' #crossing y-axis only
#' circle_area_positive(1, 0, 2)
#'
#' #crossing both axes
#' circle_area_positive(0.5, 0.5, 2)
#' circle_area_positive(-0.5, -0.5, 2)
#'
#' #fully on negative side
#' circle_area_positive(-1, -1, 2)
#' @export

circle_area_positive <- function(x, y, d) {
  #get area
  if ((y <= -0.5*d) | (x <= -0.5*d) | ((x <= 0) & (y <= 0) & ((x^2 + y^2) >= (0.5*d)^2))) {
    #fully negative
    return(0)
  } else if ((x >= 0.5*d) & (y >= 0.5*d)) {
    #full area
    return(pi/4*d^2)
  } else if (x >= 0.5*d) {
    #crosses with x-axis
    apex <- acos(2*abs(y)/d) + pi/2*(y < 0)
    return((pi - apex)/4*d^2 - sin(apex)*cos(apex)*(0.5*d)^2)
  } else if (y >= 0.5*d) {
    #crosses with y-axis
    apex <- acos(2*x/d)
    return((pi - apex)/4*d^2 - sin(apex)*cos(apex)*(0.5*d)^2)
  } else {
    #crosses both axes
    apex_x <- acos(2*y/d)
    apex_y <- acos(2*x/d)
    return(
      0.5*(1.5*pi - apex_x - apex_y)/4*d^2 +
        cos(apex_x)*cos(apex_y)*(0.5*d)^2 +
        0.5*sin(apex_x)*cos(apex_x)*(0.5*d)^2 +
        0.5*sin(apex_y)*cos(apex_y)*(0.5*d)^2
    )
  }
}


#' Distribute points evenly on log-space interval
#'
#' @description
#' Equivalent of R's `seq()` function, but then on a log scale.
#' Similar to logscape function in Matlab
#'
#' @param from start of interval
#' @param to end of interval
#' @param ... additional arguments to `seq()` function
#' @return array with values
#' @export

lseq <- function(from, to, ...){
  exp(seq(log(from), log(to), ...))
}


#' Generate coordinates of circles cropped by a polygon
#'
#' @description
#' Generates the coordinates of a number of circles based on positions of the
#' midpoint and their radii. The circles can be cropped by a polygon if
#' polygon points are provided.
#'
#' @param xm,ym x and y coordinates of midpoints of circles (arrays)
#' @param r radii of circles (array with same length as `xm` and `ym`)
#' @param xp,yp x and y coordinates of polygon, defined in clockwise order
#' @param nc number of points to use to draw a full circle
#' @importFrom magrittr `%>%`
#' @importFrom rlang .data
#' @return a tibble with fields `x` and `y` for the coordinates, and a field
#'   `id` to indicate which circle the point belongs to
#' @export
#'
#' @examples
#' #polygon coordinates
#' xp <- c(0, 1, 2, 1)
#' yp <- c(0, 1, 1, 0)
#'
#' #circle coordinates and radii
#' xm <- c(0.5, 0.8)
#' ym <- c(0.5, 0.1)
#' r <- c(0.4, 0.6)
#'
#' #get coordinates, and crop with polygon
#' dc <- circles_polygon(xm, ym, r, xp = xp, yp = yp)
#'
#' #plot
#' ggplot2::ggplot() +
#'   ggplot2::geom_polygon(
#'   ggplot2::aes(x = xp, y = yp),
#'   fill = "grey80"
#' ) +
#' ggplot2::geom_polygon(
#'   data = dc,
#'   ggplot2::aes(x = x, y = y, fill = as.factor(id), color = as.factor(id)),
#'   alpha = 0.2
#' ) +
#' ggplot2::coord_fixed(ratio = 1) +
#' ggplot2::theme(legend.position = "none")

circles_polygon <- function(
  xm = 0.1,
  ym = 0.1,
  r = 1.0,
  xp = NULL,
  yp = NULL,
  nc = 91
){
  #tibble with all circles
  dm <- tibble::tibble(
    xm = xm,
    ym = ym,
    r = r,
    id = seq(length(xm)),
  )
  #generate all points
  dc <- tidyr::expand_grid(
    dm,
    t = seq(-pi, pi, l = nc)
  ) %>%
    dplyr::mutate(
      x = .data$xm + .data$r*cos(.data$t),
      y = .data$ym + .data$r*sin(.data$t)
    )
  #if polygon is inputted
  if (!is.null(xp) & !is.null(yp)) {
    #tibble with polygon points
    dp <- tibble::tibble(
      xp = xp,
      yp = yp,
      tp = atan2(diff(c(yp, yp[1])), diff(c(xp, xp[1])))
    )
    #polygon points within circle - for each circle
    dpm <- tidyr::expand_grid(dm, dp) %>%
      dplyr::mutate(
        rmp = sqrt((.data$xp - .data$xm)^2 + (.data$yp - .data$ym)^2),
        tmp = atan2((.data$yp - .data$ym), (.data$xp - .data$xm))
      ) %>%
      dplyr::filter(.data$rmp < .data$r) %>%
      dplyr::select(.data$id, .data$xp, .data$yp, .data$tmp) %>%
      dplyr::rename(x = .data$xp, y = .data$yp, t = .data$tmp)
    #crossing points for each circle
    dcr <- tidyr::expand_grid(dm, dp) %>%
      dplyr::mutate(
        xmd = (.data$xm - .data$xp)*cos(.data$tp) + (.data$ym - .data$yp)*sin(.data$tp),
        ymd = -(.data$xm - .data$xp)*sin(.data$tp) + (.data$ym - .data$yp)*cos(.data$tp),
        rd = purrr::map2_dbl(
          .data$ymd, .data$r,
          ~ifelse(
            abs(.x) <= .y,
            sqrt(.y^2 - .x^2),
            NA
          )
        )
      ) %>%
      dplyr::filter(!is.na(.data$rd)) %>%
      dplyr::mutate(
        xd0 = .data$xmd - .data$rd,
        xd1 = .data$xmd + .data$rd
      ) %>%
      tidyr::pivot_longer(
        cols = c("xd0", "xd1"),
        names_to = "side",
        values_to = "xd"
      ) %>%
      dplyr::mutate(
        x = .data$xp + .data$xd*cos(.data$tp),
        y = .data$yp + .data$xd*sin(.data$tp),
        t = atan2((.data$y - .data$ym), (.data$x - .data$xm))
      ) %>%
      dplyr::select(.data$id, .data$x, .data$y, .data$t)
    #only keep points that are with polygon
    return(
      tidyr::expand_grid(
        dplyr::bind_rows(dplyr::select(dc, .data$id, .data$x, .data$y, .data$t), dcr),
        dp
      ) %>%
        dplyr::mutate(inside = ((-(.data$x - .data$xp)*sin(.data$tp) + (.data$y - .data$yp)*cos(.data$tp)) <= 0)) %>%
        dplyr::group_by(.data$id, .data$x, .data$y, .data$t) %>%
        dplyr::summarize(inside = all(.data$inside)) %>%
        dplyr::filter(.data$inside == TRUE) %>%
        dplyr::select(.data$id, .data$x, .data$y, .data$t) %>%
        dplyr::bind_rows(dpm) %>%
        dplyr::arrange(.data$id, .data$t)
    )
  } else {
    return(dc)
  }
}
