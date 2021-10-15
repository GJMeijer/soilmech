#' Get sparse matrix elements for 1st or 2nd order differentiation
#'
#' @description
#' Function generates rows, column and values for sparse matrix to get a
#' first order approximation of values for all real nodes in a regularly
#' spaced, rectangular grid in either the x or y direction.
#'
#' The function assumes a regular, square grid where nodes are numbered
#' row-first. The grid contains ghost nodes, and has a size 1*1.
#'
#' The grid contains 'ghost' nodes on the side of all edges. For a 4*3 grid,
#' the numbering is ('g' indicates a ghost node):
#'
#' |    | ix=0 | ix=1 | ix=2 | ix=3 | ix=4 | ix=5 |
#' | :---: | :---: | :---: | :---: | :---: | :---: | :---: |
#' | iy=4 |    | 23 (g) | 24 (g) | 25 (g) | 26 (g) |    |
#' | iy=3 | 17 (g) | 18 | 19 | 20 | 21 | 22 (g) |
#' | iy=2 | 11 (g) | 12 | 13 | 14 | 15 | 16 (g) |
#' | iy=1 | 5 (g) | 6  | 7  | 8  | 9  | 10 (g) |
#' | iy=0 |    | 1 (g) | 2 (g) | 3 (g)  | 4 (g)  |    |
#'
#' The finite difference approximation is of the second order (central
#' differences). When only real nodes are used and the node lies on an edge,
#' a second order approximation is used (forwards or backwards, depending
#' on the position)
#' @md
#'
#' @param nx,ny number of real nodes in x and y-directions. There need to be
#'   at least 3 (real or ghost) nodes in the direction of differentiation,
#'   or 4 in case 2nd order differentiation is requested with
#'   `real_only == TRUE`
#' @param direction direction of differentiation
#'   - `direction = "x"`: 1st order in x-direction (d/dx)
#'   - `direction = "y"`: 1st order in y-direction (d/dy)
#'   - `direction = "xx"`: 2nd order in x-direction (d^2/dx^2)
#'   - `direction = "xy"`: 2nd order in y-direction (d^2/dy^2)
#'   - `direction = "xy"`: 1st order in x and y-direction (d^2/dxdy)
#' @param i0 optional index offset for node offset
#' @param multiplier optional multipliers for all values in the finite
#'   difference matrix to differentiate real nodes
#' @param ... potential extra arguments
#' @return a tibble with node indices for row (`row`), column (`col`) and
#'   value (`val`) columns of non-zero entries in the matrix
#' @examples
#' nx <- 4
#' ny <- 3
#' findiff_sparse_elements(nx, ny, "x")
#' findiff_sparse_elements(nx, ny, "yy")
#' findiff_sparse_elements(nx, ny, "xy")
#' @export

findiff_sparse_elements <- function(
  nx,
  ny,
  direction,
  i0 = 0,
  multiplier = 1,
  ...
){
  #enforce at least three nodes in each direction
  nx <- max(3, nx)
  ny <- max(3, ny)
  #step sizes
  hx <- 1/(nx - 1)
  hy <- 1/(ny - 1)
  #indices of real nodes
  i <- index_real(nx, ny, i0 = i0)
  n <- length(i)
  #generate sparse elements of derivative matrix for all real nodes
  if (direction == "x") {
    ## d/dx - central differences
    return(list(
      row = rep(i, 2),
      col = c(
        index_offset(i, nx, ny, dix = -1, i0 = i0),
        index_offset(i, nx, ny, dix = 1, i0 = i0)
      ),
      val = rep(c(-0.5, 0.5)/hx, each = n)*rep(multiplier, 2)
    ))
  } else if (direction == "y") {
    ## d/dy - central differences
    return(list(
      row = rep(i, 2),
      col = c(
        index_offset(i, nx, ny, diy = -1, i0 = i0),
        index_offset(i, nx, ny, diy = 1, i0 = i0)
      ),
      val = rep(c(-0.5, 0.5)/hy, each = n)*rep(multiplier, 2)
    ))
  } else if (direction == "xx") {
    ## d^2/dx^2 - central differences
    return(list(
      row = rep(i, 3),
      col = c(
        index_offset(i, nx, ny, dix = -1, i0 = i0),
        i,
        index_offset(i, nx, ny, dix = 1, i0 = i0)
      ),
      val = rep(c(1, -2, 1)/hx^2, each = n)*rep(multiplier, 3)
    ))
  } else if (direction == "yy") {
    ## d^2/dy^2 - central differences
    return(list(
      row = rep(i, 3),
      col = c(
        index_offset(i, nx, ny, diy = -1, i0 = i0),
        i,
        index_offset(i, nx, ny, diy = 1, i0 = i0)
      ),
      val = rep(c(1, -2, 1)/hy^2, each = n)*rep(multiplier, 3)
    ))
  } else if ((direction == "xy") | (direction == "yx")) {
    ## d^2/dxdy - central differences (apart from cornerpoints)
    #indices of corner points
    ic00 <- index_grid2vector(1, 1, nx, ny, i0 = i0)
    ic10 <- index_grid2vector(nx, 1, nx, ny, i0 = i0)
    ic01 <- index_grid2vector(1, ny, nx, ny, i0 = i0)
    ic11 <- index_grid2vector(nx, ny, nx, ny, i0 = i0)
    #indices of real points
    ic00_real <- 1
    ic10_real <- nx
    ic01_real <- 1 + nx*(ny - 1)
    ic11_real <- nx*ny
    #entries for non-corner points
    cfilter <- (i %in% c(ic00, ic10, ic01, ic11))
    row_nc <- rep(i[!cfilter], 4)
    col_nc <- c(
      index_offset(i[!cfilter], nx, ny, dix = -1, diy = -1, i0 = i0),
      index_offset(i[!cfilter], nx, ny, dix = 1, diy = -1, i0 = i0),
      index_offset(i[!cfilter], nx, ny, dix = -1, diy = 1, i0 = i0),
      index_offset(i[!cfilter], nx, ny, dix = 1, diy = 1, i0 = i0)
    )
    val_nc = rep(c(1, -1, -1, 1)/(4*hx*hy), each = sum(!cfilter))*rep(multiplier[!cfilter], 4)
    #entries for cornerpoints
    row_c00 <- rep(ic00, 5)
    col_c00 <- index_offset(ic00, nx, ny, dix = c(-1, -1, -1, 1, 1), diy = c(0, 1, 2, -1, 1), i0 = i0)
    val_c00 <- multiplier[ic00_real]*c(0.75, -1, 0.25, -0.25, 0.25)/(hx*hy)
    row_c10 <- rep(ic10, 5)
    col_c10 <- index_offset(ic10, nx, ny, dix = c(1, 1, 1, -1, -1), diy = c(0, 1, 2, -1, 1), i0 = i0)
    val_c10 <- multiplier[ic10_real]*c(-0.75, 1, -0.25, 0.25, -0.25)/(hx*hy)
    row_c01 <- rep(ic01, 5)
    col_c01 <- index_offset(ic01, nx, ny, dix = c(-1, -1, -1, 1, 1), diy = c(-2, -1, 0, -1, 1), i0 = i0)
    val_c01 <- multiplier[ic01_real]*c(-0.75, 1, -0.25, -0.25, 0.25)/(hx*hy)
    row_c11 <- rep(ic11, 5)
    col_c11 <- index_offset(ic11, nx, ny, dix = c(1, 1, 1, -1, -1), diy = c(-2, -1, 0, -1, 1), i0 = i0)
    val_c11 <- multiplier[ic11_real]*c(0.75, -1, 0.25, 0.25, -0.25)/(hx*hy)
    #return all entries
    return(list(
      row = c(row_c00, row_c10, row_nc, row_c01, row_c11),
      col = c(col_c00, col_c10, col_nc, col_c01, col_c11),
      val = c(val_c00, val_c10, val_nc, val_c01, val_c11)
    ))
  }
}


#' Get sparse matrix elements for 1st order differentiation (no ghost nodes)
#'
#' @description
#' Function generates rows, column and values for sparse matrix to get a
#' first order approximation of values for all real nodes in a regularly
#' spaced, rectangular grid in either the x or y direction.
#'
#' The function assumes a regular, square grid where nodes are numbered
#' row-first. The grid does not contain any ghost nodes. The grid has a
#' size 1*1.
#'
#' For a 4*3 grid, the numbering is:
#'
#' |    | ix=1 | ix=2 | ix=3 | ix=4 |
#' | :---: | :---: | :---: | :---: | :---: |
#' | iy=3 | 9  | 10 | 11 | 12 |
#' | iy=2 | 5  | 6  | 7  | 8  |
#' | iy=1 | 1  | 2  | 3  | 4  |
#'
#' The finite difference approximation is of the second order (central
#' differences). When only real nodes are used and the node lies on an edge,
#' a second order approximation is used (forwards or backwards, depending
#' on the position)
#' @md
#'
#' @param nx,ny number of real nodes in x and y-directions. There need to be
#'   at least 3 (real or ghost) nodes in the direction of differentiation,
#'   or 4 in case 2nd order differentiation is requested with
#'   `real_only == TRUE`
#' @param direction direction of differentiation
#'   - `direction = "x"`: 1st order in x-direction (d/dx)
#'   - `direction = "y"`: 1st order in y-direction (d/dy)
#' @param i0 optional index offset for node offset
#' @param multiplier optional multipliers for all values in the finite
#'   difference matrix to differentiate real nodes
#' @param ... potential extra arguments
#' @return a tibble with node indices for row (`row`), column (`col`) and
#'   value (`val`) columns of non-zero entries in the matrix
#' @examples
#' findiff_sparse_elements_realonly(9, 5, "y")
#' @export

findiff_sparse_elements_realonly <- function(
  nx,
  ny,
  direction,
  i0 = 0,
  multiplier = 1,
  ...
){
  #enforce at least three nodes in each direction
  nx <- max(3, nx)
  ny <- max(3, ny)
  #step sizes
  hx <- 1/(nx - 1)
  hy <- 1/(ny - 1)
  #first derivative in x-direction
  if (direction == "x") {
    row_single <- c(
      rep(1, 3),
      rep(seq(2, nx - 1), each = 2),
      rep(nx, 3)
    )
    col_single <- c(
      seq(1, 3),
      rep(seq(2, nx - 1), each = 2) + rep(c(-1, 1), nx - 2),
      seq(nx - 2, nx)
    )
    val_single <- c(
      c(-1.5, 2, -0.5),
      rep(c(-0.5, 0.5), (nx - 2)),
      c(0.5, -2, 1.5)
    )/hx
    return(tibble::tibble(
      row = i0 + rep(row_single, ny) + rep(nx*(seq(ny) - 1), each = length(row_single)),
      col = i0 + rep(col_single, ny) + rep(nx*(seq(ny) - 1), each = length(row_single)),
      val = multiplier*rep(val_single, ny)
    ))
  }
  #first derivative in y-direction
  if (direction == "y") {
    row_single <- c(
      rep(1, 3),
      rep(1 + seq(1, ny - 2)*nx, each = 2),
      1 + rep(nx*(ny - 1), 3)
    )
    col_single <- c(
      1 + seq(0, 2)*nx,
      1 + (rep(seq(1, ny - 2), each = 2) + rep(c(-1, 1), ny - 2))*nx,
      1 + seq(ny - 3, ny - 1)*nx
    )
    val_single <- c(
      c(-1.5, 2, -0.5),
      rep(c(-0.5, 0.5), (ny - 2)),
      c(0.5, -2, 1.5)
    )/hy
    return(tibble::tibble(
      row = i0 + rep(row_single, nx) + rep(seq(0, nx - 1), each = length(row_single)),
      col = i0 + rep(col_single, nx) + rep(seq(0, nx - 1), each = length(row_single)),
      val = multiplier*rep(val_single, nx)
    ))
  }
}


#' Determine total number of nodes
#'
#' @description
#' Determine the total number of nodes on a grid of nodes. Ghost nodes
#' are or are not included based on `real_only` setting.
#'
#' @inheritParams index_grid2vector
#' @return number of nodes (scalar)
#' @export

nodes_total <- function(nx, ny, real_only = FALSE, ...) {
  if (real_only == TRUE) {
    return(nx*ny)
  } else {
    return(nx*ny + 2*(nx + ny))
  }
}


#' Convert node grid positions to indices
#'
#' @description
#' Function takes a vectors of x and y node indices and returns the index
#' (numbered row-wise)
#'
#' @inheritParams index_vector2grid
#' @param ix,iy arrays with index of position in x and y directions
#' @return array with node index for each node
#' @examples
#' #define grid
#' nx <- 4
#' ny <- 3
#'
#' #including ghost nodes
#' i <- seq(nodes_total(nx, ny))
#' ixy <- index_vector2grid(i, nx, ny)
#' index_grid2vector(ixy$ix, ixy$iy, nx, ny)
#'
#' #excluding ghost nodes
#' i <- seq(nodes_total(nx, ny, real_only = TRUE))
#' ixy <- index_vector2grid(i, nx, ny, real_only = TRUE)
#' index_grid2vector(ixy$ix, ixy$iy, nx, ny, real_only = TRUE)
#' @export

index_grid2vector <- function(ix, iy, nx, ny, i0 = 0, real_only = FALSE, ...) {
  if (real_only == TRUE) {
    return(i0 + ix + nx*(iy - 1))
  } else {
    return(
      i0 + (iy - 1)*(nx + 2) + 1 + nx + ix +
        as.integer(iy == 0) -
        as.integer(iy == (ny + 1))
    )
  }
}


#' Convert node position to grid x-y position
#'
#' @description
#' Function takes a vector of node indices (numbered row-wise) and
#' returns a list of their x and y indices.
#'
#' @param i array with node indices
#' @param nx,ny number of real nodes in x and y directions
#' @param i0 node index offset
#' @param real_only if `TRUE`, only real nodes are numbered and ghost nodes
#'   are ignored in numbering
#' @param ... additional arguments to pass
#' @return a list with two fields: `ix` and `iy` with indices in x and y
#'   directions
#' @examples
#' #convert to grid positions
#' nx <- 4
#' ny <- 3
#' i <- seq(1, (nx + 2)*(ny + 2) - 4)
#' index_vector2grid(i, nx, ny)
#' @export

index_vector2grid <- function(i, nx, ny, i0 = 0, real_only = FALSE, ...) {
  if (real_only == TRUE) {
    return(list(
      ix = 1 + (i - i0 - 1)%%nx,
      iy = 1 + floor((i - i0 - 1)/nx)
    ))
  } else {
    return(list(
      ix = (i - i0 - nx - 1)%%(nx + 2) -
        as.integer((i - i0) <= nx) +
        as.integer((i - i0) >= ((nx + 2)*(ny + 1) - 1)),
      iy = floor((i - i0 + 1)/(nx + 2))
    ))
  }
}


#' Get indices of all real nodes
#'
#' @description
#' Get an array of the indices of all real nodes in a rectangular grid. For
#' node numbering, see function `findiff_sparse_elements()`
#'
#' @inheritParams index_edge
#' @return array of indices for all real nodes (non-ghost nodes)
#' @export

index_real <- function(nx, ny, i0 = 0, real_only = FALSE, ...) {
  if (real_only == TRUE) {
    return(i0 + seq(nx*ny))
  } else {
    return(i0 + rep(seq(nx), ny) + rep(seq(ny), each = nx)*(nx + 2) - 1)
  }
}


#' Find indices of neighbouring nodes
#'
#' @description
#' Find indices of neighbouring nodes in a rectangular grid with or without
#' ghost nodes. Nodes are numbered row-wise. Function assumes that these
#' nodes exist, e.g. no check is conducted to see if this is so...
#'
#' @inheritParams index_vector2grid
#' @param dix,diy index offsets in x and y-directions
#' @return a list of indices of neightbouring nodes
#' @examples
#' nx <- 4
#' ny <- 3
#' i <- c(1, 11, 20)
#' index_offset(i, nx, ny, dix = 1)
#' index_offset(i, nx, ny, diy = 1)
#' @export

index_offset <- function(
  i,
  nx, ny,
  dix = 0, diy = 0,
  i0 = 0, real_only = FALSE, ...
){
  ixy <- index_vector2grid(i, nx, ny, i0 = i0, real_only = real_only)
  return(index_grid2vector(
    ixy$ix + dix, ixy$iy + diy,
    nx, ny,
    i0 = i0, real_only = real_only))
}


#' Get indices of nodes near edge of grid
#'
#' @description
#' Get an array of the indices of nodes near an edge of a rectangular grid.
#' For node numbering, see function `findiff_sparse_elements()`.
#'
#' @param nx,ny number of real nodes in x and y-direction
#' @param edge edge number, numbered clockwise. `1` for left (negative x-face),
#'   `2` for top (positive y-face), `3` for right (positive x-face) and `4`
#'   for bottom (negative y-face)
#' @param real_only if `TRUE`, only real nodes are numbered. if `FALSE`, all
#'   ghost nodes are included in the numbering as well
#' @param i0 node starting number (index of first node minus 1)
#' @param offset row (in case `edge = 2` or `edge = 4`, or column (in case
#'   `edge = 1` or `edge = 3`) offset from edge
#' @param ... additional optional arguments
#' @return array of indices of ghost nodes on selected edge
#' @examples
#' #real nodes + ghost nodes
#' index_edge(4, 3, 1)
#' index_edge(4, 3, 2, offset = 1)
#'
#' #real nodes only
#' index_edge(4, 3, 1, real_only = TRUE)
#' index_edge(4, 3, 2, offset = 1, real_only = TRUE)
#' @export

index_edge <- function(nx, ny, edge, real_only = FALSE, i0 = 0, offset = 0, ...) {
  if (real_only == TRUE) {
    if (edge == 1) {
      return(i0 + 1 + (seq(ny) - 1)*nx + offset)
    } else if (edge == 2) {
      return(i0 + (ny - 1 - offset)*nx + seq(nx))
    } else if (edge == 3) {
      return(i0 + seq(ny)*nx - offset)
    } else if (edge == 4) {
      return(i0 + seq(nx) + nx*offset)
    }
  } else {
    if (edge == 1) {
      return(i0 + seq(ny)*(nx + 2) - 1 + offset)
    } else if (edge == 2) {
      return(i0 + seq(nx) + (nx + 2)*(ny + 1) - 2 - offset*(nx + 2) + ifelse(offset > 0, 1, 0))
    } else if (edge == 3) {
      return(i0 + (seq(ny) + 1)*(nx + 2) - 2 - offset)
    } else if (edge == 4) {
      return(i0 + seq(nx) + offset*(nx + 2) - ifelse(offset > 0, 1, 0))
    }
  }
}
