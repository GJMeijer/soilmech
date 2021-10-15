#' Get sparse matrix elements for 1st or 2nd order differentiation
#'
#' @description
#' Function generates rows, column and values for sparse matrix to get a
#' first order approximation of values for all real nodes in a regularly
#' spaced, rectangular grid in either the x or y direction.
#'
#' The function assumes a regular, rectangular grid where nodes are numbered
#' row-first. The grid may or may not contain ghost nodes, depending on the
#' input settings.
#'
#' If `real_only = TRUE`, the grid contains only 'real' nodes. For a 4*3 grid,
#' these are numbered:
#'
#' |    | ix=1 | ix=2 | ix=3 | ix=4 |
#' | :---: | :---: | :---: | :---: | :---: |
#' | iy=3 | 9  | 10 | 11 | 12 |
#' | iy=2 | 5  | 6  | 7  | 8  |
#' | iy=1 | 1  | 2  | 3  | 4  |
#'
#' If `real_only = FALSE`, the grid also contains 'ghost' nodes on the side
#' of all edges. For a 4*3 grid, the numbering is:
#'
#' |    | ix=0 | ix=1 | ix=2 | ix=3 | ix=4 | ix=5 |
#' | :---: | :---: | :---: | :---: | :---: | :---: | :---: |
#' | iy=4 |    | 23 | 24 | 25 | 26 |    |
#' | iy=3 | 17 | 18 | 19 | 20 | 21 | 22 |
#' | iy=2 | 11 | 12 | 13 | 14 | 15 | 16 |
#' | iy=1 | 5  | 6  | 7  | 8  | 9  | 10 |
#' | iy=0 |    | 1  | 2  | 3  | 4  |    |
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
#' @param order the order of differentiation. `order = 1` for first-order
#'   differentiation, and `order = 2` for second order differentiation
#' @param h step size
#' @param direction direction of finite differences. `direction` should be
#'   equal to either `x` or `y`
#' @param real_only if `TRUE`, the grid is assumed to contain only real nodes,
#'   so there are `nx*ny` nodes. If `FALSE`, the grid is assumed to contain
#'   ghost nodes on the sides of each of the four edges of the grid, resulting
#'   in `nx*ny + 2*nx + 2*ny` nodes. These are numbered in row-wise order,
#'   starting with the left-most ghost node underneath the lower edge
#' @param i0_row,i0_col optional index offset for rows and columns
#' @param multiplier optional multipliers for all values in the finite
#'   difference matrix to differentiate real nodes
#' @param ... potential extra arguments
#' @return a tibble with node indices for row (`row`), column (`col`) and
#'   value (`val`) columns of non-zero entries in the matrix
#' @examples
#' findiff_sparse_entries(5, 4, order = 1, direction = "x", real_only = TRUE)
#' findiff_sparse_entries(5, 4, order = 1, direction = "y", real_only = FALSE)
#' findiff_sparse_entries(5, 4, order = 2, direction = "x", real_only = TRUE)
#' findiff_sparse_entries(5, 4, order = 2, direction = "y", real_only = FALSE)
#' @export

findiff_sparse_entries <- function(
  nx,
  ny,
  order = 1,
  h = 1,
  direction = "x",
  real_only = FALSE,
  i0_row = 0,
  i0_col = 0,
  multiplier = 1,
  ...
){
  if (order == 1) {
    ## FIRST ORDER DIFFERENTATION
    if (direction == "x") {
      #x-direction
      if (real_only == TRUE) {
        #elements for the first row
        row_first <- c(rep(1, 3), rep(seq(2, nx - 1), each = 2), rep(nx, 3))
        col_first <- c(1, 2, 3, rep(seq(2, nx - 1), each = 2) + rep(c(-1, 1), (nx - 2)), nx - c(2, 1, 0))
        val_first <- c(-1.5, 2, -0.5, rep(c(-0.5, 0.5), (nx - 2)), 0.5, -2, 1.5)/h
        #elements for all rows
        return(tibble::tibble(
          row = i0_row + rep(row_first, ny) + rep((seq(ny) - 1)*nx, each = length(row_first)),
          col = i0_col + rep(col_first, ny) + rep((seq(ny) - 1)*nx, each = length(col_first)),
          val = multiplier*rep(val_first, ny)
        ))
      } else {
        #elements for the first row
        row_first <- (nx + 1) + rep(seq(nx), each = 2)
        col_first <- row_first + rep(c(-1, 1), nx)
        val_first <- rep(c(-0.5, 0.5), nx)/h
        #elements for all rows
        return(tibble::tibble(
          row = i0_row + rep(row_first, ny) + rep((seq(ny) - 1)*(nx + 2), each = length(row_first)),
          col = i0_col + rep(col_first, ny) + rep((seq(ny) - 1)*(nx + 2), each = length(col_first)),
          val = multiplier*rep(val_first, ny)
        ))
      }
    } else if (direction == "y") {
      #y-direction
      if (real_only == TRUE) {
        #elements for the first column
        row_first <- 1 + (c(rep(1, 3), rep(seq(2, ny - 1), each = 2), rep(ny, 3)) - 1)*nx
        col_first <- 1 + (c(1, 2, 3, rep(seq(2, ny - 1), each = 2) + rep(c(-1, 1), (ny - 2)), ny - c(2, 1, 0)) - 1)*nx
        val_first <- c(-1.5, 2, -0.5, rep(c(-0.5, 0.5), (ny - 2)), 0.5, -2, 1.5)/h
        #elements for all columns
        return(tibble::tibble(
          row = i0_row + rep(row_first, nx) + rep(seq(0, nx - 1), each = length(row_first)),
          col = i0_col + rep(col_first, nx) + rep(seq(0, nx - 1), each = length(row_first)),
          val = multiplier*rep(val_first, nx)
        ))
      } else {
        #elements for the first column
        row_first <- nx + 2 + (rep(seq(ny), each = 2) - 1)*(nx + 2)
        col_first <- row_first + c(-(nx + 1), (nx + 2), rep(c(-1, 1), ny - 2)*(nx + 2), -(nx + 2), (nx + 1))
        val_first <- rep(c(-0.5, 0.5), ny)/h
        #elements for all columns
        return(tibble::tibble(
          row = i0_row + rep(row_first, nx) + rep(seq(0, nx - 1), each = length(row_first)),
          col = i0_col + rep(col_first, nx) + rep(seq(0, nx - 1), each = length(row_first)),
          val = multiplier*rep(val_first, nx)
        ))
      }
    }
  } else if (order == 2) {
    ## SECOND ORDER DIFFERENTATION
    if (direction == "x") {
      #x-direction
      if (real_only == TRUE) {
        # CHANGE #elements for the first row
        row_first <- c(rep(1, 4), rep(seq(2, nx - 1), each = 3), rep(nx, 4))
        col_first <- c(seq(4), rep(seq(2, nx - 1), each = 3) + rep(c(-1, 0, 1), (nx - 2)), nx - seq(3, 0))
        val_first <- c(2, -5, 4, -1, rep(c(1, -2, 1), (nx - 2)), -1, 4, -5, 2)/h^2
        #elements for all rows
        return(tibble::tibble(
          row = i0_row + rep(row_first, ny) + rep((seq(ny) - 1)*nx, each = length(row_first)),
          col = i0_col + rep(col_first, ny) + rep((seq(ny) - 1)*nx, each = length(col_first)),
          val = multiplier*rep(val_first, ny)
        ))
      } else {
        #elements for the first row
        row_first <- (nx + 1) + rep(seq(nx), each = 3)
        col_first <- row_first + rep(c(-1, 0, 1), nx)
        val_first <- rep(c(1, -2, 1), nx)/h^2
        #elements for all rows
        return(tibble::tibble(
          row = i0_row + rep(row_first, ny) + rep((seq(ny) - 1)*(nx + 2), each = length(row_first)),
          col = i0_col + rep(col_first, ny) + rep((seq(ny) - 1)*(nx + 2), each = length(col_first)),
          val = multiplier*rep(val_first, ny)
        ))
      }
    } else if (direction == "y") {
      #y-direction
      if (real_only == TRUE) {
        #elements for the first column
        row_first <- 1 + (c(rep(1, 4), rep(seq(2, ny - 1), each = 3), rep(ny, 4)) - 1)*nx
        col_first <- 1 + (c(seq(4), rep(seq(2, ny - 1), each = 3) + rep(c(-1, 0, 1), (ny - 2)), ny - seq(3, 0)) - 1)*nx
        val_first <- c(2, -5, 4, -1, rep(c(1, -2, 1), (ny - 2)), -1, 4, -5, 2)/h^2
        #elements for all columns
        return(tibble::tibble(
          row = i0_row + rep(row_first, nx) + rep(seq(0, nx - 1), each = length(row_first)),
          col = i0_col + rep(col_first, nx) + rep(seq(0, nx - 1), each = length(row_first)),
          val = multiplier*rep(val_first, nx)
        ))
      } else {
        #elements for the first column
        row_first <- nx + 2 + (rep(seq(ny), each = 3) - 1)*(nx + 2)
        col_first <- row_first + c(-(nx + 1), 0, (nx + 2), rep(c(-1, 0, 1), ny - 2)*(nx + 2), -(nx + 2), 0, (nx + 1))
        val_first <- rep(c(1, -2, 1), ny)/h^2
        #elements for all columns
        return(tibble::tibble(
          row = i0_row + rep(row_first, nx) + rep(seq(0, nx - 1), each = length(row_first)),
          col = i0_col + rep(col_first, nx) + rep(seq(0, nx - 1), each = length(row_first)),
          val = multiplier*rep(val_first, nx)
        ))
      }
    }
  }
}


#' Get indices of all real nodes
#'
#' @description
#' Get an array of the indices of all real nodes in a rectangular grid. For
#' node numbering, see function `findiff_sparse_entries()`
#'
#' @inheritParams index_edge
#' @return array of indices for all real nodes (non-ghost nodes)
#' @export

index_real <- function(nx, ny, i0 = 0, real_only = FALSE, ...) {
  if (real_only == FALSE) {
    return(i0 + rep(seq(nx), ny) + rep(seq(ny), each = nx)*(nx + 2) - 1)
  } else {
    return(i0 + seq(nx*ny))
  }
}


#' Determine total number of nodes
#'
#' @description
#' Determine the total number of nodes on a grid of nodes, including any
#' ghost nodes outside all four sides of the grid
#'
#' @inheritParams index_edge
#' @return number of nodes (scalar)
#' @export

nodes_total <- function(nx, ny, real_only = FALSE, ...) {
  if (real_only == TRUE) {
    return(nx*ny)
  } else {
    return(nx*ny + 2*(nx + ny))
  }
}


#' Get indices of nodes near edge of grid
#'
#' @description
#' Get an array of the indices of nodes near an edge of a rectangular grid.
#' For node numbering, see function `findiff_sparse_entries()`.
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


#' Get grid index of nodes in specific direction
#'
#' @description
#' Gets the x or y position index from a list of node indices in a rectangular
#' grid. For the node numbering system, see documentation of function
#' `findiff_sparse_entries()`
#'
#' @param i indices of nodes
#' @param nx,ny number of real nodes in x and y-direction
#' @param direction direction of index to find. Should be `x` or `y`
#' @param real_only if `TRUE`, only real nodes are numbered. if `FALSE`, all
#'   ghost nodes are included in the numbering as well
#' @param i0 node starting number (index of first node minus 1)
#' @param real_only if `TRUE`, only include real nodes in numbering,
#' @param ... additional named arguments to pass to function
#' @return vector with domain indices in x or y-direction
#' @examples
#' #number of real nodes in grid
#' nx <- 4
#' ny <- 3
#'
#' #real + ghost nodes
#' i <- seq(nodes_total(nx, ny, real_only = FALSE))
#' index_vector2grid(i, nx, ny, direction = "x", real_only = FALSE)
#' index_vector2grid(i, nx, ny, direction = "y", real_only = FALSE)
#'
#' #real nodes only
#' i <- seq(nodes_total(nx, ny, real_only = TRUE))
#' index_vector2grid(i, nx, ny, direction = "x", real_only = TRUE)
#' index_vector2grid(i, nx, ny, direction = "y", real_only = TRUE)
#' @export

index_vector2grid <- function(i, nx,ny, direction = "x", i0 = 0, real_only = FALSE, ...) {
  if (real_only == TRUE) {
    if (direction == "x") {
      return(1 + ((i - i0) - 1)%%nx)
    } else if (direction == "y") {
      return(1 + floor(((i - i0) - 1)/nx))
    }
  } else {
    if (direction == "x") {
      return(((i - i0) - 1 - nx)%%(nx + 2) - as.double((i - i0) <= nx) + as.double(((i - i0) + 1) >= (nx + 2)*(ny + 1)))
    } else if (direction == "y") {
      return(1 + floor(((i - i0) - 1 - nx)/(nx + 2)))
    }
  }
}


#' Get array index of node based on x and y-indices
#'
#' @description
#' Get the indices in the 1-D node array for rectangle with regularly spaced
#' nodes in a rectangular grid, based on x and y-position indices. For the
#' node numbering system used, see documentation of function
#' `findiff_sparse_entries()`
#'
#' @param nx,ny number of real nodes in x and y-direction
#' @param ix,iy scalar of x and y-index. If one of these is left empty
#'   (`NULL`), all array indices in that row or column are returned
#' @param i0 offset for the index of the first node
#' @param real_only if `TRUE`, the grid is assumed to contain only real nodes,
#'   so there are `nx*ny` nodes. If `FALSE`, the grid is assumed to contain
#'   ghost nodes on the sides of each of the four edges of the grid, resulting
#'   in `nx*ny + 2*nx + 2*ny` nodes. These are numbered in row-wise order,
#'   starting with the left-most ghost node underneath the lower edge
#' @param ... extra arguments
#' @return array of indices
#' @examples
#' index_grid2vector(4, 3, ix = 2, real_only = TRUE)
#' index_grid2vector(4, 3, ix = 2, real_only = FALSE)
#' @export

index_grid2vector <- function(nx, ny, ix = NULL, iy = NULL, i0 = 0, real_only = FALSE, ...) {
  if (real_only == TRUE) {
    if (is.null(ix) & is.null(iy)) {
      return(NULL)
    } else if (is.null(ix)) {
      ix <- seq(nx)
    } else if (is.null(iy)) {
      iy <- seq(ny)
    }
    return(i0 + ix + (iy - 1)*nx)
  } else {
    if (is.null(ix) & is.null(iy)) {
      return(NULL)
    } else if (is.null(ix)) {
      if (iy %in% c(0, ny + 1)) {
        ix <- seq(1, nx)
      } else {
        ix <- seq(0, nx + 1)
      }
    } else if (is.null(iy)) {
      if (ix %in% c(0, nx + 1)) {
        iy <- seq(1, ny)
      } else {
        iy <- seq(0, ny + 1)
      }
    }
    return(ix + iy*(nx + 2) - ifelse(iy > 0, 1, 0) - ifelse(iy > ny, 1, 0))
  }
}


#' Get coordinates of all real nodes in the finite difference grid
#'
#' @description
#' Function obtains the x,y positions of all real nodes in a rectangular
#' finite difference grid
#'
#' @param nx,ny number of real nodes on the grid
#' @param x0,x1,y0,y1 x and y-positions of domain edges
#' @param id array with domain identifiers
#' @param ... additional named arguments to pass to function
#' @importFrom magrittr `%>%`
#' @return a tibble with fields `x` and `y` for positions, and `id` to indicate
#'   which grid the point belongs to
#' @examples
#' df <- nodal_coordinates_real(
#'   nx = c(4, 5),
#'   ny = c(3, 4),
#'   x0 = c(0, 4),
#'   y0 = c(0, 4),
#'   x1 = c(2, 7),
#'   y1 = c(3, 5)
#' )
#' plot(df$x, df$y, "b")
#' @export

nodal_coordinates_real <- function(nx, ny, x0 = 0, y0 = 0, x1 = 1, y1 = 1, id = NULL, ...) {
  df <- tibble::tibble(nx = nx, ny = ny, x0 = x0, y0 = y0, x1 = x1, y1 = y1)
  if (is.null(id)) {
    df$id <- seq(nrow(df))
  } else {
    df$id <- id
  }
  df %>%
    dplyr::rowwise() %>%
    dplyr::summarise(
      id = id,
      ix = rep(seq(nx), ny),
      iy = rep(seq(ny), each = nx),
      x = rep(seq(x0, x1, l = nx), ny),
      y = rep(seq(y0, y1, l = ny), each = nx)
    )
}
