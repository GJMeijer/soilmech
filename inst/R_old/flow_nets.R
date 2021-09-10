#' Get indices of ghost nodes on an edge
#'
#' @description
#' Get an array of the indices of all ghost nodes on a side of a rectangular
#' grid of (real) nodes. Nodes are number per row --> per column. Ghost nodes
#' are placed outside the edges of all sides of the rectangular grid
#'
#' @param nx,ny number of real nodes in x and y-direction
#' @param edge edge number, numbered clockwise. `1` for left (negative x-face),
#'   `2` for top (positive y-face), `3` for right (positive x-face) and `4`
#'   for bottom (negative y-face)
#' @param i0 node starting number.
#' @return array of indices of ghost nodes on selected edge
#' @export

index_ghost <- function(nx, ny, edge, i0 = 0) {
  if (edge == 1) {
    return(i0 + seq(ny)*(nx + 2) - 1)
  } else if (edge == 2) {
    return(i0 + seq(nx) + (nx + 2)*(ny + 1) - 2)
  } else if (edge == 3) {
    return(i0 + (seq(ny) + 1)*(nx + 2) - 2)
  } else if (edge == 4) {
    return(i0 + seq(nx))
  }
}


#' Get indices of all real nodes
#'
#' @description
#' Get an array of the indices of all real nodes in a rectangular grid of
#' (real) nodes. Nodes are number per row --> per column. Ghost nodes
#' are placed outside the edges of all sides of the rectangular grid
#'
#' @inheritParams index_ghost
#' @param ... potential extra arguments. Required to use the fuction with the
#'   `purrr::pmap()` function
#' @return array of indices for all real nodes (non-ghost nodes)
#' @export

index_real <- function(nx, ny, i0 = 0, ...) {
  return(i0 + rep(seq(nx), ny) + rep(seq(ny), each = nx)*(nx + 2) - 1)
}


#' Get index of real nodes on an edge
#'
#' @description
#' Get indices of real nodes lying on an edge. Indices are numbered row-first
#' (x first).
#'
#' @inheritParams index_ghost
#' @param real_only if `FALSE`, ghost nodes are assumed to be included in the
#'   numbering of indices. If `TRUE`, all ghost nodes are assumed to have been
#'   removed, and only real points are numbered
#' @param ... additional arguments, in case the function is passed to
#'   `purrr::pmap()`
#' @return array with indices of real nodes on edge
#' @export

index_real_edge <- function(nx, ny, edge, i0 = 0, real_only = FALSE, ...) {
  if (real_only == TRUE) {
    if (edge == 1) {
      return(i0 + 1 + (seq(ny) - 1)*nx)
    } else if (edge == 2) {
      return(i0 + seq(nx) + nx*(ny - 1))
    } else if (edge == 3) {
      return(i0 + seq(ny)*nx)
    } else if (edge == 4) {
      return(i0 + seq(nx))
    }
  } else {
    if (edge == 1) {
      return(i0 + seq(ny)*(nx + 2))
    } else if (edge == 2) {
      return(i0 + seq(nx) + ny*(nx + 2) - 1)
    } else if (edge == 3) {
      return(i0 + (seq(ny) + 1)*(nx + 2) - 3)
    } else if (edge == 4) {
      return(i0 + seq(nx) + nx + 1)
    }
  }
}


#' Determine total number of nodes
#'
#' @description
#' Determine the total number of nodes on a grid of nodes, including any
#' ghost nodes outside all four sides of the grid
#'
#' @inheritParams index_ghost
#' @param if `real_only = FALSE`, all nodes, including ghost nodes are
#'   counted. If `real_only = TRUE`, only real nodes are counted
#' @return number of nodes (scalar)
#' @export

nodes_total <- function(nx, ny, real_only = FALSE) {
  if (real == FALSE) {
    return(nx*ny + 2*(nx + ny))
  } else {
    return(nx*ny)
  }
}


#' Get indices of nodes bordering nodes
#'
#' @description
#' Get the indices of neighbouring nodes in a grid with real and ghost nodes.
#' This function inherently assumes those nodes exists. Make sure you don't
#' request any nodes outside the grid.
#'
#' @param i array with node indices
#' @param nx,ny number of nodes in x and y-direction
#' @param edge side on which to look, numbered clockwise. `1` for left
#'   (negative x-face), `2` for top (positive y-face), `3` for right
#'   (positive x-face) and `4` for bottom (negative y-face)
#' @importFrom magrittr `%>%`
#' @return an array with node indices for neighbours
#' @export

index_neighbour <- function(i, nx, ny, edge, i0 = 0) {
  tibble::tibble(i = i, nx = nx, ny = ny, edge = edge, i0 = i0) %>%
    dplyr::mutate(
      inew = ifelse(
        edge == 1,
        i - 1,
        ifelse(
          edge == 2,
          ifelse(
            ((i - i0) >= (nx + 2)*ny) | ((i - i0) <= nx),
            i + (nx + 1),
            i + (nx + 2)
          ),
          ifelse(
            edge == 3,
            i + 1,
            ifelse(
              ((i - i0) <= (2*nx + 1)) | ((i - i0) > ((ny + 1)*(nx + 2) - 2)),
              i - (nx + 1),
              i - (nx + 2)
            )
          )
        )
      )
    ) %>%
    dplyr::pull(inew)
}


#' Get finite difference matrix entries for real nodes
#'
#' @description
#' Get row, column and value arrays for all non-zero entries in the second
#' order finite difference grid, for all real points in a grid.
#'
#' @param nx,ny number of nodes in x and y-direction
#' @param hx,hy grid size in x and y-direction
#' @param kx,ky permeability in x and y-direction
#' @param i0 starting node
#' @return tibble with non-zero matrix entries. Fields `row` for row index,
#'   `col` for column index and `val` for values. These can be used to
#'   construct full solver matrix later
#' @export

findiff_entries_real <- function(nx, ny, hx = 1, hy = 1, kx = 1, ky = 1, i0 = 0, ...) {
  #indices of real points and their neighbours
  i <- i0 + index_real(nx, ny)
  i1 <- index_neighbour(i, nx, ny, 1, i0 = i0)
  i2 <- index_neighbour(i, nx, ny, 2, i0 = i0)
  i3 <- index_neighbour(i, nx, ny, 3, i0 = i0)
  i4 <- index_neighbour(i, nx, ny, 4, i0 = i0)
  #generate lists
  tibble::tibble(
    row = c(i, i, i, i, i),
    col = c(i, i1, i3, i4, i2),
    val = c(
      rep(-2*kx/hx^2 - 2*ky/hy^2, nx*ny),
      rep(kx/hx^2, nx*ny),
      rep(kx/hx^2, nx*ny),
      rep(ky/hy^2, nx*ny),
      rep(ky/hy^2, nx*ny)
    )
  )
}


#' Get finite difference matrix entries for all boundaries
#'
#' @description
#' Get row, column and value arrays for all non-zero entries in the problem
#' finite difference grid, for all boundaries of each rectangular grid section
#'
#' @param df tibble with all properties for the rectangular grids.
#' @return tibble with non-zero matrix entries. Fields `row` for row index,
#'   `col` for column index and `val` for values. These can be used to
#'   construct full head solver matrix later
#' @export

findiff_entries_bc <- function(df) {
  #initiate list and left-hand sides
  mlist <- vector(mode = "list", length = 4*nrow(df))
  lhs <- rep(0, sum(df$ntotal))
  #loop through all blocks
  for (i in 1:nrow(df)) {
    ## left edge (1)
    edge <- 1
    ilist <- (i - 1)*4 + edge
    ig <- index_ghost(df$nx[i], df$ny[i], edge, i0 = df$i0[i])
    if (df$type1[i] == "h") {
      mlist[[ilist]] <- tibble::tibble(
        row = ig,
        col = ig + 1,
        val = 1
      )
      lhs[ig] <- df$value1[i]
    } else if (df$type1[i] == "q") {
      mlist[[ilist]] <- tibble::tibble(
        row = c(ig, ig),
        col = c(ig, ig + 2),
        val = -df$kx[i]/df$hx[i]*c(rep(-0.5, length(ig)), rep(0.5, length(ig)))
      )
      lhs[ig] <- df$value1[i]
    } else if (df$type1[i] == "c") {
      ibn <- df$value1[i]
      ign <- index_ghost(df$nx[ibn], df$ny[ibn], 3, i0 = df$i0[ibn])
      #same flow
      mlist[[ilist]] <- tibble::tibble(
        row = c(ig, ig, ig, ig),
        col = c(ig, ig + 2, ign - 2, ign),
        val = c(
          rep(df$kx[i]/df$hx[i], length(ig)),
          rep(-df$kx[i]/df$hx[i], length(ig)),
          rep(-df$kx[ibn]/df$hx[ibn], length(ign)),
          rep(df$kx[ibn]/df$hx[ibn], length(ign))
        )
      )
    }

    ## top edge (2)
    edge <- 2
    ilist <- (i - 1)*4 + edge
    ig <- index_ghost(df$nx[i], df$ny[i], edge, i0 = df$i0[i])
    if (df$type2[i] == "h") {
      mlist[[ilist]] <- tibble::tibble(
        row = ig,
        col = ig - df$nx[i] - 1,
        val = 1
      )
      lhs[ig] <- df$value2[i]
    } else if (df$type2[i] == "q") {
      mlist[[ilist]] <- tibble::tibble(
        row = c(ig, ig),
        col = c(ig - 2*df$nx[i] - 3, ig),
        val = -df$ky[i]/df$hy[i]*c(rep(-0.5, length(ig)), rep(0.5, length(ig)))
      )
      lhs[ig] <- df$value2[i]
    } else if (df$type2[i] == "c") {
      ibn <- df$value2[i]
      ign <- index_ghost(df$nx[ibn], df$ny[ibn], 4, i0 = df$i0[ibn])
      #same head
      mlist[[ilist]] <- tibble::tibble(
        row = c(ig, ig),
        col = c(ig - df$nx[i] - 1, ign + df$nx[ibn] + 1),
        val = c(rep(1, length(ig)), rep(-1, length(ign)))
      )
    }

    ## right edge (3)
    edge <- 3
    ilist <- (i - 1)*4 + edge
    ig <- index_ghost(df$nx[i], df$ny[i], edge, i0 = df$i0[i])
    if (df$type3[i] == "h") {
      mlist[[ilist]] <- tibble::tibble(
        row = ig,
        col = ig - 1,
        val = 1
      )
      lhs[ig] <- df$value3[i]
    } else if (df$type3[i] == "q") {
      mlist[[ilist]] <- tibble::tibble(
        row = c(ig, ig),
        col = c(ig - 2, ig),
        val = -df$kx[i]/df$hx[i]*c(rep(-0.5, length(ig)), rep(0.5, length(ig)))
      )
      lhs[ig] <- df$value3[i]
    } else if (df$type3[i] == "c") {
      ibn <- df$value3[i]
      ign <- index_ghost(df$nx[ibn], df$ny[ibn], 1, i0 = df$i0[ibn])
      #same head
      mlist[[ilist]] <- tibble::tibble(
        row = c(ig, ig),
        col = c(ig - 1, ign + 1),
        val = c(rep(1, length(ig)), rep(-1, length(ign)))
      )
    }

    ## bottom edge (4)
    edge <- 4
    ilist <- (i - 1)*4 + edge
    ig <- index_ghost(df$nx[i], df$ny[i], edge, i0 = df$i0[i])
    if (df$type4[i] == "h") {
      mlist[[ilist]] <- tibble::tibble(
        row = ig,
        col = ig + df$nx[i] + 1,
        val = 1
      )
      lhs[ig] <- df$value4[i]
    } else if (df$type4[i] == "q") {
      mlist[[ilist]] <- tibble::tibble(
        row = c(ig, ig),
        col = c(ig, ig + 2*df$nx[i] + 3),
        val = -df$ky[i]/df$hy[i]*c(rep(-0.5, length(ig)), rep(0.5, length(ig)))
      )
      lhs[ig] <- df$value4[i]
    } else if (df$type4[i] == "c") {
      ibn <- df$value4[i]
      ign <- index_ghost(df$nx[ibn], df$ny[ibn], 2, i0 = df$i0[ibn])
      #same flow
      mlist[[ilist]] <- tibble::tibble(
        row = c(ig, ig, ig, ig),
        col = c(ig, ig + 2*df$nx[i] + 3, ign - 2*df$nx[ibn] - 3, ign),
        val = c(
          rep(df$ky[i]/df$hy[i], length(ig)),
          rep(-df$ky[i]/df$hy[i], length(ig)),
          rep(-df$ky[ibn]/df$hy[ibn], length(ign)),
          rep(df$ky[ibn]/df$hy[ibn], length(ign))
        )
      )
    }
  }
  #return
  return(list(
    lhs = lhs,
    matrix_el = dplyr::bind_rows(mlist)
  ))
}


#' Get coordinates of all real nodes in the finite difference grid
#'
#' @description
#' Function obtains the x,y positions of all real nodes in the finite difference
#' grid
#'
#' @param df tibble with all properties for the rectangular domains
#' @importFrom magrittr `%>%`
#' @return a tibble with fields `x` and `y` for positions, and `id` to indicate
#'   which grid the point belongs to
#' @export

nodal_coordinates_real <- function(df) {
  #get x-y positions relative to bottom left corner point
  dp <- purrr::pmap_dfr(
    df %>% dplyr::mutate(id = seq(nrow(df))),
    function(id, nx, ny, hx, hy, x0, y0, ...) {
      tibble::tibble(
        id = id,
        x = x0 + hx*rep((seq(nx) - 1), ny),
        y = y0 + hy*rep((seq(ny) - 1), each = nx)
      )
    }
  )
  return(dp)
}


#' Generate tibble describing flow net problem
#'
#' @description
#' The flow net problem consists of connected rectangular domains. The
#' permabilities can be set seperately for each domain.
#'
#' @param x0 x and y position of the bottom left corner of the first domain
#' @param width,height width (x-direction) and height (y-direction) of each
#'   domain in the problem
#' @param height_dry thickness of any dry soil directly on top of this domain.
#'   This is only used to generate the soil geometry later. This extra height
#'   is only plotted when the top edge of the domain (edge 2) has a fixed
#'   head as boundary condition
#' @param kx,ky permeabilities in x and y-direction. If scalars, they are set
#'   equal for each domain
#' @param type1,type2,type3,type4 Type of boundary condition on each of the
#'   four edges of each domain. Edges are numbered clockwise:
#'   * `type1` for the left edge (negative x-face)
#'   * `type2` for the top edge (positive y-face)
#'   * `type3` for the right edge (positive x-face)
#'   * `type4` for the bottom edge (negative y-face)
#' Possible values for the boundary condition types are:
#'   * `type = "q"`: known flow (in positive direction)
#'   * `type = "h"`: known head
#'   * `type = "c"`: face is connected to a face in another domain
#' @param value1,value2,value3,value4 values of boundary conditions:
#'   * if `type = "q"`, `value` is the flow rate in positive direction
#'   * if `type = "h"`: `value` is the known head
#'   * if `type = "c"`: `value` gives the number of the domain that the face
#'     is connected to (domains are numbered 1, 2, 3, in order of appearance)
#' @param grid_size requested finite difference grid size. The grid size used
#'   may vary slighly per domain, in order to fit a integer number of nodes in
#'   the domain with a given size
#' @importFrom rlang .data
#' @return a tibble describing the entire tibble. In addition to the input,
#'   fields are added for `x0` and `y0` (the x,y position of the bottom left)
#'   point in each domain, `nx` and `ny` (the number of real nodes in each
#'   domain in each direction), and `hx`, and `hy` (the distance between
#'   nodes in each direction)
#' @export

generate_flownet_properties <- function(
  x0 = 0,
  y0 = 0,
  width = c(15, 15, 5, 20, 20),
  height = c(7.5, 10, 10, 10, 5),
  height_dry = c(0, 0, 0, 0, 0),
  kx = 1e-6,
  ky = 1e-6,
  type1 = c("q","q","c","c","q"),
  type2 = c("h","c","q","c","h"),
  type3 = c("q","c","c","q","q"),
  type4 = c("c","q","q","q","c"),
  value1 = c(0, 0, 2, 3, 0),
  value2 = c(20, 1, 0, 5, 18),
  value3 = c(0, 3, 4, 0, 0),
  value4 = c(2, 0, 0, 0, 4),
  grid_size = 0.25
){
  #create tibble with all input properties
  df <- tibble::tibble(
    width = width,
    height = height,
    height_dry = height_dry,
    kx = kx,
    ky = ky,
    type1 = type1,
    type2 = type2,
    type3 = type3,
    type4 = type4,
    value1 = value1,
    value2 = value2,
    value3 = value3,
    value4 = value4,
  )
  #add grid sizes
  # * nx,ny: the number of grid nodes in x and y-directions
  # * hx,hy: the distance between grid points in x and y-directions
  # * ntotal: the total number of grid points for each rectangle, including
  #           ghost nodes
  # * i0: the index of the first node in each rectangle, minus 1
  # * i0_real: the index of the first node in each rectangle, after removal
  #            of all ghost nodes, minus 1
  df <- dplyr::mutate(
    df,
    nx = pmax(2, ceiling(1 + .data$width / grid_size)),
    ny = pmax(2, ceiling(1 + .data$height / grid_size)),
    hx = .data$width/(.data$nx - 1),
    hy = .data$height/(.data$ny - 1),
    ntotal = nodes_total(.data$nx, .data$ny),
    i0 = cumsum(c(0, utils::head(.data$ntotal, -1))),
    i0_real = cumsum(c(0, utils::head(.data$nx, -1)*utils::head(.data$ny, -1)))
  )
  #positions of origins for each block
  df$x0 <- x0
  df$y0 <- y0
  if (nrow(df) > 1) {
    for (i in seq(nrow(df) - 1)) {
      #find rows of connecting segments
      ic1 <- which((df$type1 == "c") & (df$value1 == i))
      ic1 <- ic1[ic1 > i]
      ic2 <- which((df$type2 == "c") & (df$value2 == i))
      ic2 <- ic2[ic2 > i]
      ic3 <- which((df$type3 == "c") & (df$value3 == i))
      ic3 <- ic3[ic3 > i]
      ic4 <- which((df$type4 == "c") & (df$value4 == i))
      ic4 <- ic4[ic4 > i]
      #adjust position
      df$x0[ic1] <- df$x0[i] + df$width[i]
      df$y0[ic1] <- df$y0[i]
      df$x0[ic2] <- df$x0[i]
      df$y0[ic2] <- df$y0[i] - df$height[ic2]
      df$x0[ic3] <- df$x0[i] - df$width[ic3]
      df$y0[ic3] <- df$y0[i]
      df$x0[ic4] <- df$x0[i]
      df$y0[ic4] <- df$y0[i] + df$height[i]
    }
  }
  #return
  df
}


#' Solve flow net head problem for concatenated rectangular domains
#'
#' @description
#' Solves a flow net problem for a problem consisting of connected rectangular
#' grids. This function solves the hydroloc head using boundary conditions and
#' a finite difference approach.
#'
#' @param df tibble with problem description, generated by function
#'   `generate_flownet_properties()`
#' @return a tibble with heads (field `h`) and flow rates in x and y-direction
#'   (`qx`, `qy`) for each real point in the finite difference grid. Each point
#'   is characterised by a position (`x`, `y`) and an index `id` indicating which
#'   domain the point belongs to
#' @examples
#' df <- generate_flownet_properties()
#' dp <- solve_flownet(df)
#' @export

solve_flownet <- function(df) {
  #generate all matrix entries for all real nodes
  entries_real <- purrr::pmap_dfr(df, findiff_entries_real)
  #generate matrix entries and left-hand side of equation for all boundaries
  entries_bc <- findiff_entries_bc(df)
  #generate sparse matrix
  mat <- Matrix::spMatrix(
    nrow = sum(df$ntotal),
    ncol = sum(df$ntotal),
    i = c(entries_real$row, entries_bc$matrix_el$row),
    j = c(entries_real$col, entries_bc$matrix_el$col),
    x = c(entries_real$val, entries_bc$matrix_el$val)
  )
  #solve matrix equation for head h
  h_sol <- Matrix::solve(mat, entries_bc$lhs, sparse = TRUE)
  #only keep real values
  i_real <- unlist(purrr::pmap(df, index_real))
  #return all unique nodes and positions
  dp <- dplyr::mutate(nodal_coordinates_real(df), h = h_sol[i_real, 1])
  #add flow velocities
  layer <- rep(seq(nrow(df)), df$ntotal)[i_real]
  i1 <- index_neighbour(i_real, df$nx[layer], df$ny[layer], 1, i0 = df$i0[layer])
  i2 <- index_neighbour(i_real, df$nx[layer], df$ny[layer], 2, i0 = df$i0[layer])
  i3 <- index_neighbour(i_real, df$nx[layer], df$ny[layer], 3, i0 = df$i0[layer])
  i4 <- index_neighbour(i_real, df$nx[layer], df$ny[layer], 4, i0 = df$i0[layer])
  dp$qx <- -df$kx[layer]*(0.5*h_sol[i3] - 0.5*h_sol[i1])/df$hx[layer]
  dp$qy <- -df$ky[layer]*(0.5*h_sol[i2] - 0.5*h_sol[i4])/df$hy[layer]
  #return
  return(dp)
}





#matrix elements for a single domain
flow_potential_matrixel_single <- function(x, y, qx, qy, i0, ntotal_real) {
  #initial calcs
  xu <- unique(x)
  yu <- unique(y)
  nx <- length(xu)
  ny <- length(yu)
  hx <- (max(xu) - min(xu)) / (nx - 1)
  hy <- (max(yu) - min(yu)) / (ny - 1)
  #first order differentiation - x-direction
  row_x_single <- c(rep(1, 3), rep(seq(2, nx - 1), each = 2), rep(nx, 3))
  col_x_single <- c(1, 2, 3, rep(seq(2, nx - 1), each = 2) + rep(c(-1, 1), (nx - 2)), nx - c(2, 1, 0))
  val_x_single <- c(-1.5, 2, -0.5, rep(c(-0.5, 0.5), (nx - 2)), 0.5, -2, 1.5)/hx
  row_x <- i0[1] + rep(row_x_single, ny) + rep((seq(ny) - 1)*nx, each = length(row_x_single))
  col_x <- i0[1] + rep(col_x_single, ny) + rep((seq(ny) - 1)*nx, each = length(col_x_single))
  val_x <- rep(val_x_single, ny)
  #first order differentiation - y-direction
  row_y_single <- c(rep(1, 3), rep(seq(2, ny - 1), each = 2), rep(ny, 3))
  col_y_single <- c(1, 2, 3, rep(seq(2, ny - 1), each = 2) + rep(c(-1, 1), (ny - 2)), ny - c(2, 1, 0))
  val_y_single <- c(-1.5, 2, -0.5, rep(c(-0.5, 0.5), (ny - 2)), 0.5, -2, 1.5)/hy
  row_y <- ntotal_real + i0[1] + rep((row_y_single - 1)*nx, nx) + rep(seq(nx), each = length(row_y_single))
  col_y <- i0[1] + rep((col_y_single - 1)*nx, nx) + rep(seq(nx), each = length(col_y_single))
  val_y <- rep(val_y_single, nx)
  #return sparse non-zero elements
  return(tibble::tibble(
    row = c(row_x, row_y),
    col = c(col_x, col_y),
    val = c(val_x, val_y)
  ))
}


#' Calculate the flow potential at every point
#'
#' @description
#' Function takes the solved heads and flow from the finite difference solver
#' for the head, and calculates the flow potential for every real point in the
#' grid. This is done by solving the first-order system of equations
#' dpsi/dx = qy & dpsi/dy = -qx as well as extra conditions to ensure nodes on
#' edge of each domain match
#'
#' @param df tibble with problem prescription
#' @param dp tibble with finite difference solution. Should contain fields for
#'   `id` (domain identifier), `x` and `y` (for position of node) and `qx` and
#'   `qy` for flow rate in x and y-directions
#' @examples
#' df <- generate_flownet_properties()
#' dp <- solve_flownet(df)
#' dp2 <- flow_potential(df, dp)
#' ggplot2::ggplot(
#'   dp2,
#'   ggplot2::aes(x = x, y = y, fill = psi, z = psi)) +
#'   ggplot2::geom_tile() +
#'   ggplot2::geom_contour()
#' @export

flow_potential <- function(df, dp) {
  #number of real nodes
  ntotal_real <- sum(df$nx * df$ny)
  #get all non-zero sparse elements for matrix equations
  mat_el1 <- dp %>%
    dplyr::mutate(i0 = df$i0_real[.data$id]) %>%
    dplyr::group_by(.data$id) %>%
    dplyr::summarize(flow_potential_matrixel_single(.data$x, .data$y, .data$qx, .data$qy, .data$i0, ntotal_real))
  #get left-hand side for matrix equations
  lhs_el1 <- c(-dp$qy, dp$qx)
  #function to get matrix elements for a single connection
  flow_potential_connel_single <- function(df, id1, id2, edge, i0) {
    edge2 <- 1 + (edge + 1)%%4
    i1 <- index_real_edge(df$nx[id1], df$ny[id1], edge, i0 = df$i0_real[id1], real_only = TRUE)
    i2 <- index_real_edge(df$nx[id2], df$ny[id2], edge2, i0 = df$i0_real[id2], real_only = TRUE)
    tibble::tibble(
      row = i0 + rep(seq(length(i1)), 2),
      col = c(i1, i2),
      val = c(rep(1, length(i1)), rep(-1, length(i2)))
    )
  }
  #get all domain connections in positive direction
  dconn <- tibble::tibble(
    type = c(df$type1, df$type2, df$type3, df$type4),
    value = c(df$value1, df$value2, df$value3, df$value4),
    id = rep(seq(nrow(df)), 4),
    edge = rep(seq(4), each = nrow(df))
  ) %>%
    dplyr::filter((.data$type == "c") & (.data$edge %in% c(2, 3))) %>%
    dplyr::mutate(
      i0_temp = ifelse(.data$edge == 2, df$nx[.data$id], df$ny[.data$id]),
      i0_real = cumsum(c(0, utils::head(.data$i0_temp, -1))) + 2*sum(df$nx*df$ny)
    )
  mat_el2 <- dconn %>%
    dplyr::rowwise() %>%
    dplyr::summarise(flow_potential_connel_single(df, .data$id, .data$value, .data$edge, .data$i0_real))
  lhs_el2 <- rep(0, length(unique(mat_el2$row)))
  #bind all matrices and vectors together
  mat <- Matrix::sparseMatrix(
    i = c(mat_el1$row, mat_el2$row),
    j = c(mat_el1$col, mat_el2$col),
    x = c(mat_el1$val, mat_el2$val),
    dims = c(2*ntotal_real + length(lhs_el2), ntotal_real)
  )
  lhs <- c(lhs_el1, lhs_el2)
  #solve and assign
  dp <- dplyr::mutate(dp, psi = Matrix::solve(
    (Matrix::t(mat) %*% mat),
    (Matrix::t(mat) %*% lhs)
  )[, 1])
  #scale
  dp$psi <- dp$psi - min(dp$psi)
  #return
  return(dp)
}


#' Decide on number of equipotential intervals
#'
#' @description
#' This function takes a finite difference solution for both head `h` and flow
#' potential function `psi` to get the best number of equipotential intervals
#' `Nd` based on a given number of flow channels `Nf`. This function assumes
#' the permeabilities in all domains is the same
#'
#' @param df tibble with information about domain
#' @param dp tibble with head and flow potential results for each node
#' @param Nf number of requested flow channels
#' @examples
#' df <- generate_flownet_properties()
#' dp <- solve_flownet(df)
#' dp <- flow_potential(df, dp)
#' amount_equipotentialintervals(df, dp, Nf = 5)
#' @export

amount_equipotentialintervals <- function(df, dp, Nf = 5) {
  #total head difference
  dh <- max(dp$h) - min(dp$h)
  #total flow per second
  Q <- max(dp$psi) - min(dp$psi)
  #k_avg
  k_avg <- mean(sqrt(df$kx*df$ky))
  #Nd
  Nd <- round(k_avg*dh*Nf/Q)
  #return
  return(Nd)
}


#' Plot flow net in part of soil that is saturated
#'
#' @param df tibble with domains/problem description
#' @param dp flow net solution tibble. If not specified, it is calculated
#'   using `df`
#' @param Nf number of requested flow paths
#' @param xlab,ylab x and y-axis label
#' @param fill_soilwet fill colour for saturated soil
#' @param fill_soildry fill colour for dry soil,
#' @param colour_soil line colour for soil
#' @param colour_flowline colour of flow lines
#' @param colour_equipotentialline colour of equipotential lines
#' @param fill_water fill colour for water
#' @param colour_water line colour for water
#' @param line_water line type for water edges (defined head)
#' @param line_width line thickness for all lines
#' @param xlim,ylim user-defined axis limits
#' @param axes if `FALSE` no axes are plotted
#' @return a ggplot object
#' @importFrom magrittr `%>%`
#' @importFrom rlang .data
#' @examples
#' df <- generate_flownet_properties()
#' dp <- solve_flownet(df)
#' dp <- flow_potential(df, dp)
#' ggplot_flownet_only(df, dp, axes = FALSE)
#' @export

ggplot_flownet_only <- function(
  df,
  dp = NULL,
  Nf = 5,
  xlab = "x [m]",
  ylab = "z [m]",
  fill_soilwet = "#aebab7",
  fill_soildry = "#d3bc5f",
  colour_soil = "#65571d",
  colour_flowline = "black",
  colour_equipotentialline = "black",
  fill_water = "#2a7fff",
  colour_water = "#000080",
  line_water = 2,
  line_width = 0.5,
  xlim = c(NA, NA),
  ylim = c(NA, NA),
  axes = TRUE
){
  #if dp not provided, calculate heads and flow potential
  if (is.null(dp)) {
    dp <- flow_potential(df, solve_flownet(df))
  }
  #calculate number of equipotential intervals
  Nd <- amount_equipotentialintervals(df, dp, Nf = Nf)
  Nd_levels <- seq(min(dp$h), max(dp$h), l = Nd + 1)[2:Nd]
  Nf_levels <- seq(min(dp$psi), max(dp$psi), l = Nf + 1)
  #get head edges
  dedge <- tibble::tibble(
    type = c(df$type1, df$type2, df$type3, df$type4),
    value = c(df$value1, df$value2, df$value3, df$value4),
    id = rep(seq(nrow(df)), 4),
    edge = rep(seq(4), each = nrow(df)),
    x0 = rep(df$x0, 4),
    y0 = rep(df$y0, 4),
    width = rep(df$width, 4),
    height = rep(df$height, 4)
  ) %>%
    dplyr::filter(.data$type == "h") %>%
    dplyr::mutate(
      x = ifelse(.data$edge == 3, .data$x0 + .data$width, .data$x0),
      y = ifelse(.data$edge == 2, .data$y0 + .data$height, .data$y0),
      xend = ifelse(.data$edge %in% c(1, 3), .data$x, .data$x + .data$width),
      yend = ifelse(.data$edge %in% c(1, 3), .data$y + height, .data$y)
    )
  #plot rectangular soil elements
  plt <- ggplot2::ggplot() +
    ggplot2::annotate(
      "rect",
      xmin = df$x0,
      ymin = df$y0,
      xmax = df$x0 + df$width,
      ymax = df$y0 + df$height,
      fill = fill_soilwet,
      color = NA
    ) +
    ggplot2::coord_equal(ratio = 1, xlim = xlim, ylim = ylim, expand = FALSE)
  #switch axes off if requested
  if (axes == FALSE) {
    plt <- plt + ggplot2::theme_void()
  } else {
    plt <- plt +
      theme_soilmech() +
      ggplot2::theme(
        panel.background = ggplot2::element_blank(),
        panel.border = ggplot2::element_blank()
      ) +
      ggplot2::xlab(xlab) +
      ggplot2::ylab(ylab)
  }
  #add flow and equipotential lines using contour lines
  plt <- plt +
    ggplot2::geom_contour(
      data = dp,
      ggplot2::aes(x = .data$x, y = .data$y, z = .data$h),
      breaks = Nd_levels,
      color = colour_equipotentialline,
      size = line_width
    ) +
    ggplot2::geom_contour(
      data = dp,
      ggplot2::aes(x = .data$x, y = .data$y, z = .data$psi),
      breaks = Nf_levels,
      color = colour_flowline,
      size = line_width
    )
  #add water tables
  plt <- plt +
    ggplot2::annotate(
      "segment",
      x = dedge$x,
      xend = dedge$xend,
      y = dedge$y,
      yend = dedge$yend,
      linetype = line_water,
      color = colour_water,
      size = line_width
    )
  #return
  return(plt)
}


#' Plot flow net and soil geometry for a dam/sheet pile problem
#'
#' @param df tibble with domains/problem description, generated by function
#'   `generate_flownet_properties()`.#'
#' @param dp solution for heads and flow potential for each node. If not
#'   defined, it is calculated in this plotting function
#' @param Nf number of requested flow paths
#' @param xlab,ylab x and y-axis label
#' @param fill_soilwet fill colour for saturated soil
#' @param fill_soildry fill colour for dry soil,
#' @param colour_soil line colour for soil
#' @param colour_flowline colour of flow lines
#' @param colour_equipotentialline colour of equipotential lines
#' @param fill_water fill colour for water
#' @param colour_water line colour for water
#' @param fill_impermeable fill colour for impermeable parts of the domain
#' @param fill_air fill colour for air on top of soil/water
#' @param line_water line type for water edges (defined head)
#' @param line_width line thickness for all lines
#' @param xlim,ylim user-defined axis limits
#' @param axes if `FALSE` no axes are plotted
#' @param thickness_impermeable thickness of impermeable border around the
#'   entire problem
#' @param scale_soilmarker size of soil marker, in fraction of height
#' @param scale_watermarker size of water table marker, in fraction of height
#' @return a ggplot object
#' @importFrom magrittr `%>%`
#' @examples
#' df <- generate_flownet_properties()
#' ggplot_flownet_full(df, colour_flowline = "red")
#' @export

ggplot_flownet_full <- function(
  df,
  dp = NULL,
  Nf = 5,
  xlab = "x [m]",
  ylab = "z [m]",
  fill_soilwet = "#aebab7",
  fill_soildry = "#d3bc5f",
  colour_soil = "#65571d",
  colour_flowline = "black",
  colour_equipotentialline = "black",
  fill_water = "#2a7fff",
  colour_water = "#000080",
  fill_impermeable = "#333333",
  fill_air = "white",
  line_water = 2,
  line_width = 0.5,
  xlim = c(NA, NA),
  ylim = c(NA, NA),
  axes = TRUE,
  thickness_impermeable = 1,
  scale_soilmarker = 0.075,
  scale_watermarker = 0.050
){
  #if dp not provided, calculate heads and flow potential
  if (is.null(dp)) {
    dp <- flow_potential(df, solve_flownet(df))
  }
  #calculate number of equipotential intervals
  Nd <- amount_equipotentialintervals(df, dp, Nf = Nf)
  Nd_levels <- seq(min(dp$h), max(dp$h), l = Nd + 1)[2:Nd]
  Nf_levels <- seq(min(dp$psi), max(dp$psi), l = Nf + 1)
  #get head edges
  dedge <- tibble::tibble(
    type = c(df$type1, df$type2, df$type3, df$type4),
    value = c(df$value1, df$value2, df$value3, df$value4),
    id = rep(seq(nrow(df)), 4),
    edge = rep(seq(4), each = nrow(df)),
    x0 = rep(df$x0, 4),
    y0 = rep(df$y0, 4),
    width = rep(df$width, 4),
    height = rep(df$height, 4)
  )
  if ("height_dry" %in% colnames(df)) {
    dedge$height_dry <- rep(df$height_dry, 4)
  } else {
    dedge$height_dry <- dedge$height
  }
  dedge <- dedge %>%
    dplyr::filter(.data$type == "h") %>%
    dplyr::mutate(
      x = ifelse(.data$edge == 3, .data$x0 + .data$width, .data$x0),
      y = ifelse(.data$edge == 2, .data$y0 + .data$height, .data$y0),
      xend = ifelse(.data$edge %in% c(1, 3), .data$x, .data$x + width),
      yend = ifelse(.data$edge %in% c(1, 3), .data$y + height, .data$y)
    )
  #plot limits
  if (is.na(xlim[1])) {
    xlim[1] <- min(df$x0) - thickness_impermeable
  }
  if (is.na(xlim[2])) {
    xlim[2] <- max(df$x0 + df$width) + thickness_impermeable
  }
  if (is.na(ylim[1])) {
    ylim[1] <- min(df$y0) - thickness_impermeable
  }
  if (is.na(ylim[2])) {
    ylim[2] <- max(c(df$y0 + df$height + df$height_dry, dedge$value)) + thickness_impermeable
  }
  #plot impermeable polygon
  dimpe <- tibble::tibble(
    x = c(rep(xlim[1], 2), rep(xlim[2], 2)),
    y = c(ylim[1], rep(ylim[2], 2), ylim[1])
  )
  #initiate plot
  plt <- ggplot2::ggplot() +
    ggplot2::coord_equal(ratio = 1, xlim = xlim, ylim = ylim, expand = FALSE)
  #impermeable background
  plt <- plt +
    ggplot2::geom_polygon(
      data = dimpe,
      ggplot2::aes(x = .data$x, y = .data$y),
      color = NA,
      fill = fill_impermeable
    )
  #plot rectangular saturated soil elements
  plt <- plt + ggplot2::annotate(
      "rect",
      xmin = df$x0,
      ymin = df$y0,
      xmax = df$x0 + df$width,
      ymax = df$y0 + df$height,
      fill = fill_soilwet,
      color = NA
    )
  #plot rectangles for dry soil
  ddry <- dplyr::filter(dedge, .data$height_dry > 0)
  plt <- plt +
    ggplot2::geom_rect(
      data = ddry,
      ggplot2::aes(xmin = .data$x, xmax = .data$xend, ymin = .data$y, ymax = .data$y + .data$height_dry),
      fill = fill_soildry,
      colour = NA
    )
  #plot whitespace above dry soil water
  plt <- plt +
    ggplot2::geom_rect(
      data = dedge,
      ggplot2::aes(xmin = .data$x, xmax = .data$xend, ymin = .data$y + .data$height_dry, ymax = ylim[2]),
      fill = fill_air,
      colour = NA
    )
  #plot rectangles for ponding water
  dwat <- dplyr::filter(dedge, .data$value > .data$y + .data$height_dry)
  plt <- plt +
    ggplot2::geom_rect(
      data = dwat,
      ggplot2::aes(xmin = .data$x, xmax = .data$xend, ymin = .data$y + .data$height_dry, ymax = .data$value),
      fill = fill_water,
      color = NA
    )
  #plot soil surfaces
  plt <- plt +
    ggplot2::annotate(
      "segment",
      x = dedge$x,
      xend = dedge$xend,
      y = dedge$y + dedge$height_dry,
      yend = dedge$y + dedge$height_dry,
      color = colour_soil,
      linetype = 1,
      size = line_width
    )
  #plot soil surface marker and water table marker
  for (i in 1:nrow(dedge)) {
    if ((dedge$xend[i] > xlim[1]) & (dedge$x[i] <= xlim[2])) {
      if (!is.na(scale_soilmarker)) {
        plt <- ggplot_add_soilmarker(
          plt,
          1/3*max(xlim[1], dedge$x[i]) + 2/3*min(xlim[2], dedge$xend[i]),
          dedge$y[i] + dedge$height_dry[i],
          scale = scale_soilmarker*diff(ylim)
        )
      }
      if (!is.na(scale_watermarker)) {
        plt <- ggplot_add_watermarker(
          plt,
          2/3*max(xlim[1], dedge$x[i]) + 1/3*min(xlim[2], dedge$xend[i]),
          dedge$value[i],
          scale = scale_watermarker*diff(ylim)
        )
      }
    }
  }
  #switch axes off if requested
  if (axes == FALSE) {
    plt <- plt + ggplot2::theme_void()
  } else {
    plt <- plt +
      theme_soilmech() +
      ggplot2::theme(
        panel.background = ggplot2::element_blank(),
        panel.border = ggplot2::element_blank()
      ) +
      ggplot2::xlab(xlab) +
      ggplot2::ylab(ylab)
  }
  #add flow and equipotential lines using contour lines
  plt <- plt +
    ggplot2::geom_contour(
      data = dp,
      ggplot2::aes(x = .data$x, y = .data$y, z = .data$h),
      breaks = Nd_levels,
      color = colour_equipotentialline,
      size = line_width
    ) +
    ggplot2::geom_contour(
      data = dp,
      ggplot2::aes(x = .data$x, y = .data$y, z = .data$psi),
      breaks = Nf_levels,
      color = colour_flowline,
      size = line_width
    )
  #add water tables
  plt <- plt +
    ggplot2::annotate(
      "segment",
      x = dedge$x,
      xend = dedge$xend,
      y = dedge$y,
      yend = dedge$yend,
      linetype = line_water,
      color = colour_water,
      size = line_width
    )
  #return
  return(plt)
}
