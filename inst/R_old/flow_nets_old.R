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
#' @param if `real = FALSE`, all nodes, including ghost nodes are counted.
#'   if `real = TRUE`, only real nodes are counted
#' @return number of nodes (scalar)
#' @export

nodes_total <- function(nx, ny, real = FALSE) {
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
#' @param db tibble with all properties for the rectangular grids.
#' @return tibble with non-zero matrix entries. Fields `row` for row index,
#'   `col` for column index and `val` for values. These can be used to
#'   construct full solver matrix later
#' @export

findiff_entries_bc <- function(db) {
  #initiate list and left-hand sides
  mlist <- vector(mode = "list", length = 4*nrow(db))
  lhs <- rep(0, sum(db$ntotal))
  #loop through all blocks
  for (i in 1:nrow(db)) {
    ## left edge (1)
    edge <- 1
    ilist <- (i - 1)*4 + edge
    ig <- index_ghost(db$nx[i], db$ny[i], edge, i0 = db$i0[i])
    if (db$type1[i] == "h") {
      mlist[[ilist]] <- tibble::tibble(
        row = ig,
        col = ig + 1,
        val = 1
      )
      lhs[ig] <- db$value1[i]
    } else if (db$type1[i] == "q") {
      mlist[[ilist]] <- tibble::tibble(
        row = c(ig, ig),
        col = c(ig, ig + 2),
        val = -db$kx[i]/db$hx[i]*c(rep(-0.5, length(ig)), rep(0.5, length(ig)))
      )
      lhs[ig] <- db$value1[i]
    } else if (db$type1[i] == "c") {
      ibn <- db$value1[i]
      ign <- index_ghost(db$nx[ibn], db$ny[ibn], 3, i0 = db$i0[ibn])
      #same flow
      mlist[[ilist]] <- tibble::tibble(
        row = c(ig, ig, ig, ig),
        col = c(ig, ig + 2, ign - 2, ign),
        val = c(
          rep(db$kx[i]/db$hx[i], length(ig)),
          rep(-db$kx[i]/db$hx[i], length(ig)),
          rep(-db$kx[ibn]/db$hx[ibn], length(ign)),
          rep(db$kx[ibn]/db$hx[ibn], length(ign))
        )
      )
    }

    ## top edge (2)
    edge <- 2
    ilist <- (i - 1)*4 + edge
    ig <- index_ghost(db$nx[i], db$ny[i], edge, i0 = db$i0[i])
    if (db$type2[i] == "h") {
      mlist[[ilist]] <- tibble::tibble(
        row = ig,
        col = ig - db$nx[i] - 1,
        val = 1
      )
      lhs[ig] <- db$value2[i]
    } else if (db$type2[i] == "q") {
      mlist[[ilist]] <- tibble::tibble(
        row = c(ig, ig),
        col = c(ig - 2*db$nx[i] - 3, ig),
        val = -db$ky[i]/db$hy[i]*c(rep(-0.5, length(ig)), rep(0.5, length(ig)))
      )
      lhs[ig] <- db$value2[i]
    } else if (db$type2[i] == "c") {
      ibn <- db$value2[i]
      ign <- index_ghost(db$nx[ibn], db$ny[ibn], 4, i0 = db$i0[ibn])
      #same head
      mlist[[ilist]] <- tibble::tibble(
        row = c(ig, ig),
        col = c(ig - db$nx[i] - 1, ign + db$nx[ibn] + 1),
        val = c(rep(1, length(ig)), rep(-1, length(ign)))
      )
    }

    ## right edge (3)
    edge <- 3
    ilist <- (i - 1)*4 + edge
    ig <- index_ghost(db$nx[i], db$ny[i], edge, i0 = db$i0[i])
    if (db$type3[i] == "h") {
      mlist[[ilist]] <- tibble::tibble(
        row = ig,
        col = ig - 1,
        val = 1
      )
      lhs[ig] <- db$value3[i]
    } else if (db$type3[i] == "q") {
      mlist[[ilist]] <- tibble::tibble(
        row = c(ig, ig),
        col = c(ig - 2, ig),
        val = -db$kx[i]/db$hx[i]*c(rep(-0.5, length(ig)), rep(0.5, length(ig)))
      )
      lhs[ig] <- db$value3[i]
    } else if (db$type3[i] == "c") {
      ibn <- db$value3[i]
      ign <- index_ghost(db$nx[ibn], db$ny[ibn], 1, i0 = db$i0[ibn])
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
    ig <- index_ghost(db$nx[i], db$ny[i], edge, i0 = db$i0[i])
    if (db$type4[i] == "h") {
      mlist[[ilist]] <- tibble::tibble(
        row = ig,
        col = ig + db$nx[i] + 1,
        val = 1
      )
      lhs[ig] <- db$value4[i]
    } else if (db$type4[i] == "q") {
      mlist[[ilist]] <- tibble::tibble(
        row = c(ig, ig),
        col = c(ig, ig + 2*db$nx[i] + 3),
        val = -db$ky[i]/db$hy[i]*c(rep(-0.5, length(ig)), rep(0.5, length(ig)))
      )
      lhs[ig] <- db$value4[i]
    } else if (db$type4[i] == "c") {
      ibn <- db$value4[i]
      ign <- index_ghost(db$nx[ibn], db$ny[ibn], 2, i0 = db$i0[ibn])
      #same flow
      mlist[[ilist]] <- tibble::tibble(
        row = c(ig, ig, ig, ig),
        col = c(ig, ig + 2*db$nx[i] + 3, ign - 2*db$nx[ibn] - 3, ign),
        val = c(
          rep(db$ky[i]/db$hy[i], length(ig)),
          rep(-db$ky[i]/db$hy[i], length(ig)),
          rep(-db$ky[ibn]/db$hy[ibn], length(ign)),
          rep(db$ky[ibn]/db$hy[ibn], length(ign))
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
#' @param db tibble with all properties for the rectangular grids
#' @return a tibble with fields `x` and `y` for positions, and `id` to indicate
#'   which grid the point belongs to
#' @export

nodal_coordinates_real <- function(db) {
  #get x-y positions relative to bottom left corner point
  dp <- purrr::pmap_dfr(
    db %>% dplyr::mutate(id = seq(nrow(db))),
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
#' @param kx,ky permeabilities in x and y-direction. If scalars, they are set
#'   equal for each domain
#' @param type1,type2,type3,type4 Type of boundary condition on each of the
#'   four edges of each domain. Edges are numbered clockwise:
#'
#'   * `type1` for the left edge (negative x-face)
#'   * `type2` for the top edge (positive y-face)
#'   * `type3` for the right edge (positive x-face)
#'   * `type4` for the bottom edge (negative y-face)
#'
#' Possible values for the boundary condition types are:
#'
#'   * `type = "q"`: known flow (in positive direction)
#'   * `type = "h"`: known head
#'   * `type = "c"`: face is connected to a face in another domain
#'
#' @param value1,value2,value3,value4 values of boundary conditions:
#'
#'   * if `type = "q"`, `value` is the flow rate in positive direction
#'   * if `type = "h"`: `value` is the known head
#'   * if `type = "c"`: `value` gives the number of the domain that the face
#'     is connected to (domains are numbered 1, 2, 3, in order of appearance)
#'
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
  width = c(15, 15, 10),
  height = c(5, 10, 10),
  kx = 1e-6,
  ky = 1e-6,
  type1 = c("q","q","c"),
  type2 = c("h","c","q"),
  type3 = c("q","c","h"),
  type4 = c("c","q","q"),
  value1 = c(0, 0, 2),
  value2 = c(20, 1, 0),
  value3 = c(0, 3, 10),
  value4 = c(2, 0, 0),
  grid_size = 0.5
){
  #create tibble with all input properties
  df <- tibble::tibble(
    width = width,
    height = height,
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
    nx = ceiling(1 + .data$width / grid_size),
    ny = ceiling(1 + .data$height / grid_size),
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


#' Solve a flow net problem for concatenated rectangular grids
#'
#' @description
#' Solves a flow net problem for a problem consisting of connected rectangular
#' grids.
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
  dp$qx <- -df$kx[layer]*(0.5*h_sol[i3] - 0.5*h_sol[i1])*df$hx[layer]
  dp$qy <- -df$ky[layer]*(0.5*h_sol[i2] - 0.5*h_sol[i4])*df$hy[layer]
  #return
  return(dp)
}


#' @examples
#' df <- generate_flownet_properties()
#' dp <- solve_flownet(df)
#' dp2 <- flow_potential(df, dp)
#' ggplot2::ggplot(dp2, ggplot2::aes(x = x, y = y, fill = psi)) + ggplot2::geom_tile()

flow_potential <- function(df, dp) {
  #function for single domain
  flow_potential_single <- function(x, y, qx, qy) {
    #initial calcs
    xu <- sort(unique(x))
    yu <- sort(unique(y))
    nx <- length(xu)
    ny <- length(yu)
    hx <- (max(xu) - min(xu)) / (nx - 1)
    hy <- (max(yu) - min(yu)) / (ny - 1)
    #first order differentiation - x-direction
    row_x_single <- c(rep(1, 3), rep(seq(2, nx - 1), each = 2), rep(nx, 3))
    col_x_single <- c(1, 2, 3, rep(seq(2, nx - 1), each = 2) + rep(c(-1, 1), (nx - 2)), nx - c(2, 1, 0))
    val_x_single <- c(-1.5, 2, -0.5, rep(c(-0.5, 0.5), (nx - 2)), 0.5, -2, 1.5)/hx
    row_x <- rep(row_x_single, ny) + rep((seq(ny) - 1)*nx, each = length(row_x_single))
    col_x <- rep(col_x_single, ny) + rep((seq(ny) - 1)*nx, each = length(col_x_single))
    val_x <- rep(val_x_single, ny)
    #first order differentiation - y-direction
    row_y_single <- c(rep(1, 3), rep(seq(2, ny - 1), each = 2), rep(ny, 3))
    col_y_single <- c(1, 2, 3, rep(seq(2, ny - 1), each = 2) + rep(c(-1, 1), (ny - 2)), ny - c(2, 1, 0))
    val_y_single <- c(-1.5, 2, -0.5, rep(c(-0.5, 0.5), (ny - 2)), 0.5, -2, 1.5)/hy
    row_y <- rep((row_y_single - 1)*nx, nx) + rep(seq(nx), each = length(row_y_single))
    col_y <- rep((col_y_single - 1)*nx, nx) + rep(seq(nx), each = length(col_y_single))
    val_y <- rep(val_y_single, nx)
    #generate single sparse matrix
    mat <- Matrix::sparseMatrix(
      i = c(row_x, row_y + nx*ny),
      j = c(col_x, col_y),
      x = c(val_x, val_y),
      dims = c(2*nx*ny, nx*ny)
    )
    #lhs
    lhs <- c(-qy, qx)
    #solve inverse problem
    psi <- Matrix::solve((Matrix::t(mat) %*% mat), (Matrix::t(mat) %*% lhs))
    #return
    as.vector(psi)
  }
  #apply to each domain
  dp <- dp %>%
    dplyr::group_by(.data$id) %>%
    dplyr::mutate(psi = flow_potential_single(.data$x, .data$y, .data$qx, .data$qy))
  #get all connections
  dconn <- tibble::tibble(
    type = c(df$type1, df$type2, df$type3, df$type4),
    value = c(df$value1, df$value2, df$value3, df$value4),
    id = rep(seq(nrow(df)), 4),
    edge = rep(seq(4), each = nrow(df))
  ) %>%
    dplyr::filter((.data$type == "c") & (.data$value > .data$id)) %>%
    dplyr::arrange(id, value, edge)
  #loop through connections
  for (i in 1:nrow(dconn)) {
    iedge1 <- index_real_edge(
      df$nx[dconn$id[i]],
      df$ny[dconn$id[i]],
      dconn$edge[i],
      i0 = df$i0_real[dconn$id[i]],
      real_only = TRUE
    )
    iedge2 <- index_real_edge(
      df$nx[dconn$value[i]],
      df$ny[dconn$value[i]],
      1 + (dconn$edge[i] + 1)%%4,
      i0 = df$i0_real[dconn$value[i]],
      real_only = TRUE
    )
    psi_adjust <- dp$psi[iedge1] - dp$psi[iedge2]
    ind <- (dp$id == dconn$value[i])
    if (dconn$edge[i] %in% c(1, 3)) {
      dp$psi[ind] <- dp$psi[ind] + rep(psi_adjust, each = df$nx[dconn$value[i]])
    } else {
      dp$psi[ind] <- dp$psi[ind] + rep(psi_adjust, df$ny[dconn$value[i]])
    }
  }
  #scale to between zero and one
  dp$psi <- (dp$psi - min(dp$psi))/(max(dp$psi) - min(dp$psi))
  #return
  return(dp)
}



#' @examples
#' df <- generate_flownet_properties()
#' dp <- solve_flownet(df)
#' dp2 <- flow_potential(df, dp)
#' ggplot2::ggplot(dp2, ggplot2::aes(x = x, y = y, z = psi, fill = psi)) + ggplot2::geom_tile() + ggplot2::geom_contour()

flow_potential_old3 <- function(df, dp) {
  #integrate per domain (integral = 0 at bottom-left corner of each domain)
  dp <- dp %>%
    dplyr::group_by(id, y) %>%
    dplyr::mutate(psix = c(0, cumsum(diff(.data$x)*0.5*(utils::head(.data$qy, -1) + utils::tail(.data$qy, -1)))))
  dp <- dp %>%
    dplyr::group_by(id, x) %>%
    dplyr::mutate(psiy = -c(0, cumsum(diff(.data$y)*0.5*(utils::head(.data$qx, -1) + utils::tail(.data$qx, -1)))))
  dp <- dp %>%
    dplyr::mutate(psi = psix + psiy) %>%
    dplyr::select(id, x, y, h, qx, qy, psi)
  #get all connections
  dconn <- tibble::tibble(
    type = c(df$type1, df$type2, df$type3, df$type4),
    value = c(df$value1, df$value2, df$value3, df$value4),
    id = rep(seq(nrow(df)), 4),
    edge = rep(seq(4), each = nrow(df))
  ) %>%
    dplyr::filter((.data$type == "c") & (.data$value > .data$id)) %>%
    dplyr::arrange(id, value, edge)
  #loop through connections
  for (i in 1:nrow(dconn)) {
    iedge1 <- index_real_edge(
      df$nx[dconn$id[i]],
      df$ny[dconn$id[i]],
      dconn$edge[i],
      i0 = df$i0_real[dconn$id[i]],
      real_only = TRUE
    )
    iedge2 <- index_real_edge(
      df$nx[dconn$value[i]],
      df$ny[dconn$value[i]],
      1 + (dconn$edge[i] + 1)%%4,
      i0 = df$i0_real[dconn$value[i]],
      real_only = TRUE
    )
    psi_adjust <- dp$psi[iedge1] - dp$psi[iedge2]
    ind <- (dp$id == dconn$value[i])
    if (dconn$edge[i] %in% c(1, 3)) {
      dp$psi[ind] <- dp$psi[ind] + rep(psi_adjust, each = df$nx[dconn$value[i]])
    } else {
      dp$psi[ind] <- dp$psi[ind] + rep(psi_adjust, df$ny[dconn$value[i]])
    }
  }
  #scale
  #dp$psi <- dp$psi - min(dp$psi)
  #return
  dp
}


flow_potential_old <- function(x, y, qx, qy, ...) {
  xu <- sort(unique(x))
  yu <- sort(unique(y))
  nx <- length(xu)
  ny <- length(yu)
  hx <- (max(xu) - min(xu)) / (nx - 1)
  hy <- (max(yu) - min(yu)) / (ny - 1)
  qxm <- matrix(qx, nrow = nx, byrow = TRUE)
  qym <- matrix(qy, nrow = nx, byrow = TRUE)
  qint <- (
    apply(
      rbind(matrix(0, nrow = 1, ncol = ny), 0.5*hx*(qxm[2:nx, ] + qxm[1:(nx - 1), ])),
      2,
      cumsum
    ) +
    t(apply(
      cbind(matrix(0, nrow = nx, ncol = 1), 0.5*hy*(qym[, 2:ny] + qym[, 1:(ny - 1)])),
      1,
      cumsum
    ))
  )
  c(qint)
}


get_streamlines <- function(df, dp, Nf = 4, nt = 100) {
  #get domain and edge which highest <h> applied
  dh <- tibble::tibble(
    type = c(df$type1, df$type2, df$type3, df$type4),
    value = c(df$value1, df$value2, df$value3, df$value4),
    id = rep(seq(nrow(df)), 4),
    edge = rep(seq(4), each = nrow(df))
  ) %>%
    dplyr::filter(.data$type == "h") %>%
    dplyr::slice_max(.data$value)
  #get real nodes on this edge
  iedge <- index_real_edge(
    df$nx[dh$id],
    df$ny[dh$id],
    dh$edge,
    i0 = df$i0_real[dh$id],
    real_only = TRUE
  )
  #integrate total flow through edge (trapezoid integration)
  qedge <- sqrt(dp$qx[iedge]^2 + dp$qy[iedge]^2)
  Ledge <- sqrt(diff(dp$x[iedge])^2 + diff(dp$y[iedge])^2)
  qint <- c(0, cumsum(Ledge*0.5*(qedge[1:(length(qedge) - 1)] + qedge[2:length(qedge)])))
  #starting positions of streampaths
  xf0 <- stats::approx(qint, dp$x[iedge], xout = max(qint)*seq(1/Nf, 1 - 1/Nf, l = Nf - 1))$y
  yf0 <- stats::approx(qint, dp$y[iedge], xout = max(qint)*seq(1/Nf, 1 - 1/Nf, l = Nf - 1))$y
  #estimated number of time steps, based on average flow rate
  qavg <- mean(sqrt(dp$qx^2 + dp$qy^2))
  Lavg <- 0.5*(sum(df$width) + sum(df$width))
  dt <- Lavg / qavg / (nt)
  #initiate vectors
  nt_max <- nt*10
  #function to get gradients
  get_gradient <- function(df, x, y) {
    #check if points within domains
    in_domain <- ((x >= df$x0) & (x <= (df$x0 + df$width)) & (y >= df$y0) & (y <= (df$y0 + df$height)))
    #check if in domain
    if (all(in_domain == FALSE)) {
      #outside domain - zero gradients
      return(c(0, 0))
    } else {
      #get id of box
      id <- which.max(in_domain)
      #interpolate flow rate
      di <- interp_bilinear_longform(
        dp$x[dp$id == id],
        dp$y[dp$id == id],
        xout = x,
        yout = y,
        qx = dp$qx[dp$id == id],
        qy = dp$qy[dp$id == id]
      )
      #return
      return(c(dplyr::pull(di, "qx"), dplyr::pull(di, "qy")))
    }
  }
  #function for Runge-Kutta integration
  path_rk4 <- function(df, x, y, dt) {
    k1 <- get_gradient(df, x, y)
    k2 <- get_gradient(df, x + k1[1]*dt/2, y + k1[2]*dt/2)
    k3 <- get_gradient(df, x + k2[1]*dt/2, y + k2[2]*dt/2)
    k4 <- get_gradient(df, x + k3[1]*dt, y + k3[2]*dt)
    c(x, y) + dt/6*(k1 + 2*k2 + 2*k3 + k4)
  }
  #loop trough all flow channels
  for (j in seq(Nf - 1)) {
    xf <- rep(NA, nt_max)
    yf <- rep(NA, nt_max)
    xf[1] <- xf0[j]
    yf[1] <- yf0[j]
    #loop
    i <- 1
    while(i <= nt_max) {
      #find domain
      in_domain <- ((xf[i] >= df$x0) & (xf[i] <= (df$x0 + df$width)) & (yf[i] >= df$y0) & (yf[i] <= (df$y0 + df$height)))
      #break if point has left the domain
      if (all(in_domain == FALSE)) {
        break
      } else {
        #runge kutta integration to get new position
        xynew <- path_rk4(df, xf[i], yf[i], dt)
        xf[i + 1] <- xynew[1]
        yf[i + 1] <- xynew[2]
        #update index
        i <- i + 1
      }
    }
    #generate tibble
    if (j == 1) {
      ds <- tibble::tibble(id = j, t = (seq(nt_max) - 1)*dt, x = xf, y = yf) %>% tidyr::drop_na()
    } else {
      ds <- dplyr::bind_rows(
        ds,
        tibble::tibble(id = j, t = (seq(nt_max) - 1)*dt, x = xf, y = yf) %>% tidyr::drop_na()
      )
    }
  }
  #return
  return(ds)
}



#' Plot flow net
#' @param dp flow net solution tibble
#' @examples
#' df <- generate_flownet_properties(grid_size = 0.5)
#' dp <- solve_flownet(df)
#' ds <- get_streamlines(df, dp)
#' ggplot_flownet(dp, ds = ds)
#' ggplot_flownet(dp)

ggplot_flownet <- function(dp, ds = NULL, Nd = 5){
  plt <- ggplot2::ggplot() +
    ggplot2::geom_tile(
      data = dp,
      ggplot2::aes(x = x, y = y, fill = psi)
    ) +
    ggplot2::geom_contour(
      data = dp,
      ggplot2::aes(x = x, y = y, z = psi),
      bins = Nd + 1
    ) +
    ggplot2::coord_equal(ratio = 1)
  if (!is.null(ds)) {
    plt <- plt +
      ggplot2::geom_path(
        data = ds,
        ggplot2::aes(x = x, y = y, group = as.factor(id))
      )
  }
  return(plt)
}


#' Bilinear interpolation of data in long format
#'
#' @description
#' Bilinearly interpolates values in 2-dimensional fields.
#'
#' @param x array with x-data for each point. Should be sorted in row-first
#'   order, e.g. 1, 2, 3, 1, 2, 3, 1, 2, 3 for a 3x3 matrix
#' @param x array with x-data for each point. Should be sorted in row-first
#'   order, e.g. 1, 1, 1, 2, 2, 2, 3, 3, 3 for a 3x3 matrix
#' @param xout,yout arrays with x and y positions for output
#' @param ... named fields to be interpolated. fields should have same length
#'   as `x` and `y`
#' @return a tibble with fields `x` and `y` for position and columns
#'   `names(...)` with interpolated values for each fields in input `...`
#' @export

interp_bilinear_longform <- function(x, y, xout, yout, ...) {
  #convert input to be interpolated to tibble
  z <- tibble::as_tibble(list(...))
  #unique x and y values
  xu <- sort(unique(x))
  yu <- sort(unique(y))
  #length of vectors
  nx <- length(xu)
  ny <- length(yu)
  #relative index of requested values in a grid position
  ix <- 1 + (xout - min(xu)) / (max(xu) - min(xu)) * (nx - 1)
  iy <- 1 + (yout - min(yu)) / (max(yu) - min(yu)) * (ny - 1)
  #interpolate and return
  tibble::tibble(
    x = xout,
    y = yout
  ) %>%
    dplyr::mutate(
      z[floor(ix) + nx*(floor(iy) - 1), ]*(1 - ix%%1)*(1 - iy%%1) +
      z[ceiling(ix) + nx*(floor(iy) - 1), ]*(ix%%1)*(1 - iy%%1) +
      z[floor(ix) + nx*(ceiling(iy) - 1), ]*(1 - ix%%1)*(iy%%1) +
      z[ceiling(ix) + nx*(ceiling(iy) - 1), ]*(ix%%1)*(iy%%1)
    )
}







#' Bilinear interpolation on a regular grid
#'
#' @description
#' Interpolate on a 2d grid using bilinear interpolation
#'
#' @param x,y arrays of grid x and y-values in increasing order
#' @param z matrix of values to interpolate. Matrix has `length(x)` number
#'   of rows and `length(y)` number of columns
#' @param xout,yout arrays of x and y values at which to interpolate
#' @examples
#' x = c(1, 2, 3)
#' y = c(2, 3)
#' z = matrix(c(1, 2, 2, 4, 4, 7), nrow = 3, byrow = TRUE)
#' xout = c(1.0, 2.4)
#' yout = c(2.6, 3.0)
#' interp_bilinear(x, y, z, xout, yout)
#' @export

interp_bilinear <- function(x, y, z, xout, yout) {
  #relative index of requested values
  ix <- 1 + (xout - min(x)) / (max(x) - min(x)) * (length(x) - 1)
  iy <- 1 + (yout - min(y)) / (max(y) - min(y)) * (length(y) - 1)
  #interpolation
  return(
    z[matrix(c(floor(ix), floor(iy)), ncol = 2)]*(1 - ix%%1)*(1 - iy%%1) +
    z[matrix(c(floor(ix), ceiling(iy)), ncol = 2)]*(1 - ix%%1)*(iy%%1) +
    z[matrix(c(ceiling(ix), floor(iy)), ncol = 2)]*(ix%%1)*(1 - iy%%1) +
    z[matrix(c(ceiling(ix), ceiling(iy)), ncol = 2)]*(ix%%1)*(iy%%1)
  )
}
