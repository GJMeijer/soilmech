#' Define flow net problem using grided domain system
#'
#' @description
#' Defines a geometry for a flow net problem. rectangular domains are defined
#' using a grid system. For each domain, horizontal and vertical permeabilities
#' may be defined seperately.
#'
#' Domains that are next to each other in the grid are automatically connected.
#' All boundary conditions are assumed impermeable, unless otherwise defined
#' using the inputs `bc_domain`, `bc_edge`, `bc_type` and `bc_value`.
#'
#' @param x x-coordinates of boundary lines in the grid
#' @param y y-coordinates of boundary lines in the grid
#' @param ix x-index in grid of soil rectangle (length n)
#' @param iy y-index in grid of soil rectangle (length n)
#' @param kx,ky horizontal and vertical permeabilities. Either define as
#'   arrays (with length n), or as scalars if the permeability is the same
#'   in each domain
#' @param bc_domain domain at which the boundary condition applies (array with
#'   length m, the integers relate to the position in `ix`,`iy`)
#' @param bc_edge edge in domain at which the boundary condition applies (
#'   array with length m). Edges are numbered clockwise: `1` for left
#'   (negative x-face), `2` for top (positive y-face), `3` for right (
#'   positive x-face) or `4` for bottom (negative y-face)
#' @param bc_type type of boundary condition (array with length m). Use `h`
#'   for a known head, and `q` for a known flow (positive when flowing into
#'   the domain)
#' @param bc_value Values for boundary conditions (array with length m). If
#'   `bc_type == "h"`, it defines the known hydraulic head on the boundary.
#'   If `bc_type == "q"`, it defines the fixed discharge into the domain on
#'   this boundary.
#' @param grid_size the requested grid size. Actual sizes may be slightly
#'   larger to ensure a regular number of nodes fits in each domain. Thus
#'   grid sizes may be slightly different in each domain
#' @param node_min minimum number of nodes in each direction in each
#'   domain. Set to at least 3 to ensure finite difference approximation
#'   works properly
#' @importFrom rlang .data
#' @examples
#' flownet_geometry_rectangular()
#' @export

flownet_geometry_rectangular <- function(
  x = c(0, 15, 20, 40),
  y = c(0, 10, 12, 17.5, 20.0),
  ix = c(1, 1, 1, 1, 2, 3, 3),
  iy = c(4, 3, 2, 1, 1, 1, 2),
  kx = 1e-6,
  ky = 1e-6,
  bc_domain = c(1, 7),
  bc_edge = c(2, 2),
  bc_type  = c("h", "h"),
  bc_value = c(17.5, 15),
  grid_size = 0.25,
  node_min = 3
){
  #generate tibble with all soil domains
  dom <- tibble::tibble(
    ix = ix,
    iy = iy,
    x0 = x[ix],
    y0 = y[iy],
    x1 = x[ix + 1],
    y1 = y[iy + 1],
    kx = kx,
    ky = ky
  )
  #generate numbers of nodes
  dom <- dplyr::mutate(
    dom,
    domain = seq(dplyr::n()),
    nx = pmax(node_min, 1 + ceiling((.data$x1 - .data$x0)/grid_size)),
    ny = pmax(node_min, 1 + ceiling((.data$y1 - .data$y0)/grid_size)),
    hx = (.data$x1 - .data$x0)/(.data$nx - 1),
    hy = (.data$y1 - .data$y0)/(.data$ny - 1)
  )
  #add offsets for node numbering (both with and without ghost nodes)
  dom <- dplyr::mutate(
    dom,
    n = nodes_total(.data$nx, .data$ny, real_only = FALSE),
    n_real = nodes_total(.data$nx, .data$ny, real_only = TRUE),
    i0 = cumsum(c(0, utils::head(.data$n, -1))),
    i0_real = cumsum(c(0, utils::head(.data$n_real, -1))),
  )
  #generate tibble with all default boundary conditions (impermeable boundaries)
  bc <- tibble::tibble(
    domain = rep(seq(nrow(dom)), each = 4),
    edge = rep(seq(4), nrow(dom)),
    type = "q",
    value = 0
  )
  #add connecting boundary conditions
  for (id in 1:nrow(dom)) {
    #other domain connected at the right (positive x-side, edge 3)
    id_next <- which((dom$ix == (dom$ix[id] + 1)) & (dom$iy == dom$iy[id]))
    if (length(id_next) == 1) {
      bc$edge[4*(id - 1) + 3] <- 3
      bc$type[4*(id - 1) + 3] <- "c"
      bc$value[4*(id - 1) + 3] <- id_next
      bc$edge[4*(id_next - 1) + 1] <- 1
      bc$type[4*(id_next - 1) + 1] <- "c"
      bc$value[4*(id_next - 1) + 1] <- id
    }
    #other domain connected at the bottom (positive y-side, edge 2)
    id_next <- which((dom$ix == dom$ix[id]) & (dom$iy == (dom$iy[id] + 1)))
    if (length(id_next) == 1) {
      bc$edge[4*(id - 1) + 2] <- 2
      bc$type[4*(id - 1) + 2] <- "c"
      bc$value[4*(id - 1) + 2] <- id_next
      bc$edge[4*(id_next - 1) + 4] <- 4
      bc$type[4*(id_next - 1) + 4] <- "c"
      bc$value[4*(id_next - 1) + 4] <- id
    }
  }
  #add boundary conditions for head or flow
  bc$domain[4*(bc_domain - 1) + bc_edge] <- bc_domain
  bc$edge[4*(bc_domain - 1) + bc_edge] <- bc_edge
  bc$type[4*(bc_domain - 1) + bc_edge] <- bc_type
  bc$value[4*(bc_domain - 1) + bc_edge] <- bc_value
  #return
  return(list(dom = dom, bc = bc))
}


#' get all sparse entries and lhs for boundary conditions
#'
#' @description
#' function to get sparse matrix element positions for flow net finite
#' difference boundary conditions (known head, known flow or connected to
#' another domain)
#'
#' @param df tibble describing flow problem. Generated by function
#'   `generate_flownet_properties()`
#' @importFrom magrittr `%>%`
#' @importFrom rlang .data
#' @return a list with two fields: `mat` containing a tibble with sparse
#'   row (`row`), column (`col`) and values (`val`) of finite difference
#'   matrix elements related to boundary conditions; and `lhs` which is a
#'   vector with all left-hand side of equation values
#' @examples
#' df <- flownet_geometry_rectangular()
#' ds <- findiff_sparse_entries_bc_rectangular(df)
#' @export

findiff_sparse_entries_bc_rectangular <- function(df) {
  #join domain properties onto boundary conditions
  dbc <- dplyr::left_join(df$bc, df$dom, by = "domain") %>%
    dplyr::mutate(
      k = ifelse((.data$edge %in% c(1, 3)), .data$kx, .data$ky),
      h = ifelse((.data$edge %in% c(1, 3)), .data$hx, .data$hy)
    )
  #initiate list and left-hand sides of bx's
  mlist <- vector(mode = "list", length = nrow(dbc))
  lhs <- rep(0, sum(df$dom$n))
  #initiate list of connecting nodes
  clist <- vector(mode = "list", length = nrow(dbc))
  #loop through all bcs
  for (i in 1:nrow(dbc)) {
    if (dbc$type[i] == "h") {
      #Defined head
      i1 <- index_edge(dbc$nx[i], dbc$ny[i], dbc$edge[i], i0 = dbc$i0[i], offset = 0)
      i2 <- index_edge(dbc$nx[i], dbc$ny[i], dbc$edge[i], i0 = dbc$i0[i], offset = 1)
      mlist[[i]] <- tibble::tibble(
        row = i1,
        col = i2,
        val = 1
      )
      lhs[i1] <- dbc$value[i]
    } else if (dbc$type[i] == "q") {
      #Defined flow rate (positive into domain)
      i1 <- index_edge(dbc$nx[i], dbc$ny[i], dbc$edge[i], i0 = dbc$i0[i], offset = 0)
      i2 <- index_edge(dbc$nx[i], dbc$ny[i], dbc$edge[i], i0 = dbc$i0[i], offset = 2)
      mlist[[i]] <- tibble::tibble(
        row = rep(i1, 2),
        col = c(i1, i2),
        val = rep(c(0.5, -0.5)*dbc$k[i]/dbc$h[i], each = length(i1))
      )
      lhs[i1] <- dbc$value[i]
    } else if (dbc$type[i] == "c") {
      #connected to other domain (b)
      #index in <dbc> of connecting edge
      ib <- which((dbc$domain == dbc$value[i]) & (dbc$edge == (1 + (dbc$edge[i] + 1)%%4)))
      if (dbc$edge[i] %in% c(2, 3)) {
        #other domain connected on positive side - same head
        ia1 <- index_edge(dbc$nx[i], dbc$ny[i], dbc$edge[i], i0 = dbc$i0[i], offset = 0)
        ia2 <- index_edge(dbc$nx[i], dbc$ny[i], dbc$edge[i], i0 = dbc$i0[i], offset = 1)
        ib2 <- index_edge(dbc$nx[ib], dbc$ny[ib], dbc$edge[ib], i0 = dbc$i0[ib], offset = 1)
        mlist[[i]] <- tibble::tibble(
          row = rep(ia1, 2),
          col = c(ia2, ib2),
          val = rep(c(-1, 1), each = length(ia1))
        )
        #add connecting real nodes - real node indexing - used later for solving flow potential
        clist[[i]] <- tibble::tibble(
          i1 = index_edge(dbc$nx[i], dbc$ny[i], dbc$edge[i], i0 = dbc$i0_real[i], offset = 0, real_only = TRUE),
          i2 = index_edge(dbc$nx[ib], dbc$ny[ib], dbc$edge[ib], i0 = dbc$i0_real[ib], offset = 0, real_only = TRUE)
        )
      } else {
        #other domain on negative side - same flow rate
        ia1 <- index_edge(dbc$nx[i], dbc$ny[i], dbc$edge[i], i0 = dbc$i0[i], offset = 0)
        ia2 <- index_edge(dbc$nx[i], dbc$ny[i], dbc$edge[i], i0 = dbc$i0[i], offset = 2)
        ib1 <- index_edge(dbc$nx[ib], dbc$ny[ib], dbc$edge[ib], i0 = dbc$i0[ib], offset = 0)
        ib2 <- index_edge(dbc$nx[ib], dbc$ny[ib], dbc$edge[ib], i0 = dbc$i0[ib], offset = 2)
        mlist[[i]] <- tibble::tibble(
          row = rep(ia1, 4),
          col = c(ia1, ia2, ib1, ib2),
          val = c(
            rep(c(-0.5, 0.5)*dbc$k[i]/dbc$h[i], each = length(ia1)),
            rep(c(-0.5, 0.5)*dbc$k[ib]/dbc$h[ib], each = length(ib1))
          )
        )
      }
    }
  }
  #return
  return(list(
    mat = dplyr::bind_rows(mlist),
    lhs = lhs,
    conn = dplyr::bind_rows(clist)
  ))
}


#' Get coordinates of all real nodes in the finite difference grid
#'
#' @description
#' Function obtains the x,y positions of all real nodes in a rectangular
#' finite difference grid
#'
#' @param nx,ny number of real nodes on the grid
#' @param x0,x1,y0,y1 x and y-positions of domain edges
#' @param domain array with domain identifiers
#' @param ... additional named arguments to pass to function
#' @importFrom magrittr `%>%`
#' @return a tibble with fields `x` and `y` for positions, and `id` to indicate
#'   which grid the point belongs to
#' @examples
#' df <- nodal_coordinates_real_rectangular(
#'   nx = c(4, 5),
#'   ny = c(3, 4),
#'   x0 = c(0, 4),
#'   y0 = c(0, 4),
#'   x1 = c(2, 7),
#'   y1 = c(3, 5)
#' )
#' plot(df$x, df$y, "b")
#' @export

nodal_coordinates_real_rectangular <- function(nx, ny, x0 = 0, y0 = 0, x1 = 1, y1 = 1, domain = NULL, ...) {
  df <- tibble::tibble(nx = nx, ny = ny, x0 = x0, y0 = y0, x1 = x1, y1 = y1)
  if (is.null(domain)) {
    df$domain <- seq(nrow(df))
  } else {
    df$domain <- domain
  }
  df %>%
    dplyr::rowwise() %>%
    dplyr::summarise(
      domain = domain,
      ix = rep(seq(nx), ny),
      iy = rep(seq(ny), each = nx),
      x = rep(seq(x0, x1, l = nx), ny),
      y = rep(seq(y0, y1, l = ny), each = nx)
    )
}


#' Solve flow net head problem for concatenated rectangular domains
#'
#' @description
#' Solves a flow net problem for a problem consisting of connected rectangular
#' grids. This function solves the hydraulic head using boundary conditions and
#' a finite difference approach, and also calculates the flow potential psi
#'
#' @param df tibble with problem description, generated by function
#'   `flownet_geometry_rectangular()`
#' @return a tibble with heads (field `h`), flow rates in x and y-direction
#'   (`qx`, `qy`) and flow potential (`psi`) for each real point in the finite
#'   difference grid. Each point is characterised by a position (`x`, `y`)
#'   and an index `id` indicating which domain the point belongs to
#' @importFrom magrittr `%>%`
#' @importFrom rlang .data
#' @examples
#' df <- flownet_geometry_rectangular()
#' dp <- flownet_solve_rectangular(df)
#' @export

flownet_solve_rectangular <- function(df) {
  ## SOLVE FOR HEAD
  #get sparse matrix elements and left-hand side for boundary conditions
  Mbc <- findiff_sparse_entries_bc_rectangular(df)
  #get sparse matrix elements for Poissons equation for each real node
  M2x <- df$dom %>%
    dplyr::select(.data$nx, .data$ny, .data$kx, .data$x0, .data$x1, .data$i0) %>%
    purrr::pmap_dfr(
      function(nx, ny, kx, x0, x1, i0) findiff_sparse_elements(nx, ny, direction = "xx", i0 = i0, multiplier = kx/(x1 - x0)^2)
    )
  M2y <- df$dom %>%
    dplyr::select(.data$nx, .data$ny, .data$ky, .data$y0, .data$y1, .data$i0) %>%
    purrr::pmap_dfr(
      function(nx, ny, ky, y0, y1, i0) findiff_sparse_elements(nx, ny, direction = "yy", i0 = i0, multiplier = ky/(y1 - y0)^2)
    )
  #create sparse matrix
  mat <- Matrix::sparseMatrix(
    i = c(Mbc$mat$row, M2x$row, M2y$row),
    j = c(Mbc$mat$col, M2x$col, M2y$col),
    x = c(Mbc$mat$val, M2x$val, M2y$val),
    dims = rep(sum(df$dom$n), 2)
  )
  #bind together, and solve system
  #h = as.vector(Matrix::solve(mat, Mbc$lhs))
  h <- as.vector(
    Matrix::solve(
      Matrix::t(mat) %*% mat,
      Matrix::t(mat) %*% Mbc$lhs
    )
  )

  ## CALCULATE FLOW FROM HEAD
  #flow in x-direction
  M1qx <- df$dom %>%
    dplyr::select(.data$nx, .data$ny, .data$kx, .data$x0, .data$x1, .data$i0) %>%
    purrr::pmap_dfr(
      function(nx, ny, kx, x0, x1, i0) findiff_sparse_elements(nx, ny, direction = "x", i0 = i0, multiplier = -kx/(x1 - x0))
    )
  qx <- as.vector(
    Matrix::sparseMatrix(
      i = M1qx$row,
      j = M1qx$col,
      x = M1qx$val,
      dims = rep(sum(df$dom$n), 2)
    ) %*% h
  )
  #flow in y-direction
  M1qy <- df$dom %>%
    dplyr::select(.data$nx, .data$ny, .data$ky, .data$y0, .data$y1, .data$i0) %>%
    purrr::pmap_dfr(
      function(nx, ny, ky, y0, y1, i0) findiff_sparse_elements(nx, ny, direction = "y", i0 = i0, multiplier = -ky/(y1 - y0))
    )
  qy <- as.vector(
    Matrix::sparseMatrix(
      i = M1qy$row,
      j = M1qy$col,
      x = M1qy$val,
      dims = rep(sum(df$dom$n), 2)
    ) %*% h
  )

  ## ASSIGN SOLUTIONS AND POSITIONS (real nodes only)
  #indices of real nodes
  i_real <- purrr::pmap(df$dom, index_real) %>% unlist()
  #generate positions for all solution object for all nodes
  dp <- df$dom %>%
    dplyr::summarise(
      nodal_coordinates_real_rectangular(
        .data$nx,
        .data$ny,
        x0 = .data$x0,
        x1 = .data$x1,
        y0 = .data$y0,
        y1 = .data$y1,
        domain = .data$domain
      )
    )
  #add head and flow for all real points
  dp <- dplyr::mutate(
    dp,
      h = h[i_real],
      qx = qx[i_real],
      qy = qy[i_real]
    )

  ## GET FLOW POTENTIAL FROM FLOW
  #x-direction elements - first order differentiation - real nodes only
  M1x_el <- df$dom %>%
    dplyr::select(.data$nx, .data$ny, .data$x0, .data$x1, .data$i0_real) %>%
    purrr::pmap_dfr(
      function(nx, ny, x0, x1, i0_real) findiff_sparse_elements_realonly(nx, ny, direction = "x", i0 = i0_real, multiplier = 1/(x1 - x0))
    )
  #y-direction elements - first order differentiation - real nodes only
  M1y_el <- df$dom %>%
    dplyr::select(.data$nx, .data$ny, .data$y0, .data$y1, .data$i0_real) %>%
    purrr::pmap_dfr(
      function(nx, ny, y0, y1, i0_real) findiff_sparse_elements_realonly(nx, ny, direction = "y", i0 = i0_real, multiplier = 1/(y1 - y0))
    )
  #flow potential must match on each of the connecting segments
  Mc <- Mbc$conn
  Mc_el <- tibble::tibble(
    row = rep(seq(nrow(Mc)), 2),
    col = c(Mc$i1, Mc$i2),
    val = rep(c(-1, 1), each = nrow(Mc)) / mean(sqrt(df$dom$kx*df$dom$ky))
  )
  #notal number of real nodes
  ntreal <- sum(df$dom$n_real)
  #crate single matrix and left-hand side for flow potential solving
  mat_fp <- Matrix::sparseMatrix(
    i = c(M1x_el$row, M1y_el$row + ntreal, Mc_el$row + 2*ntreal),
    j = c(M1x_el$col, M1y_el$col, Mc_el$col),
    x = c(M1x_el$val, M1y_el$val, Mc_el$val),
    dims = c(2*ntreal + nrow(Mc), ntreal)
  )
  lhs_fp <- c(dp$qy, -dp$qx, rep(0, nrow(Mc)))
  #solve for flow potential
  psi <- as.vector(
    Matrix::solve(
      (Matrix::t(mat_fp) %*% mat_fp),
      (Matrix::t(mat_fp) %*% lhs_fp)
    )
  )
  #scale so flow potential ranges between 0 and Q_total, and assign
  dp$psi <- psi - min(psi)

  ## RETURN
  #return
  return(dp)
}
