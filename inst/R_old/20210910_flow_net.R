#' Grid system
#'
#' @examples
#' flownet_geometry_new()

flownet_geometry_new <- function(
  x = c(0, 15, 20, 40),
  y = c(0, 10, 15, 17.5),
  ix = c(1, 1, 1, 2, 3, 3),
  iy = c(3, 2, 1, 1, 1, 2),
  bc_id = c(1, 6),
  bc_edge = c(2, 2),
  bc_type  = c("h", "h"),
  bc_value = c(20, 15),
  dry = FALSE,
  grid_size = 0.25
){
  #generate tibble with all soil domains
  dom <- tibble::tibble(
    ix = ix,
    iy = iy,
    x0 = x[ix],
    y0 = y[iy],
    x1 = x[ix + 1],
    y1 = x[iy + 1],
    dry = dry
  )
  #generate tibble with all default boundary conditions (impermeable boundaries)
  bc <- tibble::tibble(
    id = rep(which(dom$dry == FALSE), each = 4),
    edge = rep(seq(4), sum(dom$dry == FALSE)),
    type = "q",
    value = 0
  )
  #add boundary conditions for head or flow
  bc$id[4*(bc_id - 1) + bc_edge] <- bc_id
  bc$edge[4*(bc_id - 1) + bc_edge] <- bc_edge
  bc$type[4*(bc_id - 1) + bc_edge] <- bc_type
  bc$value[4*(bc_id - 1) + bc_edge] <- bc_value
  #add connecting boundary conditions
  for (id in 1:nrow(dom)) {
    #other domain connected at the top (edge 3)
    id_next <- which((dom$ix == (dom$ix[id] + 1)) & (dom$iy == dom$iy[id]))
    if (length(id_next) == 1) {
      bc$edge[4*(id - 1) + 3] <- 3
      bc$type[4*(id - 1) + 3] <- "c"
      bc$value[4*(id_next - 1) + 3] <- id_next
      bc$edge[4*(id_next - 1) + 1] <- 1
      bc$type[4*(id_next - 1) + 1] <- "c"
      bc$value[4*(id_next - 1) + 1] <- id
    }
    #other domain connected at the bottom (edge 2)
    id_next <- which((dom$ix == dom$ix[id]) & (dom$iy == (dom$iy[id] + 1)))
    if (length(id_next) == 1) {
      bc$edge[4*(id - 1) + 2] <- 2
      bc$type[4*(id - 1) + 2] <- "c"
      bc$value[4*(id_next - 1) + 2] <- id_next
      bc$edge[4*(id_next - 1) + 4] <- 4
      bc$type[4*(id_next - 1) + 4] <- "c"
      bc$value[4*(id_next - 1) + 4] <- id
    }
  }
  #return
  return(list(dom = dom, bc = bc))
}
a <- flownet_geometry_new()



#' Generate tibble describing flow net problem
#'
#' @description
#' The flow net problem consists of connected rectangular domains. The
#' permabilities can be set seperately for each domain.
#'
#' @md
#' @param x0,y0 x and y position of the bottom left corner of the first domain
#' @param Lx,Ly size in x and y direction of each rectangular domain
#' @param Ly_dry thickness of any dry soil directly on top of this domain.
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
#'   * `type = "q"`: known flow (positive in inflow)
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
#' @param nodes_min the minimum number of real nodes in any direction in any
#'   domain. Set to at least 3 to allow all finite difference approximations
#'   to work
#' @importFrom rlang .data
#' @return a tibble describing the entire tibble. In addition to the input,
#'   fields are added for `x0` and `y0` (the x,y position of the bottom left)
#'   point in each domain, `nx` and `ny` (the number of real nodes in each
#'   domain in each direction), and `hx`, and `hy` (the distance between
#'   nodes in each direction)
#' @examples
#' #create the default flow net
#' generate_flownet_properties()
#' @export

generate_flownet_properties <- function(
  x0 = 0,
  y0 = 0,
  Lx = c(15, 15, 5, 20, 20),
  Ly = c(7.5, 10, 10, 10, 5),
  Ly_dry = c(0, 0, 0, 0, 0),
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
  grid_size = 0.25,
  nodes_min = 3
){
  #create tibble with all input properties
  df <- tibble::tibble(
    Lx = Lx,
    Ly = Ly,
    Ly_dry = Ly_dry,
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
    nx = pmax(nodes_min, ceiling(1 + .data$Lx / grid_size)),
    ny = pmax(nodes_min, ceiling(1 + .data$Ly / grid_size)),
    hx = .data$Lx/(.data$nx - 1),
    hy = .data$Ly/(.data$ny - 1),
    ntotal = nodes_total(.data$nx, .data$ny, real_only = FALSE),
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
      df$x0[ic1] <- df$x0[i] + df$Lx[i]
      df$y0[ic1] <- df$y0[i]
      df$x0[ic2] <- df$x0[i]
      df$y0[ic2] <- df$y0[i] - df$Ly[ic2]
      df$x0[ic3] <- df$x0[i] - df$Lx[ic3]
      df$y0[ic3] <- df$y0[i]
      df$x0[ic4] <- df$x0[i]
      df$y0[ic4] <- df$y0[i] + df$Ly[i]
    }
  }
  #return
  df
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
#' @export

findiff_sparse_entries_bc <- function(df) {
  #get long form for all boundary conditions
  dbc <- tibble::tibble(
    type = c(df$type1, df$type2, df$type3, df$type4),
    value = c(df$value1, df$value2, df$value3, df$value4),
    edge = rep(seq(4), each = nrow(df)),
    id = rep(seq(nrow(df)), 4)
  ) %>%
    dplyr::left_join(
      df %>%
        dplyr::select(dplyr::all_of(c("kx", "ky", "nx", "ny", "hx", "hy", "Lx", "Ly", "x0", "y0", "i0", "i0_real"))) %>%
        dplyr::mutate(id = seq(nrow(df))),
      by = "id"
    ) %>%
    dplyr::mutate(
      k = ifelse((.data$edge %in% c(1, 3)), .data$kx, .data$ky),
      h = ifelse((.data$edge %in% c(1, 3)), .data$hx, .data$hy),
    )
  #initiate list and left-hand sides of bx's
  mlist <- vector(mode = "list", length = nrow(dbc))
  lhs <- rep(0, sum(df$ntotal))
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
      ib <- which((dbc$id == dbc$value[i]) & (dbc$edge == (1 + (dbc$edge[i] + 1)%%4)))
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


#' Solve flow net head problem for concatenated rectangular domains
#'
#' @description
#' Solves a flow net problem for a problem consisting of connected rectangular
#' grids. This function solves the hydraulic head using boundary conditions and
#' a finite difference approach, and also calculates the flow potential psi
#'
#' @param df tibble with problem description, generated by function
#'   `generate_flownet_properties()`
#' @return a tibble with heads (field `h`), flow rates in x and y-direction
#'   (`qx`, `qy`) and flow potential (`psi`) for each real point in the finite
#'   difference grid. Each point is characterised by a position (`x`, `y`)
#'   and an index `id` indicating which domain the point belongs to
#' @importFrom magrittr `%>%`
#' @importFrom rlang .data
#' @examples
#' df <- generate_flownet_properties()
#' dp <- solve_flownet(df)
#' @export

solve_flownet <- function(df) {
  ## SOLVE FOR HEAD
  #get sparse matrix elements and left-hand side for boundary conditions
  Mbc <- findiff_sparse_entries_bc(df)
  #get sparse matrix elements for Poissons equation for each real node
  M2x <- df %>%
    dplyr::select(.data$nx, .data$ny, .data$kx, .data$hx, .data$i0) %>%
    purrr::pmap_dfr(
      function(nx, ny, kx, hx, i0) findiff_sparse_entries(nx, ny, order = 2, h = hx, direction = "x", multiplier = kx, i0_row = i0, i0_col = i0)
    )
  M2y <- df %>%
    dplyr::select(.data$nx, .data$ny, .data$ky, .data$hy, .data$i0) %>%
    purrr::pmap_dfr(
      function(nx, ny, ky, hy, i0) findiff_sparse_entries(nx, ny, order = 2, h = hy, direction = "y", multiplier = ky, i0_row = i0, i0_col = i0)
    )
  #create sparse matrix
  mat <- Matrix::sparseMatrix(
    i = c(Mbc$mat$row, M2x$row, M2y$row),
    j = c(Mbc$mat$col, M2x$col, M2y$col),
    x = c(Mbc$mat$val, M2x$val, M2y$val),
    dims = rep(sum(df$ntotal), 2)
  )
  #bind together, and solve system
  h = Matrix::solve(mat, Mbc$lhs)

  ## CALCULATE FLOW FROM HEAD
  #flow in x-direction
  M1qx <- df %>%
    dplyr::select(.data$nx, .data$ny, .data$kx, .data$hx, .data$i0) %>%
    purrr::pmap_dfr(
      function(nx, ny, kx, hx, i0) findiff_sparse_entries(nx, ny, order = 1, h = hx, direction = "x", multiplier = -kx, i0_row = i0, i0_col = i0)
    )
  qx <- as.vector(
    Matrix::sparseMatrix(
      i = M1qx$row,
      j = M1qx$col,
      x = M1qx$val,
      dims = rep(sum(df$ntotal), 2)
    ) %*% h
  )
  #flow in y-direction
  M1qy <- df %>%
    dplyr::select(.data$nx, .data$ny, .data$ky, .data$hy, .data$i0) %>%
    purrr::pmap_dfr(
      function(nx, ny, ky, hy, i0) findiff_sparse_entries(nx, ny, order = 1, h = hy, direction = "y", multiplier = -ky, i0_row = i0, i0_col = i0)
    )
  qy <- as.vector(
    Matrix::sparseMatrix(
      i = M1qy$row,
      j = M1qy$col,
      x = M1qy$val,
      dims = rep(sum(df$ntotal), 2)
    ) %*% h
  )

  ## ASSIGN SOLUTIONS AND POSITIONS (real nodes only)
  #indices of real nodes
  i_real <- purrr::pmap(df, index_real) %>% unlist()
  #generate solution object for all nodes
  dp <- purrr::pmap_dfr(df, nodal_coordinates_real) %>%
    dplyr::mutate(
      h = h[i_real],
      qx = qx[i_real],
      qy = qy[i_real]
    )

  ## GET FLOW POTENTIAL FROM FLOW
  #x-direction elements - first order differentiation - real nodes only
  M1x_el <- df %>%
    dplyr::select(.data$nx, .data$ny, .data$kx, .data$hx, .data$i0_real) %>%
    purrr::pmap_dfr(
      function(nx, ny, kx, hx, i0_real) findiff_sparse_entries(nx, ny, order = 1, h = hx, direction = "x", i0_row = i0_real, i0_col = i0_real, real_only = TRUE)
    )
  #x-direction elements - first order differentiation - real nodes only
  M1y_el <- df %>%
    dplyr::select(.data$nx, .data$ny, .data$ky, .data$hy, .data$i0_real) %>%
    purrr::pmap_dfr(
      function(nx, ny, ky, hy, i0_real) findiff_sparse_entries(nx, ny, order = 1, h = hy, direction = "y", i0_row = i0_real, i0_col = i0_real, real_only = TRUE)
    )
  #flow potential must match on each of the connecting segments
  Mc <- Mbc$conn
  Mc_el <- tibble::tibble(
    row = rep(seq(nrow(Mc)), 2),
    col = c(Mc$i1, Mc$i2),
    val = rep(c(-1, 1), each = nrow(Mc))
  )
  #notal number of real nodes
  ntreal <- sum(df$nx*df$ny)
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
#' df <- generate_flownet_properties(grid_size = 0.25)
#' ggplot_flownet_only(df, axes = FALSE)
#' ggplot_flownet_only(df, axes = TRUE)
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
    dp <- solve_flownet(df)
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
    Lx = rep(df$Lx, 4),
    Ly = rep(df$Ly, 4)
  ) %>%
    dplyr::filter(.data$type == "h") %>%
    dplyr::mutate(
      x = ifelse(.data$edge == 3, .data$x0 + .data$Lx, .data$x0),
      y = ifelse(.data$edge == 2, .data$y0 + .data$Ly, .data$y0),
      xend = ifelse(.data$edge %in% c(1, 3), .data$x, .data$x + .data$Lx),
      yend = ifelse(.data$edge %in% c(1, 3), .data$y + .data$Ly, .data$y)
    )
  #plot rectangular soil elements
  plt <- ggplot2::ggplot() +
    ggplot2::annotate(
      "rect",
      xmin = df$x0,
      ymin = df$y0,
      xmax = df$x0 + df$Lx,
      ymax = df$y0 + df$Ly,
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
#'   entire problem. Defined as an array for all edges (left, top, right,
#'   bottom), or single value for all edges
#' @param scale_soilmarker size of soil marker, in fraction of height
#' @param scale_watermarker size of water table marker, in fraction of height
#' @return a ggplot object
#' @importFrom magrittr `%>%`
#' @examples
#' #default example
#' df <- generate_flownet_properties()
#' ggplot_flownet_dam(df)
#'
#' #anisotropic permeability - 6 flow channels
#' df <- generate_flownet_properties(kx = 4e-6, ky = 1e-6)
#' ggplot_flownet_dam(df, Nf = 6)
#'
#' #excavation next to river
#' df <- generate_flownet_properties(
#'   Lx = c(20, 20, 0.2, 5, 5),
#'   Ly = c(10, 5, 5, 5, 2),
#'   Ly_dry = c(2, 0, 0, 0, 0),
#'   type1 = c("q", "q", "c", "c", "q"),
#'   type2 = c("h", "c", "q", "c", "h"),
#'   type3 = c("q", "c", "c", "q", "q"),
#'   type4 = c("c", "q", "q", "q", "c"),
#'   value1 = c(0, 0, 2, 3, 0),
#'   value2 = c(12, 1, 0, 5, 6),
#'   value3 = c(0, 3, 4, 0, 0),
#'   value4 = c(2, 0, 0, 0, 4)
#' )
#' ggplot_flownet_dam(df)
#' @export

ggplot_flownet_dam <- function(
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
  thickness_impermeable = c(0, 1, 0, 1),
  scale_soilmarker = 0.075,
  scale_watermarker = 0.050
){
  #if dp not provided, calculate heads and flow potential
  if (is.null(dp)) {
    dp <- solve_flownet(df)
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
    Lx = rep(df$Lx, 4),
    Ly = rep(df$Ly, 4)
  )
  if ("Ly_dry" %in% colnames(df)) {
    dedge$Ly_dry <- rep(df$Ly_dry, 4)
  } else {
    dedge$Ly_dry <- dedge$Ly
  }
  dedge <- dedge %>%
    dplyr::filter(.data$type == "h") %>%
    dplyr::mutate(
      x = ifelse(.data$edge == 3, .data$x0 + .data$Lx, .data$x0),
      y = ifelse(.data$edge == 2, .data$y0 + .data$Ly, .data$y0),
      xend = ifelse(.data$edge %in% c(1, 3), .data$x, .data$x + .data$Lx),
      yend = ifelse(.data$edge %in% c(1, 3), .data$y + .data$Ly, .data$y)
    )
  #impermeable thickness
  if (length(thickness_impermeable) == 1) {
    thickness_impermeable <- rep(thickness_impermeable, 4)
  }
  #plot limits
  if (is.na(xlim[1])) {
    xlim[1] <- min(df$x0) - thickness_impermeable[1]
  }
  if (is.na(xlim[2])) {
    xlim[2] <- max(df$x0 + df$Lx) + thickness_impermeable[3]
  }
  if (is.na(ylim[1])) {
    ylim[1] <- min(df$y0) - thickness_impermeable[4]
  }
  if (is.na(ylim[2])) {
    ylim[2] <- max(c(df$y0 + df$Ly + df$Ly_dry, dedge$value)) + thickness_impermeable[2]
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
      xmax = df$x0 + df$Lx,
      ymax = df$y0 + df$Ly,
      fill = fill_soilwet,
      color = NA
    )
  #plot rectangles for dry soil
  ddry <- dplyr::filter(dedge, .data$Ly_dry > 0)
  plt <- plt +
    ggplot2::geom_rect(
      data = ddry,
      ggplot2::aes(xmin = .data$x, xmax = .data$xend, ymin = .data$y, ymax = .data$y + .data$Ly_dry),
      fill = fill_soildry,
      colour = NA
    )
  #plot whitespace above dry soil water
  plt <- plt +
    ggplot2::geom_rect(
      data = dedge,
      ggplot2::aes(xmin = .data$x, xmax = .data$xend, ymin = .data$y + .data$Ly_dry, ymax = ylim[2]),
      fill = fill_air,
      colour = NA
    )
  #plot rectangles for ponding water
  dwat <- dplyr::filter(dedge, .data$value > .data$y + .data$Ly_dry)
  plt <- plt +
    ggplot2::geom_rect(
      data = dwat,
      ggplot2::aes(xmin = .data$x, xmax = .data$xend, ymin = .data$y + .data$Ly_dry, ymax = .data$value),
      fill = fill_water,
      color = NA
    )
  #plot soil surfaces
  plt <- plt +
    ggplot2::annotate(
      "segment",
      x = dedge$x,
      xend = dedge$xend,
      y = dedge$y + dedge$Ly_dry,
      yend = dedge$y + dedge$Ly_dry,
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
          dedge$y[i] + dedge$Ly_dry[i],
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

