#' Define flow net problem using grided domain system
#'
#' @description
#' Defines a geometry for a flow net problem. rectangular domains are defined
#' using a grid system. For each domain, horizontal and vertical permeabilities
#' may be defined seperately.
#'
#' Domains that are next to each other in the grid are automatically connected.
#' All boundary conditions are assumed impermeable, unless otherwise defined
#' using the inputs `bc_id`, `bc_edge`, `bc_type` and `bc_value`.
#'
#' @param x x-coordinates of boundary lines in the grid
#' @param y y-coordinates of boundary lines in the grid
#' @param ix x-index in grid of soil rectangle (length n)
#' @param iy y-index in grid of soil rectangle (length n)
#' @param kx,ky horizontal and vertical permeabilities. Either define as
#'   arrays (with length n), or as scalars if the permeability is the same
#'   in each domain
#' @param bc_id domain at which the boundary condition applies (array with
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
#' @param dry if `TRUE`, the soil domain is dry (above the water table) and
#'   will not be included in the calculations (but may be plotted). Array
#'   size n
#' @param grid_size the requested grid size. Actual sizes may be slightly
#'   larger to ensure a regular number of nodes fits in each domain. Thus
#'   grid sizes may be slightly different in each domain
#' @param node_min minimum number of nodes in each direction in each
#'   domain. Set to at least 3 to ensure finite difference approximation
#'   works properly
#' @importFrom rlang .data
#' @examples
#' flownet_geometry()
#' @export

flownet_geometry <- function(
  x = c(0, 15, 20, 40),
  y = c(0, 10, 12, 17.5, 20.0),
  ix = c(1, 1, 1, 1, 2, 3, 3),
  iy = c(4, 3, 2, 1, 1, 1, 2),
  kx = 1e-6,
  ky = 1e-6,
  bc_id = c(2, 7),
  bc_edge = c(2, 2),
  bc_type  = c("h", "h"),
  bc_value = c(17.5, 15),
  dry = c(TRUE, rep(FALSE, 6)),
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
    ky = ky,
    dry = dry
  )
  #generate numbers of nodes
  dom <- dplyr::mutate(
    dom,
    id = seq(dplyr::n()),
    nx = pmax(node_min, 1 + ceiling((.data$x1 - .data$x0)/grid_size)),
    ny = pmax(node_min, 1 + ceiling((.data$y1 - .data$y0)/grid_size)),
    hx = (.data$x1 - .data$x0)/(.data$nx - 1),
    hy = (.data$y1 - .data$y0)/(.data$ny - 1)
  )
  #add offsets for node numbering (both with and without ghost nodes)
  dom <- dplyr::mutate(
    dom,
    n = ifelse((.data$dry == FALSE), nodes_total(.data$nx, .data$ny, real_only = FALSE), 0),
    n_real = ifelse((.data$dry == FALSE), nodes_total(.data$nx, .data$ny, real_only = TRUE), 0),
    i0 = cumsum(c(0, utils::head(.data$n, -1))),
    i0_real = cumsum(c(0, utils::head(.data$n_real, -1))),
  )
  #generate tibble with all default boundary conditions (impermeable boundaries)
  bc <- tibble::tibble(
    id = rep(seq(nrow(dom)), each = 4),
    edge = rep(seq(4), nrow(dom)),
    type = "q",
    value = 0
  )
  #add connecting boundary conditions
  for (id in which(dom$dry == FALSE)) {
    #other domain connected at the right (positive x-side, edge 3)
    id_next <- which((dom$ix == (dom$ix[id] + 1)) & (dom$iy == dom$iy[id]))
    if (length(id_next) == 1) {
      if (dom$dry[id_next] == FALSE) {
        bc$edge[4*(id - 1) + 3] <- 3
        bc$type[4*(id - 1) + 3] <- "c"
        bc$value[4*(id - 1) + 3] <- id_next
        bc$edge[4*(id_next - 1) + 1] <- 1
        bc$type[4*(id_next - 1) + 1] <- "c"
        bc$value[4*(id_next - 1) + 1] <- id
      }
    }
    #other domain connected at the bottom (positive y-side, edge 2)
    id_next <- which((dom$ix == dom$ix[id]) & (dom$iy == (dom$iy[id] + 1)))
    if (length(id_next) == 1) {
      if (dom$dry[id_next] == FALSE) {
        bc$edge[4*(id - 1) + 2] <- 2
        bc$type[4*(id - 1) + 2] <- "c"
        bc$value[4*(id - 1) + 2] <- id_next
        bc$edge[4*(id_next - 1) + 4] <- 4
        bc$type[4*(id_next - 1) + 4] <- "c"
        bc$value[4*(id_next - 1) + 4] <- id
      }
    }
  }
  #add boundary conditions for head or flow
  bc$id[4*(bc_id - 1) + bc_edge] <- bc_id
  bc$edge[4*(bc_id - 1) + bc_edge] <- bc_edge
  bc$type[4*(bc_id - 1) + bc_edge] <- bc_type
  bc$value[4*(bc_id - 1) + bc_edge] <- bc_value
  #exclude boundary conditions for dry domains
  bc <- bc[bc$id %in% which(dom$dry == FALSE), ]
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
#' df <- flownet_geometry()
#' ds <- findiff_sparse_entries_bc(df)
#' @export

findiff_sparse_entries_bc <- function(df) {
  #join domain properties onto boundary conditions
  dbc <- dplyr::left_join(df$bc, df$dom, by = "id") %>%
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
#'   `flownet_geometry()`
#' @return a tibble with heads (field `h`), flow rates in x and y-direction
#'   (`qx`, `qy`) and flow potential (`psi`) for each real point in the finite
#'   difference grid. Each point is characterised by a position (`x`, `y`)
#'   and an index `id` indicating which domain the point belongs to
#' @importFrom magrittr `%>%`
#' @importFrom rlang .data
#' @examples
#' df <- flownet_geometry()
#' dp <- solve_flownet(df)
#' @export

solve_flownet <- function(df) {
  ## SOLVE FOR HEAD
  #get sparse matrix elements and left-hand side for boundary conditions
  Mbc <- findiff_sparse_entries_bc(df)
  #get sparse matrix elements for Poissons equation for each real node
  M2x <- df$dom %>%
    dplyr::filter(.data$dry == FALSE) %>%
    dplyr::select(.data$nx, .data$ny, .data$kx, .data$hx, .data$i0) %>%
    purrr::pmap_dfr(
      function(nx, ny, kx, hx, i0) findiff_sparse_entries(nx, ny, order = 2, h = hx, direction = "x", multiplier = kx, i0_row = i0, i0_col = i0)
    )
  M2y <- df$dom %>%
    dplyr::filter(.data$dry == FALSE) %>%
    dplyr::select(.data$nx, .data$ny, .data$ky, .data$hy, .data$i0) %>%
    purrr::pmap_dfr(
      function(nx, ny, ky, hy, i0) findiff_sparse_entries(nx, ny, order = 2, h = hy, direction = "y", multiplier = ky, i0_row = i0, i0_col = i0)
    )
  #create sparse matrix
  mat <- Matrix::sparseMatrix(
    i = c(Mbc$mat$row, M2x$row, M2y$row),
    j = c(Mbc$mat$col, M2x$col, M2y$col),
    x = c(Mbc$mat$val, M2x$val, M2y$val),
    dims = rep(sum(df$dom$n), 2)
  )
  #bind together, and solve system
  h = as.vector(Matrix::solve(mat, Mbc$lhs))

  ## CALCULATE FLOW FROM HEAD
  #flow in x-direction
  M1qx <- df$dom %>%
    dplyr::filter(.data$dry == FALSE) %>%
    dplyr::select(.data$nx, .data$ny, .data$kx, .data$hx, .data$i0) %>%
    purrr::pmap_dfr(
      function(nx, ny, kx, hx, i0) findiff_sparse_entries(nx, ny, order = 1, h = hx, direction = "x", multiplier = -kx, i0_row = i0, i0_col = i0)
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
    dplyr::filter(.data$dry == FALSE) %>%
    dplyr::select(.data$nx, .data$ny, .data$ky, .data$hy, .data$i0) %>%
    purrr::pmap_dfr(
      function(nx, ny, ky, hy, i0) findiff_sparse_entries(nx, ny, order = 1, h = hy, direction = "y", multiplier = -ky, i0_row = i0, i0_col = i0)
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
  i_real <- purrr::pmap(df$dom[df$dom$dry == FALSE, ], index_real) %>% unlist()
  #generate positions for all solution object for all nodes
  dp <- df$dom %>%
    dplyr::filter(.data$dry == FALSE) %>%
    dplyr::summarise(
      nodal_coordinates_real(
        .data$nx,
        .data$ny,
        x0 = .data$x0,
        x1 = .data$x1,
        y0 = .data$y0,
        y1 = .data$y1,
        id = .data$id
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
    dplyr::filter(.data$dry == FALSE) %>%
    dplyr::select(.data$nx, .data$ny, .data$kx, .data$hx, .data$i0_real) %>%
    purrr::pmap_dfr(
      function(nx, ny, kx, hx, i0_real) findiff_sparse_entries(nx, ny, order = 1, h = hx, direction = "x", i0_row = i0_real, i0_col = i0_real, real_only = TRUE)
    )
  #x-direction elements - first order differentiation - real nodes only
  M1y_el <- df$dom %>%
    dplyr::filter(.data$dry == FALSE) %>%
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
#' @return the optimal number of equipotential intervals `Nd`
#' @importFrom magrittr `%>%`
#' @importFrom rlang .data
#' @examples
#' df <- flownet_geometry()
#' dp <- solve_flownet(df)
#' amount_equipotentialintervals(df, dp, Nf = 5)
#' @export

amount_equipotentialintervals <- function(df, dp, Nf = 5) {
  #total head difference
  dh <- max(dp$h) - min(dp$h)
  #total flow per second
  Q <- max(dp$psi) - min(dp$psi)
  #k_avg
  k_avg <- df$dom %>%
    dplyr::filter(.data$dry == FALSE) %>%
    dplyr::summarize(k = mean(sqrt(.data$kx*.data$ky))) %>%
    dplyr::pull(.data$k)
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
#' df <- flownet_geometry()
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
  Nd_levels <- seq(min(dp$h), max(dp$h), l = (Nd + 1))[2:Nd]
  Nf_levels <- seq(min(dp$psi), max(dp$psi), l = (Nf + 1))
  #get head edges
  dedge <- df$bc %>%
    dplyr::filter(.data$type == "h") %>%
    dplyr::left_join(df$dom, by = "id") %>%
    dplyr::mutate(
      x = ifelse(.data$edge == 3, .data$x1, .data$x0),
      xend = ifelse(.data$edge == 1, .data$x0, .data$x1),
      y = ifelse(.data$edge == 2, .data$y1, .data$y0),
      yend = ifelse(.data$edge == 4, .data$y0, .data$y1)
    )
  #plot rectangular soil elements
  plt <- ggplot2::ggplot() +
    ggplot2::geom_rect(
      data = df$dom %>% dplyr::filter(.data$dry == FALSE),
      ggplot2::aes(xmin = .data$x0, xmax = .data$x1, ymin = .data$y0, ymax = .data$y1),
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
#' df <- flownet_geometry()
#' ggplot_flownet_dam(df)
#'
#' #anisotropic permeability - with 6 flow channels
#' df <- flownet_geometry(kx = 4e-6, ky = 1e-6)
#' ggplot_flownet_dam(df, Nf = 6)
#'
#' #deep excavation example
#' df <- flownet_geometry(
#'   x = c(0, 5, 5.25, 30),
#'   y = c(0, 5, 7.5, 15, 20),
#'   ix = c(1, 1, 2, 3, 3, 3, 3),
#'   iy = c(2, 1, 1, 1, 2, 3, 4),
#'   bc_id = c(1, 6),
#'   bc_edge = c(2, 2),
#'   bc_type = c("h", "h"),
#'   bc_value = c(7.5, 15),
#'   dry = c(rep(FALSE, 6), TRUE)
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
  Nd_levels <- seq(min(dp$h), max(dp$h), l = (Nd + 1))[2:Nd]
  Nf_levels <- seq(min(dp$psi), max(dp$psi), l = (Nf + 1))
  #impermeable thickness
  if (length(thickness_impermeable) == 1) {
    thickness_impermeable <- rep(thickness_impermeable, 4)
  }
  #get head edges
  dedge <- df$bc %>%
    dplyr::filter(.data$type == "h") %>%
    dplyr::left_join(df$dom, by = "id") %>%
    dplyr::mutate(
      x = ifelse(.data$edge == 3, .data$x1, .data$x0),
      xend = ifelse(.data$edge == 1, .data$x0, .data$x1),
      y = ifelse(.data$edge == 2, .data$y1, .data$y0),
      yend = ifelse(.data$edge == 4, .data$y0, .data$y1)
    )
  #plot limits
  if (is.na(xlim[1])) {
    xlim[1] <- min(df$dom$x0) - thickness_impermeable[1]
  }
  if (is.na(xlim[2])) {
    xlim[2] <- max(df$dom$x1) + thickness_impermeable[3]
  }
  if (is.na(ylim[1])) {
    ylim[1] <- min(df$dom$y0) - thickness_impermeable[4]
  }
  if (is.na(ylim[2])) {
    ylim[2] <- max(c(df$dom$y1, dedge$yend)) + thickness_impermeable[2]
  }
  #initiate plot
  plt <- ggplot2::ggplot() +
    ggplot2::coord_equal(ratio = 1, xlim = xlim, ylim = ylim, expand = FALSE)
  #impermeable background
  plt <- plt +
    ggplot2::annotate(
      "polygon",
      x = c(rep(xlim[1], 2), rep(xlim[2], 2)),
      y = c(ylim[1], rep(ylim[2], 2), ylim[1]),
      color = NA,
      fill = fill_impermeable
    )
  #plot rectangles for ponding water
  plt <- plt +
    ggplot2::geom_rect(
      data = dedge,
      ggplot2::aes(xmin = .data$x, xmax = .data$xend, ymin = .data$y0, ymax = .data$value),
      fill = fill_water,
      color = NA
    )
  #plot air
  plt <- plt +
    ggplot2::geom_rect(
      data = dedge,
      ggplot2::aes(xmin = .data$x, xmax = .data$xend, ymin = .data$value, ymax = ylim[2]),
      fill = fill_air,
      colour = NA
    )
  #plot rectangles for dry soil
  plt <- plt +
    ggplot2::geom_rect(
      data = dplyr::filter(df$dom, .data$dry == TRUE),
      ggplot2::aes(xmin = .data$x0, xmax = .data$x1, ymin = .data$y0, ymax = .data$y1),
      fill = fill_soildry,
      colour = NA
    )
  #plot rectangles for saturated soil
  plt <- plt +
    ggplot2::geom_rect(
      data = dplyr::filter(df$dom, .data$dry == FALSE),
      ggplot2::aes(xmin = .data$x0, xmax = .data$x1, ymin = .data$y0, ymax = .data$y1),
      fill = fill_soilwet,
      colour = NA
    )
  #plot soil surface lines
  dsurf <- dplyr::left_join(
    dedge %>%
      dplyr::select(.data$ix, .data$x, .data$xend, .data$y, .data$yend),
    df$dom %>%
      dplyr::filter(.data$dry == TRUE) %>%
      dplyr::select(.data$ix, .data$y1),
    by = "ix"
  ) %>%
    dplyr::mutate(
      y = pmax(.data$y, .data$y1, na.rm = TRUE),
      yend = pmax(.data$yend, .data$y1, na.rm = TRUE)
    )
  plt <- plt +
    ggplot2::geom_segment(
      data = dsurf,
      ggplot2::aes(x = .data$x, xend = .data$xend, y = .data$y, yend = .data$yend),
      color = colour_soil,
      size = line_width
    )
  #plot water table lines
  plt <- plt +
    ggplot2::geom_segment(
      data = dedge,
      ggplot2::aes(x = .data$x, xend = .data$xend, y = .data$value, yend = .data$value),
      color = colour_water,
      linetype = line_water,
      size = line_width
    )
  #plot soil surface marker
  if (!is.na(scale_soilmarker)) {
    for (i in 1:nrow(dsurf)) {
      if ((dsurf$xend[i] > xlim[1]) & (dsurf$x[i] <= xlim[2])) {
        plt <- ggplot_add_soilmarker(
          plt,
          1/3*max(xlim[1], dsurf$x[i]) + 2/3*min(xlim[2], dsurf$xend[i]),
          dsurf$y[i],
          scale = scale_soilmarker*diff(ylim)
        )
      }
    }
  }
  #plot water table markers
  if (!is.na(scale_watermarker)) {
    for (i in 1:nrow(dedge)) {
      if ((dedge$xend[i] > xlim[1]) & (dedge$x[i] <= xlim[2])) {
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
  #return
  return(plt)
}


#' ggplot 2D map of hydraulic heads and pressures
#'
#' @description
#' Plots a 2D map of the pressure or head distributions in a flow net problem.
#' The plot shows both a heat map and contour lines.
#'
#' @inheritParams ggplot_flownet_dam
#' @param type parameter to plot. `type = "h"` for hydraulic heads,
#'   `type = "hz"` for elevation head, `type = `hb` for pressure head
#'   or `type = "u"` for pore water pressures
#' @param binwidth contour line distance
#' @param label_size contour line text label size
#' @param label_n number of times contour label appears on each line
#' @param gamma_w unit weight of water, in kN/m3 (required for calculating
#'   pore water pressures)
#' @param palette RColorBrewer color map for shading
#' @param palette_direction RColorBrewer color map direction (`1` or `-1`)
#' @param colour_contour colour of contour lines
#' @return a ggplot object
#' @examples
#' #generate flow net geometry
#' df <- flownet_geometry()
#'
#' #solve flow net problem
#' dp <- solve_flownet(df)
#'
#' #plot heads
#' ggplot_hydraulichead_dam(df, dp = dp, type = "h", binwidth = 0.25)
#'
#' #plot elevation head
#' ggplot_hydraulichead_dam(df, dp = dp, type = "hz", binwidth = 2.5)
#'
#' #plot pressure head
#' ggplot_hydraulichead_dam(df, dp = dp, type = "hb", binwidth = 2.5)
#'
#' #plot pore water pressure
#' ggplot_hydraulichead_dam(df, dp = dp, type = "u", binwidth = 25)
#' @export

ggplot_hydraulichead_dam <- function(
  df,
  dp = NULL,
  type = "h",  #h, hz, hb, u
  binwidth = 0.25,
  gamma_w = 10,
  label_size = 3,
  label_n = 1,
  xlab = "x [m]",
  ylab = "z [m]",
  palette = "Blues",  #YlGnBu
  palette_direction = 1,
  fill_soilwet = "white",
  fill_soildry = "white",
  colour_soil = "#65571d",
  colour_contour = "black",
  fill_water = "white",
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
  #impermeable thickness
  if (length(thickness_impermeable) == 1) {
    thickness_impermeable <- rep(thickness_impermeable, 4)
  }
  #get head edges
  dedge <- df$bc %>%
    dplyr::filter(.data$type == "h") %>%
    dplyr::left_join(df$dom, by = "id") %>%
    dplyr::mutate(
      x = ifelse(.data$edge == 3, .data$x1, .data$x0),
      xend = ifelse(.data$edge == 1, .data$x0, .data$x1),
      y = ifelse(.data$edge == 2, .data$y1, .data$y0),
      yend = ifelse(.data$edge == 4, .data$y0, .data$y1)
    )
  #plot limits
  if (is.na(xlim[1])) {
    xlim[1] <- min(df$dom$x0) - thickness_impermeable[1]
  }
  if (is.na(xlim[2])) {
    xlim[2] <- max(df$dom$x1) + thickness_impermeable[3]
  }
  if (is.na(ylim[1])) {
    ylim[1] <- min(df$dom$y0) - thickness_impermeable[4]
  }
  if (is.na(ylim[2])) {
    ylim[2] <- max(c(df$dom$y1, dedge$yend)) + thickness_impermeable[2]
  }
  #initiate plot
  plt <- ggplot2::ggplot() +
    ggplot2::coord_equal(ratio = 1, xlim = xlim, ylim = ylim, expand = FALSE)
  #impermeable background
  plt <- plt +
    ggplot2::annotate(
      "polygon",
      x = c(rep(xlim[1], 2), rep(xlim[2], 2)),
      y = c(ylim[1], rep(ylim[2], 2), ylim[1]),
      color = NA,
      fill = fill_impermeable
    )
  #plot rectangles for ponding water
  plt <- plt +
    ggplot2::geom_rect(
      data = dedge,
      ggplot2::aes(xmin = .data$x, xmax = .data$xend, ymin = .data$y0, ymax = .data$value),
      fill = fill_water,
      color = NA
    )
  #plot air
  plt <- plt +
    ggplot2::geom_rect(
      data = dedge,
      ggplot2::aes(xmin = .data$x, xmax = .data$xend, ymin = .data$value, ymax = ylim[2]),
      fill = fill_air,
      colour = NA
    )
  #plot rectangles for dry soil
  plt <- plt +
    ggplot2::geom_rect(
      data = dplyr::filter(df$dom, .data$dry == TRUE),
      ggplot2::aes(xmin = .data$x0, xmax = .data$x1, ymin = .data$y0, ymax = .data$y1),
      fill = fill_soildry,
      colour = NA
    )
  #plot rectangles for saturated soil
  plt <- plt +
    ggplot2::geom_rect(
      data = dplyr::filter(df$dom, .data$dry == FALSE),
      ggplot2::aes(xmin = .data$x0, xmax = .data$x1, ymin = .data$y0, ymax = .data$y1),
      fill = fill_soilwet,
      colour = NA
    )
  #plot hydraulic head profiles
  if (type == "h") {
    dp$val <- dp$h
    label <- expression(h~"[m]")
  } else if (type == "hz") {
    dp$val <- dp$y
    label <- expression(h[z]~"[m]")
  } else if (type == "hb") {
    dp$val <- dp$h - dp$y
    label <- expression(h[b]~"[m]")
  } else if (type == "u") {
    dp$val <- (dp$h - dp$y)*gamma_w
    label <- expression(u~"[kPa]")
  }
  plt <- plt +
    ggplot2::stat_summary_2d(
      data = dp,
      ggplot2::aes(x = .data$x, y = .data$y, z = .data$val),
      binwidth = min(c(df$dom$hx, df$dom$hy))
    ) +
    ggplot2::geom_contour(
      data = dp,
      ggplot2::aes(x = .data$x, y = .data$y, z = .data$val),
      binwidth = binwidth,
      color = colour_contour,
      size = line_width
    ) +
    ggplot2::scale_fill_distiller(
      name = label,
      palette = palette,
      direction = palette_direction
    ) +
    metR::geom_text_contour(
      data = dp,
      binwidth = binwidth,
      ggplot2::aes(x = .data$x, y = .data$y, z = .data$val),
      label.placement = metR::label_placement_n(label_n),
      size = label_size
    )
  #plot soil surface lines
  dsurf <- dplyr::left_join(
    dedge %>%
      dplyr::select(.data$ix, .data$x, .data$xend, .data$y, .data$yend),
    df$dom %>%
      dplyr::filter(.data$dry == TRUE) %>%
      dplyr::select(.data$ix, .data$y1),
    by = "ix"
  ) %>%
    dplyr::mutate(
      y = pmax(.data$y, .data$y1, na.rm = TRUE),
      yend = pmax(.data$yend, .data$y1, na.rm = TRUE)
    )
  plt <- plt +
    ggplot2::geom_segment(
      data = dsurf,
      ggplot2::aes(x = .data$x, xend = .data$xend, y = .data$y, yend = .data$yend),
      color = colour_soil,
      size = line_width
    )
  #plot water table lines
  plt <- plt +
    ggplot2::geom_segment(
      data = dedge,
      ggplot2::aes(x = .data$x, xend = .data$xend, y = .data$value, yend = .data$value),
      color = colour_water,
      linetype = line_water,
      size = line_width
    )
  #plot soil surface marker
  if (!is.na(scale_soilmarker)) {
    for (i in 1:nrow(dsurf)) {
      if ((dsurf$xend[i] > xlim[1]) & (dsurf$x[i] <= xlim[2])) {
        plt <- ggplot_add_soilmarker(
          plt,
          1/3*max(xlim[1], dsurf$x[i]) + 2/3*min(xlim[2], dsurf$xend[i]),
          dsurf$y[i],
          scale = scale_soilmarker*diff(ylim)
        )
      }
    }
  }
  #plot water table markers
  if (!is.na(scale_watermarker)) {
    for (i in 1:nrow(dedge)) {
      if ((dedge$xend[i] > xlim[1]) & (dedge$x[i] <= xlim[2])) {
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
  #return
  return(plt)
}

