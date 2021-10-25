#' Define flow net problem using gridded quadrilateral domain system
#'
#' @description
#' Defines a geometry for a flow net problem. quadrilateral domains are defined
#' using a grid system. For each domain, horizontal and vertical permeabilities
#' may be defined seperately. All quadrilateral domains should be strictly
#' convex in shape.
#'
#' Domains that are next to each other in the grid are automatically connected.
#' All boundary conditions are assumed impermeable, unless otherwise defined
#' using the inputs `bc_id`, `bc_edge`, `bc_type` and `bc_value`.
#'
#' @param x matrix of x-coordinates of corner points in domain
#' @param y matrix of y-coordinates of corner points in domain
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
#' flownet_geometry_quadrilateral()
#' @export

flownet_geometry_quadrilateral <- function(
  x = matrix(c(1,3,6, 1,3,5, 1,3,5), ncol = 3, byrow = FALSE),
  y = matrix(c(0,0,0, 1,1,1.5, 2.5,2,2), ncol = 3, byrow = FALSE),
  ix = c(1, 1, 2),
  iy = c(2, 1, 1),
  kx = 1e-6,
  ky = 1e-6,
  grid_size = 0.05,
  bc_domain = c(1, 3),
  bc_edge = c(2, 3),
  bc_type = c("h", "h"),
  bc_value = c(10, 5),
  node_min = 3
){
  #estimate max grid size for each row/column
  Lx <- apply(sqrt(diff(x)^2 + diff(t(x))^2), 1, max)
  Ly <- apply(sqrt(diff(t(y))^2 + diff(y)^2), 1, max)
  #create dataframe
  dom <- tibble::tibble(ix = ix, iy = iy, kx = kx, ky = ky) %>%
    dplyr::mutate(
      #domain identifier
      domain = seq(dplyr::n()),
      #positions of corner points, numbered clockwise
      x1 = x[matrix(c(ix, iy), ncol = 2)],
      x2 = x[matrix(c(ix, iy + 1), ncol = 2)],
      x3 = x[matrix(c(ix + 1, iy + 1), ncol = 2)],
      x4 = x[matrix(c(ix + 1, iy), ncol = 2)],
      y1 = y[matrix(c(ix, iy), ncol = 2)],
      y2 = y[matrix(c(ix, iy + 1), ncol = 2)],
      y3 = y[matrix(c(ix + 1, iy + 1), ncol = 2)],
      y4 = y[matrix(c(ix + 1, iy), ncol = 2)],
      #number of nodes in each direction
      nx = pmax(node_min, (1 + ceiling(Lx/grid_size))[.data$ix]),
      ny = pmax(node_min, (1 + ceiling(Ly/grid_size))[.data$iy]),
      #total nodes in each domain and starting number
      n = nodes_total(.data$nx, .data$ny),
      i0 = cumsum(c(0, utils::head(.data$n, -1)))
    )
  #generate tibble with all default boundary conditions (impermeable boundaries)
  bc <- tibble::tibble(
    domain = rep(seq(nrow(dom)), each = 4),
    edge = rep(seq(4), nrow(dom)),
    type = "q",
    value = 0
  )
  #add connecting boundary conditions
  for (i in 1:nrow(dom)) {
    #other domain connected at the right (positive x-side, edge 3)
    i_next <- which((dom$ix == (dom$ix[i] + 1)) & (dom$iy == dom$iy[i]))
    if (length(i_next) == 1) {
      bc$edge[4*(i - 1) + 3] <- 3
      bc$type[4*(i - 1) + 3] <- "c"
      bc$value[4*(i - 1) + 3] <- i_next
      bc$edge[4*(i_next - 1) + 1] <- 1
      bc$type[4*(i_next - 1) + 1] <- "c"
      bc$value[4*(i_next - 1) + 1] <- i
    }
    #other domain connected at the bottom (positive y-side, edge 2)
    i_next <- which((dom$ix == dom$ix[i]) & (dom$iy == (dom$iy[i] + 1)))
    if (length(i_next) == 1) {
      bc$edge[4*(i - 1) + 2] <- 2
      bc$type[4*(i - 1) + 2] <- "c"
      bc$value[4*(i - 1) + 2] <- i_next
      bc$edge[4*(i_next - 1) + 4] <- 4
      bc$type[4*(i_next - 1) + 4] <- "c"
      bc$value[4*(i_next - 1) + 4] <- i
    }
  }
  #add user-defined boundary conditions for head or flow
  bc$domain[4*(bc_domain - 1) + bc_edge] <- bc_domain
  bc$edge[4*(bc_domain - 1) + bc_edge] <- bc_edge
  bc$type[4*(bc_domain - 1) + bc_edge] <- bc_type
  bc$value[4*(bc_domain - 1) + bc_edge] <- bc_value
  #return
  return(list(dom = dom, bc = bc))
}


#' Get properties of real edge nodes in quadrilateral domain
#'
#' @description
#' Function returns all real nodes on the edges of a quadrilateral domain,
#' as well as some key properties such as positions and a normal vector
#'
#' @param nx,ny number of real nodes in a and b-directions
#' @param x1,x2,x3,x4,y1,y2,y3,y4 x and y positions of the corner points of the
#'   quadrilateral. Corners are numbered in clockwise order.
#' @param i0 node offset
#' @param domain domain identifier
#' @param ... additional arguments to pass
#' @return tibble with all nodes on edges. Contains fields for the domain
#'   identifier (`domain`), edge (`edge`, numbered clockwise), node index
#'   (`i`), row and column indices (`ix`, `iy`), the unit vector normal
#'   to the edge (positing towards the domain) (`vx`, `vy`) and first
#'   order derivatives of position (`x_a`, `x_b`, `y_a`, `y_b`)
#' @examples
#' nx <- 20
#' ny <- 10
#' x1 <- 0
#' x2 <- 0
#' x3 <- 1
#' x4 <- 1
#' y1 <- 0
#' y2 <- 1
#' y3 <- 2
#' y4 <- 0
#' edges_quadrilateral(nx, ny, x1, x2, x3, x4, y1, y2, y3, y4)
#' @export

# Get nodes on edges
#get nodes on edges + normal vectors
edges_quadrilateral <- function(
  nx, ny,
  x1, x2, x3, x4, y1, y2, y3, y4,
  i0 = 0,
  domain = 1,
  ...
){
  #get details for each node on an edge
  return(tibble::tibble(
    #get properties for each node on each edge
    domain = domain,
    edge = c(rep(1, ny), rep(2, nx), rep(3, ny), rep(4, nx)),
    ix = c(rep(1, ny), seq(nx), rep(nx, ny), seq(nx)),
    iy = c(seq(ny), rep(ny, nx), seq(ny), rep(1, nx)),
    a = c(rep(0, ny), seq(0, 1, l = nx), rep(1, ny), seq(0, 1, l = nx)),
    b = c(seq(0, 1, l = ny), rep(1, nx), seq(0, 1, l = ny), rep(0, nx))
  ) %>%
    dplyr::mutate(
      #get node index
      i = index_grid2vector(.data$ix, .data$iy, nx, ny, i0 = i0),
      #get positions and normal vector pointing towards the element
      position_xy_quadrilateral(.data$a, .data$b, x1, x2, x3, x4, y1, y2, y3, y4),
      vx_temp = ifelse(
        (.data$edge %in% c(1, 3)),
        -.data$y_b/sqrt(.data$x_b^2 + .data$y_b^2),
        .data$y_a/sqrt(.data$x_a^2 + .data$y_a^2)
      ),
      vy_temp = ifelse(
        (.data$edge %in% c(1, 3)),
        .data$x_b/sqrt(.data$x_b^2 + .data$y_b^2),
        -.data$x_a/sqrt(.data$x_a^2 + .data$y_a^2)
      ),
      vx = ifelse(.data$edge %in% c(1, 4), -.data$vx_temp, .data$vx_temp),
      vy = ifelse(.data$edge %in% c(1, 4), -.data$vy_temp, .data$vy_temp)
    ) %>%
    dplyr::select(.data$domain, .data$edge, .data$i,
                  .data$ix, .data$iy,
                  .data$vx, .data$vy,
                  .data$x_a, .data$x_b, .data$y_a, .data$y_b)
  )
}


#' Get positions of all real nodes in quadrilateral domains
#'
#' @description
#' Generate positions of all real nodes in a quadrilateral domain, both in
#' terms of normalised positions a and b, and actual positions x and y
#'
#' @param nx,ny number of real nodes in x and y directions
#' @param x1,x2,x3,x4,y1,y2,y3,y4 x and y positions of the corner points of the
#'   quadrilateral. Corners are numbered in clockwise order.
#' @param domain domain identifier
#' @param i0 optional index numbering offset
#' @param ... additional arguments to pass
#' @return a tibble with positions, e.g. `a`, `b`, `x`, `y` as well as
#'   1st and 2nd order derivatives of positions (e.g. `x_a`, `y_bb`)
#' @examples
#' nx <- 4
#' ny <- 3
#' x1 <- 0
#' x2 <- 0
#' x3 <- 1
#' x4 <- 1
#' y1 <- 0
#' y2 <- 1
#' y3 <- 2
#' y4 <- 0
#' positions_real_quadrilateral(nx, ny, x1, x2, x3, x4, y1, y2, y3, y4)
#' @export

positions_real_quadrilateral <- function(
  nx, ny, x1, x2, x3, x4, y1, y2, y3, y4, domain = 1, i0 = 0, ...
){
  return(
    tidyr::expand_grid(
      b = seq(0, 1, l = ny),
      a = seq(0, 1, l = nx)
    ) %>%
      dplyr::mutate(
        domain = domain,
        x1 = x1, x2 = x2, x3 = x3, x4 = x4,
        y1 = y1, y2 = y2, y3 = y3, y4 = y4,
        i0 = i0,
        i = index_real(nx, ny, i0 = i0),
        position_xy_quadrilateral(
          .data$a, .data$b,
          .data$x1, .data$x2, .data$x3, .data$x4,
          .data$y1, .data$y2, .data$y3, .data$y4
        )
      )
  )
}


#' Get x,y positions and derivatives of real nodes in quadrilateral grid
#'
#' @description
#' Function calculates the x and y positions of nodes in a quadrilateral grid
#' based on their a and b positions. The derivatives of x and y with respect
#' to a and b are also returned
#'
#' @param a,b arrays with quadrilateral positions
#' @param x1,x2,x3,x4,y1,y2,y3,y4 x and y positions of the corner points of the
#'   quadrilateral. Corners are numbered in clockwise order.
#' @param ... additional arguments
#' @return a tibble with fields for position `x` and `y`, as well as 1st order
#'   derivatives (`x_a`, `x_b`, `y_a`, `y_b`) and second order derivatives
#'   (`x_aa`, `x_bb`, `x_ab`, `y_aa`, `y_bb`, `y_ab`)
#' @examples
#' a <- c(0, 0.5, 1)
#' b <- c(0, 0.5, 1)
#' x1 <- 0
#' x2 <- 0
#' x3 <- 1
#' x4 <- 1
#' y1 <- 0
#' y2 <- 1
#' y3 <- 2
#' y4 <- 0
#' position_xy_quadrilateral(a, b, x1, x2, x3, x4, y1, y2, y3, y4)
#' @export

position_xy_quadrilateral <- function(a, b, x1, x2, x3, x4, y1, y2, y3, y4, ...) {
  return(tibble::tibble(
    #positions
    x = x1*(1 - a)*(1 - b) + x4*a*(1 - b) + x2*(1 - a)*b + x3*a*b,
    y = y1*(1 - a)*(1 - b) + y4*a*(1 - b) + y2*(1 - a)*b + y3*a*b,
    #derivatives of x and y with respect to a and b
    x_a = -x1*(1 - b) + x4*(1 - b) - x2*b + x3*b,
    y_a = -y1*(1 - b) + y4*(1 - b) - y2*b + y3*b,
    x_b = -x1*(1 - a) - x4*a + x2*(1 - a) + x3*a,
    y_b = -y1*(1 - a) - y4*a + y2*(1 - a) + y3*a,
    #second derivatives
    x_aa = 0,
    x_ab = x1 - x4 - x2 + x3,
    x_bb = 0,
    y_aa = 0,
    y_ab = y1 - y4 - y2 + y3,
    y_bb = 0
  ))
}


#' Get sparse Poisson's matrix elements for quadrilateral domain
#'
#' @description
#' Function generates the sparse matrix elements for the Poisson's
#' equation in all real nodes, for a single quadrilateral domain
#' @param nx,ny number of real nodes in a and b-directions
#' @param kx,ky horizontal and vertical permeability in domain
#' @param x1,x2,x3,x4,y1,y2,y3,y4 x and y positions of the corner points of the
#'   quadrilateral. Corners are numbered in clockwise order.
#' @param i0 index offset for nodes
#' @param ... extra arguments to pass
#' @return a tibble with rows (`row`), columns (`col`) and values (`val`)
#'   of non-zero elements in the sparse matrix
#' @examples
#' nx <- 20
#' ny <- 10
#' kx <- 1
#' ky <- 1
#' x1 <- 0
#' x2 <- 0
#' x4 <- 1
#' x3 <- 1
#' y1 <- 0
#' y2 <- 1
#' y3 <- 2
#' y4 <- 0
#' poissons_eq_quadrilateral(nx, ny, kx, ky, x1, x2, x3, x4, y1, y2, y3, y4)
#' @export

poissons_eq_quadrilateral <- function(
  nx, ny,
  kx, ky,
  x1, x2, x3, x4, y1, y2, y3, y4,
  i0 = 0,
  ...
) {
  #generate tibble
  #get a and b for all real nodes
  a <- rep(seq(0, 1, l = nx), ny)
  b <- rep(seq(0, 1, l = ny), each = nx)
  #get positions etc
  dxy <- position_xy_quadrilateral(a, b, x1, x2, x3, x4, y1, y2, y3, y4)
  # matrix elements
  #   [h_a, h_b] = A * [h_x, h_y]
  #   [h_aa, h_bb, h_ab] = B * [h_xx, h_yy, h_xy] + C * [h_x, h_y]
  # so:
  #   [h_xx, h_yy, h_xy] = B^-1 * [h_aa, h_bb, h_ab] - B^-1 * C * A^-1 * [h_a, h_b]
  # calculate C * A^-1
  temp <- with(dxy, x_a*y_b - x_b*y_a)
  CAinv_11 <- with(dxy, (x_aa*y_b - x_b*y_aa)/temp)
  CAinv_12 <- with(dxy, (x_a*y_aa - x_aa*y_a)/temp)
  CAinv_21 <- with(dxy, -(x_b*y_bb - x_bb*y_b)/temp)
  CAinv_22 <- with(dxy, (x_a*y_bb - x_bb*y_a)/temp)
  CAinv_31 <- with(dxy, (x_ab*y_b - x_b*y_ab)/temp)
  CAinv_32 <- with(dxy, (x_a*y_ab - x_ab*y_a)/temp)
  #calculate D = B^-1 (only calculate first two rows)
  temp2 <- with(dxy, x_a^2*y_b^2 - 2*x_a*x_b*y_a*y_b + x_b^2*y_a^2)
  D11 <- with(dxy, y_b^2/temp2)
  D12 <- with(dxy, y_a^2/temp2)
  D13 <- with(dxy, -(2*y_a*y_b)/temp2)
  D21 <- with(dxy, x_b^2/temp2)
  D22 <- with(dxy, x_a^2/temp2)
  D23 <- with(dxy, -(2*x_a*x_b)/temp2)
  #calculate E = B^-1 * (A * C^-1)
  E11 <- D11*CAinv_11 + D12*CAinv_21 + D13*CAinv_31
  E12 <- D11*CAinv_12 + D12*CAinv_22 + D13*CAinv_32
  E21 <- D21*CAinv_11 + D22*CAinv_21 + D23*CAinv_31
  E22 <- D21*CAinv_12 + D22*CAinv_22 + D23*CAinv_32
  #get finite difference elements
  el_aa <- findiff_sparse_elements(nx, ny, "xx", i0 = i0, multiplier = kx*D11 + ky*D21)
  el_bb <- findiff_sparse_elements(nx, ny, "yy", i0 = i0, multiplier = kx*D12 + ky*D22)
  el_ab <- findiff_sparse_elements(nx, ny, "xy", i0 = i0, multiplier = kx*D13 + ky*D23)
  el_a <- findiff_sparse_elements(nx, ny, "x", i0 = i0, multiplier = -kx*E11 - ky*E21)
  el_b <- findiff_sparse_elements(nx, ny, "y", i0 = i0, multiplier = -kx*E12 - ky*E22)
  #return combined list
  return(tibble::tibble(
    row = c(el_aa$row, el_bb$row, el_ab$row, el_a$row, el_b$row),
    col = c(el_aa$col, el_bb$col, el_ab$col, el_a$col, el_b$col),
    val = c(el_aa$val, el_bb$val, el_ab$val, el_a$val, el_b$val)
  ))
}


#' Obtain sparse linear system for solving hydraulic head (quadrilateral)
#'
#' @description
#' Take the geometry and boundary conditions of a flownet problem consisting
#' of connected quadrilateral domains, and returns the matrix and vector of
#' the linear finite difference problem
#'
#' @param df list with tibbles describing the problem, see function
#'   `flownet_geometry_quadrilateral()` for more information
#' @return list with two elements. `mat` contains a tibble with all sparse
#'   matrix elements (`row`, `col`, `val`), and `lhs` contains a vector
#'   with all left-hand side of equation values
#' @importFrom magrittr `%>%`
#' @importFrom rlang .data
#' @examples
#' df <- flownet_geometry_quadrilateral()
#' matrix_sparse_head_quadrilateral(df)
#' @export

matrix_sparse_head_quadrilateral <- function(df) {
  #vector with solutions
  lhs <- rep(0, sum(df$dom$n))

  #matrix elements for each domain
  M_el <- purrr::pmap_dfr(df$dom, poissons_eq_quadrilateral)

  #get all data for all edges + normal vectors
  dedge <- purrr::pmap_dfr(df$dom, edges_quadrilateral) %>%
    #merge nx, ny and i0, and boundary conditions
    dplyr::left_join(
      df$dom %>% dplyr::select(.data$domain, .data$nx, .data$ny, .data$i0, .data$kx, .data$ky),
      by = "domain"
    ) %>%
    dplyr::left_join(
      df$bc,
      by = c("domain", "edge")
    ) %>%
    dplyr::mutate(
      dix = (.data$edge %% 2)*(.data$edge - 2),
      diy = ((.data$edge - 1) %% 2)*(3 - .data$edge),
      i_ghost = index_offset(.data$i, .data$nx, .data$ny, dix = .data$dix, diy = .data$diy, i0 = .data$i0)
    )

  #fixed heads
  dedge_h <- dplyr::filter(dedge, .data$type == "h")
  lhs[dedge_h$i_ghost] <- dedge_h$value
  M_bc_h <- dedge_h %>%
    dplyr::select(.data$i_ghost, .data$i) %>%
    dplyr::rename(row = .data$i_ghost, col = .data$i) %>%
    dplyr::mutate(val = 1)

  #connections - elements
  dedge_c <- dedge %>%
    dplyr::filter(.data$type == "c") %>%
    dplyr::left_join(
      df$dom %>%
        dplyr::select(.data$domain, .data$nx, .data$ny, .data$i0) %>%
        dplyr::rename("domain_2" = "domain", "nx_2" = "nx", "ny_2" = "ny", "i0_2" = "i0"),
      by = c("value" = "domain_2")
    ) %>%
    dplyr::mutate(
      edge_2 = 1 + (1 + .data$edge)%%4,
      ix_2 = ifelse(.data$edge_2 %in% c(2, 4), .data$ix, ifelse(.data$edge_2 == 1, 1, .data$nx_2)),
      iy_2 = ifelse(.data$edge_2 %in% c(1, 3), .data$iy, ifelse(.data$edge_2 == 4, 1, .data$ny_2)),
      i_2 = index_grid2vector(.data$ix_2, .data$iy_2, .data$nx_2, .data$ny_2, i0 = .data$i0_2)
    )

  #connections - same head (edges 2 and 3 - positive)
  dedge_ch <- dplyr::filter(dedge_c, .data$edge %in% c(2, 3))
  M_bc_ch <- tibble::tibble(
    row = rep(dedge_ch$i_ghost, 2),
    col = c(dedge_ch$i, dedge_ch$i_2),
    val = rep(c(-1, 1), each = nrow(dedge_ch))
  )

  #get elements and multipliers for flow
  dedge_q_c <- dedge %>%
    dplyr::filter(.data$type %in% c("q", "c")) %>%
    dplyr::mutate(
      #multipliers: q = mult_ha * dh/da + mult_hb * dh/db
      mult_ha = (-.data$vx*.data$kx*.data$y_b + .data$vy*.data$ky*.data$x_b)/(.data$x_a*.data$y_b - .data$x_b*.data$y_a),
      mult_hb = (.data$vx*.data$kx*.data$y_a - .data$vy*.data$ky*.data$x_a)/(.data$x_a*.data$y_b - .data$x_b*.data$y_a),
      #step sizes
      hx = 1/(.data$nx + 1),
      hy = 1/(.data$ny + 1),
      #indices of neighbouring nodes
      i1 = index_offset(.data$i, .data$nx, .data$ny, dix = -1, i0 = .data$i0),
      i2 = index_offset(.data$i, .data$nx, .data$ny, diy = 1, i0 = .data$i0),
      i3 = index_offset(.data$i, .data$nx, .data$ny, dix = 1, i0 = .data$i0),
      i4 = index_offset(.data$i, .data$nx, .data$ny, diy = -1, i0 = .data$i0),
      #matrix values of all neighbouring nodes
      val1 = -0.5/.data$hx*.data$mult_ha,
      val2 = 0.5/.data$hy*.data$mult_hb,
      val3 = 0.5/.data$hx*.data$mult_ha,
      val4 = -0.5/.data$hy*.data$mult_hb
    )

  #fixed flow
  dedge_q <- dplyr::filter(dedge_q_c, .data$type == "q")
  lhs[dedge_q$i_ghost] <- dedge_q$value
  M_bc_q <- tibble::tibble(
    row = rep(dedge_q$i_ghost, 4),
    col = c(dedge_q$i1, dedge_q$i2, dedge_q$i3, dedge_q$i4),
    val = c(dedge_q$val1, dedge_q$val2, dedge_q$val3, dedge_q$val4)
  )

  #connections - same flow (edges 1 and 4 - negative)
  dedge_cq <- dedge_q_c %>%
    dplyr::filter((.data$type == "c")) %>%
    dplyr::left_join(
      dedge_c %>% dplyr::select(.data$i, .data$edge, .data$i_2, .data$edge_2),
      by = c("i", "edge")
    )
  dedge_cq2 <- dplyr::left_join(
    dedge_cq %>%
      dplyr::filter(.data$edge %in% c(1, 4)),
    dedge_cq %>%
      dplyr::filter(.data$edge %in% c(2, 3)) %>%
      dplyr::select(
        .data$i, .data$edge,
        .data$i1, .data$i2, .data$i3, .data$i4,
        .data$val1, .data$val2, .data$val3, .data$val4
      ),
    by = c("i_2" = "i", "edge_2" = "edge"),
    suffix = c("_a", "_b")
  )
  M_bc_cq <- tibble::tibble(
    row = rep(dedge_cq2$i_ghost, 8),
    col = c(
      dedge_cq2$i1_a, dedge_cq2$i2_a, dedge_cq2$i3_a, dedge_cq2$i4_a,
      dedge_cq2$i1_b, dedge_cq2$i2_b, dedge_cq2$i3_b, dedge_cq2$i4_b
    ),
    val = c(
      dedge_cq2$val1_a, dedge_cq2$val2_a, dedge_cq2$val3_a, dedge_cq2$val4_a,
      dedge_cq2$val1_b, dedge_cq2$val2_b, dedge_cq2$val3_b, dedge_cq2$val4_b
    )
  )

  #return matrix elements and lhs
  return(list(
    mat = dplyr::bind_rows(M_el, M_bc_h, M_bc_q, M_bc_ch, M_bc_cq),
    lhs = lhs
  ))
}


#' Solve a flow net problem with quadrilateral domains
#'
#' @description
#' Take a quadrilateral flow net problem and calculate results for the
#' hydraulic head, flow rates and flow potential, using finite differences
#'
#' @param df list with geometry and boundary conditions of the problem. See
#'   function `flownet_geometry_quadrilateral()` for more information
#' @return a tibble with domains (`domain`) and positions of all nodes
#'   (`a`, `b`, `x`, `y`) and solutions for head (`h`), flow rates (`qx`,
#'   `qy`) and flow potential (`psi`)
#' @importFrom magrittr `%>%`
#' @importFrom rlang .data
#' @examples
#' df <- flownet_geometry_quadrilateral(grid_size = 0.1)
#' flownet_solve_quadrilateral(df)
#' @export

flownet_solve_quadrilateral <- function(df) {
  #positions of real nodes
  dp <- purrr::pmap_dfr(df$dom, positions_real_quadrilateral)
  #solve for flow potential in all real nodes
  n_all <- sum(df$dom$n)
  n_real <- sum(df$dom$nx * df$dom$ny)
  #filter for real nodes
  i_real <- unlist(purrr::pmap(df$dom, index_real))
  #get position in all node from real node
  f_real <- rep(FALSE, n_all)
  f_real[i_real] <- TRUE
  i_real2 <- cumsum(f_real)

  ## SOLVE HEAD
  #create matrix elements for Poissons equation
  M_el <- matrix_sparse_head_quadrilateral(df)
  #solve matrix for head h
   #h <- as.vector(Matrix::solve(
   #  Matrix::sparseMatrix(
   #    i = M_el$mat$row,
   #    j = M_el$mat$col,
   #    x = M_el$mat$val,
   #    dims = rep(sum(df$dom$n), 2)
   #  ),
   #  M_el$lhs)
   #)
  mat <- Matrix::sparseMatrix(
    i = M_el$mat$row,
    j = M_el$mat$col,
    x = M_el$mat$val,
    dims = rep(sum(df$dom$n), 2)
  )
  h <- as.vector(
    Matrix::solve(
      Matrix::t(mat) %*% mat,
      Matrix::t(mat) %*% M_el$lhs
    )
  )
  #add solution of <h> to real points
  dp$h <- h[i_real]

  ## CALCULATE FLOW RATES
  #first order derivative of head in a-direction
  D1a_el <- purrr::pmap_dfr(
    df$dom %>% dplyr::select(.data$nx, .data$ny, .data$i0),
    function(nx, ny, i0) findiff_sparse_elements(nx, ny, "x", i0 = i0)
  )
  D1a <- Matrix::sparseMatrix(
    i = D1a_el$row,
    j = D1a_el$col,
    x = D1a_el$val,
    dims = rep(sum(df$dom$n), 2)
  )
  h_a <- as.vector(D1a %*% h)
  #first order derivative of head in b-direction
  D1b_el <- purrr::pmap_dfr(
    df$dom %>% dplyr::select(.data$nx, .data$ny, .data$i0),
    function(nx, ny, i0) findiff_sparse_elements(nx, ny, "y", i0 = i0)
  )
  D1b <- Matrix::sparseMatrix(
    i = D1b_el$row,
    j = D1b_el$col,
    x = D1b_el$val,
    dims = rep(sum(df$dom$n), 2)
  )
  h_b <- as.vector(D1b %*% h)
  #flow rates in x and y-directions
  temp <- with(dp, x_a*y_b - x_b*y_a)
  dp$qx <- -h_a[i_real]*df$dom$kx[dp$domain]*dp$y_b/temp + h_b[i_real]*df$dom$kx[dp$domain]*dp$y_a/temp
  dp$qy <- h_a[i_real]*df$dom$ky[dp$domain]*dp$x_b/temp - h_b[i_real]*df$dom$ky[dp$domain]*dp$x_a/temp

  ## CALCULATE FLOW POTENTIAL psi
  #offset to go from real nodes to all nodes
  i_offset <- i_real - seq(n_real)
  #lhs - q
  lhs <- with(dp, c(-x_a*qy + y_a*qx, -x_b*qy + y_b*qx))
  #real starting node
  df$dom$i0_real <- utils::head(cumsum(c(0, df$dom$nx*df$dom$ny)), -1)
  #rhs F * psi
  D1a_F <- purrr::pmap_dfr(
    df$dom %>% dplyr::select(.data$nx, .data$ny, .data$i0_real),
    function(nx, ny, i0_real) findiff_sparse_elements_realonly(nx, ny, "x", i0 = i0_real)
  )
  D1b_F <- purrr::pmap_dfr(
    df$dom %>% dplyr::select(.data$nx, .data$ny, .data$i0_real),
    function(nx, ny, i0_real) findiff_sparse_elements_realonly(nx, ny, "y", i0 = i0_real)
  )
  #ensure same flow potential on edges
  fun_temp_edge <- function(nx, ny, edge, i0_real = 0, ...){
    if (edge == 1) {
      return(i0_real + 1 + seq(0, ny - 1)*nx)
    } else if (edge == 2) {
      return(i0_real + nx*(ny - 1) + seq(nx))
    } else if (edge == 3) {
      return(i0_real + nx + seq(0, ny - 1)*nx)
    } else {
      return(i0_real + seq(nx))
    }
  }
  edge_real_1 <- df$bc %>%
    dplyr::filter((.data$type == "c") & (.data$edge %in% c(2, 3))) %>%
    dplyr::left_join(
      df$dom %>% dplyr::select(.data$domain, .data$nx, .data$ny, .data$i0_real),
      by = "domain"
    )
  i_edge_1 <- purrr::pmap(edge_real_1, fun_temp_edge) %>% unlist()
  edge_real_2 <- edge_real_1 %>%
    dplyr::left_join(
      df$dom %>% dplyr::select(.data$domain, .data$nx, .data$ny, .data$i0_real),
      by = c("value" = "domain"),
      suffix = c("_1", "_2")
    ) %>%
    dplyr::mutate(edge_2 = 1 + (1 + .data$edge)%%4) %>%
    dplyr::select(.data$nx_2, .data$ny_2, .data$i0_real_2, .data$edge_2) %>%
    dplyr::rename("nx" = "nx_2", "ny" = "ny_2", "edge" = "edge_2", "i0_real" = "i0_real_2")
  i_edge_2 <- purrr::pmap(edge_real_2, fun_temp_edge) %>% unlist()
  #generate sparse matrix
  Dpsi <- Matrix::sparseMatrix(
    i = c(
      D1a_F$row,
      D1b_F$row + n_real,
      rep(seq(length(i_edge_1)), 2) + 2*n_real
    ),
    j = c(
      D1a_F$col,
      D1b_F$col,
      i_edge_1, i_edge_2
    ),
    x = c(
      D1a_F$val,
      D1b_F$val,
      rep(c(-1, 1) / mean(sqrt(df$dom$kx*df$dom$ky)), each = length(i_edge_1))
    ),
    dims = c(2*n_real + length(i_edge_1), n_real)
  )
  lhs2 <- c(lhs, rep(0, length(i_edge_1)))
  #solve flow potential
  psi <- as.vector(Matrix::solve(
    (Matrix::t(Dpsi) %*% Dpsi),
    (Matrix::t(Dpsi) %*% lhs2)
  ))
  dp$psi <- psi - min(psi)

  #return
  return(dplyr::select(
    dp,
    .data$domain, .data$a, .data$b, .data$x, .data$y,
    .data$h, .data$qx, .data$qy, .data$psi
  ))
}


#' Interpolate quadrilateral flownet solution on rectangular grid
#'
#' @description
#' Take a flow net solution calculated on a quadrilateral grid, and
#' interpolate the results on a regular grid. This grid can be used
#' to calculate contour lines (flow lines etc)
#'
#' @param df list with geometry and boundary conditions of the problem. See
#'   function `flownet_geometry_quadrilateral()` for more information
#' @param dp tibble with solution, see function `flownet_solve_quadrilateral()`
#' @param grid_size grid size to use for rectangular interpolation grid. If
#'   not defined, it is chosen roughly in line with the smallest quadrilateral
#'   grid cell
#' @return a tibble with flow net results at an rectangular grid
#' @importFrom magrittr `%>%`
#' @importFrom rlang .data
#' @examples
#' df <- flownet_geometry_quadrilateral(grid_size = 0.25)
#' dp <- flownet_solve_quadrilateral(df)
#' di <- flownet_interpolate_quadrilateral(df, dp)
#'
#' ggplot2::ggplot() +
#'   ggplot2::geom_tile(
#'     data = di,
#'     ggplot2::aes(x = x, y = y, fill = h)
#'   ) +
#'   ggplot2::geom_point(
#'     data = dp,
#'     ggplot2::aes(x = x, y = y, fill = h),
#'     color = "black",
#'     pch = 21
#'   ) +
#'   ggplot2::scale_fill_distiller(palette = "RdYlBu")
#' @export

flownet_interpolate_quadrilateral <- function(
  df,
  dp,
  grid_size = NULL
){
  ## if not specifically defined, get grid size from data
  if (is.null(grid_size)) {
    h1 <- with(df$dom, sqrt((x2 - x1)^2 + (y2 - y1)^2)/ny)
    h2 <- with(df$dom, sqrt((x3 - x2)^2 + (y3 - y2)^2)/ny)
    h3 <- with(df$dom, sqrt((x4 - x3)^2 + (y4 - y3)^2)/ny)
    h4 <- with(df$dom, sqrt((x1 - x4)^2 + (y1 - y4)^2)/ny)
    grid_size <- min(c(h1, h2, h3, h4))
  }
  #number of interpolating points
  nx <- ceiling((max(dp$x) - min(dp$x))/grid_size)
  ny <- ceiling((max(dp$y) - min(dp$y))/grid_size)
  ## generate all points
  di <- tidyr::expand_grid(
    x = seq(min(dp$x), max(dp$x), l = nx),
    y = seq(min(dp$y), max(dp$y), l = ny)
  ) %>%
    #find which domain each point belongs to
    tidyr::expand_grid(
      df$dom %>% dplyr::select(
        .data$domain,
        .data$x1, .data$x2, .data$x3, .data$x4,
        .data$y1, .data$y2, .data$y3, .data$y4
      )
    ) %>%
    dplyr::mutate(
      b1 = ((.data$x - .data$x1)*(.data$y2 - .data$y1)) >= ((.data$y - .data$y1)*(.data$x2 - .data$x1)),
      b2 = ((.data$x - .data$x2)*(.data$y3 - .data$y2)) >= ((.data$y - .data$y2)*(.data$x3 - .data$x2)),
      b3 = ((.data$x - .data$x3)*(.data$y4 - .data$y3)) >= ((.data$y - .data$y3)*(.data$x4 - .data$x3)),
      b4 = ((.data$x - .data$x4)*(.data$y1 - .data$y4)) >= ((.data$y - .data$y4)*(.data$x1 - .data$x4))
    ) %>%
    dplyr::filter(.data$b1*.data$b2*.data$b3*.data$b4 == 1) %>%
    dplyr::distinct(.data$x, .data$y, .keep_all = TRUE)
  ## get a and b values for each point
  #initiate
  di$a <- NA
  di$b <- NA
  #constants
  di <- dplyr::mutate(
    di,
    c1 = -.data$x1 + .data$x4,
    c2 = -.data$x1 + .data$x2,
    c3 = .data$x1 - .data$x2 + .data$x3 - .data$x4,
    c4 = -.data$y1 + .data$y4,
    c5 = -.data$y1 + .data$y2,
    c6 = .data$y1 - .data$y2 + .data$y3 - .data$y4
  )
  #c3 = 0, c6 = 0 - all a-lines parallel, all b-lines parallel
  i36 <- with(di, (c3 == 0) & (c6 == 0))
  di$a[i36] <- with(di[i36, ], (c5*(x - x1) - c2*(y - y1))/(c1*c5 - c2*c4))
  di$b[i36] <- with(di[i36, ], (-c4*(x - x1) + c1*(y - y1))/(c1*c5 - c2*c4))
  #c1 = 0, c3 = 0 - all a-lines vertical
  i13 <- with(di, (c2 == 0 ) & (c2 != 0) & (c3 == 0) & (c6 != 0))
  di$b[i13] <- with(di[i13, ], (x - x1)/c2)
  di$a[i13] <- with(di[i13, ], (y - y1 - c5*b)/(c4 + c6*b))
  #c2 = 0, c3 = 0 - all a-lines horizontal
  i23 <- with(di, (c1 != 0 ) & (c2 == 0) & (c3 == 0) & (c6 != 0))
  di$a[i23] <- with(di[i23, ], (x - x1)/c1)
  di$b[i23] <- with(di[i23, ], (y - y1 - c4*a)/(c5 + c6*a))
  #c4 = 0, c6 = 0 - all b-lines vertical
  i46 <- with(di, (c4 == 0) & (c5 != 0) & (c6 == 0) & (c3 != 0))
  di$b[i46] <- with(di[i46, ], (y - y1)/c5)
  di$a[i46] <- with(di[i46, ], (x - x1 - c2*b)/(c1 + c3*b))
  #c5 = 0, c6 = 0 - all b-lines vertical
  i56 <- with(di, (c4 != 0) & (c5 == 0) & (c6 == 0) & (c3 != 0))
  di$a[i56] <- with(di[i56, ], (y - y1)/c4)
  di$b[i56] <- with(di[i56, ], (x - x1 - c1*a)/(c2 + c3*a))
  #c1 != 0, c2 != 0, c3 = 0, c6 != 0
  i3 <- with(di, (c1 != 0) & (c2 != 0) & (c3 == 0) & (c6 != 0))
  aa <- with(di[i3, ], -c1*c6)
  bb <- with(di[i3, ], c2*c4 + c6*(x - x1) - c1*c5)
  cc <- with(di[i3, ], c5*(x - x1) - c2*(y - y1))
  a1 <- (-bb + sqrt(bb^2 - 4*aa*cc))/(2*aa)
  a2 <- (-bb - sqrt(bb^2 - 4*aa*cc))/(2*aa)
  di$a[i3] <- a1*(a1>=0 & a1<=1) + a2*(a2>=0 & a2<=1)
  di$b[i3] <- with(di[i3, ], (x - x1 - c1*a)/c2)
  #c3 != 0, c4 != 0, c5 != 0, c6 = 0
  i6 <- with(di, (c3 != 0) & (c4 != 0) & (c5 == 0) & (c3 != 0))
  aa <- with(di[i6, ], -c3*c5)
  bb <- with(di[i6, ], c2*c4 + c3*(y - y1) - c1*c5)
  cc <- with(di[i6, ], c1*(y - y1) - c4*(x - x1))
  b1 <- (-bb + sqrt(bb^2 - 4*aa*cc))/(2*aa)
  b2 <- (-bb - sqrt(bb^2 - 4*aa*cc))/(2*aa)
  di$b[i6] <- b1*(b1>=0 & b1<=1) + b2*(b2>=0 & b2<=1)
  di$a[i6] <- with(di[i6, ], (y - y1 - c5*b)/c4)
  #c3 != 0, c6 != 0
  i0 <- with(di, (c3 != 0) & (c6 != 0))
  aa <- with(di[i0, ], -c3*(c3*c5 - c2*c6))
  bb <- with(di[i0, ], c2*(c3*c4 - c1*c6) + c3*(c3*(y - y1) - c6*(x - x1)) - c1*(c3*c5 - c2*c6))
  cc <- with(di[i0, ], c1*c3*(y - y1) - c3*c4*(x - x1))
  b1 <- (-bb + sqrt(bb^2 - 4*aa*cc))/(2*aa)
  b2 <- (-bb - sqrt(bb^2 - 4*aa*cc))/(2*aa)
  di$b[i0] <- b1*(b1>=0 & b1<=1) + b2*(b2>=0 & b2<=1)
  di$a[i0] <- with(di[i0, ], (x - x1 - c2*b)/(c1 + c3*b))
  #within limits (sort potential rounding issues)
  di$a = pmax(0, pmin(di$a, 1))
  di$b = pmax(0, pmin(di$b, 1))
  ## join domain properties
  di <- di %>%
    dplyr::left_join(
      df$dom %>%
        dplyr::mutate(i0_real = utils::head(cumsum(c(0, .data$nx*.data$ny)), -1)) %>%
        dplyr::select(.data$domain, .data$nx, .data$ny, .data$i0_real),
      by = "domain"
    ) %>%
    dplyr::mutate(
      #get indices and fractions for each node in grid
      ix0 = 1 + floor(.data$a*(.data$nx - 1)),
      ix1 = 1 + ceiling(.data$a*(.data$nx - 1)),
      iy0 = 1 + floor(.data$b*(.data$ny - 1)),
      iy1 = 1 + ceiling(.data$b*(.data$ny - 1)),
      fx0 = ifelse(
        (.data$ix0 == .data$ix1),
        1,
        (.data$ix1 - (1 + .data$a*(.data$nx - 1)))/(.data$ix1 - .data$ix0)
      ),
      fx1 = 1 - .data$fx0,
      fy0 = ifelse(
        (.data$iy0 == .data$iy1),
        1,
        (.data$iy1 - (1 + .data$b*(.data$ny - 1)))/(.data$iy1 - .data$iy0)),
      fy1 = 1 - .data$fy0,
      i00 = .data$i0_real + .data$ix0 + .data$nx*(.data$iy0 - 1),
      i10 = .data$i0_real + .data$ix1 + .data$nx*(.data$iy0 - 1),
      i01 = .data$i0_real + .data$ix0 + .data$nx*(.data$iy1 - 1),
      i11 = .data$i0_real + .data$ix1 + .data$nx*(.data$iy1 - 1),
      v00 = .data$fx0*.data$fy0,
      v10 = .data$fx1*.data$fy0,
      v01 = .data$fx0*.data$fy1,
      v11 = .data$fx1*.data$fy1
    )
  #get interpolation sparse matrix
  Mint <- Matrix::sparseMatrix(
    i = rep(seq(nrow(di)), 4),
    j = c(di$i00, di$i10, di$i01, di$i11),
    x = c(di$v00, di$v10, di$v01, di$v11),
    dims = c(nrow(di), sum(df$dom$nx*df$dom$ny))
  )
  #calculate
  if ("h" %in% colnames(dp)) {
    di$h <- as.vector((Mint %*% dp$h))
  }
  if ("qx" %in% colnames(dp)) {
    di$qx <- as.vector((Mint %*% dp$qx))
  }
  if ("qy" %in% colnames(dp)) {
    di$qy <- as.vector((Mint %*% dp$qy))
  }
  if ("psi" %in% colnames(dp)) {
    di$psi <- as.vector((Mint %*% dp$psi))
  }
  if ("val" %in% colnames(dp)) {
    di$val <- as.vector((Mint %*% dp$val))
  }
  #return
  return(di)
}
