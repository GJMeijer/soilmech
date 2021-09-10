###########
### OLD ###
###########




#' Get indices of nodes bordering nodes
#'
#' @description
#' Get the indices of neighbouring nodes in a rectangular grid with nodes.
#' For node numbering, see function `findiff_sparse_entries()`
#'
#' @inheritParams index_edge
#' @param i array with node indices
#' @param offset number of nodes to move in specific direction
#' @importFrom magrittr `%>%`
#' @return an array with node indices for neighbours
#' @export

index_neighbour <- function(i, nx, ny, edge, offset = 1, i0 = 0, real_only = FALSE, ...) {
  #create tibble
  df <- tibble::tibble(i = i, nx = nx, ny = ny, edge = edge, i0 = i0, offset = offset, ioffset = 0)
  #grid with only real nodes
  if (real_only == TRUE) {
    #edge 1 (left)
    ind <- ((1 + (df$i-1)%%df$nx) <= (df$nx - df$offset)) & (df$edge == 1)
    df$ioffset[ind] <- df$i[ind] + 1
    #edge 2 (top)
    ind <- ((1 + (df$i-1)%%df$nx) <= (df$nx - df$offset)) & (df$edge == 1)
    df$ioffset[ind] <- df$i[ind] + 1
    #edge 3 (right)
    ind <- ((1 + (df$i-1)%%df$nx) > (df$offset)) & (df$edge == 3)
    df$ioffset[ind] <- df$i[ind] - 1
    #edge 4 (bottom)
  } else


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

