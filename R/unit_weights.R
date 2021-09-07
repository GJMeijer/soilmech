#' Function to obtain unit weight from soil description
#'
#' @description
#' Get soil unit weight, above and below the water table, from a soil
#' description
#' @param description array with descriptions
#' @round round unit weights to be nearest multiple of this number
#' @examples
#' description <- c("very dense sand", "loose gravel")
#' get_unitweight(description)
#' @export

get_unitweight <- function(
  description,
  round = 0.5
){
  #convert to lower vase
  do <- tibble::tibble(
    description = stringr::str_to_lower(description),
    gamma_b_above = NA,
    gamma_b_below = NA
  )
  #determine relative density - if not provided, take value halfway range
  ID <- rep(0.5, length(description))
  ID[stringr::str_detect(do$description, "loose")] <- 0.5*(0.15 + 0.35)
  ID[stringr::str_detect(do$description, "very loose")] <- 0.5*(0.00 + 0.15)
  ID[stringr::str_detect(do$description, "dense")] <- 0.5*(0.65 + 0.85)
  ID[stringr::str_detect(do$description, "medium dense")] <- 0.5*(0.35 + 0.65)
  ID[stringr::str_detect(do$description, "very dense")] <- 0.5*(0.85 + 1.00)
  #find primary soil type - last word
  primary <- purrr::map_chr(
    stringr::str_split(description, " "),
    ~utils::tail(.x, 1)
  )
  #assign - gravel
  ind <- (primary == "gravel")
  do$gamma_b_above[ind] <- 18.5 + (22.0 - 18.5)*ID[ind]
  do$gamma_b_below[ind] <- 20.5 + (23.5 - 20.5)*ID[ind]
  #assign - sand
  ind <- (primary == "sand")
  do$gamma_b_above[ind] <- 16.5 + (20.5 - 16.5)*ID[ind]
  do$gamma_b_below[ind] <- 18.0 + (22.0 - 18.0)*ID[ind]
  #assign - silt
  ind <- (primary == "silt")
  do$gamma_b_above[ind] <- 15.5 + (21.0 - 15.5)*ID[ind]
  do$gamma_b_below[ind] <- 18.0 + (21.5 - 18.0)*ID[ind]
  #assign - clay
  ind <- (primary == "clay")
  do$gamma_b_above[ind] <- 14.0 + (17.0 - 14.0)*ID[ind]
  do$gamma_b_below[ind] <- 16.0 + (18.5 - 16.0)*ID[ind]
  #assign - peat
  ind <- (primary == "peat")
  do$gamma_b_above[ind] <- 12.0
  do$gamma_b_below[ind] <- 12.0
  #round
  do$gamma_b_above <- round(do$gamma_b_above/round)*round
  do$gamma_b_below <- round(do$gamma_b_below/round)*round
  #return
  do %>% dplyr::select(gamma_b_above, gamma_b_below)
}
