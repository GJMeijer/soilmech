#assign unit weights based on soil description
assign_unitweight <- function(description){
  #COURSE SOILS - relative density terms for sands - based on ISO 14688
  relative_density <- tibble::tibble(
    description = c("very loose", "loose", "medium dense", "dense", "very dense"),
    ID0 = c(0.00, 0.15, 0.35, 0.65, 0.85),
    ID1 = c(0.15, 0.35, 0.65, 0.85, 1.00)
  )
  #COURSE SOILS - void ratio and rounded
  void_ratio_uniformity <- tibble::tibble(
    Cu = rep(c(1, 2, 5, 10), 2),
    shape = rep(c("rounded", "angular"), each = 4),
    emin = c(0.50, 0.40, 0.30, 0.25, 0.80, 0.65, 0.50, 0.40),
    emax = c(0.80, 0.60, 0.50, 0.40, 1.30, 1.20, 0.90, 0.70),
  )
  10*(2.65+void_ratio_uniformity$emin)/(1+void_ratio_uniformity$emin)
  10*(2.65+void_ratio_uniformity$emax)/(1+void_ratio_uniformity$emax)


  #clay softness
  clay_density <- tibble::tibble(
    label = c("very soft", "soft", "firm", "stiff", "very stiff"),
    factor0 = c(0.0, 0.2, 0.4, 0.6, 0.8),
    factor1 = c(0.2, 0.4, 0.6, 0.8, 1.0)
  )
}


#stress_profile <- function(
  thickness <- c(20, 21, 18)
  unitweight <- 20
  surcharge <- 10
  zsurface <- 0
  zw = NULL

#create tibble with stresses
ds <- tibble::tibble(
  z0 = zsurface + c(0, cumsum(utils::head(thickness, -1))),
  z1 = zsurface + cumsum(thickness),
  gammab = surcharge,
  q = c(surcharge, rep(0, length(thickness) - length(surcharge))),
  sigv0 = cumsum(q) + cumsum((z1 - z0)*gammab),
  sigv1 = 1
)
