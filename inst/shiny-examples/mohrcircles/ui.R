
################
### DEFINE UI ###
#################

#make UI
ui <- fluidPage(
  titlePanel("Mohr circles"),
  sidebarPanel(
    sliderInput(
      "sigz",
      HTML(paste0("Vertical stress \u03C3", tags$sub("z"), " [kPa]")),
      #"Vertical stress \u03C3<sub>z<\sub> [kPa]",
      min = 0,
      max = 100,
      value = 40,
      step = 1
    ),
    sliderInput(
      "sigx",
      HTML(paste0("Horizontal stress \u03C3", tags$sub("x"), " [kPa]")),
      min = 0,
      max = 100,
      value = 20,
      step = 1
    ),
    sliderInput(
      "tau",
      HTML(paste0("Shear stress \u03c4", tags$sub("xz"), " [kPa]")),
      min = 0,
      max = 100,
      value = 15,
      step = 1
    ),
    checkboxInput(
      "counterclockwise_shear",
      "Counter-clockwise sheear",
      value = TRUE
    ),
    sliderInput(
      "theta",
      "Plane rotation [\u00B0]",
      min = -180,
      max = 180,
      value = 0,
      step = 1
    )
  ),
  mainPanel(
    plotOutput("plot_all")
  )
)
