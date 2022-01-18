#v0.1 - 20210416 - first working version

#####################
### DEFINE SERVER ###
#####################

# Define server
server <- function(input, output) {
  #plot initial element
  p1 <- reactive({
    soilmech::ggplot_stresselement(
      sigz = input$sigz,
      sigx = input$sigx,
      tau = input$tau,
      theta = 0,
      clockwise_shear = !input$counterclockwise_shear,
      coordinate_system = TRUE
    )
  })
  #plot rotated element
  p2 <- reactive({
    soilmech::ggplot_stresselement(
      sigz = input$sigz,
      sigx = input$sigx,
      tau = input$tau,
      theta = input$theta/180*pi,
      clockwise_shear = !input$counterclockwise_shear,
      coordinate_system = FALSE
    )
  })
  #plot Mohr circle
  p3 <- reactive({
    soilmech::ggplot_mohrcircle(
      sigz = input$sigz,
      sigx = input$sigx,
      tau = input$tau,
      theta = input$theta/180*pi,
      clockwise_shear = !input$counterclockwise_shear
    )
  })
  #join
  output$plot_all <- renderPlot({
    gridExtra::grid.arrange(
      p1(), p2(), p3(),
      layout_matrix = rbind(c(1,2), c(3,3))
    )
  })
}
