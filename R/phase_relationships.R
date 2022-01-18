#' Plotly linked volumetric relationships
#'
#' @description
#' Plot a plotly chart that shows how various phase relationships are
#' interconnected.
#'
#' @param path_param,path_linking,path_equations paths to raw data
#' @param width,height plot width and height
#' @importFrom magrittr `%>%`
#' @examples
#' plotly_phase_relationships()
#' @export

plotly_phase_relationships <- function(
  path_param = "./inst/extdata/phaserelationships_parameters.csv",
  path_linking = "./inst/extdata/phaserelationships_parameter_relationships.csv",
  path_equations = "./inst/extdata/phaserelationships_relationships.csv",
  width = 400,
  height = 400
) {
  #load data
  df1 <- readr::read_csv(
    path_param,
    col_types = readr::cols()
  )
  df2 <- readr::read_csv(
    path_linking,
    col_types = readr::cols()
  )
  df3 <- readr::read_csv(
    path_equations,
    col_types = readr::cols()
  )
  #merge positions
  df2b <- df2 %>%
    dplyr::group_by(relation_id, .groups = "drop") %>%
    dplyr::summarize(parameter_id = c(parameter_id, parameter_id[1])) %>%
    dplyr::ungroup() %>%
    dplyr::left_join(
      df1 %>% dplyr::select(parameter_id, x, y),
      by = "parameter_id"
    ) %>%
    dplyr::left_join(
      df3 %>% dplyr::select(relation_id, equation_html),
      by = "relation_id"
    )
  #hoverlabels
  df1$hoverlabel <- paste0(df1$parameter_html, "<br>", df1$name)
  #plot
  plotly::plot_ly(
    height = height,
    width = width
  ) %>%
    plotly::add_polygons(
      data = df2b,
      x = ~x,
      y = ~y,
      split = ~relation_id,
      hoverinfo = "text",
      text = ~equation_html,
      hoveron = "fills",
      showlegend = FALSE
    ) %>%
    plotly::add_trace(
      type = "scatter",
      mode = "markers",
      data = df1,
      x = ~x,
      y = ~y,
      text = ~hoverlabel,
      showlegend = FALSE,
      marker = list(
        size = 20,
        color = "white",
        line = list(
          color = "black",
          width = 1
        )
      ),
      #name = ~hoverlabel,
      hoverinfo = "text"
    ) %>%
    plotly::add_trace(
      type = "scatter",
      mode = "text",
      data = df1,
      x = ~x,
      y = ~y,
      text = ~parameter_html,
      showlegend = FALSE,
      hoverinfo = "skip"
    ) %>%
    plotly::layout(
      xaxis = list(
        showgrid = FALSE,
        zeroline = FALSE,
        visible = FALSE
      ),
      yaxis = list(
        showgrid = FALSE,
        zeroline = FALSE,
        visible = FALSE,
        scaleanchor = "x",
        scaleratio = 1
      )
    )
}

plotly_phase_relationships()
