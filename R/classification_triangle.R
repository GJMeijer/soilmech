#' Convert fractions to coordinates - plot 1
#'
#' @description
#' Function takes masses of gravel, sand and fines, and converts their
#' quantities to cartesian x-y positions in the triangular part of the
#' classification chart. x=0 and y=0 conincides with the lower left
#' corner, and x=1 and y=0 with the lower right corner
#'
#' @param gravel mass of gravel (array)
#' @param sand mass of sand (array)
#' @param fines mass of fines (array)
#' @return dataframe with x and y positions
#' @export

convert_to_grid_part1 <- function(gravel, sand, fines){
  #get axes and normalise
  fines_normalised <- fines/(gravel + sand + fines)
  sand_normalised <- sand/(gravel + sand + fines)
  #convert to x-y - in triangle plot
  x <- 1 - fines_normalised - sand_normalised*cos(pi/3)
  y <- sand_normalised*cos(pi/6)
  #return
  tibble::tibble(x = x, y = y)
}


#' Convert fractions to coordinates - plot 2
#'
#' @description
#' Function takes masses of fines and clay (already normalised by total mass),
#' and converts their quantities to cartesian x-y positions in the
#' bottom square part of the classification chart. x=0 and y=0 conincides with
#' the lower left corner of the graph, and x=1 and y=0 with the lower right
#' corner.
#'
#' @param fines_fraction normalised mass of gravel (array)
#' @param clay_fraction normalised mass of sand (array)
#' @return dataframe with x and y positions
#' @export

convert_to_grid_part2 <- function(fines_fraction, clay_fraction){
  #convert to x-y - in square plot
  x <- 1 - fines_fraction
  y <- 1 - clay_fraction
  #return
  tibble::tibble(x = x, y = y)
}


#' Create soil classification chart - part 1
#'
#' @description
#' Function returns ggplot of the first parts of the BS soil classification
#' chart: the triagular upper part
#'
#' @param gravel fraction of gravel for crosshair annotation
#' @param sand fraction of sand for crosshair annotation
#' @param fines fraction of fines for crosshair annotation
#' @param group grouping for crosshair annotation (label)
#' @param grid_ticks distance between axis ticks, in \%
#' @param grid_label distance between axis tick labels, in \%
#' @param tick_length length of tick mark, in \%
#' @param label_length distance of tick label to chart axis, in \%
#' @param label_size font size of tick label
#' @param title_length distance to axis title to chart, in \%
#' @param title_size font size of axis title
#' @param soil_size font size of soil labels
#' @param soil_size2 font size of soil labels, for the small ones in the
#'   middle of the triangular plot
#' @param soil_spacing line spacing for soil labels
#' @param polygon_line thickness of soil polygon lines
#' @param polygon_color color of soil polygon lines
#' @param xlim plot x-limits
#' @param ylim plot y-limits
#' @importFrom rlang .data
#' @return a ggplot object
#' @examples
#' #empty chart
#' ggplot_classificationtriangle_part1()
#'
#' #chart with some annotations
#' ggplot_classificationtriangle_part1(
#'   gravel = 0.25,
#'   sand = 0.30,
#'   fines = 0.45
#' )
#' @export

ggplot_classificationtriangle_part1 <- function(
  gravel = NULL,
  sand = NULL,
  fines = NULL,
  group = NULL,
  grid_ticks = 10,
  grid_label = 20,
  tick_length = 2,
  label_length = 7,
  label_size = 2,
  title_length = 15,
  title_size = 3,
  soil_size = 2.0,
  soil_size2 = 1.4,
  soil_spacing = 0.75,
  polygon_line = 0.4,
  polygon_color = "grey60",
  xlim = c(-0.2, 1.2),
  ylim = c(-0.2, 1.0)
){

  #################
  ### LOAD DATA ###
  #################

  #file path
  file_path <- system.file(
    "extdata",
    "classificationtriangle_polygons_part1.csv",
    package = "soilmech"
  )
  #read polygon data
  dp1f <- readr::read_csv(file_path, col_types = readr::cols())
  #grid polygons
  dp1 <- dplyr::mutate(dp1f, convert_to_grid_part1(.data$gravel, .data$sand, .data$fines))


  ##############################
  ### PLOT AXIS ETC - PLOT 1 ###
  ##############################

  #generate ticks for plot 1
  nt <- 1 + floor(100/grid_ticks)
  dt1 <- tibble::tibble(
    fines = c(
      seq(0, 100, grid_ticks), seq(0, 100, grid_ticks),
      seq(100, 0, -grid_ticks), seq(100, 0, -grid_ticks) + tick_length,
      rep(0, nt), rep(0, nt) - tick_length
    ),
    sand = c(
      rep(0, nt), rep(0,nt) - tick_length,
      seq(0, 100, grid_ticks), seq(0, 100, grid_ticks),
      seq(100, 0, -grid_ticks), seq(100, 0, -grid_ticks) + tick_length
    ),
    gravel = c(
      seq(100, 0, -grid_ticks), seq(100, 0, -grid_ticks) + tick_length,
      rep(0, nt), rep(0, nt) - tick_length,
      seq(0, 100, grid_ticks), seq(0, 100, grid_ticks)
    )
  )
  dt1$tick_group <- rep(rep(seq(nt), 2), 3) + nt*rep(seq(3) - 1, each = nt*2)
  dt1 <- dplyr::mutate(dt1, convert_to_grid_part1(.data$gravel, .data$sand, .data$fines))

  #generate labels for plot1
  nl <- 1 + floor(100/grid_label)
  dl1 <- tibble::tibble(
    fines = c(
      seq(0, 100, grid_label),
      seq(100, 0, -grid_label) + label_length,
      rep(0, nl) - label_length
    ),
    sand = c(
      rep(0, nl) - label_length,
      seq(0, 100, grid_label),
      seq(100, 0, -grid_label) + label_length
    ),
    gravel = c(
      seq(100, 0, -grid_label) + label_length,
      rep(0, nl) - label_length,
      seq(0, 100, grid_label)
    ),
    rotation = c(rep(-60, nl), rep(0, nl), rep(60, nl))
  )
  dl1$label <- rep(paste0(seq(0, 100, grid_label), "%"), 3)
  dl1 <- dplyr::mutate(dl1, convert_to_grid_part1(.data$gravel, .data$sand, .data$fines))

  #generate axis titles for plot 1
  da1 <- tibble::tibble(
    fines = c(50, 50 + title_length, 0 - title_length),
    sand = c(0 - title_length, 50, 50 + title_length),
    gravel = c(50 + title_length, 0 - title_length, 50),
    label = c("Fines", "Sand", "Gravel"),
    rotation = c(0, 60, -60)
  )
  da1 <- dplyr::mutate(da1, convert_to_grid_part1(.data$gravel, .data$sand, .data$fines))

  #positions and labels for soil
  dc1 <- dp1 %>%
    dplyr::group_by(.data$soiltype) %>%
    dplyr::summarize(
      n = length(.data$x),
      polygon_centroid(.data$x, .data$y)
    ) %>%
    dplyr::mutate(label = gsub(",", "\n", .data$soiltype))


  ############
  ### PLOT ###
  ############

  #plot 1 - triangle
  plt <- ggplot2::ggplot() +
    ggplot2::theme_void() +
    ggplot2::geom_polygon(
      data = dp1,
      ggplot2::aes(x = .data$x, y = .data$y, group = .data$soiltype),
      fill = NA,
      color = polygon_color,
      size = polygon_line
    ) +
    ggplot2::geom_polygon(
      ggplot2::aes(x = c(0,1, 0.5), y = c(0,0, cos(pi/6))),
      fill = NA,
      color = "black"
    ) +
    ggplot2::geom_line(
      data = dt1,
      ggplot2::aes(x = .data$x, y = .data$y, group = .data$tick_group)
    ) +
    ggplot2::geom_text(
      data = dl1,
      ggplot2::aes(x = .data$x, y = .data$y, label = .data$label, angle = .data$rotation),
      hjust = 0.5,
      vjust = 0.5,
      size = label_size
    ) +
    ggplot2::geom_text(
      data = da1,
      ggplot2::aes(x = .data$x, y = .data$y, label = .data$label, angle = .data$rotation),
      hjust = 0.5,
      vjust = 0.5,
      size = title_size
    ) +
    ggplot2::geom_text(
      data = dc1[dc1$n >= 4, ],
      ggplot2::aes(x = .data$x, y = .data$y, label = .data$label),
      size = soil_size,
      hjust = 0.5,
      vjust = 0.5,
      lineheight = soil_spacing
    ) +
    ggplot2::geom_text(
      data = dc1[dc1$n < 4, ],
      ggplot2::aes(x = .data$x, y = .data$y, label = .data$label),
      size = soil_size2,
      hjust = 0.5,
      vjust = 0.5,
      lineheight = soil_spacing
    ) +
    ggplot2::coord_fixed(ratio = 1, xlim = xlim, ylim = ylim, expand = FALSE)
  #add annotations
  if (!is.null(gravel) & !is.null(sand) & !is.null(fines)) {
    plt <- ggplot_classificationtriangle_addcrosshairs_part1(
      plt,
      gravel,
      sand,
      fines,
      group = group
    )
  }

  ##############
  ### RETURN ###
  ##############

  #return
  plt
}


#' Create soil classification chart - part 2
#'
#' @description
#' Function returns ggplot of the second part of the BS soil classification
#' chart: the square bottom part
#'
#' @inheritParams ggplot_classificationtriangle_part1
#' @param clay fraction of clay for crosshair annotation
#' @importFrom rlang .data
#' @return a ggplot object
#' @export

ggplot_classificationtriangle_part2 <- function(
  fines = NULL,
  clay = NULL,
  group = NULL,
  grid_ticks = 10,
  grid_label = 20,
  tick_length = 2,
  label_length = 7,
  label_size = 2,
  title_length = 15,
  title_size = 3,
  soil_size = 2.0,
  soil_size2 = 1.4,
  soil_spacing = 0.75,
  polygon_line = 0.4,
  polygon_color = "grey60",
  xlim = c(-0.2, 1.2),
  ylim = c(-0.2, 1.2)
){

  #################
  ### LOAD DATA ###
  #################

  #file path
  file_path <- system.file(
    "extdata",
    "classificationtriangle_polygons_part2.csv",
    package = "soilmech"
  )
  #read polygon data
  dp2f <- readr::read_csv(file_path, col_types = readr::cols())
  #grid polygons
  dp2 <- dplyr::mutate(dp2f, convert_to_grid_part2(.data$fines/100, .data$clay/100))


  ##############################
  ### PLOT AXIS ETC - PLOT 2 ###
  ##############################

  #generate ticks for plot 2
  nt <- 1 + floor(100/grid_ticks)
  dt2 <- tibble::tibble(
    fines = c(
      seq(0, 100, grid_ticks), seq(0, 100, grid_ticks),
      rep(100, nt), rep(100, nt) + tick_length
    ),
    clay = c(
      rep(0, nt), rep(0, nt) - tick_length,
      seq(0, 100, grid_ticks), seq(0, 100, grid_ticks)
    )
  )
  dt2$tick_group <- rep(rep(seq(nt), 2), 2) + nt*rep(seq(2) - 1, each = nt*2)
  dt2 <- dplyr::mutate(dt2, convert_to_grid_part2(.data$fines/100, .data$clay/100))

  #generate labels for plot2
  nl <- 1 + floor(100/grid_label)
  dl2 <- tibble::tibble(
    fines = c(
      seq(0, 100, grid_label),
      rep(100, nl) + label_length
    ),
    clay = c(
      rep(0, nl) - label_length,
      seq(0, 100, grid_label)
    ),
    rotation = c(rep(90, nl), rep(0, nl))
  )
  dl2$label <- rep(paste0(seq(0, 100, grid_label), "%"), 2)
  dl2 <- dplyr::mutate(dl2, convert_to_grid_part2(.data$fines/100, .data$clay/100))

  #generate axis titles for plot 1
  da2 <- tibble::tibble(
    fines = c(50, 100 + title_length),
    clay = c(0 - title_length, 50),
    label = c("Fines", "Clay"),
    rotation = c(0, 90)
  )
  da2 <- dplyr::mutate(da2, convert_to_grid_part2(.data$fines/100, .data$clay/100))

  #positions and labels for soil
  dc2 <- dp2 %>%
    dplyr::group_by(.data$soiltype) %>%
    dplyr::summarize(polygon_centroid(.data$x, .data$y)) %>%
    dplyr::mutate(label = gsub(",", "\n", .data$soiltype))


  ############
  ### PLOT ###
  ############

  #plot 2 - square
  plt <- ggplot2::ggplot() +
    ggplot2::theme_void() +
    ggplot2::geom_polygon(
      data = dp2,
      ggplot2::aes(x = .data$x, y = .data$y, group = .data$soiltype),
      fill = NA,
      color = polygon_color,
      size = polygon_line
    ) +
    ggplot2::geom_polygon(
      ggplot2::aes(x = c(0, 1, 1, 0), y = c(0, 0, 1, 1)),
      fill = NA,
      color = "black"
    ) +
    ggplot2::geom_line(
      data = dt2,
      ggplot2::aes(x = .data$x, y = .data$y, group = .data$tick_group)
    ) +
    ggplot2::geom_text(
      data = dl2,
      ggplot2::aes(x = .data$x, y = .data$y, label = .data$label, angle = .data$rotation),
      hjust = 0.5,
      vjust = 0.5,
      size = label_size
    ) +
    ggplot2::geom_text(
      data = da2,
      ggplot2::aes(x = .data$x, y = .data$y, label = .data$label, angle = .data$rotation),
      hjust = 0.5,
      vjust = 0.5,
      size = title_size
    ) +
    ggplot2::geom_text(
      data = dc2,
      ggplot2::aes(x = .data$x, y = .data$y, label = .data$label),
      size = soil_size,
      hjust = 0.5,
      vjust = 0.5,
      lineheight = soil_spacing,
      show.legend = FALSE
    ) +
    ggplot2::coord_fixed(ratio = 1, xlim = xlim, ylim = ylim, expand = FALSE)
  #add annotations
  if (!is.null(fines) & !is.null(clay)) {
    plt <- ggplot_classificationtriangle_addcrosshairs_part2(
      plt,
      fines_fraction = fines,
      clay_fraction = clay,
      group = group
    )
  }

  ##############
  ### RETURN ###
  ##############

  #return
  plt
}


#' Add data to classification chart - part 1
#'
#' @description
#' Add measured data to part 1 of the soil classification chart
#'
#' @param plt a ggplot object with chart
#' @param gravel fraction of gravel (array)
#' @param sand fraction of sand (array)
#' @param fines fraction of fines (array)
#' @param group name of group
#' @param crosshairs if `TRUE`, crosshairs are plotted
#' @param legend_position position of legend. If `none`, no legend is shown
#' @param legend_title title of legend
#' @param palette RColorBrewer palette used for discrete colors
#' @param label_parse if `TRUE`, parse `group` labels
#' @importFrom rlang .data
#' @returns ggplot object
#' @export

ggplot_classificationtriangle_addcrosshairs_part1 <- function(
  plt,
  gravel,
  sand,
  fines,
  group = NULL,
  crosshairs = TRUE,
  legend_position = "none",
  legend_title = "",
  palette = "Set1",
  label_parse = TRUE
){
  #calculate positions of points
  dp <- tibble::tibble(gravel = gravel, sand = sand, fines = fines)
  dp <- dplyr::mutate(dp, convert_to_grid_part1(.data$gravel, .data$sand, .data$fines))
  #add group
  if (is.null(group)){
    dp$group <- as.character(seq(nrow(dp)))
  } else {
    dp$group <- group
  }
  #add crosshairs
  if (crosshairs == TRUE){
    #positions
    dc <- tibble::tibble(
      gravel2 = c(rep(0, nrow(dp)), gravel, gravel + sand, gravel, gravel, gravel),
      sand2 = c(sand, sand, rep(0, nrow(dp)), sand, sand + fines, sand),
      fines2 = c(fines + gravel, fines, fines, fines, rep(0, nrow(dp)), fines),
      group = rep(dp$group, 6),
      group2 = rep(seq(3), each = 2*nrow(dp))
    )
    dc <- dplyr::mutate(dc, convert_to_grid_part1(.data$gravel2, .data$sand2, .data$fines2))
    dc$grouping <- paste0(dc$group, "-", dc$group2)
    #add to plot
    plt <- plt +
      ggplot2::geom_path(
        data = dc,
        ggplot2::aes(x = .data$x, y = .data$y, group = .data$grouping, color = as.factor(.data$group))
      )
  }
  #add points to plot
  plt <- plt +
    ggplot2::geom_point(
      data = dp,
      ggplot2::aes(x = .data$x, y = .data$y, color = as.factor(.data$group), shape = as.factor(.data$group))
    ) +
    ggplot2::scale_color_brewer(name = legend_title, palette = palette) +
    ggplot2::scale_shape_discrete(name = legend_title) +
    ggplot2::theme(legend.position = legend_position)
  # add labels to plot
  if (!is.null(group)) {
    plt <- plt +
      ggplot2::geom_label(
        data = dp,
        ggplot2::aes(x = .data$x, y = .data$y, label = .data$group, color = as.factor(.data$group)),
        hjust = 0,
        vjust = 0.5,
        parse = label_parse,
        label.padding = ggplot2::unit(0.1, "lines"),
        nudge_x = 0.02,
        nudge_y = 0
      )
  }
  #return
  plt
}


#' Add data to classification chart - part 2
#'
#' @description
#' Add measured data to part 2 of the soil classification chart
#'
#' @param plt a ggplot object with chart
#' @param fines_fraction mass fraction of fines (array)
#' @param clay_fraction mass fraction of clay (array)
#' @param group name of group
#' @param crosshairs if `TRUE`, crosshairs are plotted
#' @param legend_position position of legend. If `none`, no legend is shown
#' @param legend_title title of legend
#' @param palette RColorBrewer palette used for discrete colors
#' @param label_parse if `TRUE`, parse `group` labels
#' @importFrom rlang .data
#' @returns ggplot object
#' @export

ggplot_classificationtriangle_addcrosshairs_part2 <- function(
  plt,
  fines_fraction,
  clay_fraction,
  group = NULL,
  crosshairs = TRUE,
  legend_position = "none",
  legend_title = "",
  palette = "Set1",
  label_parse = FALSE
){
  #calculate positions of points
  dp <- tibble::tibble(fines = fines_fraction, clay = clay_fraction)
  dp <- dplyr::mutate(dp, convert_to_grid_part2(.data$fines, .data$clay))
  #add group
  if (is.null(group)){
    dp$group <- as.character(seq(nrow(dp)))
  } else {
    dp$group <- group
  }
  #add crosshairs
  if (crosshairs == TRUE){
    #positions
    dc <- tibble::tibble(
      fines2 = c(rep(1, nrow(dp)), dp$fines, dp$fines, dp$fines),
      clay2 = c(dp$clay, dp$clay, rep(0, nrow(dp)), dp$clay),
      group = rep(dp$group, 4),
      group2 = rep(seq(2), each = 2*nrow(dp))
    )
    dc <- dplyr::mutate(dc, convert_to_grid_part2(.data$fines2, .data$clay2))
    dc$grouping <- paste0(dc$group, "-", dc$group2)
    #add to plot
    plt <- plt +
      ggplot2::geom_path(
        data = dc,
        ggplot2::aes(x = .data$x, y = .data$y, group = .data$grouping, color = as.factor(.data$group))
      )
  }
  #add points to plot
  plt <- plt +
    ggplot2::geom_point(
      data = dp,
      ggplot2::aes(x = .data$x, y = .data$y, color = as.factor(.data$group), shape = as.factor(.data$group))
    ) +
    ggplot2::scale_color_brewer(name = legend_title, palette = palette) +
    ggplot2::scale_shape_discrete(name = legend_title) +
    ggplot2::theme(legend.position = legend_position)
  # add labels to plot
  if (!is.null(group)) {
    plt <- plt +
      ggplot2::geom_label(
        data = dp,
        ggplot2::aes(x = .data$x, y = .data$y, label = .data$group, color = as.factor(.data$group)),
        hjust = 0,
        vjust = 1,
        parse = label_parse,
        label.padding = ggplot2::unit(0.1, "lines"),
        nudge_x = 0.01,
        nudge_y = -0.01
      )
  }
  #return
  plt
}


#' Classify soil using BS classification triangle
#'
#' @description
#' Add classification to soil using BS classification chart, based on measured
#' quantities of gravel, .data$ssand, silt and clay
#'
#' @param gravel mass of gravel
#' @param sand mass of sand
#' @param silt mass of silt
#' @param clay mass of clay
#' @param multiple if `TRUE`, soils that are on boundary of classifications may
#'   belong to either side. As a result, for a single entry multiple
#'   classifications may be returned. if `FALSE`, only one soil type is
#'   returned, based on best interpretation of code: a) primary type only if
#'   >40\%. Secondary type if >=15 and <=40\% for fines and >20\% and <=40\% for
#'   sand and gravel. If equal portions, gravel takes precedence over sand,
#'   and clay takes precendence over silt.
#' @returns tibble object with fractions and classification. Multiple
#'   classifications may be returned per input if classification lies
#'   on a boundary. The field `index` identifies the index of the input.
#' @examples
#' gravel <- c(10, 30, 50)
#' sand <- c(40, 20, 15)
#' silt <- c(30, 50, 10)
#' clay <- c(20, 0, 25)
#' classify_classificationtriangle(gravel, sand, silt, clay, multiple = FALSE)
#' classify_classificationtriangle(gravel, sand, silt, clay, multiple = TRUE)
#' @importFrom magrittr `%>%`
#' @importFrom rlang .data
#' @export

classify_classificationtriangle <- function(
  gravel,
  sand,
  silt,
  clay,
  multiple = FALSE
){
  #normalise & express in percentages for easier comparion
  total <- gravel + sand + silt + clay
  #dataframe
  df <- tibble::tibble(
    gravel = 100*gravel/total,
    sand = 100*sand/total,
    fines = 100*(silt + clay)/total,
    silt = 100*silt/total,
    clay = 100*clay/total,
  ) %>%
    dplyr::mutate(index = seq(length(total)))
  #Options
  if (multiple == FALSE){
    #soils can only belong to one side:
    #- Prioritises gravel over .data$sand when equal amounts
    #- needs >40 for major type
    #- needs >=20 to <=40 for secondary fraction (.data$sand/gravel) (or 15 for fines)
    #- prioritises clay over silt when equal amounts
    df1 <- dplyr::mutate(df,
      Sa = ((.data$fines < 15) & (.data$gravel < 20)),
      grSa = ((.data$fines < 15) & (.data$gravel >= 20) & (.data$sand > .data$gravel)),
      saGr = ((.data$fines < 15) & (.data$sand >= 20) & (.data$sand <= .data$gravel)),
      Gr = ((.data$fines < 15) & (.data$sand < 20)),
      siSa = ((.data$gravel < 20) & (.data$fines <= 40) & (.data$fines >= 15) & (.data$clay < 0.2*.data$fines)),
      clSa = ((.data$gravel < 20) & (.data$fines <= 40) & (.data$fines >= 15) & (.data$clay >= 0.2*.data$fines)),
      grsiSa = ((.data$sand > 40) & (.data$gravel >= 20) & (.data$fines >= 15) & (.data$sand > .data$gravel) & (.data$clay < 0.2*.data$fines)),
      grclSa = ((.data$sand > 40) & (.data$gravel >= 20) & (.data$fines >= 15) & (.data$sand > .data$gravel) & (.data$clay >= 0.2*.data$fines)),
      grsasiS = ((.data$sand <= 40) & (.data$sand > .data$gravel) & (.data$fines <= 40) & (.data$clay < 0.2*.data$fines)),
      grsaclS = ((.data$sand <= 40) & (.data$sand > .data$gravel) & (.data$fines <= 40) & (.data$clay >= 0.2*.data$fines)),
      sagrsiS = ((.data$gravel <= 40) & (.data$gravel >= .data$sand) & (.data$fines <= 40) & (.data$clay < 0.2*.data$fines)),
      sagrclS = ((.data$gravel <= 40) & (.data$gravel >= .data$sand) & (.data$fines <= 40) & (.data$clay >= 0.2*.data$fines)),
      sasiGr = ((.data$sand >= 20) & (.data$gravel > 40) & (.data$fines >= 15) & (.data$sand <= .data$gravel) & (.data$clay < 0.2*.data$fines)),
      saclGr = ((.data$sand >= 20) & (.data$gravel > 40) & (.data$fines >= 15) & (.data$sand <= .data$gravel) & (.data$clay >= 0.2*.data$fines)),
      siGr = ((.data$sand < 20) & (.data$fines >= 15) & (.data$fines <= 40) & (.data$clay < 0.2*.data$fines)),
      clGr = ((.data$sand < 20) & (.data$fines >= 15) & (.data$fines <= 40) & (.data$clay >= 0.2*.data$fines)),
      saSi = ((.data$sand >= 20) & (.data$gravel < 20) & (.data$fines > 40) & (.data$clay < 0.1*.data$fines)),
      saclSi = ((.data$sand >= 20) & (.data$gravel < 20) & (.data$fines > 40) & (.data$clay >= 0.1*.data$fines) & (.data$clay < 0.2*.data$fines)),
      sasiCl = ((.data$sand >= 20) & (.data$gravel < 20) & (.data$fines > 40) & (.data$clay >= 0.2*.data$fines) & (.data$clay < 0.4*.data$fines)),
      saCl = ((.data$sand >= 20) & (.data$gravel < 20) & (.data$fines > 40) & (.data$clay >= 0.4*.data$fines)),
      grsaSi = ((.data$gravel >= 20) & (.data$fines > 40) & (.data$sand > .data$gravel) & (.data$clay < 0.2*.data$fines)),
      grsaCl = ((.data$gravel >= 20) & (.data$fines > 40) & (.data$sand > .data$gravel) & (.data$clay >= 0.2*.data$fines)),
      sagrSi = ((.data$sand >= 20) & (.data$fines > 40) & (.data$gravel >= .data$sand) & (.data$clay < 0.2*.data$fines)),
      sagrCl = ((.data$sand >= 20) & (.data$fines > 40) & (.data$gravel >= .data$sand) & (.data$clay >= 0.2*.data$fines)),
      grSi = ((.data$sand < 20) & (.data$gravel >= 20) & (.data$fines > 40) & (.data$clay < 0.1*.data$fines)),
      grclSi = ((.data$sand < 20) & (.data$gravel >= 20) & (.data$fines > 40) & (.data$clay >= 0.1*.data$fines) & (.data$clay < 0.2*.data$fines)),
      grsiCl = ((.data$sand < 20) & (.data$gravel >= 20) & (.data$fines > 40) & (.data$clay >= 0.2*.data$fines) & (.data$clay < 0.4*.data$fines)),
      grCl = ((.data$sand < 20) & (.data$gravel >= 20) & (.data$fines > 40) & (.data$clay >= 0.4*.data$fines)),
      Si = ((.data$sand < 20) & (.data$gravel < 20) & (.data$clay < 0.1*.data$fines)),
      clSi = ((.data$sand < 20) & (.data$gravel < 20) & (.data$clay >= 0.1*.data$fines) & (.data$clay < 0.2*.data$fines)),
      siCl = ((.data$sand < 20) & (.data$gravel < 20) & (.data$clay >= 0.2*.data$fines) & (.data$clay < 0.4*.data$fines)),
      Cl = ((.data$sand < 20) & (.data$gravel < 20) & (.data$clay >= 0.4*.data$fines))
    )
  } else {
    #soils on boundary may belong to either side
    df1 <- dplyr::mutate(df,
      Sa = ((.data$fines <= 15) & (.data$gravel <= 20)),
      grSa = ((.data$fines <= 15) & (.data$gravel >= 20) & (.data$sand >= .data$gravel)),
      saGr = ((.data$fines <= 15) & (.data$sand >= 20) & (.data$sand <= .data$gravel)),
      Gr = ((.data$fines <= 15) & (.data$sand <= 20)),
      siSa = ((.data$gravel <= 20) & (.data$fines <= 40) & (.data$fines >= 15) & (.data$clay <= 0.2*.data$fines)),
      clSa = ((.data$gravel <= 20) & (.data$fines <= 40) & (.data$fines >= 15) & (.data$clay >= 0.2*.data$fines)),
      grsiSa = ((.data$sand >= 40) & (.data$gravel >= 20) & (.data$fines >= 15) & (.data$sand >= .data$gravel) & (.data$clay <= 0.2*.data$fines)),
      grclSa = ((.data$sand >= 40) & (.data$gravel >= 20) & (.data$fines >= 15) & (.data$sand >= .data$gravel) & (.data$clay >= 0.2*.data$fines)),
      grsasiS = ((.data$sand <= 40) & (.data$sand >= .data$gravel) & (.data$fines <= 40) & (.data$clay <= 0.2*.data$fines)),
      grsaclS = ((.data$sand <= 40) & (.data$sand >= .data$gravel) & (.data$fines <= 40) & (.data$clay >= 0.2*.data$fines)),
      sagrsiS = ((.data$gravel <= 40) & (.data$gravel >= .data$sand) & (.data$fines <= 40) & (.data$clay <= 0.2*.data$fines)),
      sagrclS = ((.data$gravel <= 40) & (.data$gravel >= .data$sand) & (.data$fines <= 40) & (.data$clay >= 0.2*.data$fines)),
      sasiGr = ((.data$sand >= 20) & (.data$gravel >= 40) & (.data$fines >= 15) & (.data$sand <= .data$gravel) & (.data$clay <= 0.2*.data$fines)),
      saclGr = ((.data$sand >= 20) & (.data$gravel >= 40) & (.data$fines >= 15) & (.data$sand <= .data$gravel) & (.data$clay >= 0.2*.data$fines)),
      siGr = ((.data$sand <= 20) & (.data$fines >= 15) & (.data$fines <= 40) & (.data$clay <= 0.2*.data$fines)),
      clGr = ((.data$sand <= 20) & (.data$fines >= 15) & (.data$fines <= 40) & (.data$clay >= 0.2*.data$fines)),
      saSi = ((.data$sand >= 20) & (.data$gravel <= 20) & (.data$fines >= 40) & (.data$clay <= 0.1*.data$fines)),
      saclSi = ((.data$sand >= 20) & (.data$gravel <= 20) & (.data$fines >= 40) & (.data$clay >= 0.1*.data$fines) & (.data$clay <= 0.2*.data$fines)),
      sasiCl = ((.data$sand >= 20) & (.data$gravel <= 20) & (.data$fines >= 40) & (.data$clay >= 0.2*.data$fines) & (.data$clay <= 0.4*.data$fines)),
      saCl = ((.data$sand >= 20) & (.data$gravel <= 20) & (.data$fines >= 40) & (.data$clay >= 0.4*.data$fines)),
      grsaSi = ((.data$gravel >= 20) & (.data$fines >= 40) & (.data$sand >= .data$gravel) & (.data$clay <= 0.2*.data$fines)),
      grsaCl = ((.data$gravel >= 20) & (.data$fines >= 40) & (.data$sand >= .data$gravel) & (.data$clay >= 0.2*.data$fines)),
      sagrSi = ((.data$sand >= 20) & (.data$fines >= 40) & (.data$gravel >= .data$sand) & (.data$clay <= 0.2*.data$fines)),
      sagrCl = ((.data$sand >= 20) & (.data$fines >= 40) & (.data$gravel >= .data$sand) & (.data$clay >= 0.2*.data$fines)),
      grSi = ((.data$sand <= 20) & (.data$gravel >= 20) & (.data$fines >= 40) & (.data$clay <= 0.1*.data$fines)),
      grclSi = ((.data$sand <= 20) & (.data$gravel >= 20) & (.data$fines >= 40) & (.data$clay >= 0.1*.data$fines) & (.data$clay <= 0.2*.data$fines)),
      grsiCl = ((.data$sand <= 20) & (.data$gravel >= 20) & (.data$fines >= 40) & (.data$clay >= 0.2*.data$fines) & (.data$clay <= 0.4*.data$fines)),
      grCl = ((.data$sand <= 20) & (.data$gravel >= 20) & (.data$fines >= 40) & (.data$clay >= 0.4*.data$fines)),
      Si = ((.data$sand <= 20) & (.data$gravel <= 20) & (.data$clay <= 0.1*.data$fines)),
      clSi = ((.data$sand <= 20) & (.data$gravel <= 20) & (.data$clay >= 0.1*.data$fines) & (.data$clay <= 0.2*.data$fines)),
      siCl = ((.data$sand <= 20) & (.data$gravel <= 20) & (.data$clay >= 0.2*.data$fines) & (.data$clay <= 0.4*.data$fines)),
      Cl = ((.data$sand <= 20) & (.data$gravel <= 20) & (.data$clay >= 0.4*.data$fines))
    )
  }
  #to long form, and only keep those that are true
  dfl <- tidyr::pivot_longer(df1,
    cols = c(-.data$gravel, -.data$sand, -.data$silt, -.data$clay, -.data$fines, -.data$index),
    values_to = "boolean",
    names_to = "soiltype"
  ) %>%
    dplyr::filter(.data$boolean == TRUE) %>%
    dplyr::select(-.data$boolean)
  #add long name, and return
  dfl %>%
    dplyr::mutate(
      soiltype_full = .data$soiltype,
      soiltype_full = stringr::str_replace(.data$soiltype_full, "gr", "gravelly "),
      soiltype_full = stringr::str_replace(.data$soiltype_full, "sa", "sandy "),
      soiltype_full = stringr::str_replace(.data$soiltype_full, "si", "silty "),
      soiltype_full = stringr::str_replace(.data$soiltype_full, "cl", "clayey "),
      soiltype_full = stringr::str_replace(.data$soiltype_full, "Gr$", "GRAVEL"),
      soiltype_full = stringr::str_replace(.data$soiltype_full, "Sa$", "CLAY"),
      soiltype_full = stringr::str_replace(.data$soiltype_full, "Si$", "CLAY"),
      soiltype_full = stringr::str_replace(.data$soiltype_full, "Cl$", "CLAY"),
      soiltype_full = stringr::str_replace(.data$soiltype_full, "S$", "soil")
    )
}
