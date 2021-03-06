#' Function to rotate coordinates
#'
#' @description
#' Simple rotation function used for plotting rotated soil element
#' and stress arrows
#'
#' @param x x-coordinate (array)
#' @param y y-coordinate (array)
#' @param theta rotation angle (array, in radians)
#' @return a tibble with rotated positions in fields `x2` and `y2`
#' @export

rotate_elements <- function(x, y, theta = 0) {
  tibble::tibble(
    x2 = x*cos(theta) - y*sin(theta),
    y2 = x*sin(theta) + y*cos(theta)
  )
}


#' Function to calculate rotated stresses
#'
#' @description
#' Function to apply stress rotation
#'
#' @param sigx normal stress on x-plane (array)
#' @param sigz normal stress on z-plane (array)
#' @param tau shear stress (array)
#' @param theta rotation angle (array, in radians)
#' @return a tibble with rotated stresses in fields `sigx2`, `sigy2` and `tau2`
#' @export

rotate_stresses <- function(sigx, sigz, tau, theta) {
  tibble::tibble(
    sigx2 = sigx*cos(theta)^2 + sigz*sin(theta)^2 + 2*tau*sin(theta)*cos(theta),
    sigz2 = sigz*cos(theta)^2 + sigx*sin(theta)^2 - 2*tau*sin(theta)*cos(theta),
    tau2 = (sigz - sigx)*sin(2*theta)/2 + tau*cos(2*theta)
  )
}


#' ggplot to plot rotated stress element
#'
#' @description
#' Function plots a rotated stress elements and stress arrows. Magnitude of
#' the arrows scales with magnitude of stress. Normal stress arrows scale by
#' the major principle stress, shear arrows scale by radius of Mohr circle.
#'
#' Unit cube of soil has side 1 by 1, and its centre is centred on x=0, y=0
#'
#' @param sigz normal stress on z-plane (scalar)
#' @param sigx normal stress on x-plane (scalar)
#' @param tau shear stress
#' @param theta rotation, in radians (scalar)
#' @param rotation_label label to use for rotation
#' @param face_label labels for x- and z-faces
#' @param arrow_offset distance between shear arrow and cube
#' @param arrow_length maximum length of arrow
#' @param palette RColorBrewer color palette to use
#' @param fill_soil color of fill of soil cube
#' @param color_soil color of outline of soil cube
#' @param coordinate_system if `TRUE`, show arrows with the coordinate system
#'   used
#' @param coordinate_arrow_length length or coordinate system arrows. Defined
#'   as fraction of `arrow_length`
#' @param stress_label if `TRUE`, plot magnitude of (rotated) stresses in top
#'   right-hand corner
#' @param stress_unit unit for stresses, e.g. `stress_unit = "kPa"`
#' @param stress_label_size size of stress labels
#' @param stress_nround number of decimals in stress labels
#' @param effective_stress if `TRUE`, stresses are assumed effective and
#'   'primes' are added to the stress labels to indicate so
#' @param clockwise_shear if `TRUE`, shear stresses are positive on
#'   positive faces of the element and when pointing in the positive
#'   direction. When `FALSE`, the default structural mechanics notation
#'   (right-hand rule for moments) is used instead, resulting in positive
#'   shear stresses pointing in the opposite direction. Function calculations
#'   assume (`clockwise_shear = TRUE`)
#' @importFrom magrittr `%>%`
#' @importFrom rlang .data
#' @return ggplot object
#' @examples
#' ggplot_stresselement(sigz = 40, sigx = 20, tau = 5, theta = pi/8)
#' @export

ggplot_stresselement <- function(
  sigz = 40,
  sigx = 20,
  tau = 10,
  theta = 0/180*pi,
  rotation_label = "theta",
  face_label = c("X", "Z"),
  arrow_offset = 0.1,
  arrow_length = 0.5,
  palette = "Set1",
  fill_soil = "#d3bc5f",
  color_soil = "#65571d",
  coordinate_system = TRUE,
  coordinate_arrow_length = 0.5,
  stress_label = TRUE,
  stress_unit = "kPa",
  stress_label_size = 3,
  stress_nround = 1,
  effective_stress = FALSE,
  clockwise_shear = FALSE
){
  #shear direction
  if (clockwise_shear == TRUE) {
    dir <- 1
  } else {
    dir <- -1
    tau <- -tau
  }
  #stress invariants
  p <- 0.5*(sigx + sigz)
  q <- sqrt((0.5*sigx - 0.5*sigz)^2 + tau^2)
  #max principal stresses (used to scale arrows)
  sig1 <- p + q
  #values of stresses in rotated state
  df <- rotate_stresses(sigx, sigz, tau, theta)
  #coordinates of unrotated soil element
  ds <- tibble::tibble(
    x = 0.5*c(-1, -1, 1, 1),
    y = 0.5*c(-1, 1, 1, -1)
  ) %>%
    dplyr::mutate(rotate_elements(.data$x, .data$y, theta))
  #dataframe with arrows
  da1 <- tibble::tibble(
    type = c(
      rep("sig", 2),
      rep("sig", 2),
      rep("tau", 2),
      rep("tau", 2)
    ),
    side = c(
      rep("x", 2),
      rep("z", 2),
      rep("x", 2),
      rep("z", 2)
    ),
    x = c(
      0.5 + 2*arrow_offset + 0.5*arrow_length*(abs(df$sigx2/sig1) + c(1, -1)*df$sigx2/abs(sig1)) ,
      c(0, 0),
      0.5 + arrow_offset + c(0, 0),
      c(0.5, -0.5)*df$tau2/q*arrow_length
    ),
    y = c(
      c(0, 0),
      0.5 + 2*arrow_offset + 0.5*arrow_length*(abs(df$sigz2/sig1) + c(1, -1)*df$sigz2/abs(sig1)),
      c(0.5, -0.5)*df$tau2/q*arrow_length,
      0.5 + arrow_offset + c(0, 0)
    ),
    plane = "positive"
  )
  #remove shear arrows if shear stress = 0
  if (dplyr::near(df$tau2, 0)){
    da1 <- dplyr::filter(da1, .data$type == "sig")
  }
  #all arrows
  da <- rbind(
    da1,
    dplyr::mutate(da1,
      plane = "negative",
      x = -.data$x,
      y = -.data$y
    )
  ) %>% dplyr::mutate(
    plotgroup = paste0(.data$type, .data$side, "-", .data$plane),
    rotate_elements(.data$x, .data$y, theta)
  )
  #letters indicating faces
  dlet <- tibble::tibble(
    x = c(0.4, 0),
    y = c(0, 0.4),
    label = face_label,
    side = c("x", "z")
  ) %>%
    dplyr::mutate(rotate_elements(.data$x, .data$y, theta))
  if (!dplyr::near(theta, 0)){
    dlet$label <- paste(dlet$label, "'")
  }
  #ggplot
  plt <- ggplot2::ggplot() +
    ggplot2::theme_void() +
    ggplot2::geom_polygon(
      data = ds,
      ggplot2::aes(x = .data$x2, y = .data$y2),
      color = color_soil,
      fill = fill_soil
    ) +
    ggplot2::geom_path(
      data = da,
      ggplot2::aes(x = .data$x2, y = .data$y2, group = .data$plotgroup, color = .data$side),
      arrow = ggplot2::arrow(length = ggplot2::unit(0.10, "inches"))
    ) +
    ggplot2::geom_text(
      data = dlet,
      ggplot2::aes(x = .data$x2, y = .data$y2, label = .data$label, color = .data$side),
      hjust = 0.5,
      vjust = 0.5
    ) +
    ggplot2::coord_fixed(
      ratio = 1,
      xlim = c(-1, 1)*(0.5 + 2*arrow_offset + arrow_length),
      ylim = c(-1, 1)*(0.5 + 2*arrow_offset + arrow_length),
      expand = TRUE
    ) +
    ggplot2::scale_color_brewer(palette = palette) +
    ggplot2::theme(legend.position = "none")
  #add rotation arrow
  if (!dplyr::near(theta, 0)){
    darr <- tibble::tibble(
      x = (0.5 + 2*arrow_offset + 0.5*arrow_length)*cos(seq(0, theta, l = 91)),
      y = (0.5 + 2*arrow_offset + 0.5*arrow_length)*sin(seq(0, theta, l = 91)),
      side = "zz"
    )
    plt <- plt +
      ggplot2::geom_text(
        ggplot2::aes(
          x = (0.5 + 2*arrow_offset + 0.7*arrow_length)*cos(0.5*theta),
          y = (0.5 + 2*arrow_offset + 0.7*arrow_length)*sin(0.5*theta),
          label = rotation_label,
          color = "zz"
        ),
        hjust = 0.5,
        vjust = 0.5,
        parse = TRUE
      ) +
      ggplot2::geom_path(
        data = darr,
        ggplot2::aes(x = .data$x, y = .data$y, color = .data$side),
        arrow = ggplot2::arrow(length = ggplot2::unit(0.5, "lines"))
      )
  }
  #add coordinate system
  if (coordinate_system == TRUE){
    plt <- plt +
      ggplot2::annotate(
        "segment",
        x = -(0.5 + arrow_length + arrow_offset),
        y = (0.5 + arrow_length + arrow_offset),
        xend = -(0.5 + arrow_length + arrow_offset) + coordinate_arrow_length*arrow_length,
        yend = (0.5 + arrow_length + arrow_offset),
        arrow = ggplot2::arrow(length = ggplot2::unit(0.5, "lines")),
        color = "black"
      ) +
      ggplot2::annotate(
        "segment",
        x = -(0.5 + arrow_length + arrow_offset),
        y = (0.5 + arrow_length + arrow_offset),
        xend = -(0.5 + arrow_length + arrow_offset),
        yend = (0.5 + arrow_length + arrow_offset) - coordinate_arrow_length*arrow_length,
        arrow = ggplot2::arrow(length = ggplot2::unit(0.5, "lines")),
        color = "black"
      ) +
      ggplot2::annotate(
        "text",
        x = -(0.5 + arrow_length + arrow_offset) + coordinate_arrow_length*(arrow_length + arrow_offset),
        y = (0.5 + arrow_length + arrow_offset),
        label = "x",
        color = "black",
        hjust = 0,
        vjust = 0.5
      ) +
      ggplot2::annotate(
        "text",
        x = -(0.5 + arrow_length + arrow_offset),
        y = (0.5 + arrow_length + arrow_offset) - coordinate_arrow_length*(arrow_length + arrow_offset),
        label = "z",
        color = "black",
        hjust = 0.5,
        vjust = 1
      )
  }
  #plot magnitude of stresses in top right-hand corner
  if (stress_label == TRUE) {
    #labels
    if (effective_stress == TRUE) {
      dlab <- tibble::tibble(label = c(
        paste0("sigma*minute[z]==", round(df$sigz2, stress_nround), "~", stress_unit),
        paste0("sigma*minute[x]==", round(df$sigx2, stress_nround), "~", stress_unit),
        paste0("tau*minute[z]==", round(dir*df$tau2, stress_nround), "~", stress_unit)
      ))
    } else {
      dlab <- tibble::tibble(label = c(
        paste0("sigma[z]==", round(df$sigz2, stress_nround), "~", stress_unit),
        paste0("sigma[x]==", round(df$sigx2, stress_nround), "~", stress_unit),
        paste0("tau[xz]==", round(dir*df$tau2, stress_nround), "~", stress_unit)
      ))
    }
    #positions
    dlab$x <- 0.99*(0.5 + 2*arrow_offset + arrow_length)
    dlab$y <- c(0.99, 0.84, 0.69)*(0.5 + 2*arrow_offset + arrow_length)
    #add to plot
    plt <- plt + ggplot2::geom_text(
      data = dlab,
      ggplot2::aes(x = .data$x, y = .data$y, label = .data$label),
      hjust = 1,
      vjust = 1,
      size = stress_label_size,
      parse = TRUE
    )
  }
  #return
  return(plt)
}


#' ggplot to plot Mohr circle of rotated stress element
#'
#' @description
#' Function plots a Mohr circle for a rotated stress element and indicates the
#' pole and the rotation angles
#'
#' @inheritParams ggplot_stresselement
#' @param pole plot the pole point, label and lines if `pole = TRUE`. If
#'   `pole = FALSE`, do not plot anything related to the pole
#' @param pole_label label to plot at pole
#' @param n_circle number of points to use for drawing Mohr circle
#' @param color_circle color of the Mohr circle
#' @param color_lines color of lines crossing midpoint of circle
#' @param double_angle if `TRUE`, an indicator for twice the applied rotation in
#'   is plotted in the Mohr circle. If `FALSE`, this is omitted
#' @param effective_stress if `TRUE`, effective stress label is plotted on the
#'   x-axis. if `FALSE`, total stress
#' @param xlim,ylim user defined min and max x and y axis limits. If not defined,
#'   the are automatically chosen
#' @importFrom magrittr `%>%`
#' @return ggplot object
#' @examples
#' ggplot_mohrcircle(sigz = 40, sigx = 20, tau = 5, theta = pi/8)
#' @export

ggplot_mohrcircle <- function(
  sigz = 40,
  sigx = 20,
  tau = 10,
  theta = 0,
  rotation_label = "theta",
  face_label = c("X", "Z"),
  palette = "Set1",
  pole = TRUE,
  pole_label = "Pole",
  n_circle = 181,
  color_circle = "black",
  color_lines = "grey50",
  double_angle = FALSE,
  effective_stress = FALSE,
  clockwise_shear = FALSE,
  xlim = c(0, NA),
  ylim = c(NA, NA)
){
  #shear direction
  if (clockwise_shear == TRUE) {
    dir <- 1
  } else {
    dir <- -1
    tau <- -tau
  }
  #stress invariants
  p <- 0.5*(sigx + sigz)
  q <- sqrt((0.5*sigx - 0.5*sigz)^2 + tau^2)
  #rotated stresses
  df <- rotate_stresses(sigx, sigz, tau, theta)
  #circle coordinates
  dc <- tibble::tibble(
    x = p + q*cos(seq(0, 2*pi, l = n_circle)),
    y = q*sin(seq(0, 2*pi, l = n_circle))
  )
  #point coordinates - undisplaced
  dp0 <- tibble::tibble(
    x = c(sigx, sigz),
    y = c(-tau, tau),
    side = c("x", "z"),
    label = face_label
  )
  #point coordinates - displaced
  dp <- tibble::tibble(
    x = c(df$sigx2, df$sigz2),
    y = c(-df$tau2, df$tau2),
    side = c("x", "z"),
    label = paste0(face_label, "'")
  )
  #lines of unrotated faces
  dl0 <- tibble::tibble(
    x = c(sigx, sigx, sigx, sigz),
    y = c(-tau, tau, tau, tau),
    side = c("x", "x", "z", "z")
  )
  #lines of rotated faces
  dl <- tibble::tibble(
    x = c(sigx, df$sigx2, sigx, df$sigz2),
    y = c(tau, -df$tau2, tau, df$tau2),
    side = c("x", "x", "z", "z")
  )
  #rotation arrows
  theta0 <- atan2(tau, sigz - 0.5*(sigz + sigx))
  drot1 <- tibble::tibble(
    x = sigx + 0.6*(df$sigz2 - sigx)*cos(seq(0, theta, l = 91)),
    y = tau + 0.6*(df$sigz2 - sigx)*sin(seq(0, theta, l = 91)),
    side = "zz"
  )
  drot2 <- tibble::tibble(
    x = p + 0.8*q*cos(seq(theta0, theta0 + 2*theta, l = 91)),
    y = 0.8*q*sin(seq(theta0, theta0 + 2*theta, l = 91)),
    side = "zz"
  )
  drotlabel <- tibble::tibble(
    x = c(
      sigx + 0.5*(df$sigz2 - sigx)*cos(0.5*theta),
      p + 0.7*q*cos(theta0 + theta)
    ),
    y = c(
      tau + 0.5*(df$sigz2 - sigx)*sin(0.5*theta),
      0.7*q*sin(theta0 + theta)
    ),
    label = c(rotation_label, paste0("2*", rotation_label)),
    side = "zz"
  )
  #plot
  plt <- ggplot2::ggplot() +
    ggplot2::theme_bw() +
    ggplot2::geom_hline(yintercept = 0) +
    ggplot2::annotate("segment", x = sigx, xend = sigz, y = -tau, yend = tau, color = color_lines, linetype = 2) +
    ggplot2::annotate("segment", x = df$sigx2, xend = df$sigz2, y = -df$tau2, yend = df$tau2, color = color_lines) +
    ggplot2::geom_path(
      data = dc,
      ggplot2::aes(x = .data$x, y = .data$y),
      color = color_circle
    )
  if ((pole == TRUE) | !dplyr::near(theta, 0)) {
    plt <- plt + ggplot2::geom_path(
      data = dl0,
      ggplot2::aes(x = .data$x, y = .data$y, color = .data$side),
      linetype = 2
    )
  }
  plt <- plt +
    ggplot2::geom_point(
      data = dp0,
      ggplot2::aes(x = .data$x, y = .data$y, color = .data$side)
    ) +
    ggplot2::geom_text(
      data = dp0,
      ggplot2::aes(x = .data$x, y = .data$y, color = .data$side, label = .data$label, hjust = as.double(p >= .data$x), vjust = as.double(0 >= .data$y)),
      parse = FALSE
    )
  if (!dplyr::near(theta, 0)) {
    plt <- plt +
      ggplot2::geom_path(
        data = dl,
        ggplot2::aes(x = .data$x, y = .data$y, color = .data$side)
      ) +
      ggplot2::geom_point(
        data = dp,
        ggplot2::aes(x = .data$x, y = .data$y, color = .data$side)
      ) +
      ggplot2::geom_text(
        data = dp,
        ggplot2::aes(x = .data$x, y = .data$y, color = .data$side, label = .data$label, hjust = as.double(p >= .data$x), vjust = as.double(0 >= .data$y)),
        parse = FALSE
      ) +
      ggplot2::geom_path(
        data = drot1,
        ggplot2::aes(x = .data$x, y = .data$y, color = .data$side),
        arrow = ggplot2::arrow(length = ggplot2::unit(0.5, "lines"))
      )
    if (double_angle == TRUE) {
      plt <- plt +
        ggplot2::geom_path(
          data = drot2,
          ggplot2::aes(x = .data$x, y = .data$y, color = .data$side),
          arrow = ggplot2::arrow(length = ggplot2::unit(0.5, "lines"))
        ) +
        ggplot2::geom_text(
          data = drotlabel,
          ggplot2::aes(x = .data$x, y = .data$y, label = .data$label, color = .data$side),
          hjust = 0.5,
          vjust = 0.5,
          parse = TRUE
        )
    } else {
      plt <- plt + ggplot2::geom_text(
        data = drotlabel[1,],
        ggplot2::aes(x = .data$x, y = .data$y, label = .data$label, color = .data$side),
        hjust = 0.5,
        vjust = 0.5,
        parse = TRUE
      )
    }
  }
  plt <- plt +
    ggplot2::annotate("point", x = p, y = 0) +
    ggplot2::coord_fixed(
      xlim = round_limits(p + q, lower = xlim[1], upper = xlim[2]),
      ylim = round_limits(c(-1.2*q, 1.2*q), lower = ylim[1], upper = ylim[2]),
      ratio = 1,
      expand = FALSE
    ) +
    ggplot2::scale_color_brewer(palette = palette) +
    ggplot2::theme(legend.position = "none") +
    ggplot2::ylab(expression(tau~"[kPa]"))
  #add pole point and text
  if (pole == TRUE) {
    plt <- plt +
      ggplot2::annotate("point", x = sigx, y = tau) +
      ggplot2::annotate(
        "text",
        x = sigx,
        y = tau,
        label = pole_label,
        hjust = as.double(p >= sigx),
        vjust = as.double(0 >= tau)
      )
  }
  #x-label: effective stress or not
  if (effective_stress == TRUE) {
    plt <- plt + ggplot2::xlab(expression(sigma*"'"~"[kPa]"))
  } else {
    plt <- plt + ggplot2::xlab(expression(sigma~"[kPa]"))
  }
  #return
  return(plt)
}
