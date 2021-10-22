#' Generate MCC elasto-plastic matrix elements
#'
#' @description
#' Returns the matrix elements of the elasto-plastic stiffness matrix
#' according to the modified cam-clay model, in terms of isotropic and
#' deviatoric stress invariants p and q, and volumetric and deviatoric
#' strain invariants epsv and epq (epsv is assumed positive in
#' compression)
#'
#' dp = M11*depsv + M12*depsq
#' dq = M21*depsv + M22*depsq
#'
#' @param p,q current stress invariants
#' @param pc current preconsolidation pressure
#' @param v0 initial specific volume
#' @param M M-parameter
#' @param kappa,lambda MCC compression parameters
#' @param nu Poisson's ratio
#' @return array with four stiffness elements [M11, M12, M21, M22]
#' @export

mcc_matrix_elastoplastic <- function(p, q, pc, v0, M, kappa, lambda, nu) {
  #elastic stiffnesses
  p_ev <- v0*p/kappa
  q_eq <- p_ev*(3*(1 - 2*nu))/(2*(1 + nu))
  #derivatives of yield surface
  f_p <- M^2*(2*p - pc)
  f_q <- 2*q
  f_pc <- -M^2*p
  #flow potential
  g_p <- f_p
  g_q <- f_q
  #derivative of preconsolidation pressure
  pc_evp <- v0*pc/(lambda - kappa)
  #temp - plastic multiplier - term underneath division
  tmp <- f_p*p_ev*g_p + f_q*q_eq*g_q - f_pc*pc_evp*g_p
  #matrix elements
  M11 <- p_ev - p_ev*f_p*p_ev*g_p/tmp
  M12 <- 0 - p_ev*f_q*q_eq*g_p/tmp
  M21 <- 0 - q_eq*f_p*p_ev*g_q/tmp
  M22 <- q_eq - q_eq*f_q*q_eq*g_q/tmp
  #return
  return(c(M11, M12, M21, M22))
}


#' Solve incremental plastic MCC response - drained triaxial compression
#'
#' @description
#' Solve the incremental response for a drained triaxial compression
#' test according to MCC theory, when behaviour is plastic
#'
#' @inheritParams mcc_matrix_elastoplastic
#' @param de1 axial strain increment
#' @return array with increment in radial strain and increment in axial stress
#' @export

mcc_triaxcomp_drained_plastic <- function(de1, p, q, pc, v0, M, kappa, lambda, nu) {
  #matrix elements
  M <- mcc_matrix_elastoplastic(p, q, pc, v0, M, kappa, lambda, nu)
  #solve system of equations
  de3 <- -(9*M[1] + 6*M[2] - 3*M[3] - 2*M[4])/(2*(9*M[1] - 3*M[2] - 3*M[3] + M[4]))*de1
  ds1 <- 9*(M[1]*M[4] - M[2]*M[3])/(9*M[1] - 3*M[2] - 3*M[3] + M[4])*de1
  #return
  return(c(de3, ds1))
}


#' Solve incremental plastic MCC response - undrained triaxial compression
#'
#' @description
#' Solve the incremental response for a undrained triaxial compression
#' test according to MCC theory, when behaviour is plastic
#'
#' @inheritParams mcc_matrix_elastoplastic
#' @param de1 axial strain increment
#' @return array with increment in p and q
#' @export

mcc_triaxcomp_undrained_plastic <- function(de1, p, q, pc, v0, M, kappa, lambda, nu) {
  #matrix elements
  M <- mcc_matrix_elastoplastic(p, q, pc, v0, M, kappa, lambda, nu)
  #solve system of equations (dev = 0, deq = deq)
  dp <- M[2]*de1
  dq <- M[4]*de1
  #return
  return(c(dp, dq))
}


#' Predict full drained triaxial compression test - Modified Cam-Clay
#'
#' @description
#' Predict the full MCC behaviour for a drained triaxial compression
#' test
#'
#' @param p0 initial isotropic effective stress
#' @param pc initial preconsolidation pressure
#' @param M M-parameter
#' @param Gamma CSL specific volume intercept
#' @param kappa,lambda MCC compression parameters
#' @param nu Poisson's ratio
#' @param e1 values of axial strain at which to calculate response
#' @return a tibble with stress and strain properties at each level of
#'   axial strain `e1`
#' @importFrom magrittr `%>%`
#' @examples
#' mcc_solve_triax_drained()
#' @export

mcc_solve_triax_drained <- function(
  p0 = 10,
  pc = 100,
  M = 1.35,
  Gamma = 2.00,
  lambda = 0.115,
  kappa = 0.015,
  nu = 0.3,
  e1 = seq(0, 0.20, l = 101)
){
  ## INITIALISE
  #assign output
  df <- tibble::tibble(
    e1 = e1,
    e3 = 0,
    p = p0,
    q = 0,
    pc = pc,
    plastic = FALSE
  )
  #cell pressure
  s3 <- p0
  #initial specific volume
  v0 <- Gamma + log(2)*(lambda - kappa) - lambda*log(pc) + kappa*log(pc/p0)
  #ratio K and G
  alpha <- 3*(1 - 2*nu)/(2*(1 + nu))

  ## ELASTIC
  #strain to plastic yield
  a <- (9 + M^2)*p0^2
  b <- -(18*s3 + M^2*pc)*p0
  c <- 9*s3^2
  zeta <- (-b + sqrt(b^2 - 4*a*c))/(2*a)
  e1y <- log(zeta)*(3 + alpha)*kappa/(3*v0*alpha)
  #parameters at yielding
  py <- p0*zeta
  qy <- 3*p0*zeta - 3*s3
  e3y <- -e1y*0.5*nu
  #assign elastic solutions
  ie <- (e1 <= e1y)
  zeta <- exp(3*v0*alpha*e1[ie]/((3 + alpha)*kappa))
  df$p[ie] <- p0*zeta
  df$q[ie] <- 3*p0*zeta - 3*s3
  df$e3[ie] <- -e1[ie]*0.5*nu
  df$plastic[!ie] <- TRUE

  ## PLASTIC
  #plastic steps
  ip <- which(ie == FALSE)
  #loop through plastic steps
  for (i in ip) {
    if (i == ip[1]) {
      #first plastic step - continue from known yield surface
      sol <- mcc_triaxcomp_drained_plastic(
        e1[i] - e1y, py, qy, pc,
        v0, M, kappa, lambda, nu
      )
      df$e3[i] <- e3y + sol[1]
      df$p[i] <- py + 1/3*sol[2]
      df$q[i] <- qy + sol[2]
    } else {
      #continue from previous plastic step
      sol <- mcc_triaxcomp_drained_plastic(
        e1[i] - e1[i-1], df$p[i - 1], df$q[i - 1], df$pc[i - 1],
        v0, M, kappa, lambda, nu
      )
      df$e3[i] <- df$e3[i - 1] + sol[1]
      df$p[i] <- df$p[i - 1] + 1/3*sol[2]
      df$q[i] <- df$q[i - 1] + sol[2]
    }
    df$pc[i] <- (df$q[i]^2 + M^2*df$p[i]^2)/(M^2*df$p[i])
  }

  ## RETURN
  #add yield point and sort
  if (max(e1) >= e1y) {
    df <- dplyr::bind_rows(
      df,
      tibble::tibble(e1 = e1y, e3 = e3y, p = py, q = qy, pc = pc, plastic = TRUE)
    ) %>%
      dplyr::arrange(e1)
  }
  #calculate volumetric strain (negative in compression)
  df$ev <- -(df$e1 + 2*df$e3)
  #specific volume
  df$v <- v0*(1 + df$ev)
  #return
  return(df)
}


#' Predict full undrained triaxial compression test - Modified Cam-Clay
#'
#' @description
#' Predict the full MCC behaviour for a undrained triaxial compression
#' test
#'
#' @inheritParams mcc_solve_triax_drained
#' @return a tibble with stress and strain properties at each level of
#'   axial strain `e1`
#' @importFrom magrittr `%>%`
#' @examples
#' mcc_solve_triax_undrained()
#' @export

mcc_solve_triax_undrained <- function(
  p0 = 10,
  pc = 100,
  M = 1.35,
  Gamma = 2.00,
  lambda = 0.115,
  kappa = 0.015,
  nu = 0.3,
  e1 = seq(0, 0.20, l = 101)
){
  ## INITIALISE
  #assign output
  df <- tibble::tibble(
    e1 = e1,
    e3 = 0,
    p = p0,
    q = 0,
    pc = pc,
    plastic = FALSE
  )
  #cell pressure
  s3 <- p0
  #initial specific volume
  v0 <- Gamma + log(2)*(lambda - kappa) - lambda*log(pc) + kappa*log(pc/p0)
  #ratio G over K
  alpha <- 3*(1 - 2*nu)/(2*(1 + nu))

  ## ELASTIC
  #strain to plastic yield
  e1y <- sqrt(M^2*(pc - p0)*kappa^2/(9*alpha^2*v0^2*p0))
  py <- p0
  qy <- 3*alpha*v0*p0/kappa*e1y
  e3y <- -0.5*e1y
  #assign elastic solutions
  ie <- (e1 <= e1y)
  df$p[ie] <- p0
  df$q[ie] <- 3*alpha*v0*p0/kappa*df$e1[ie]
  df$plastic[!ie] <- TRUE
  #radial strains
  df$e3 <- -0.5*e1

  ## PLASTIC
  #plastic steps
  ip <- which(ie == FALSE)
  #loop through plastic steps
  for (i in ip) {
    if (i == ip[1]) {
      #first plastic step - continue from known yield surface
      sol <- mcc_triaxcomp_undrained_plastic(
        e1[i] - e1y, py, qy, pc,
        v0, M, kappa, lambda, nu
      )
      df$p[i] <- py + sol[1]
      df$q[i] <- qy + sol[2]
    } else {
      #continue from previous plastic step
      sol <- mcc_triaxcomp_undrained_plastic(
        e1[i] - e1[i-1], df$p[i - 1], df$q[i - 1], df$pc[i - 1],
        v0, M, kappa, lambda, nu
      )
      df$p[i] <- df$p[i - 1] + sol[1]
      df$q[i] <- df$q[i - 1] + sol[2]
    }
    df$pc[i] <- (df$q[i]^2 + M^2*df$p[i]^2)/(M^2*df$p[i])
  }

  ## RETURN
  #add yield point and sort
  if (max(e1) >= e1y) {
    df <- dplyr::bind_rows(
      df,
      tibble::tibble(e1 = e1y, e3 = e3y, p = py, q = qy, pc = pc, plastic = TRUE)
    ) %>%
      dplyr::arrange(e1)
  }
  #calculate volumetric strain (negative in compression)
  df$ev <- -(df$e1 + 2*df$e3)
  #specific volume
  df$v <- v0*(1 + df$ev)
  #return
  return(df)
}
