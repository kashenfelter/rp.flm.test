

#' @title fdata.cfs.2003.a
#' 
#' @description The function fdata.cfs.2003.a() generates simulated data of the functional linear linear model with
# scalar response, Cardot, Ferraty and Sarda (2003).
#'  The function data.cfs.2003.b() generates simulated data of the functional linear linear model with
#' scalar response, Cardot, Ferraty and Sarda (2003).
#' The function data.cfs.2003.c() generates simulated data of the functional linear model with scalar
#' response, Cardot, Ferraty and Sarda (2003), but with the inclusion of different eigenfunctions.
#' 
#' @param n TODO
#' @param t TODO
#' @param b TODO
#' @param beta TODO
#' @param nc TODO
#' @param R2 TODO
#' @return TODO
#' @example 
#' # TODO
#' @author Manuel Febrero-Bande (\email{manuel.febrero@@usc.es}) and Eduardo Garcia-Portugues (\email{edgarcia@@est-econ.uc3m.es}).
#' @export 
fdata.cfs.2003.a <- function(n = 100, t = seq(0, 1, len = 101), 
                             b = c(2, 4, 5) / sqrt(2), nc = length(b), R2 = 0.95) {

  x <- fda.usc::rproc2fdata(n = n, t = t, sigma = "brownian")
  rtt <- x$rangeval
  X <- x$data
  J <- length(t)
  if (is.null(nc))
    stop("No vector b")
  v <- fda.usc::fdata(sqrt(2) * sin(t * pi * 0.5), argvals = t)
  for (i in 2:length(b)) {
    v <- c(v, fda.usc::fdata(sqrt(2) * sin(t * pi * (i - 0.5)), argvals = t))
  }
  lamb <- 1/(pi * (1:length(b) - 0.5))^2

  # bet <- b1 * v1 + b2 * v2 + b3 * v3
  bet <- fda.usc::fdata(matrix(b, nrow = 1) %*% v$data, argvals = t)
  varXg <- sum(lamb * b^2)
  # varXg <- b1^2 * lamb1 + b2^2 * lamb2 + b3^2 * lamb3
  s2eps <- varXg * (1/R2 - 1)
  # Xmean <- fdata.cen(x) yp <- (Xmean$Xcen$data %*% bet) * (1 / J)
  yp <- fda.usc::inprod.fdata(x, bet)
  e <- matrix(rnorm(n, 0, sqrt(s2eps)), ncol = 1)
  y <- yp + e

  return(list(x = x, y = y, v = v, beta = bet, b = b, s2eps = s2eps))

}


#' @rdname fdata.cfs.2003.a
#' @export
fdata.cfs.2003.b <- function(n = 100, t = seq(0, 1, len = 101), beta = NULL, R2 = 0.95) {

  x <- fda.usc::rproc2fdata(n = n, t = t, sigma = "brownian")
  rtt <- x$rangeval
  X <- x$data
  J <- length(t)
  maxK <- qr(X)$rank
  varXg <- 0
  if (is.null(beta) | !fda.usc::is.fdata(beta)) {
    bet <- fda.usc::fdata(log(15 * t^2 + 10) + cos(4 * pi * t), argvals = t)
  }

  for (k in 1:maxK) {
    integrand <- function(t) {
      sqrt(2) * (log(15 * t^2 + 10) + cos(4 * pi * t)) * sin((k - 0.5) * pi * t)
    }
    bk <- integrate(integrand, lower = rtt[1], upper = rtt[2])$value
    varXg <- varXg + bk^2/((k - 0.5) * pi)^2
  }
  s2eps <- varXg * (1/R2 - 1)
  yp <- fda.usc::inprod.fdata(x, bet)
  e <- matrix(rnorm(n, 0, sqrt(s2eps)), ncol = 1)
  y <- yp + e

  return(list(x = x, y = y, beta = bet, s2eps = s2eps))

}


#' @rdname fdata.cfs.2003.c
#' @export
fdata.cfs.2003.c <- function(n = 100, t = seq(0, 1, len = 101), b = c(0, 0, 2, 0, 4, 0, 5) / sqrt(2), 
                             nc = length(b), R2 = 0.95) {

  x <- fda.usc::rproc2fdata(n = n, t = t, sigma = "brownian")
  rtt <- x$rangeval
  X <- x$data
  J <- length(t)
  if (is.null(nc))
    stop("No vector b")
  v <- fda.usc::fdata(sqrt(2) * sin(t * pi * 0.5), argvals = t)
  for (i in 2:length(b)) {
    v <- c(v, fda.usc::fdata(sqrt(2) * sin(t * pi * (i - 0.5)), argvals = t))
  }
  lamb <- 1/(pi * (1:length(b) - 0.5))^2

  bet <- fda.usc::fdata(matrix(b, nrow = 1) %*% v$data, argvals = t)
  varXg <- sum(lamb * b^2)
  s2eps <- varXg * (1/R2 - 1)
  yp <- fda.usc::inprod.fdata(x, bet)
  e <- matrix(rnorm(n, 0, sqrt(s2eps)), ncol = 1)
  y <- yp + e

  return(list(x = x, y = y, v = v, beta = bet, b = b, s2eps = s2eps))

}


#' @title fdata.cfs.2003.a
#' @description The function data.hh.2006 generates simulated data of the functional linear model with scalar
# response from Hall and Housseini-Nasab (2006).  imod=1,2,3 --> model (i), model (ii), model (iii)
#' @param n TODO
#' @param t TODO
#' @param imod TODO
#' @param ncb TODO
#' @param R2 TODO
#' @return TODO
#' @example 
#' # TODO
#' @author Manuel Febrero-Bande (\email{manuel.febrero@@usc.es}) and Eduardo Garcia-Portugues (\email{edgarcia@@est-econ.uc3m.es}).
#' @export 
fdata.hh.2006 <- function(n = 100, t = seq(0, 1, len = 101), imod = 1, ncb = 20, R2 = 0.95) {

  v <- fda.usc::fdata(sqrt(2) * cos(pi * t), argvals = t)
  for (i in 2:ncb) {
    v <- c(v, fda.usc::fdata(sqrt(2) * cos(i * pi * t), argvals = t))
  }
  bet <- fda.usc::fdata(pi^2 * (t^2 - 1/3), argvals = t)
  sdcoef <- (1:ncb)^(-imod)
  bk <- (-1)^(1:ncb) * (1:ncb)^(-2) * 2^(3/2)
  varXg <- sum(bk^2/((1:ncb)^(2 * imod)))
  coefsim <- matrix(NA, ncol = ncb, nrow = n)
  for (i in 1:ncb) {
    coefsim[, i] <- rnorm(n, mean = 0, sd = sdcoef[i])
  }
  s2eps <- varXg * (1/R2 - 1)
  x <- fda.usc::fdata(coefsim %*% v$data, argvals = t)
  yp <- fda.usc::inprod.fdata(x, bet)
  e <- matrix(rnorm(n, 0, sqrt(s2eps)), ncol = 1)
  y <- yp + e

  return(list(x = x, y = y, v = v, beta = bet, s2eps = s2eps))

}


#' @title Geometric brownian motion
#' @export
GBM <- function(n = 100, t = seq(0, 1, len = 101), S0 = 1) {
  
  S0 * exp(fda.usc::rproc2fdata(n = n, t = t, sigma = "brownian"))

}


#' @title fdata.bridge
#' @description
#' Brownian Motion: eigenvectors: sqrt(2)*sin((k-.5)*pi*t) --- Lambda: 1/((k-.5)^2*pi^2) Brownian
#' Bridge: eigenvectors: sqrt(2)*sin(k*pi*t) --- Lambda: 1/(k^2*pi^2)
#' # Make Bridgebrownian to have a certain R2
#' @param n TODO
#' @param t TODO
#' @param imod TODO
#' @param ncb TODO
#' @param R2 TODO
#' @return TODO
#' @example 
#' # TODO
#' @author Manuel Febrero-Bande (\email{manuel.febrero@@usc.es}) and Eduardo Garcia-Portugues (\email{edgarcia@@est-econ.uc3m.es}).
#' @export 
fdata.bridge <- function(n = 100, t = seq(0, 1, len = 101), b = c(2, 4, 5) / sqrt(2), 
                         nc = length(b), R2 = 0.95) {

  x <- Bridgebrownian(n = n, t = t)
  rtt <- x$rangeval
  X <- x$data
  J <- length(t)
  if (is.null(nc))
    stop("No vector b")
  v <- fda.usc::fdata(sqrt(2) * sin(t * pi), argvals = t)
  for (i in 2:length(b)) {
    v <- c(v, fda.usc::fdata(sqrt(2) * sin(t * pi * i), argvals = t))
  }
  lamb <- 1/(pi * (1:length(b)))^2

  # bet <- b1 * v1 + b2 * v2 + b3 * v3
  bet <- fda.usc::fdata(matrix(b, nrow = 1) %*% v$data, argvals = t)
  varXg <- sum(lamb * b^2)
  # varXg <- b1^2 * lamb1 + b2^2 * lamb2 + b3^2 * lamb3
  s2eps <- varXg * (1/R2 - 1)
  # Xmean <- fdata.cen(x) yp <- (Xmean$Xcen$data %*% bet) * (1 / J)
  yp <- fda.usc::inprod.fdata(x, bet)
  e <- matrix(rnorm(n, 0, sqrt(s2eps)), ncol = 1)
  y <- yp + e

  return(list(x = x, y = y, v = v, beta = bet, b = b, s2eps = s2eps))

}


#' @rdname fdata.bridge
#' @export
Bridgebrownian <- function(n = 100, t = seq(0, 1, len = 101)) {
  
  x <- fda.usc::rproc2fdata(n = n, t = t, sigma = "brownian")
  for (i in 1:n) {
    
    x$data[i, ] <- x$data[i, ] - t * x$data[i, length(t)]
    
  }
  
  return(x)
  
}


#' @title Quadratic model
#' @export
quadratic <- function(x, b) {

  if (is.matrix(x)) {
    nc <- ncol(x)
    nr <- nrow(x)
    rtt <- 1/nr
  }
  if (fda.usc::is.fdata(x)) {
    nc <- ncol(x$data)
    nr <- nrow(x$data)
    t <- x$argvals
    x <- x$data
    rtt <- t[2] - t[1]
  }

  if (is.null(nc) | is.null(nr))
    stop("x is not a fdata or a matrix")
  if (ncol(b) != nc | nrow(b) != nc)
    stop("Revised the dimensions of b")
  mat <- diag(x %*% b %*% t(x)) * rtt^2
  return(mat)
  # return(apply(mat,1,mean))

}


