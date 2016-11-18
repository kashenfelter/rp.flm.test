

#' @title fdata.cfs.2003
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
#' @param type TODO
#' @return TODO
#' @examples 
#' plot(r.cfs.2003(n = 100, type = "a")$X.fdata)
#' @author Manuel Febrero-Bande (\email{manuel.febrero@@usc.es}) and Eduardo Garcia-Portugues (\email{edgarcia@@est-econ.uc3m.es}).
#' @export 
r.cfs.2003 <- function(n = 100, t = seq(0, 1, len = 101), b = c(2, 4, 5) / sqrt(2), 
                       type = "a") {

  # X.fdata
  X.fdata <- fda.usc::rproc2fdata(n = n, t = t, sigma = "brownian")
  
  if (type == "a") {
    
    # beta = b1 * v1 + b2 * v2 + b3 * v3
    v <- fda.usc::fdata(sqrt(2) * sin(t * pi * 0.5), argvals = t)
    for (i in 2:length(b)) {
      
      v <- c(v, fda.usc::fdata(mdata = sqrt(2) * sin(t * pi * (i - 0.5)), 
                               argvals = t))
      
    }
    beta.fdata <- fda.usc::fdata(mdata = matrix(b, nrow = 1) %*% v$data, argvals = t)
    
  } else if (type == "b") {
    
    beta.fdata <- fda.usc::fdata(mdata = log(15 * t^2 + 10) + cos(4 * pi * t), 
                                 argvals = t)
    
  } else {
    
    stop("Wrong type")
    
  }
  
  return(list("X.fdata" = X.fdata, "beta.fdata" = beta.fdata))
    
}


#' @title fdata.hh.2006
#' 
#' @description The function data.hh.2006 generates simulated data of the functional linear model with scalar
# response from Hall and Housseini-Nasab (2006).  imod=1,2,3 --> model (i), model (ii), model (iii)
# 
#' @inheritParams r.cfs.2003
#' @param imod TODO
#' @param ncb TODO
#' @return TODO
#' @examples 
#' plot(r.hh.2006(n = 100)$X.fdata)
#' @author Manuel Febrero-Bande (\email{manuel.febrero@@usc.es}) and Eduardo Garcia-Portugues (\email{edgarcia@@est-econ.uc3m.es}).
#' @export 
r.hh.2006 <- function(n = 100, t = seq(0, 1, len = 101), imod = 1, ncb = 20) {

  # beta
  beta.fdata <- fda.usc::fdata(pi^2 * (t^2 - 1/3), argvals = t)
  
  # X.fdata
  v <- fda.usc::fdata(sqrt(2) * cos(pi * t), argvals = t)
  for (i in 2:ncb) {
    v <- c(v, fda.usc::fdata(sqrt(2) * cos(i * pi * t), argvals = t))
  }
  sdcoef <- (1:ncb)^(-imod)
  coefsim <- matrix(NA, ncol = ncb, nrow = n)
  for (i in 1:ncb) {
    
    coefsim[, i] <- rnorm(n, mean = 0, sd = sdcoef[i])
    
  }
  X.fdata <- fda.usc::fdata(coefsim %*% v$data, argvals = t)

  return(list("X.fdata" = X.fdata, "beta.fdata" = beta.fdata))
  
}


#' @title Geometric brownian motion
#' 
#' @description TODO
#' 
#' @param n TODO
#' @param t TODO
#' @param S0 TODO
#' @return TODO
#' @examples 
#' plot(r.gbm(n = 100))
#' @export
r.gbm <- function(n = 100, t = seq(0, 1, len = 101), S0 = 1) {
  
  S0 * exp(fda.usc::rproc2fdata(n = n, t = t, sigma = "brownian"))

}


#' @title fdata.bridge
#' 
#' @description
#' Brownian Motion: eigenvectors: sqrt(2)*sin((k-.5)*pi*t) --- Lambda: 1/((k-.5)^2*pi^2) Brownian
#' Bridge: eigenvectors: sqrt(2)*sin(k*pi*t) --- Lambda: 1/(k^2*pi^2)
#' # Make Bridgebrownian to have a certain R2
#' 
#' @inheritParams r.cfs.2003
#' @return TODO
#' @examples 
#' plot(r.bridge(n = 100)$X.fdata)
#' @author Manuel Febrero-Bande (\email{manuel.febrero@@usc.es}) and Eduardo Garcia-Portugues (\email{edgarcia@@est-econ.uc3m.es}).
#' @export 
r.bridge <- function(n = 100, t = seq(0, 1, len = 101), b = c(2, 4, 5) / sqrt(2)) {

  # X.fdata
  X.fdata <- fda.usc::rproc2fdata(n = n, t = t, sigma = "brownian")
  lt <- length(t)
  for (i in 1:n) {
    
    X.fdata$data[i, ] <- X.fdata$data[i, ] - t * X.fdata$data[i, lt]
    
  }
  
  # beta = b1 * v1 + b2 * v2 + b3 * v3
  v <- fda.usc::fdata(mdata = sqrt(2) * sin(t * pi), argvals = t)
  for (i in 2:length(b)) {
    
    v <- c(v, fda.usc::fdata(sqrt(2) * sin(t * pi * i), argvals = t))
    
  }
  beta.fdata <- fda.usc::fdata(mdata = matrix(b, nrow = 1) %*% v$data, argvals = t)

  return(list("X.fdata" = X.fdata, "beta.fdata" = beta.fdata))
  
}


#' @title Quadratic model
#' 
#' @description TODO
#' 
#' @param X.fdata TODO
#' @param B TODO
#' @return TODO
#' @examples 
#' t <- seq(0, 1, len = 101)
#' X.fdata <- r.bridge(n = 100, t = t)$X.fdata
#' quadratic(X.fdata = X.fdata, B = outer(t, t, function(x, y) sin(2 * pi * x) + 
#'                                                             cos(2 * pi * y)))
#' @export
quadratic <- function(X.fdata, B) {

  x <- X.fdata$data
  diag(x %*% B %*% t(x)) * diff(X.fdata$argvals[1:2])^2

}


