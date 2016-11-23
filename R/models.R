

#' @title Functional covariates and coefficients for functional linear models
#' 
#' @description The functions \code{r.cfs.2003}, \code{r.hh.2006} and \code{r.bridge} sample the functional covariate and construct functional coefficients for their use in functional linear models:
#' \itemize{
#'   \item{r.cfs.2003}{implements examples (a) and (b) in Cardot et al. (2003).}
#'   \item{r.hh.2006}{gives models (i), (ii) and (iii) in Hall and Housseini-Nasab (2006).}
#'   \item{r.bridge}{samples a brownian motion and creates a functional coefficient made of the eigenfunctions \eqn{\sqrt(2) * \sin(k t \pi)}{sqrt(2) * sin(k*t*\pi)}.}
#' }
#' @param n sample size.
#' @param t locations where the functional data is observed.
#' @param b coefficients of the functional coefficient in the theoretical basis of principal components of the brownian motion.
#' @param type either example \code{"a"} or \code{"b"} from Cardot et al. (2003).
#' @param imod either \code{1}, \code{2} or \code{3} for denoting models (i), (ii) and (iii) in Hall and Hosseini-Nasab (2006).
#' @param ncb size of basis expansion used in Hall and Hosseini-Nasab (2006) for simulating the functional process.
#' @return A list with the following elements:
#' \itemize{
#'   \item{X.fdata}{the sample of functional data, an \code{\link[fda.usc]{fdata}} object of length \code{n}.}
#'   \item{beta.fdata}{the functional coefficient, an \code{\link[fda.usc]{fdata}} object.}
#' }
#' @examples 
#' # Cardot et al. (2003)
#' plot(r.cfs.2003(n = 100, type = "a")$X.fdata)
#' plot(r.cfs.2003(n = 100, type = "b")$X.fdata)
#' 
#' # Hall and Hosseini-Nasab (2006)
#' plot(r.hh.2006(n = 100, imod = 1)$X.fdata)
#' plot(r.hh.2006(n = 100, imod = 2)$X.fdata)
#' plot(r.hh.2006(n = 100, imod = 3)$X.fdata)
#' 
#' # Sample bridge
#' plot(r.bridge(n = 100)$X.fdata)
#' @author Manuel Febrero-Bande (\email{manuel.febrero@@usc.es}) and Eduardo Garcia-Portugues (\email{edgarcia@@est-econ.uc3m.es}).
#' @references 
#' Cardot, H., Ferraty, F., Sarda, P. (2003) Spline estimators for the functional linear model. Statistica Sinica, 13(3), 571--592. \url{http://www3.stat.sinica.edu.tw/statistica/oldpdf/a13n31.pdf}
#' Hall, P. and Hosseini-Nasab, M. (2006) On properties of functional principal components analysis. Journal of the Royal Statistical Society: Series B (Statistical Methodology), 68(1), 109--126. \url{http://dx.doi.org/10.1111/j.1467-9868.2005.00535.x}
#' @export 
r.cfs.2003 <- function(n = 100, t = seq(0, 1, len = 201), b = c(2, 4, 5) / sqrt(2), 
                       type = c("a", "b")[1]) {

  # X.fdata
  X.fdata <- fda.usc::rproc2fdata(n = n, t = t, sigma = "brownian")
  
  if (type == "a") {
    
    # beta = b1 * v1 + b2 * v2 + b3 * v3
    v <- fda.usc::fdata(sqrt(2) * sin(t * pi * 0.5), argvals = t)
    for (k in 2:length(b)) {
      
      v <- c(v, fda.usc::fdata(mdata = sqrt(2) * sin(t * pi * (k - 0.5)), 
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


#' @rdname r.cfs.2003
#' @export 
r.hh.2006 <- function(n = 100, t = seq(0, 1, len = 201), imod = 1, ncb = 20) {

  # beta
  beta.fdata <- fda.usc::fdata(pi^2 * (t^2 - 1/3), argvals = t)
  
  # X.fdata
  v <- fda.usc::fdata(sqrt(2) * cos(pi * t), argvals = t)
  for (k in 2:ncb) {
    
    v <- c(v, fda.usc::fdata(sqrt(2) * cos(k * pi * t), argvals = t))
    
  }
  sdcoef <- (1:ncb)^(-imod)
  coefsim <- matrix(NA, ncol = ncb, nrow = n)
  for (i in 1:ncb) {
    
    coefsim[, i] <- rnorm(n, mean = 0, sd = sdcoef[i])
    
  }
  X.fdata <- fda.usc::fdata(coefsim %*% v$data, argvals = t)

  return(list("X.fdata" = X.fdata, "beta.fdata" = beta.fdata))
  
}


#' @rdname r.cfs.2003
#' @export 
r.bridge <- function(n = 100, t = seq(0, 1, len = 201), b = c(2, 4, 5) / sqrt(2)) {

  # X.fdata
  X.fdata <- fda.usc::rproc2fdata(n = n, t = t, sigma = "brownian")
  lt <- length(t)
  for (i in 1:n) {
    
    X.fdata$data[i, ] <- X.fdata$data[i, ] - t * X.fdata$data[i, lt]
    
  }
  
  # beta = b1 * v1 + b2 * v2 + b3 * v3
  v <- fda.usc::fdata(mdata = sqrt(2) * sin(t * pi), argvals = t)
  for (k in 2:length(b)) {
    
    v <- c(v, fda.usc::fdata(sqrt(2) * sin(t * pi * k), argvals = t))
    
  }
  beta.fdata <- fda.usc::fdata(mdata = matrix(b, nrow = 1) %*% v$data, argvals = t)

  return(list("X.fdata" = X.fdata, "beta.fdata" = beta.fdata))
  
}


#' @title Geometric brownian motion
#' 
#' @description Sampling of paths of the geometric brownian motion.
#' 
#' @inheritParams r.cfs.2003
#' @param S0 initial value of the geometric brownian motion.
#' @return Functional sample, an \code{\link[fda.usc]{fdata}} object of length \code{n}.
#' @details The absolute value of the process will be always larger or equal to the absolute value of \code{S0}.
#' @examples 
#' plot(r.gbm(n = 100, S0 = rnorm(100)))
#' @author Eduardo Garcia-Portugues (\email{edgarcia@@est-econ.uc3m.es}).
#' @export
r.gbm <- function(n = 100, t = seq(0, 1, len = 201), S0 = 1) {
  
  S0 * exp(fda.usc::rproc2fdata(n = n, t = t, sigma = "brownian"))
  
}

