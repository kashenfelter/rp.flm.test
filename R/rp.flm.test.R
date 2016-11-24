

#' @title Data-driven sampling of random directions guided by sample of functional data
#'
#' @description Generation of random directions based on the principal components \eqn{\hat e_1,\ldots,\hat e_k}{\hat e_1,...,\hat e_k} of a sample of functional data \eqn{X_1,\ldots,X_n}{X_1,...,X_n}. The random directions are sampled as
#' \deqn{h=\sum_{j=1}^kh_j\hat e_j,}{h=\sum_{j=1}^kh_j\hat e_j,}
#' with \eqn{h_j\sim\mathcal{N}(0, \sigma_j^2)}{h_j~N(0, \sigma_j^2)}, \eqn{j=1,\ldots,k}{j=1,...,k}. Useful for sampling non-orthogonal random directions \eqn{h}{h} such that they are non-orthogonal for the random sample.
#'
#' @param n number of curves to be generated.
#' @param X.fdata an \code{\link[fda.usc]{fdata}} object used to compute the functional principal components.
#' @param ncomp if an integer vector is provided, the index for the principal components to be considered. If a threshold between \code{0} and \code{1} is given, the number of components \eqn{k}{k} is determined automatically as the minimum number that explains at least the \code{ncomp} proportion of the total variance of \code{X.fdata}.
#' @param sd whether the variances \eqn{\sigma_j} are estimated by the variances of the scores for \eqn{e_j}. If not, the \eqn{\sigma_j}'s are set to one.
#' @return A \code{\link[fda.usc]{fdata}} object with the sampled directions. 
#' @examples
#' # Simulate some data
#' set.seed(34567)
#' X.fdata <- rproc2fdata(n = 50, t = seq(0, 1, l = 201), sigma = "OU")
#' 
#' # Comparison for the variance type
#' par(mfrow = c(1, 1))
#' plot(X.fdata, col = gray(0.75), lty = 1)
#' lines(rdir.pc(n = 10, X.fdata = X.fdata), col = 2, lty = 1)
#' lines(rdir.pc(n = 10, X.fdata = X.fdata, sd = FALSE), col = 4, lty = 1)
#' legend("topleft", legend = c("Data", "Different variances", "Equal variances"), 
#'        col = c(gray(0.5), 2, 4), lwd = 2)
#'        
#' # Comparison for the threshold
#' samp1 <- rdir.pc(n = 100, X.fdata = X.fdata, ncomp = 0.5)
#' samp2 <- rdir.pc(n = 100, X.fdata = X.fdata, ncomp = 0.9)
#' samp3 <- rdir.pc(n = 100, X.fdata = X.fdata, ncomp = 0.95)
#' samp4 <- rdir.pc(n = 100, X.fdata = X.fdata, ncomp = 0.99)
#' samp5 <- rdir.pc(n = 100, X.fdata = X.fdata, ncomp = 0.999)
#' cols <- rainbow(5, alpha = 0.75)
#' par(mfrow = c(3, 2))
#' plot(X.fdata, col = gray(0.75), lty = 1, main = "Data")
#' plot(samp1, col = cols[1], lty = 1, main = "Threshold = 0.5")
#' plot(samp2, col = cols[2], lty = 1, main = "Threshold = 0.95")
#' plot(samp3, col = cols[3], lty = 1, main = "Threshold = 0.90")
#' plot(samp4, col = cols[4], lty = 1, main = "Threshold = 0.99")
#' plot(samp5, col = cols[5], lty = 1, main = "Threshold = 0.999")
#' @author Eduardo Garcia-Portugues (\email{edgarcia@@est-econ.uc3m.es}) and Manuel Febrero-Bande (\email{manuel.febrero@@usc.es}).
#' @export
rdir.pc <- function(n = 10, X.fdata, ncomp = 0.9, sd = TRUE) {
  
  # Check fdata
  if (class(X.fdata) != "fdata") {
    
    stop("X.fdata must be of class fdata")
    
  } 

  # Compute PCs: fda.usc::fdata2pc computes all the PCs and then returns 
  # the eigenvectors for the ncomp components (but returns eigenvalues of all the
  # components)
  if (ncomp < 1) {
    
    # Up to a threshold of the explained variance
    ej <- fda.usc::fdata2pc(fdataobj = X.fdata, ncomp = min(length(X.fdata$argvals),
                                                            nrow(X.fdata)))
    m <- max(2, min(which(cumsum(ej$d^2) / sum(ej$d^2) > ncomp)))
    ncomp <- 1:m
    
  } else {
    
    # Up to max(ncomp)
    m <- max(ncomp)
    ej <- fda.usc::fdata2pc(fdataobj = X.fdata, ncomp = m)
    
  }
  
  # Standard deviations of the scores of X.fdata on the eigenvectors
  if (sd) {
    
    sdarg <- apply(ej$x[, ncomp], 2, sd)
    
  } else {
    
    sdarg <- rep(1, length(ncomp))
    
  }
  
  # Eigenvectors
  eigv <- ej$rotation[ncomp]
  
  # Compute linear combinations of the eigenvectors with coefficients sampled 
  # from a centred normal with standard deviations sdarg
  x <- matrix(rnorm(n * m), ncol = m)
  x <- sweep(x, 2, sdarg, "*")
  rprojs <- x %*% eigv$data
  
  # Add mean
  rprojs <- fda.usc::fdata(sweep(rprojs, 2, ej$mean$data, "+"), argvals = fda.usc::argvals(X.fdata))
  
  return(rprojs)
  
}


#' @title Statistics for testing the functional linear model using random projections
#'
#' @description Computes the Cramer-von Mises (CvM) and Kolmogorv-Smirnov (kS) statistics on the projected process
#' \deqn{T_{n,h}(u)=\frac{1}{n}\sum_{i=1}^n (Y_i-\langle X_i,\hat\beta\rangle)1_{\{\langle X_i, h\rangle\leq u\}},}{T_{n, h}(u)=1/n\sum_{i = 1}^n (Y_i - <X_i, \hat\beta>)1_{<X_i, h> \le u},}
#' designed to test the goodness-of-fit of a functional linear model with scalar response.
#'
#' @param proj.X matrix of size \code{c(n, n.proj)} containing, for each column, the projections of the functional data \eqn{X_1,\ldots,X_n} into a random direction \eqn{h}. Not required if \code{proj.X.ord} is provided.
#' @param residuals the residuals of the fitted funtional linear model, \eqn{Y_i-\langle X_i,\hat\beta\rangle}{Y_i - <X_i, \hat\beta, Y_i>}. Either a vector of length \code{n} (same residuals for all projections) or a matrix of size \code{c(n.proj, n)} (each projection has an associated set residuals).
#' @param proj.X.ord matrix containing the row permutations of \code{proj.X} which rearranges them increasingly, for each column. So, for example \code{proj.X[proj.X.ord[, 1], 1]} equals \code{sort(proj.X[, 1])}. If not provided, it is computed internally.
#' @param F.code whether to use faster \code{FORTRAN} code or \code{R} code.
#' @return A list containing:
#' \itemize{
#'   \item{statistic}{a matrix of size \code{c(n.proj, 2)} with the the CvM (first column) and KS (second) statistics, for the \code{n.proj} different projections.}
#'   \item{proj.X.ord}{the computed row permutations of \code{proj.X}, useful for recycling in subsequent calls to \code{rp.flm.statistic} with the same projections but different residuals.}
#' }
#' @details \code{NA}'s are not allowed neither in the functional covariate nor in the scalar response.
#' @examples
#' # Simulated example
#' set.seed(345678)
#' t <- seq(0, 1, l = 101)
#' n <- 100
#' X <- rproc2fdata(n = n, t = t, sigma = "OU")
#' beta0 <- fdata(mdata = cos(2 * pi * t) - (t - 0.5)^2, argvals = t,
#'                rangeval = c(0,1))
#' Y <- inprod.fdata(X, beta0) + rnorm(n, sd = 0.1)
#'
#' # Linear model
#' mod <- fregre.pc(fdataobj = X, y = Y, l = 1:3)
#'
#' # Projections
#' proj.X1 <- inprod.fdata(X, rproc2fdata(n = 1, t = t))
#' proj.X2 <- inprod.fdata(X, rproc2fdata(n = 1, t = t))
#' proj.X12 <- cbind(proj.X1, proj.X2)
#'
#' # Statistics
#' t1 <- rp.flm.statistic(proj.X = proj.X1, residuals = mod$residuals)
#' t2 <- rp.flm.statistic(proj.X = proj.X2, residuals = mod$residuals)
#' t12 <- rp.flm.statistic(proj.X = proj.X12, residuals = mod$residuals)
#' t1$statistic
#' t2$statistic
#' t12$statistic
#'
#' # Recycling proj.X.ord
#' rp.flm.statistic(proj.X.ord = t1$proj.X.ord, residuals = mod$residuals)$statistic
#' t1$statistic
#'
#' # Sort in the columns
#' cbind(proj.X12[t12$proj.X.ord[, 1], 1], proj.X12[t12$proj.X.ord[, 2], 2]) -
#' apply(proj.X12, 2, sort)
#'
#' # FORTRAN and R code
#' rp.flm.statistic(proj.X = proj.X1, residuals = mod$residuals)$statistic -
#' rp.flm.statistic(proj.X = proj.X1, residuals = mod$residuals, 
#'                  F.code = FALSE)$statistic
#'
#' # Matrix and vector residuals
#' rp.flm.statistic(proj.X = proj.X12, residuals = mod$residuals)$statistic
#' rp.flm.statistic(proj.X = proj.X12, 
#'                  residuals = rbind(mod$residuals, mod$residuals * 2))$statistic
#' @author Eduardo Garcia-Portugues (\email{edgarcia@@est-econ.uc3m.es}) and Manuel Febrero-Bande (\email{manuel.febrero@@usc.es}).
#' @references
#' Cuesta-Albertos, J.A., Garcia-Portugues, E., Gonzalez-Manteiga, W. and Febrero-Bande, M. (2016). Goodness-of-fit tests for the functional linear model based on randomly projected empirical processes. arXiv XXXX:XXXX. \url{https://arxiv.org/abs/XXXX.XXXX}
#' @useDynLib rp.flm.test
#' @export
rp.flm.statistic <- function(proj.X, residuals, proj.X.ord = NULL, F.code = TRUE) {
  
  # Number of projections
  n.proj <- ifelse(is.null(proj.X.ord), ncol(proj.X), ncol(proj.X.ord))
  
  # Residuals as a matrix
  if (!is.matrix(residuals)) {
    
    residuals <- matrix(residuals, nrow = n.proj, ncol = length(residuals),
                        byrow = TRUE)
    
  }
  n <- ncol(residuals)
  if (nrow(residuals) != n.proj) {
    
    stop("The number of rows in residuals must be the number of projections")
    
  }
  
  # Matrix of statistics (columns) projected in n.proj projections (rows)
  rp.stat <- matrix(0, nrow = n.proj, ncol = 2)
  
  # Order projections if not provided
  if (is.null(proj.X.ord)) {
    
    proj.X.ord <- apply(proj.X, 2, order)
    
  }
  
  # Compute statistics
  if (F.code) {
    
    # Statistic
    rp.stat <- .Fortran("rp_stat", proj_X_ord = proj.X.ord, residuals = residuals,
                        n_proj = n.proj, n = n, rp_stat_proj = rp.stat,
                        PACKAGE = "rp.flm.test")$rp_stat_proj
    
  } else {
    
    # R implementation
    for (i in 1:n.proj) {
      
      # Empirical process
      y <- cumsum(residuals[i, proj.X.ord[, i]])
      
      # Statistics (CvM and KS, rows)
      CvM <- sum(y^2)
      KS <- max(abs(y))
      rp.stat[i, ] <- c(CvM, KS)
      
    }
    
    # Standardize
    rp.stat[, 1] <- rp.stat[, 1] / (n^2)
    rp.stat[, 2] <- rp.stat[, 2] / sqrt(n)
    
  }
  
  # Return both statistics
  colnames(rp.stat) <- c("CvM", "KS")
  return(list(statistic = rp.stat, proj.X.ord = proj.X.ord))
  
}


#' @title Goodness-of fit test for the functional linear model using random projections
#'
#' @description Tests the composite null hypothesis of a Functional Linear Model with scalar response (FLM),
#' \deqn{H_0:\,Y=\langle X,\beta\rangle+\epsilon\quad\mathrm{vs}\quad H_1:\,Y\neq\langle X,\beta\rangle+\epsilon.}{H_0: Y=<X,\beta>+\epsilon vs H_1: Y!=<X,\beta>+\epsilon.}
#' If \eqn{\beta=\beta_0}{\beta=\beta0} is provided, then the simple hypothesis \eqn{H_0:\,Y=\langle X,\beta_0\rangle+\epsilon}{H_0: Y=<X,\beta0>+\epsilon} is tested. The way of testing the null hypothesis is via a norm (Cramer-von Mises or Kolmogorov-Smirnov) in the empirical process indexed by the projections.
#'
#' @param X.fdata functional observations in the class \code{\link[fda.usc]{fdata}}.
#' @param Y scalar responses for the FLM. Must be a vector with the same number of elements as functions are in \code{X.fdata}.
#' @param beta0.fdata functional parameter for the simple null hypothesis, in the \code{\link[fda.usc]{fdata}} class. The \code{argvals} and \code{rangeval} arguments of \code{beta0.fdata} must be the same of \code{X.fdata}. If \code{beta0.fdata=NULL} (default), the function will test for the composite null hypothesis.
#' @param est.method estimation method for \eqn{\beta}{\beta}, only used in the composite case. There are three methods:
#' \itemize{
#'   \item{"pc"}{if \code{p} is given, then \eqn{\beta}{\beta} is estimated by \code{\link[fda.usc]{fregre.pc}}. Otherwise, \code{p} is chosen using \code{\link[fda.usc]{fregre.pc.cv}} and the \code{"SICc"} criterion.}
#'   \item{"pls"}{if \code{p} is given, \eqn{\beta}{\beta} is estimated by \code{\link[fda.usc]{fregre.pls}}. Otherwise, an optimum \code{p} is chosen using \code{\link[fda.usc]{fregre.pls.cv}} and the \code{"SICc"} criterion.}
#'   \item{"basis"}{if \code{p} is given, \eqn{\beta}{\beta} is estimated by \code{\link[fda.usc]{fregre.basis}}. Otherwise, an optimum \code{p} is chosen using \code{\link[fda.usc]{fregre.basis.cv}} and the \code{"GCV.S"} criterion. Both in \code{\link[fda.usc]{fregre.basis}} and \code{\link[fda.usc]{fregre.basis.cv}}, the same basis for \code{basis.x} and \code{basis.b} is considered.}
#' }
#' @param p number of elements for the basis representation of \code{beta0.fdata} and \code{X.fdata} with the \code{est.method} (only composite hypothesis). If not supplied, it is estimated from the data.
#' @param pmax maximum size of the basis expansion to consider in \code{\link{fregre.pc.cv}}, \code{\link{fregre.pls.cv}} and \code{\link{fregre.basis.cv}} for the optimal estimation of \eqn{\beta}{\beta}.
#' @param type.basis type of basis if \code{est.method = "basis"}.
#' @param B number of bootstrap replicates to calibrate the distribution of the test statistic.
#' @param n.proj vector with the number of projections to consider.
#' @param verbose whether to show or not information about the testing progress.
#' @param projs a \code{\link[fda.usc]{fdata}} object containing the random directions employed to project \code{X.fdata}. If numeric, the convenient value for \code{ncomp} in \code{\link{rdir.pc}}.
#' @param same.rwild wether to employ the same wild bootstrap residuals for different projections or not.
#' @param ... further arguments passed to \code{\link[fda]{create.basis}} (not \code{rangeval} that is taken as the \code{rangeval} of \code{X.fdata}).
#' @return An object with class \code{"htest"} whose underlying structure is a list containing the following components:
#' \itemize{
#'   \item{p.values.fdr}{A matrix of size \code{c(n.proj, 2)}, containing in each row the FDR p-values of the CvM and KS tests up to that projection.}
#'   \item{proj.statistics}{A matrix of size \code{c(n.proj, 2)} with the value of the test statistic on each projection.}
#'   \item{boot.proj.statistics}{An array of size \code{c(n.proj, 2, B)} with the values of the bootstrap test statistics for each projection.}
#'   \item{proj.p.values}{A matrix of size \code{c(n.proj, 2)}}
#'   \item{method}{Information about the test performed and the kind of estimation performed.}
#'   \item{B}{Number of bootstrap replicates used.}
#'   \item{n.proj}{Number of projections considered.}
#'   \item{projs}{Random directions employed to project \code{X.fdata}.}
#'   \item{type.basis}{Type of basis for \code{est.method = "basis"}.}
#'   \item{beta.est}{Estimated functional parameter \eqn{\hat\beta}{\hat\beta} in the composite hypothesis. For the simple hypothesis, \code{beta0.fdata}.}
#'   \item{p}{Number of basis elements considered for estimation of \eqn{\beta}{\beta}.}
#'   \item{data.name}{The character string "Y = <X, b> + e"}
#' }
#' @details
#' No NA's are allowed neither in the functional covariate nor in the scalar response.
#' @examples
#' # Simulated example
#'
#' set.seed(345678)
#' t <- seq(0, 1, l = 101)
#' n <- 100
#' X <- rproc2fdata(n = n, t = t, sigma = "OU")
#' beta0 <- fdata(mdata = cos(2 * pi * t) - (t - 0.5)^2, argvals = t,
#'                rangeval = c(0,1))
#' Y <- inprod.fdata(X, beta0) + rnorm(n, sd = 0.1)
#'
#' # Test all cases
#' rp.flm.test(X.fdata = X, Y = Y, est.method = "pc", B = 1000) #  p-value = 0.83, p-value = 0.75
#' rp.flm.test(X.fdata = X, Y = Y, est.method = "pls", B = 1000) #  p-value = 0.5000, p-value = 0.5889
#' rp.flm.test(X.fdata = X, Y = Y, est.method = "basis", B = 1000) #  p-value = 0.2, p-value = <2e-16
#'
#' rp.flm.test(X.fdata = X, Y = Y, est.method = "pc", p = 5, B = 1000)
#' rp.flm.test(X.fdata = X, Y = Y, est.method = "pls", p = 5, B = 1000)
#' rp.flm.test(X.fdata = X, Y = Y, est.method = "basis", p = 5, B = 1000)
#'
#' rp.flm.test(X.fdata = X, Y = Y, beta0.fdata = beta0, est.method = "pc", B = 1000)
#' rp.flm.test(X.fdata = X, Y = Y, beta0.fdata = beta0, est.method = "pls", B = 1000)
#' rp.flm.test(X.fdata = X, Y = Y, beta0.fdata = beta0, est.method = "basis", B = 1000)
#'
#' rp.flm.test(X.fdata = X, Y = Y, est.method = "pc", p = 1, B = 1000)
#' rp.flm.test(X.fdata = X, Y = Y, est.method = "pls", p = 1, B = 1000)
#' rp.flm.test(X.fdata = X, Y = Y, est.method = "basis", p = 5, B = 1000)
#' \dontrun{
#' # Composite hypothesis: do not reject FLM
#' rp.test <- rp.flm.test(X.fdata = X, Y = Y, est.method = "pc", B = 1000)
#' pcvm.test <- flm.test(X.fdata = X, Y = Y, est.method = "pc", B = 1000,
#'                       plot.it = FALSE)
#' rp.test
#' pcvm.test
#'
#' # Estimation of beta
#' par(mfrow = c(1, 3))
#' plot(X, main = "X")
#' plot(beta0, main = "beta")
#' lines(rp.test$beta.est, col = 2)
#' lines(pcvm.test$beta.est, col = 3)
#' plot(density(Y), main = "Density of Y", xlab = "Y", ylab = "Density")
#' rug(Y)
#'
#' # Simple hypothesis: do not reject beta = beta0
#' rp.flm.test(X.fdata = X, Y = Y, beta0.fdata = beta0, B = 1000)
#' flm.test(X.fdata = X, Y = Y, beta0.fdata = beta0, B = 1000, plot.it = FALSE)
#'
#' # Simple hypothesis: reject beta = beta0^2
#' rp.flm.test(X.fdata = X, Y = Y, beta0.fdata = beta0^2, B = 1000)
#' flm.test(X.fdata = X, Y = Y, beta0.fdata = beta0^2, B = 1000, plot.it = FALSE)
#'
#' # AEMET dataset
#'
#' # Load data
#' data(aemet)
#' wind.speed <- apply(aemet$wind.speed$data, 1, mean)
#' temp <- aemet$temp
#'
#' # Remove the 5% of the curves with less depth (i.e. 4 curves)
#' par(mfrow = c(1, 1))
#' res.FM <- depth.FM(temp, draw = TRUE)
#' qu <- quantile(res.FM$dep, prob = 0.05)
#' l <- which(res.FM$dep <= qu)
#' lines(aemet$temp[l], col = 3)
#'
#' # Data without outliers
#' wind.speed <- wind.speed[-l]
#' temp <- temp[-l]
#'
#' # Exploratory analysis: do not reject the FLM
#' rp.aemet <- rp.flm.test(X.fdata = temp, Y = wind.speed, B = 1000,
#'                         est.method = "pc")
#' pcvm.aemet <- flm.test(X.fdata = temp, Y = wind.speed, B = 1000,
#'                        est.method = "pc", plot.it = FALSE)
#' rp.aemet
#' pcvm.aemet
#'
#' # Estimated betas
#' plot(rp.aemet$beta.est, col = 2)
#' lines(pcvm.aemet$beta.est, col = 3)
#'
#' # B = 5000 for more precision on the calibration: do not reject the FLM
#' rp.flm.test(X.fdata = temp, Y = wind.speed, B = 5000, est.method = "pc")
#' flm.test(X.fdata = temp, Y = wind.speed, B = 5000, est.method = "pc",
#'          plot.it = FALSE)
#'
#' # Simple hypothesis: rejection of beta0 = 0
#' zero <- fdata(mdata = rep(0, length(temp$argvals)), argvals = temp$argvals,
#'               rangeval = temp$rangeval)
#' flm.test(X.fdata = temp, Y = wind.speed, beta0.fdata = zero, B = 1000,
#'          plot.it = FALSE)
#' rp.flm.test(X.fdata = temp, Y = wind.speed, beta0.fdata = zero, B = 1000)
#'
#' # Tecator dataset
#'
#' # Load data
#' data(tecator)
#' absorp <- tecator$absorp.fdata
#' ind <- 1:129 # or ind <- 1:215
#' x <- absorp[ind, ]
#' y <- tecator$y$Fat[ind]
#'
#' # Exploratory analysis for composite hypothesis with automatic choose of p
#' rp.tecat <- rp.flm.test(X.fdata = x, Y = y, est.method = "pc", B = 5000)
#' pcvm.tecat <- flm.test(X.fdata = x, Y = y, est.method = "pc", B = 5000,
#'                        plot.it = FALSE)
#' rp.tecat
#' pcvm.tecat
#'
#' # Distribution of the CvM and KS p-values
#' hist(rp.tecat$proj.p.values[, 1], lwd = 2, breaks = seq(0, 1, l = 20),
#'      freq = FALSE, main = "CvM and KS projected p-values")
#' rug(rp.tecat$proj.p.values)
#'
#' # Distribution of the PCvM statistic
#' plot(density(pcvm.tecat$boot.statistics), lwd = 2, main = "PCvM distribution",
#'      xlab = "PCvM*", ylab = "Density under the null")
#' rug(pcvm.tecat$boot.statistics)
#' abline(v = pcvm.tecat$statistic, col = 2, lwd = 2)
#' legend("topright", legend = "PCvM observed", lwd = 2, col = 2)
#'
#' # Simple hypothesis: fixed p
#' zero <- fdata(mdata = rep(0, length(x$argvals)), argvals = x$argvals,
#'               rangeval = x$rangeval)
#'
#' # Fixed p
#' rp.flm.test(X.fdata = x, Y = y, beta0.fdata = zero, B = 1000, p = 11)
#' flm.test(X.fdata = x, Y = y, beta0.fdata = zero, B = 1000, p = 11)
#'
#' # Simple hypothesis, automatic election of p
#' rp.flm.test(X.fdata = x, Y = y, beta0.fdata = zero, B = 1000)
#' flm.test(X.fdata = x, Y = y, beta0.fdata = zero, B = 1000)
#' }
#' @author Eduardo Garcia-Portugues (\email{edgarcia@@est-econ.uc3m.es}) and Manuel Febrero-Bande (\email{manuel.febrero@@usc.es}).
#' @references
#' Cuesta-Albertos, J.A., Garcia-Portugues, E., Gonzalez-Manteiga, W. and Febrero-Bande, M. (2016). Goodness-of-fit tests for the functional linear model based on randomly projected empirical processes. arXiv XXXX:XXXX. \url{https://arxiv.org/abs/XXXX.XXXX}
#' Garcia-Portugues, E., Gonzalez-Manteiga, W. and Febrero-Bande, M. (2014). A goodness-of-fit test for the functional linear model with scalar response. Journal of Computational and Graphical Statistics, 23(3), 761--778. \url{http://dx.doi.org/10.1080/10618600.2013.812519}
#' @export
rp.flm.test <- function(X.fdata, Y, beta0.fdata = NULL, est.method = "pc",
                        p = NULL, pmax = 10, type.basis = "bspline", B = 5000,
                        n.proj = 10, verbose = TRUE, projs = 0.9,
                        same.rwild = FALSE, ...) {
  
  # Sample size
  n <- dim(X.fdata)[1]
  
  # p data driven flag
  p.data.driven <- is.null(p)
  
  # Number of projections
  if (length(n.proj) > 1) {
    
    vec.nproj <- sort(n.proj)
    n.proj <- max(n.proj)
    
  } else {
    
    vec.nproj <- 1:n.proj
    
  }
  
  # Display progress
  if (verbose) {
    
    cat("Computing estimation of beta... ")
    
  }
  
  # Truncate maximum basis expansion
  pmax <- min(pmax, n)
  
  ## Estimation of beta
  
  # Composite hypothesis: optimal estimation of beta and the basis expansion
  if (is.null(beta0.fdata)) {
    
    # Center the data first
    X.fdata <- fda.usc::fdata.cen(X.fdata)$Xcen
    Y <- Y - mean(Y)
    
    # Method
    meth <- "Random projection based test for the functional linear model using"
    
    # PC
    if (est.method == "pc") {
      
      # Optimal p by SICc criterion
      if (p.data.driven) {
        
        # Method
        meth <- paste(meth, "optimal PC basis representation")
        
        # Choose the number of basis elements
        mod.pc <- fda.usc::fregre.pc.cv(fdataobj = X.fdata, y = Y, kmax = 1:pmax,
                                        criteria = "SICc")
        p.opt <- length(mod.pc$pc.opt)
        ord.opt <- mod.pc$pc.opt
        
        # PC components to be passed to the bootstrap
        pc.comp <- mod.pc$fregre.pc$fdata.comp
        pc.comp$l <- mod.pc$pc.opt
        
        # Express X.fdata and beta.est in the PC basis
        mdata <- mod.pc$fregre.pc$fdata.comp$x[, mod.pc$fregre.pc$l] %*%
          mod.pc$fregre.pc$fdata.comp$rotation$data[mod.pc$fregre.pc$l, , drop = FALSE]
        X.est <- fda.usc::fdata(mdata = mdata, argvals = X.fdata$argvals,
                                rangeval = X.fdata$rangeval)
        beta.est <- mod.pc$fregre.pc$beta.est
        norm.beta.est <- fda.usc::norm.fdata(beta.est)
        
        # Compute the residuals
        e <- mod.pc$fregre.pc$residuals
        
      # Fixed p
      } else {
        
        # Method
        meth <- paste(meth, " a representation in a PC basis of ", p, "elements")
        
        # Estimation of beta on the given fixed basis
        mod.pc <- fda.usc::fregre.pc(fdataobj = X.fdata, y = Y, l = 1:p)
        p.opt <- p
        ord.opt <- mod.pc$l
        
        # PC components to be passed to the bootstrap
        pc.comp <- mod.pc$pc
        pc.comp$l <- mod.pc$l
        
        # Express X.fdata and beta.est in the basis
        mdata <- mod.pc$fdata.comp$x[, mod.pc$l] %*%
          mod.pc$fdata.comp$rotation$data[mod.pc$l, , drop = FALSE]
        X.est <- fda.usc::fdata(mdata = mdata, argvals = X.fdata$argvals,
                                rangeval = X.fdata$rangeval)
        beta.est <- mod.pc$beta.est
        norm.beta.est <- fda.usc::norm.fdata(beta.est)
        
        # Compute the residuals
        e <- mod.pc$residuals
        
      }
      
    # PLS
    } else if (est.method == "pls") {
      
      # Optimal p by SICc criterion
      if (p.data.driven) {
        
        # Method
        meth <- paste(meth, "optimal PLS basis representation")
        
        # Choose the number of the basis: SIC is probably the best criteria
        mod.pls <- fda.usc::fregre.pls.cv(fdataobj = X.fdata, y = Y, kmax = pmax,
                                          criteria = "SICc")
        p.opt <- length(mod.pls$pls.opt)
        ord.opt <- mod.pls$pls.opt
        
        # PLS components to be passed to the bootstrap
        pls.comp <- mod.pls$fregre.pls$fdata.comp
        pls.comp$l <- mod.pls$pls.opt
        
        # Express X.fdata and beta.est in the PLS basis
        mdata <- mod.pls$fregre.pls$fdata.comp$x[, mod.pls$fregre.pls$l] %*%
          mod.pls$fregre.pls$fdata.comp$rotation$data[mod.pls$fregre.pls$l, , drop = FALSE]
        X.est <- fda.usc::fdata(mdata = mdata, argvals = X.fdata$argvals,
                                rangeval = X.fdata$rangeval)
        beta.est <- mod.pls$fregre.pls$beta.est
        norm.beta.est <- fda.usc::norm.fdata(beta.est)
        
        # Compute the residuals
        e <- mod.pls$fregre.pls$residuals
        
      # Fixed p
      } else {
        
        # Method
        meth <- paste(meth, "a representation in a PLS basis of ", p, "elements")
        
        # Estimation of beta on the given fixed basis
        mod.pls <- fda.usc::fregre.pc(fdataobj = X.fdata, y = Y, l = 1:p)
        p.opt <- p
        ord.opt <- mod.pls$l
        
        # PLS components to be passed to the bootstrap
        pls.comp <- mod.pls$fdata.comp
        pls.comp$l <- mod.pls$l
        
        # Express X.fdata and beta.est in the basis
        mdata <- mod.pls$fdata.comp$x[, mod.pls$l] %*%
          mod.pls$fdata.comp$rotation$data[mod.pls$l, , drop = FALSE]
        X.est <- fda.usc::fdata(mdata = mdata, argvals = X.fdata$argvals,
                                rangeval = X.fdata$rangeval)
        beta.est <- mod.pls$beta.est
        norm.beta.est <- fda.usc::norm.fdata(beta.est)
        
        # Compute the residuals
        e <- mod.pls$residuals
        
      }
      
    # Deterministic basis
    } else if (est.method == "basis") {
      
      # Optimal p by GCV criterion
      if (p.data.driven) {
        
        # Method
        meth <- paste(meth, "optimal", type.basis, "basis representation")
        
        # Choose the number of the bspline basis with GCV.S
        mod.basis <- fda.usc::fregre.basis.cv(fdataobj = X.fdata, y = Y,
                                              basis.x = 5:max(pmax, 5),
                                              basis.b = NULL,
                                              type.basis = type.basis,
                                              type.CV = fda.usc::GCV.S,
                                              verbose = FALSE, ...)
        p.opt <- mod.basis$basis.x.opt$nbasis
        ord.opt <- 1:p.opt
        
        # Express X.fdata and beta.est in the optimal basis
        X.est <- mod.basis$x.fd
        beta.est <- mod.basis$beta.est
        norm.beta.est <- fda.usc::norm.fd(beta.est)
        
        # Compute the residuals
        e <- mod.basis$residuals
        
      # Fixed p
      } else {
        
        # Method
        meth <- paste(meth, "a representation in a", type.basis, "basis of ",
                      p, "elements")
        
        # Estimation of beta on the given fixed basis
        basis.opt <- do.call(what = paste("create.", type.basis,
                                          ".basis", sep = ""),
                             args = list(rangeval = X.fdata$rangeval,
                                         nbasis = p, ...))
        mod.basis <- fda.usc::fregre.basis(fdataobj = X.fdata, y = Y,
                                           basis.x = basis.opt,
                                           basis.b = basis.opt)
        p.opt <- p
        ord.opt <- 1:p.opt
        
        # Express X.fdata and beta.est in the basis
        X.est <- mod.basis$x.fd
        beta.est <- mod.basis$beta.est
        norm.beta.est <- fda.usc::norm.fd(beta.est)
        
        # Compute the residuals
        e <- mod.basis$residuals
        
      }
      
    } else {
      
      stop(paste("Estimation method", est.method, "not implemented."))
      
    }
    
  # Simple hypothesis
  } else {
    
    # Method
    meth <- "Random projection based test for the simple hypothesis in a functional linear model"
    
    # Do not need to estimate beta
    beta.est <- beta0.fdata
    p.opt <- NA
    
    # Compute the residuals
    e <- drop(Y - fda.usc::inprod.fdata(X.fdata, beta.est))
    
  }
  
  ## Computation of the statistic
  
  # Sample random directions
  if (verbose) {
    
    cat("Done.\nComputing projections... ")
    
  }
  if (is.numeric(projs)) {
    
    projs <- rdir.pc(n = n.proj, X.fdata = X.fdata, ncomp = projs, sd = TRUE)
    
  }
  
  # Compute projections for the statistic and the bootstrap replicates
  proj.X <- fda.usc::inprod.fdata(X.fdata, projs) # A matrix n x n.proj
  
  # Statistic
  rp.stat <- rp.flm.statistic(proj.X = proj.X, residuals = e, F.code = TRUE)
  
  ## Bootstrap calibration
  
  # Define required objects
  rp.stat.star <- array(NA, dim = c(n.proj, 2, B))
  e.hat.star <- array(NA, dim = c(B, n.proj, n))
  if (verbose) {
    
    cat("Done.\nBootstrap calibration...\n ")
    pb <- txtProgressBar(style = 3)
    
  }
  
  # Composite hypothesis
  if (is.null(beta0.fdata)) {
    
    # Calculate the design matrix of the linear model depending on the chosen basis
    # This allows to resample efficiently the residuals without re-estimating
    # again the beta
    
    # PC
    if (est.method == "pc") {
      
      # Design matrix for the PC estimation
      X.matrix <- switch(p.data.driven + 1L,
                         mod.pc$lm$x,
                         mod.pc$fregre.pc$lm$x)
      
    # PLS
    } else if (est.method == "pls") {
      
      # Design matrix for the PLS estimation
      X.matrix <- switch(p.data.driven + 1L,
                         mod.pls$lm$x,
                         mod.pls$fregre.pls$lm$x)
      
    # Deterministic basis
    } else if (est.method == "basis") {
      
      # Design matrix for the basis estimation
      X.matrix <- mod.basis$lm$x
      
    }
    
    # Projection matrix of the linear model
    lm.proj.X.matrix <- diag(rep(1, n)) -
      X.matrix %*% solve(crossprod(X.matrix)) %*% t(X.matrix)
    
    # Bootstrap resampling
    for (i in 1:B) {
      
      # Generate bootstrap errors
      if (same.rwild) {
        
        e.hat <- matrix(fda.usc::rwild(e, "golden"), nrow = n.proj, ncol = n,
                        byrow = TRUE)
        
      } else {
        
        e.hat <- matrix(fda.usc::rwild(rep(e, n.proj), "golden"), nrow = n.proj,
                        ncol = n, byrow = TRUE)
        
      }
      
      # Calculate Y.star
      Y.star <- sweep(e.hat, 2, Y - e, "+")
      
      # Residuals from the bootstrap estimated model
      e.hat.star[i, , ] <- Y.star %*% lm.proj.X.matrix
      
      # Calculate the bootstrap statistics
      rp.stat.star[, , i] <- rp.flm.statistic(residuals = e.hat.star[i, , ],
                                              proj.X.ord = rp.stat$proj.X.ord,
                                              F.code = TRUE)$statistic
      
      # Display progress
      if (verbose) {
        
        setTxtProgressBar(pb, i / B)
        
      }
      
    }
    
  # Simple hypothesis
  } else {
    
    # Bootstrap resampling
    for (i in 1:B) {
      
      # Generate bootstrap errors
      if (same.rwild) {
        
        e.hat.star[i, , ] <- matrix(fda.usc::rwild(e, "golden"), nrow = n.proj, ncol = n,
                                    byrow = TRUE)
        
      } else {
        
        e.hat.star[i, , ] <- matrix(fda.usc::rwild(rep(e, n.proj), "golden"),
                                    nrow = n.proj, ncol = n, byrow = TRUE)
        
      }
      
      # Calculate the bootstrap statistics
      rp.stat.star[, , i] <- rp.flm.statistic(residuals = e.hat.star[i, , ],
                                              proj.X.ord = rp.stat$proj.X.ord,
                                              F.code = TRUE)$statistic
      
      # Display progress
      if (verbose) {
        
        setTxtProgressBar(pb, i/B)
        
      }
      
    }
    
  }
  
  # Compute the p-values
  pval <- matrix(nrow = n.proj, ncol = 2)
  for (i in 1:n.proj) {
    
    pval[i, 1] <- mean(rp.stat.star[i, 1, ] > rp.stat$statistic[i, 1]) # CvM
    pval[i, 2] <- mean(rp.stat.star[i, 2, ] > rp.stat$statistic[i, 2]) # KS
    
  }
  
  # Compute p-values depending for the vector of projections
  rp.pvalue <- matrix(nrow = length(vec.nproj), ncol = 2)
  for (k in seq_along(vec.nproj)) {
    
    rp.pvalue[k, ] <- apply(pval[1:vec.nproj[k], , drop = FALSE], 2, function(x) {
      
      l <- length(x)
      return(min(l / (1:l) * sort(x)))
      
    })
    
  }
  colnames(rp.pvalue) <- colnames(pval) <- c("CvM", "KS")
  rownames(rp.pvalue) <- vec.nproj
  
  # Return result
  if (verbose) {
    
    cat("\nDone.\n")
    
  }
  options(warn = -1)
  result <- structure(list(statistics.mean = colMeans(rp.stat$statistic),
                           p.values.fdr = rp.pvalue,
                           proj.statistics = rp.stat$statistic,
                           boot.proj.statistics = rp.stat.star,
                           proj.p.values = pval, method = meth, B = B,
                           n.proj = vec.nproj, projs = projs,
                           type.basis = type.basis, beta.est = beta.est, 
                           p = p.opt, data.name = "Y = <X, b> + e"))
  class(result) <- "htest"
  return(result)
  
}

