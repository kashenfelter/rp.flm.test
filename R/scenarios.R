

#' @title Functional regression models for simple and composite hypothesis
#' 
#' @description TODO
#' 
#' @param n TODO
#' @param model TODO
#' @param delta TODO
#' @param R2 TODO
#' @param comp TODO
#' @return TODO
#' @example
#' # TODO
#' @author Eduardo Garcia-Portugues (\email{edgarcia@@est-econ.uc3m.es}).
#' @export 
r.mod <- function(n, model, delta, R2 = 0.95, comp = TRUE) {
  
  # Common argvals
  t <- seq(0, 1, l = 301)
  
  # Check for delta
  if (!(delta %in% 0:3)) {
    
    stop("delta must be 0, 1, 2 or 3")
  
  }
  
  # Different models with deviations
  if (model == 1) {
    
    # M1: example (a) of CFS 2003c, R2=0.95, three PC (1,2,3)
    b <- c(2/sqrt(2), 4/sqrt(2), 5/sqrt(2))
    fdatos <- fdata.cfs.2003.a(n = n, t = t, b = b, R2 = R2)
    X.fdata <- fdatos$x
    beta0 <- fdatos$beta
    
    # Deviations
    delta <- c(0, 0.25, 0.75, 1.25)[delta + 1]
    if (delta == 0) {
      
      mdev <- 0
      
    } else {
      
      mdev <- m.dev(X.fdata, type = 1)
      
    }
    
    # Response
    if (comp) {
      
      Y <- fdatos$y + delta * mdev
      
    } else {
      
      Y <- rnorm(n, 0, sqrt(fdatos$s2eps)) + delta * mdev
      
    }
    
  } else if (model == 2) {
    
    # M2: example (a) of CFS 2003c, R2=0.95, three PC (3,5,7)
    b <- c(0, 0, 2/sqrt(2), 0, 4/sqrt(2), 0, 5/sqrt(2))
    fdatos <- fdata.cfs.2003.a(n = n, t = t, b = b, R2 = R2)
    X.fdata <- fdatos$x
    beta0 <- fdatos$beta
    
    # Deviations
    delta <- c(0, 0.05, 0.2, 0.5)[delta + 1]
    if (delta == 0) {
      
      mdev <- 0
      
    } else {
      
      mdev <- m.dev(X.fdata, type = 1)
      
    }
    
    # Response
    if (comp) {
      
      Y <- fdatos$y + delta * mdev
      
    } else {
      
      Y <- rnorm(n, 0, sqrt(fdatos$s2eps)) + delta * mdev
      
    }
    
  } else if (model == 3) {
    
    # M3: example (a) of CFS 2003c, R2=0.95, three PC (2,3,7)
    b <- c(0, 2/sqrt(2), 4/sqrt(2), 0, 0, 0, 5/sqrt(2))
    fdatos <- fdata.cfs.2003.a(n = n, t = t, b = b, R2 = R2)
    X.fdata <- fdatos$x
    beta0 <- fdatos$beta
    
    # Deviations
    delta <- -c(0, 0.2, 0.5, 1)[delta + 1]
    if (delta == 0) {
      
      mdev <- 0
      
    } else {
      
      mdev <- m.dev(X.fdata, type = 1)
      
    }
    
    # Response
    if (comp) {
      Y <- fdatos$y + delta * mdev
    } else {
      Y <- rnorm(n, 0, sqrt(fdatos$s2eps)) + delta * mdev
    }
    
  } else if (model == 4) {
    
    # M4: example b of CFS 2003, Beta= log(15*t^2 + 10) + cos(4*pi*t)
    fdatos <- fdata.cfs.2003.b(n = n, t = t, R2 = R2)
    X.fdata <- fdatos$x
    beta0 <- fdatos$beta
    
    # Deviations
    delta <- c(0, 0.2, 1, 2)[delta + 1]
    if (delta == 0) {
      
      mdev <- 0
      
    } else {
      
      mdev <- m.dev(X.fdata, type = 1)
      
    }
    
    # Response
    if (comp) {
      
      Y <- fdatos$y + delta * mdev
      
    } else {
      
      Y <- rnorm(n, 0, sqrt(fdatos$s2eps)) + delta * mdev
      
    }
    
  } else if (model == 5) {
    
    # M5: example Hall and Hoseini, imod=1
    fdatos <- fdata.hh.2006(n = n, t = t, imod = 1, R2 = R2)
    X.fdata <- fdatos$x
    beta0 <- fdatos$beta
    
    # Deviations
    delta <- -c(0, 1, 3, 7)[delta + 1]
    if (delta == 0) {
      
      mdev <- 0
      
    } else {
      
      mdev <- m.dev(X.fdata, type = 2)
      
    }
    
    # Response
    if (comp) {
      
      Y <- fdatos$y + delta * mdev
      
    } else {
      
      Y <- rnorm(n, 0, sqrt(fdatos$s2eps)) + delta * mdev
      
    }
    
  } else if (model == 6) {
    
    # M6: example Hall & Hoseini, imod=2
    fdatos <- fdata.hh.2006(n = n, t = t, imod = 2, R2 = R2)
    X.fdata <- fdatos$x
    beta0 <- fdatos$beta
    
    # Deviations
    delta <- -c(0, 1, 3, 7)[delta + 1]
    if (delta == 0) {
      
      mdev <- 0
      
    } else {
      
      mdev <- m.dev(X.fdata, type = 2)
      
    }
    
    # Response
    if (comp) {
      
      Y <- fdatos$y + delta * mdev
      
    } else {
      
      Y <- rnorm(n, 0, sqrt(fdatos$s2eps)) + delta * mdev
      
    }
    
  } else if (model == 7) {
    
    # M7: Brownian bridge with exact representation in the three first PC, R2=0.95
    fdatos <- fdata.bridge(n = n, t = t, R2 = R2)
    X.fdata <- fdatos$x
    beta0 <- fdatos$beta
    
    # Deviations
    delta <- -c(0, 2, 7.5, 15)[delta + 1]
    if (delta == 0) {
      
      mdev <- 0
      
    } else {
      
      mdev <- m.dev(X.fdata, type = 2)
      
    }
    
    # Response
    if (comp) {
      
      Y <- fdatos$y + delta * mdev
      
    } else {
      
      Y <- rnorm(n, 0, sqrt(fdatos$s2eps)) + delta * mdev
      
    }
    
  } else if (model == 8) {
    
    # M8: Ornstein-Uhlenbeck and quadratic deviation I as in G-P,F-B and G-M (2013)
    X.fdata <- fda.usc::rproc2fdata(n = n, t = t, sigma = "OrnsteinUhlenbeck")
    beta0 <- fda.usc::fdata(mdata = sin(2 * pi * t) - cos(2 * pi * t), argvals = t, rangeval = c(0, 1))
    
    # Noise
    varXg <- var(fda.usc::inprod.fdata(X.fdata, beta0))
    s2eps <- varXg * (1/R2 - 1)
    noise <- rnorm(n, mean = 0, sd = sqrt(s2eps))
    
    # Deviations
    delta <- -c(0, 0.25, 1, 2)[delta + 1]
    if (delta == 0) {
      
      mdev <- 0
      
    } else {
      
      mdev <- m.dev(X.fdata, type = 2)
      
    }
    
    # Response
    if (comp) {
      
      Y <- drop(fda.usc::inprod.fdata(X.fdata, beta0)) + delta * mdev + noise
      
    } else {
      
      Y <- noise + delta * mdev
      
    }
    
  } else if (model == 9) {
    
    # M9: Ornstein-Uhlenbeck and quadratic deviation II as in G-P,F-B and G-M (2013)
    X.fdata <- fda.usc::rproc2fdata(n = n, t = t, sigma = "OrnsteinUhlenbeck")
    beta0 <- fda.usc::fdata(mdata = t - (t - 0.75)^2, argvals = t, rangeval = c(0, 1))
    
    # Noise
    varXg <- var(fda.usc::inprod.fdata(X.fdata, beta0))
    s2eps <- varXg * (1/R2 - 1)
    noise <- rnorm(n, mean = 0, sd = sqrt(s2eps))
    
    # Deviations
    delta <- -c(0, 0.01, 0.1, 0.3)[delta + 1]
    if (delta == 0) {
      
      mdev <- 0
      
    } else {
      
      mdev <- m.dev(X.fdata, type = 3)
      
    }
    
    # Response
    if (comp) {
      
      Y <- drop(fda.usc::inprod.fdata(X.fdata, beta0)) + delta * mdev + noise
      
    } else {
      
      Y <- noise + delta * mdev
      
    }
    
  } else if (model == 10) {
    
    # M10: Ornstein-Uhlenbeck and quadratic deviation III as in G-P,F-B and G-M (2013)
    X.fdata <- fda.usc::rproc2fdata(n = n, t = t, sigma = "OrnsteinUhlenbeck")
    beta0 <- fda.usc::fdata(mdata = t + cos(2 * pi * t), argvals = t, rangeval = c(0, 1))
    
    # Noise
    varXg <- var(fda.usc::inprod.fdata(X.fdata, beta0))
    s2eps <- varXg * (1/R2 - 1)
    noise <- rnorm(n, mean = 0, sd = sqrt(s2eps))
    
    # Deviations
    delta <- c(0, 0.01, 0.05, 0.25)[delta + 1]
    if (delta == 0) {
      
      mdev <- 0
      
    } else {
      
      mdev <- m.dev(X.fdata, type = 3)
      
    }
    
    # Response
    if (comp) {
      
      Y <- drop(fda.usc::inprod.fdata(X.fdata, beta0)) + delta * mdev + noise
      
    } else {
      
      Y <- noise + delta * mdev
      
    }
    
  } else if (model == 11) {
    
    # M11: Geometrical Brownian Motion I
    beta0 <- fda.usc::fdata(log(15 * t^2 + 10) + cos(4 * pi * t), argvals = t)
    X.fdata <- GBM(n = n, t = t, S0 = 1)
    
    # Noise
    varXg <- var(fda.usc::inprod.fdata(X.fdata, beta0))
    s2eps <- varXg * (1/R2 - 1)
    noise <- rnorm(n, mean = 0, sd = sqrt(s2eps))
    
    # Deviations
    delta <- c(0, 0.5, 2, 5)[delta + 1]
    if (delta == 0) {
      
      mdev <- 0
      
    } else {
      
      mdev <- m.dev(X.fdata, type = 3)
      
    }
    
    # Response
    if (comp) {
      
      Y <- drop(fda.usc::inprod.fdata(X.fdata, beta0)) + delta * mdev + noise
      
    } else {
      
      Y <- noise + delta * mdev
      
    }
    
  } else if (model == 12) {
    
    # M12: Geometrical Brownian Motion II
    beta0 <- fda.usc::fdata(pi^2 * (t^2 - 1/3), argvals = t)
    X.fdata <- GBM(n = n, t = t, S0 = 2)
    
    # Noise
    varXg <- var(fda.usc::inprod.fdata(X.fdata, beta0))
    s2eps <- varXg * (1/R2 - 1)
    noise <- rnorm(n, mean = 0, sd = sqrt(s2eps))
    
    # Deviations
    delta <- c(0, 0.5, 2.5, 7)[delta + 1]
    if (delta == 0) {
      
      mdev <- 0
      
    } else {
      
      mdev <- m.dev(X.fdata, type = 3)
      
    }
    
    # Response
    if (comp) {
      
      Y <- drop(fda.usc::inprod.fdata(X.fdata, beta0)) + delta * mdev + noise
      
    } else {
      
      Y <- noise + delta * mdev
      
    }
    
  } else if (model == 13) {
    
    # Toy example (Edu)
    X.fdata <- fda.usc::rproc2fdata(n = n, t = t, sigma = "OrnsteinUhlenbeck")
    for (i in 1:n) X.fdata$data[i, ] <- rnorm(1, mean = 0, sd = 2) * exp(-t) + rnorm(1, mean = 1,
                                                                                     sd = 1) * exp(t) + rnorm(1, mean = -1, sd = 0.1) * exp(2 * t)  #X.fdata$data[i, ]=rnorm(1, mean = 0, sd = 2)*sqrt(t)+rnorm(1, mean = 1, sd = 1)*t+rnorm(1, mean = -1, sd = 0.1)*t^2
    beta0 <- fda.usc::fdata(mdata = 1 - 10 * exp(-t) - 5 * exp(t) + 3 * exp(2 * t), argvals = t,
                            rangeval = c(0, 1))
    
    # Noise
    varXg <- var(fda.usc::inprod.fdata(X.fdata, beta0))
    s2eps <- varXg * (1/R2 - 1)
    noise <- rnorm(n, mean = 0, sd = sqrt(s2eps))
    
    # Deviations
    delta <- -c(0, 0.25, 1, 2)[delta + 1]
    if (delta == 0) {
      
      mdev <- 0
      
    } else {
      
      mdev <- m.dev(X.fdata, type = 2)
      
    }
    
    # Response
    if (comp) {
      
      Y <- drop(fda.usc::inprod.fdata(X.fdata, beta0)) + delta * mdev + noise
      
    } else {
      
      Y <- noise + delta * mdev
    }
    
  } else if (model == 14) {
    
    # Toy example (Manolo)
    X.fdata <- fda.usc::rproc2fdata(n = n, t = t, sigma = "vexponential")
    b5 <- fda::create.bspline.basis(rangeval = range(t), nbasis = 5)
    b5 <- fda.usc::fdata(t(fda::eval.basis(b5, t)), t)
    beta0 <- fda.usc::fdata(apply(sweep(b5$data, 1, c(1, 2, 3, 2, 1), "*"), 2, sum), t)
    
    # Noise
    varXg <- var(fda.usc::inprod.fdata(X.fdata, beta0))
    s2eps <- varXg * (1/R2 - 1)
    noise <- rnorm(n, mean = 0, sd = sqrt(s2eps))
    
    # Deviations
    delta <- -c(0, 0.25, 1, 2)[delta + 1]
    if (delta == 0) {
      
      mdev <- 0
      
    } else {
      
      mdev <- m.dev(X.fdata, type = 2, comp = comp)
      
    }
    
    # Response
    if (comp) {
      
      Y <- drop(fda.usc::inprod.fdata(X.fdata, beta0)) + delta * mdev + noise
      
    } else {
      
      Y <- delta * mdev + noise
      
    }
    
  } else {
    
    stop("Incorrect model")
    
  }
  
  return(list(X.fdata = X.fdata, Y = Y, beta.fdata = beta0))
  
}


#' @title Deviations from the composite hypothesis
#' 
#' @description TODO
#' 
#' @param X.fdata TODO
#' @param type TODO
#' @param beta0 TODO
#' @param pars TODO
#' @inheritParams r.mod
#' @return TODO
#' @example
#' # TODO
#' @author Eduardo Garcia-Portugues (\email{edgarcia@@est-econ.uc3m.es}).
#' @export
m.dev <- function(X.fdata, type, beta0, pars = NULL, comp = TRUE) {

  if (comp) {

    if (type == 1) {

      # Norm
      Y <- as.numeric(fda.usc::norm.fdata(X.fdata))

    } else if (type == 2) {

      # Quadratic model
      B <- outer(fda.usc::argvals(X.fdata), fda.usc::argvals(X.fdata), function(s, t) {
        sin(2 * pi * t * s) * s * (1 - s) * t * (1 - t)
      }) * 25
      Y <- quadratic(X.fdata, B)

    } else if (type == 3) {

      # Interior product of X.fdata non-linear transformations
      Y <- sapply(1:dim(X.fdata$data)[1], function(i) fda.usc::inprod.fdata(exp((-1) * X.fdata[i]),
        X.fdata[i]^2))

    } else if (type == 4) {

      # Non-linear transformation of the interior product
      eta <- beta0
      xx <- drop(fda.usc::inprod.fdata(X.fdata, eta))
      Y <- xx * log(xx^2)

    } else if (type == 5) {

      # Difference between ranges in absolute value
      Y <- apply(apply(abs(X.fdata$data), 1, range), 2, diff)

    } else {

      stop("Incorrect type")

    }

  } else {

    stop("Not implemented")

  }

  return(Y)

}


#' @title Check density of the responses and betas and functional process
#' 
#' @description TODO
#' 
#' @param models TODO
#' @inheritParams r.mod
#' @param times TODO
#' @return TODO
#' @example 
#' # TODO
#' @author Eduardo Garcia-Portugues (\email{edgarcia@@est-econ.uc3m.es}).
#' @export
test.models <- function(models = 1:12, comp = TRUE, times = TRUE) {

  par(mfrow = c(3, 4), mar = c(3, 2, 2, 0) + 0.1)

  for (k in models) {
    
    t <- proc.time()[3]

    set.seed(12456789)
    d0 <- density(rmod(n = 1000, model = k, delta = 0, comp = comp)$Y, bw = "SJ")
    set.seed(12456789)
    d1 <- density(rmod(n = 1000, model = k, delta = 1, comp = comp)$Y, bw = "SJ")
    set.seed(12456789)
    d2 <- density(rmod(n = 1000, model = k, delta = 2, comp = comp)$Y, bw = "SJ")
    # set.seed(12456789);d3=density(rmod(n=1000,model=k,delta=3,comp=comp)$Y,bw='SJ')
    t <- proc.time()[3] - t
    cat("Model", k, "finish in", t, "secs\n")
    plot(d0, type = "l", lwd = 2, main = paste("Model", k), xlab = "", ylab = "", ylim = c(0, max(d0$y,
      d1$y, d2$y) * 1.2))
    title(xlab = expression(Y == paste(symbol("\xe1"), list(X, rho), symbol("\xf1")) + delta[k] *
      m(X) + epsilon), ylab = "Density", cex.lab = 1.55)
    lines(d1, col = 2)
    lines(d2, col = 4)
    # lines(d3,col=4)
    legend("topright", legend = expression(delta[0] == 0, delta[1], delta[2]), lwd = 2, col = c(1:2,
      4))

  }

}


#' @rdname test.models
#' @export
test.betas <- function(models = 1:12, comp = TRUE) {

  par(mfrow = c(3, 4), mar = c(5, 4, 4, 2) + 0.1)

  for (k in models) {

    samp <- rmod(n = 20, model = k, delta = 0, comp = comp)
    lylim <- range(samp$beta.fdata$data) + c(-1, 1)
    rylim <- range(samp$X.fdata$data) + c(-1, 1)
    plotrix::twoord.plot(lx = samp$X.fdata$argvals, ly = samp$beta.fdata$data, rx = samp$X.fdata$argvals,
      ry = samp$X.fdata$data[1, ], xlab = "", rylab = "", ylab = "", lylim = lylim, rylim = rylim,
      type = "n", lcol = 1, rcol = gray(0.5), axis.main = 3, mar = c(3, 2.5, 2, 2.5) + 0.1)
    lines(samp$beta.fdata, ylim = lylim, lwd = 2, col = 1)
    
    for (i in 1:20) {
      
      lines(diff(lylim)/diff(rylim) * (samp$X.fdata[i, ] - rylim[1]) + lylim[1], col = gray(0.5))
      
    }

  }

}
