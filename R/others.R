
## Bug using fourier basis, currently unsolved ##
donotrun <- TRUE
if (!donotrun) {
  
  ## CRASH: possible bug in fregre.basis, the intercept appears TWO times on the design matrix (one from
  ## the Fourier basis and other from the linear model specification). So we have a 4 x 4 matrix even if
  ## the number of parameters is 3.
  samp <- rmod(n = 100, model = 13, delta = 0, R2 = 0.95, comp = TRUE)
  fda.usc::fregre.basis(fdataobj = samp$X.fdata, y = samp$Y, basis.x = create.power.basis(rangeval = c(0,
                                                                                                       1), exponents = c(0.5, 1, 2)))
  
  # Also crashes
  samp <- rmod(n = 100, model = 13, delta = 1, R2 = 0.95, comp = TRUE)
  fda.usc::fregre.basis(fdataobj = samp$X.fdata, y = samp$Y, basis.x = create.fourier.basis(nbasis = 5))
  
}

## With exponential basis it does not work neither ##
donotrun <- TRUE
if (!donotrun) {
  
  n <- 300
  R2 <- 0.95
  t <- seq(0, 1, l = 301)
  X.fdata <- rproc2fdata(n = n, t = t, sigma = "OrnsteinUhlenbeck")  # Just to create the fdata structure quickly
  for (i in 1:n) {
    X.fdata$data[i, ] <- rnorm(1, mean = 0, sd = 2) * exp(-t) + rnorm(1, mean = 1, sd = 1) * exp(t) +
      rnorm(1, mean = -1, sd = 0.1) * exp(2 * t)
  }
  beta0 <- fda.usc::fdata(mdata = -10 * exp(-t) - 5 * exp(t) + 3 * exp(2 * t), argvals = t, rangeval = c(0,
                                                                                                         1))
  
  # Noise
  varXg <- var(fda.usc::inprod.fdata(X.fdata, beta0))
  s2eps <- varXg * (1/R2 - 1)
  noise <- rnorm(n, mean = 0, sd = sqrt(s2eps))
  
  # Response
  Y <- drop(fda.usc::inprod.fdata(X.fdata, beta0)) + noise
  
  # Model
  basis1 <- fda::create.exponential.basis(rangeval = c(0, 1), ratevec = c(-1, 1, 2))
  plot(basis1)
  mod <- fda.usc::fregre.basis(fdataobj = X.fdata, y = Y, basis.x = basis1, basis.b = basis1)
  mod
  # True coefficients: 0, -10, -5, 3
  plot(beta0, ylim = c(-20, 20))
  lines(mod$beta.est, col = 2)
  # -10*exp(-t)-5*exp(t)+3*exp(2*t)
  
  ## Toy example from Manolo
  
  set.seed(15235)
  samp <- rmod(n = 1000, model = 14, delta = 0, R2 = 0.95, comp = TRUE)
  X <- samp$X.fdata
  Y <- samp$Y
  basis.b <- fda::create.bspline.basis(rangeval = c(0, 1), nbasis = 5)
  res <- fda.usc::fregre.basis(fdataobj = X, y = Y, basis.x = fda::create.bspline.basis(c(0, 1), 51),
                               basis.b = basis.b)
  res  # Should be c(1, 2, 3, 2, 1)
  plot(samp$beta.fdata, ylim = c(-1, 5))
  lines(res$beta.est, col = 2)
  
  
}
