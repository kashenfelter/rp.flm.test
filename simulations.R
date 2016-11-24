
# Install from GitHub - auth token should work for anyone
# pkgName <- "rp.flm.test"; devtools::install_github(paste("egarpor/", pkgName, sep = ""), username = NULL, ref = "master", subdir = NULL, auth_token = "5d85d98bdd4ad74f3897192760c4bb20e311173e", host = "api.github.com")

# Load packages
library(fda.usc)
library(rp.flm.test)
library(simTool)
library(doSNOW)

# Skip chunks of code
donotrun <- TRUE

##########################
## New simulation study ##
##########################

# For general simulation
user.simulation.function <- function(scenario = 1, n = 100, delta = 0, est.method = "pc", 
                                     p.fixed = NA, n.proj = 50, type.basis = "bspline",
                                     composite = TRUE, flm = FALSE, p.criterion = "SICc",
                                     ...) {

  # Seed
  seed <- .Random.seed

  # Set p.fixed to NULL if NA
  if (is.na(p.fixed)) {
    
    p.fixed <- NULL

  }
  
  # Sample
  samp <- tryCatch(r.mod(n = n, scenario = scenario, delta = delta, R2 = 0.95, 
                         composite = composite), error = function(e) NA)

  # Composite or simple hypothesis?
  if (composite) {
    
    beta0 <- NULL
    
  } else {
    
    argvals <- argvals(samp$X.fdata)
    beta0 <- fdata(mdata = rep(0, length(argvals)), argvals = argvals)
    
  }
  
  # Random projection test
  test.rp.flm <- tryCatch(rp.flm.test(X.fdata = samp$X.fdata, Y = samp$Y, 
                                      beta0.fdata = beta0, B = B, pmax = 15,
                                      est.method = est.method, n.proj = n.proj, 
                                      p.criterion = p.criterion, p = p.fixed, 
                                      F.code = TRUE, verbose = FALSE, 
                                      same.rwild = FALSE, type.basis = type.basis,
                                      projs = 0.90), error = function(e) e)

  # Flm test
  if (flm) { 
    
    test.flm <- tryCatch(flm.test(X.fdata = samp$X.fdata, Y = samp$Y, B = B,
                                  est.method = est.method, p = test.rp.flm$p,
                                  verbose = FALSE, plot.it = FALSE),
                         error = function(e) NA)
    
  } else {
    
    test.flm <- list()
    test.flm$p.value <- NA
    
  }
  
  
  # Attach call
  result <- c("CvM" = test.rp.flm$p.values.fdr[, 1],
               "KS" = test.rp.flm$p.values.fdr[, 2],
               "PCvM" = test.flm$p.value,
               "p" = test.rp.flm$p,
               "beta.est.data." = drop(test.rp.flm$beta.est$data))

  # Attach call
  # call <- list("scenario" = scenario, "n" = n, "delta" = delta, "est.method" = est.method, 
  #              "p.fixed" = p.fixed, "n.proj" = n.proj, "composite" = composite, 
  #              "type.basis" = type.basis, "seed" = )
  # attr(result, "call") <- call
  
  return(result)

}

# Bootstrap replicates
B <- 5e3

# Monte Carlo replicates
M <- 1e3

# Models
scenarios <- c(1:12)

# Sample sizes
n <- 50 #c(50, 100, 250)

# Deviations
delta <- 0 #:2

# Test
system.time(res <- user.simulation.function(scenario = 10, n = 50, n.proj = 50))

# Fake data generation
fg <- function(scenario) scenario

# Create a cluster
cl <- makeCluster(40, type = "SOCK")

# To save fallback.R in other directory
# setwd("~/Desktop/")

# Simulation study
# Caution! If a NA is present in any of the parameters passed to expandGrid() 
# the default values will be used
time <- system.time(eg <- evalGrids(
  
  dataGrid = expandGrid("fg", scenario = scenarios),
  procGrid = expandGrid(proc = "user.simulation.function", n = n, delta = delta),
  replications = M,
  progress = TRUE,
  discardGeneratedData = TRUE,
  cluster = cl,
  clusterLibraries = c("fda.usc", "rp.flm.test"),
  clusterGlobalObjects = c("B"),
  fallback = "new-simus",
  clusterSeed = rep(12345678, 6)

))

# Stop cluster
stopCluster(cl)

# Alert by email
system(paste("echo 'Sent from R. Time: ", time[3], " secs. Check output attached.' | mutt -s 'Simulations FLM finished!' -a /home/egarcia/ScriptsFLM/new-simus.R -a /home/egarcia/ScriptsFLM/output.txt -- edugarpor@gmail.com", sep = ""))

# Alert by email (with potential large attachment)
system(paste("echo 'Sent from R. Time: ", time[3], " secs. Check RData and output attached.' | mutt -s 'Simulations FLM finished! - with attached results' -a /home/egarcia/ScriptsFLM/new-simus.R -a /home/egarcia/ScriptsFLM/output.txt -a /home/egarcia/ScriptsFLM/new-simus.Rdata -- edugarpor@gmail.com", sep = ""))


##################
## New analysis ##
##################

if (!donotrun) {

# Tried with p.fixed = 3 and only model = 1, which has a beta coefficient made of 3 PCs

# With M = 1000, n = 100, est.method = 1 and 10 clusters
#"Estimated replications per hour:  590"
# user   system  elapsed
# 95.706   21.176 6101.581

## Read results
  
# Set wd
setwd("~/Dropbox/10 Linear model/New code")

# Load brute data
load("new-simus.Rdata")
eg <- fallBackObj

# As data frame
df <- simTool:::as.data.frame.evalGrid(eg)
df <- df[, -c(1:3, 5)]

# Powers
power <- function(x, alpha = 0.05) mean(x < alpha, na.rm = TRUE)
by <- list(df$delta, df$n, df$scenario)
df2 <- df[, -c(1:4, (ncol(df) - 201 - 1):ncol(df))]
pow.10 <- aggregate.data.frame(x = df2, by = by, FUN = power, alpha = 0.10)
pow.05 <- aggregate.data.frame(x = df2, by = by, FUN = power, alpha = 0.05)
pow.01 <- aggregate.data.frame(x = df2, by = by, FUN = power, alpha = 0.01)
colnames(pow.10)[1:3] <- colnames(pow.05)[1:3] <- colnames(pow.01)[1:3] <- 
  c("delta", "n", "scenario")

# CvM and KS indexes
ind.CvM <- grep(pattern = glob2rx("CvM.*"), x = colnames(pow.05))
ind.KS <- grep(pattern = glob2rx("KS.*"), x = colnames(pow.05))
ind.PCvM <- grep(pattern = "PCvM", x = colnames(pow.05))

# 95% confidence interval for the proportion p from a sample of size M
ci <- function(p, M) {
  
  p + c(-1, 1) * qnorm(0.025) * sqrt(p * (1 - p) / M)
  
}

# Powers of the FDR p-values indexed in the number of projections
trajs.power <- function(n, delta = 0, alpha = 0.05, stat = "CvM") {
    
  # Which level?
  pow <- switch(as.character(alpha),
                "0.1" = pow.10,
                "0.05" = pow.05,
                "0.01" = pow.01)
  
  # Which statistic?
  ind <- switch(stat,
                "CvM" = ind.CvM,
                "KS" = ind.KS,
                "PCvM" = ind.PCvM)
  
  # Transposed trajectory for matplot
  t(pow[pow$delta == delta & pow$n == n, ind])
  
}

# Trajectories plot
trajs.plot <- function(n, scenarios = 1:12, alpha = 0.05, deltas = 0, 
                       stat = "CvM", n.proj = 1:50, M = 1e3, main = "", 
                       cex.main = 1) {
  
  # Colors
  col.05 <- "darkorchid4"
  col.10 <- "darkorange3"
  col.01 <- "firebrick3"
  
  # H0
  if (identical(deltas, 1)) {
    
    par(mfrow = c(1, length(stat)))
    
  } else {
    
    par(mfrow = c(1, length(stat) * length(deltas)))
    
  }
  for (tt in rev(stat)) {
    
    for (delta in deltas) {
      
      if (delta == 0) {
        
        A1 <- trajs.power(n = n, delta = delta, alpha = 0.05, stat = stat)
        matplot(n.proj, A1[, scenarios], type = "l", pch = 1, cex = 0.5, col = col.05, 
                ylim = c(0, 0.15), lty = 1, xlim = c(1, max(n.proj)), 
                xlab = "Number of projections", ylab = "Empirical rejection rate",
                main = main, cex.main = cex.main)
        A2 <- trajs.power(n = n, delta = delta, alpha = 0.10, stat = stat)
        matlines(n.proj, A2[, scenarios], type = "l", pch = 1, cex = 0.5, lty = 1, col = col.10)
        A3 <- trajs.power(n = n, delta = delta, alpha = 0.01, stat = stat)
        matlines(n.proj, A3[, scenarios], type = "l", pch = 1, cex = 0.5, lty = 1, col = col.01)
        abline(h = ci(p = 0.05, M = M), lwd = 3, col = col.05, lty = 2)
        abline(h = ci(p = 0.10, M = M), lwd = 3, col = col.10, lty = 2)
        abline(h = ci(p = 0.01, M = M), lwd = 3, col = col.01, lty = 2)
        
      } else {
        
        A <- trajs.power(n = n, delta = delta, alpha = alpha, stat = stat)
        matplot(n.proj, A[, scenarios], type = "l", 
                pch = 1, ylim = c(0, 1), lty = 1, xlim = c(1, max(n.proj)), 
                xlab = "Number of projections", ylab = "Empirical rejection rate",
                main = main, cex.main = cex.main)
        
      }
      
    }
    
  }
  
}

## Graphs FDR p-value vs projections for level

M <- 5e2
trajs.plot(n = 50, M = M, stat = "CvM")
trajs.plot(n = 100, M = M, stat = "CvM")
trajs.plot(n = 250, M = M, stat = "CvM")
trajs.plot(n = 50, M = M, stat = "KS")
trajs.plot(n = 100, M = M, stat = "KS")
trajs.plot(n = 250, M = M, stat = "KS")

## Graphs FDR p-values vs projections for power 

M <- 1e3
trajs.plot(n = 50, M = M, deltas = 1:2, stat = "CvM")
trajs.plot(n = 100, M = M, deltas = 1:2, stat = "CvM")
trajs.plot(n = 250, M = M, deltas = 1:2, stat = "CvM")
trajs.plot(n = 50, M = M, deltas = 1:2, stat = "KS")
trajs.plot(n = 100, M = M, deltas = 1:2, stat = "KS")
trajs.plot(n = 250, M = M, deltas = 1:2, stat = "KS")

## Individual plots for sizes and powers

save <- FALSE
for (n in c(50, 100, 250)) {
  
  ## H0
  
  # CvM
  if (save) pdf(paste("sizes_cvm_", n, ".pdf", sep = ""))
  trajs.plot(n = n, M = M, deltas = 0, stat = "CvM", 
             main = substitute(list(H[0], CvM, n == nn), list(nn = n)), 
             cex.main = 1)
  if (save) dev.off()
  
  # KS
  if (save) pdf(paste("sizes_ks_", n, ".pdf", sep = ""))
  trajs.plot(n = n, M = M, deltas = 0, stat = "KS", 
             main = substitute(list(H[0], KS, n == nn), list(nn = n)), 
             cex.main = 1)
  if (save) dev.off()
  
  ## H1, H2
  
  for (delta in 1:2) {
    
    # CvM
    if (save) pdf(paste("powers_cvm_", n, "_", delta, ".pdf", sep = ""))
    trajs.plot(n = n, M = M, deltas = delta, stat = "CvM", 
               main = substitute(list(H[0], CvM, n == nn), list(nn = n)), 
               cex.main = 1)
    if (save) dev.off()
    
    # KS
    if (save) pdf(paste("powers_ks_", n, "_", delta, ".pdf", sep = ""))
    trajs.plot(n = n, M = M, deltas = delta, stat = "KS", 
               main = substitute(list(H[0], KS, n == nn), list(nn = n)), 
               cex.main = 1)
    if (save) dev.off()
    
  }
  
}

## Table average dn

ind.p <- match(x = "p", table = colnames(df))
summary.p <- aggregate.data.frame(x = df[, ind.p], by = by, FUN = mean, 
                                  na.rm = TRUE)
summary.p <- cbind(summary.p, 
                   aggregate.data.frame(x = df[, ind.p], by = by, FUN = sd, 
                                        na.rm = TRUE)[, 4])
names(summary.p) <- c("delta", "n", "scenario", "mean.p", "sd.p")
summary.p

## Estimated vs real betas for n = 100

par(mfrow = c(3, 4), mar = c(5, 4, 4, 2) + 0.1)
t <- seq(0, 1, l = 201)
ind1 <- df$n == 100 & df$delta == 0
ind2 <- grep(pattern = glob2rx("beta.est.data.*"), x = colnames(df))
for (scenario in 1:12) {
  
  # Subset
  df3 <- df[df$scenario == scenario & ind1, ind2]
  
  # Plot
  beta <- r.mod(n = 1, scenario = scenario, delta = 0)$beta.fdata
  plot(beta, col = 2, lwd = 3, main = scenario)
  sapply(1:M, function(i) {
    lines(t, df3[i, ], col = gray(0.5, alpha = 0.1))
    })
  lines(beta, col = 2, lwd = 3)
  
}

## Tables n = 50, 100, 250

# Create table for alpha = 0.05
ind.k <- c(10, 15, 25)
block <- function(d = 0, a = 0.05) {

  t(rbind(
  trajs.power(n = 100, delta = d, alpha = a, stat = "KS")[ind.k, ],
  trajs.power(n = 100, delta = d, alpha = a, stat = "CvM")[ind.k, ],
  trajs.power(n = 100, delta = d, alpha = a, stat = "PCvM"),
  trajs.power(n = 250, delta = d, alpha = a, stat = "KS")[ind.k, ],
  trajs.power(n = 250, delta = d, alpha = a, stat = "CvM")[ind.k, ],
  trajs.power(n = 250, delta = d, alpha = a, stat = "PCvM")
  ))

}
tab <- rbind(block(d = 0), block(d = 1), block(d = 2))

# Print table
k <- 1
for (delta in 0:2) {
  
  for (scenario in 1:12) {
    
    H <- paste("$H_{", scenario, ",", delta, "}$", sep = "")
    cat(paste(H, paste("$", paste(sprintf(fmt = "%.3f", tab[k, ]), 
                                  collapse = "$ & $"), "$", sep = ""), 
              sep = " & "), "\\\\\n")
    k <- k + 1
    
  }
  if (delta == 2) {
    
    cat("\\bottomrule\\bottomrule\n")
    
  } else {
    
    cat("\\midrule\n")
    
  }
  
}

## Time simulation

# For general simulation
user.simulation.function.time <- function(scenario = 1, n = 100, delta = 0, 
                                          est.method = "pc", n.proj = 5, B = 1e3, 
                                          seed  = 2345678, ...) {
  
  # Fix different seeds
  set.seed(seed)
  
  # Sample
  samp <- tryCatch(r.mod(n = n, scenario = scenario, delta = delta, R2 = 0.95, 
                         composite = TRUE), error = function(e) NA)
  
  # Get the p
  aux <- tryCatch(rp.flm.test(X.fdata = samp$X.fdata, Y = samp$Y, 
                              beta0.fdata = NULL, B = 2, pmax = 10,
                              est.method = est.method, n.proj = n.proj, 
                              p = NULL, F.code = TRUE, verbose = FALSE, 
                              same.rwild = FALSE, projs = 0.9),
                  error = function(e) e)
  p <- tryCatch(aux$p, error = function(e) NULL)
  
  # Random projection test
  t.rp <- proc.time()
  test.rp.flm <- tryCatch(rp.flm.test(X.fdata = samp$X.fdata, Y = samp$Y, 
                                      beta0.fdata = NULL, B = B, pmax = 10,
                                      est.method = est.method, n.proj = n.proj, 
                                      p = p, F.code = TRUE, verbose = FALSE, 
                                      same.rwild = FALSE, projs = 0.9), 
                          error = function(e) e)
  t.rp <- proc.time() - t.rp
  
  # Flm test 
  t.flm <- proc.time()
  test.flm <- tryCatch(flm.test(X.fdata = samp$X.fdata, Y = samp$Y, B = B,
                                est.method = est.method, p = p,
                                verbose = FALSE, plot.it = FALSE),
                       error = function(e) NA)
  t.flm <- proc.time() - t.flm
  
  # Return results
  return(c(t.rp[3], t.flm[3]))
  
}

# Simulation
M <- 10
nn <- 2^(3:8)
res <- array(dim = c(length(nn), M, 2))
pb <- txtProgressBar(style = 3)
for (ni in seq_along(nn)) {
  
  cat("n =", ni, "\n")
  
  for (i in 1:M) {
  
    # Randomize models and deviations
    s <- (i - 1) %% 12 + 1
    d <- i %% 2
    
    # TIme
    res[ni, i, ] <- user.simulation.function.time(scenario = s, n = nn[ni], delta = d,
                                                  B = 1e3, est.method = "pc",
                                                  seed = round(log(i) * ni + 4335))
    
    # Display 
    setTxtProgressBar(pb = pb, value = i/M)
    
  }
  
  cat("Done\n")
  
}

# 1: RP 
# 2: FLM
means <- apply(res, c(1, 3), mean, na.rm = TRUE)
means

# Times
if (save) pdf("times.pdf")
plot(nn, means[, 2], type = "o", pch = 19, cex = 0.5, col = 3, 
     xlab = "Sample size", ylab = "Time (seconds)")
lines(nn, means[, 1], type = "o", pch = 19, cex = 0.5, col = 4)
legend("topleft", lwd = 2, col = 3:4, legend = c("CvM and KS", "PCvM"))
if (save) dev.off()

}


# Important points
# v Use SICc to penalize variability
# v Compute flm.test with the same kn
# v Remove non numeric values to be able to convert eg to a df
# v The larger variability on the estimate, the worse calibration of power
# v Some scenarios have fundamental problems when estimating beta from the PCs
