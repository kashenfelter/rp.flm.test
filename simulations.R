
# Install from GitHub - auth token should work for anyone
# pkgName <- "rp.flm.test"; devtools::install_github(paste("egarpor/", pkgName, sep = ""), username = NULL, ref = "master", subdir = NULL, auth_token = "5d85d98bdd4ad74f3897192760c4bb20e311173e", host = "api.github.com")

# Load packages
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
                                     p.fixed = FALSE, n.proj = 50, type.basis = "bspline",
                                     composite = TRUE, flm = FALSE, p.criterion = "SICc",
                                     ...) {

  # Seed
  seed <- .Random.seed

  # Set p.fixed to the true/sensible number of PCs or to NULL if FALSE
  if (p.fixed) {

    p.fixed <- switch(scenario, 3, 7, 7, 2, 2, 2, 4, 4, 2, 3, 2, 4)

  } else {

    p.fixed <- NULL

  }

  # Sample
  samp <- tryCatch(r.mod(n = n, scenario = scenario, delta = delta,
                         t = seq(0, 1, len = 201), R2 = 0.95,
                         composite = composite), error = function(e) NA)

  # Composite or simple hypothesis?
  if (composite) {

    beta0 <- NULL

  } else {

    argvals <- argvals(samp$X.fdata)
    beta0 <- fda.usc::fdata(mdata = rep(0, length(argvals)), argvals = argvals)

  }
  
  # Projections
  # projs <- 0.50
  # projs <- rdir.pc(n = n.proj, X.fdata = samp$X.fdata, ncomp = 0.95, norm = TRUE, 
  #                  zero.mean = FALSE)
  projs <- r.ou(n = n.proj, t = seq(0, 1, l = 201), mu = 0, alpha = 0.5, sigma = 1)
  
  # Random projection test
  test.rp.flm <- tryCatch(rp.flm.test(X.fdata = samp$X.fdata, Y = samp$Y,
                                      beta0.fdata = beta0, B = B, pmax = 13,
                                      est.method = est.method, n.proj = n.proj,
                                      p.criterion = p.criterion, p = p.fixed,
                                      F.code = TRUE, verbose = FALSE,
                                      same.rwild = FALSE, type.basis = type.basis,
                                      projs = projs), error = function(e) e)

  # Flm test
  if (flm) {

    test.flm <- tryCatch(fda.usc::flm.test(X.fdata = samp$X.fdata, Y = samp$Y,
                                           B = B, est.method = est.method,
                                           p = test.rp.flm$p, verbose = FALSE,
                                           plot.it = FALSE),
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
  #              "type.basis" = type.basis, "seed" = seed)
  # attr(result, "call") <- call

  return(result)

}

# Bootstrap replicates
B <- 1e3

# Monte Carlo replicates
M <- 1e3

# Models
scenarios <- c(1, 7, 3, 5, 6, 4, 8, 9, 12) # c(1:12)

# Sample sizes
n <- c(50, 100, 250)

# Deviations
delta <- 0 #:2

# Fixed projections
proj <- r.ou(n = 50, t = seq(0, 1, l = 201))

# Test
system.time(res <- user.simulation.function(scenario = 10, n = 50, n.proj = 50))

# Fake data generation
fg <- function(scenario) scenario

# Create a cluster
cl <- makeCluster(45, type = "SOCK")

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
  clusterGlobalObjects = c("B", "proj"),
  fallback = "simus.projindep.b1e3.m1e3.s",
  clusterSeed = rep(12345678, 6)

))

# Stop cluster
stopCluster(cl)

# Alert by email
system(paste("echo 'Sent from R. Time: ", time[3], " secs. Check output attached.' | mutt -s 'Simulations FLM finished!' -- edugarpor@gmail.com", sep = ""))

# Alert by email (with potential large attachment)
system(paste("echo 'Sent from R. Time: ", time[3], " secs. Check RData and output attached.' | mutt -s 'Simulations FLM finished! - with attached results' -a /home/egarcia/ScriptsFLM/new-simus.R -a /home/egarcia/ScriptsFLM/output.txt -a /home/egarcia/ScriptsFLM/new-simus.Rdata -- edugarpor@gmail.com", sep = ""))


##################
## New analysis ##
##################

if (!donotrun) {

## Read results

# Set wd
setwd("~/Dropbox/10 Linear model/New code")

# Load brute data
load("simus.test.b5e2.m1e3.s.Rdata")
load("simus.test.b1e3.m1e3.s.Rdata")
load("simus.test.b5e3.m1e3.s.Rdata")
load("simus.test.b1e4.m1e3.s.Rdata")
load("simus.test.fixproj.b1e4.m1e3.s.Rdata")
load("simus.test.095.b1e3.m1e3.s.Rdata")
load("simus.test.simp.b1e3.m1e3.s.Rdata")
load("simus.test.norm.b1e3.m1e3.s.Rdata")
load("simus.norm095mean.simp.b1e3.m1e3.s.Rdata")
load("simus.projindep.b1e3.m1e3.s.Rdata")
eg <- fallBackObj

# As data frame
df <- simTool:::as.data.frame.evalGrid(eg)
df <- df[, -c(1:3, 5)]

# Powers
power <- function(x, alpha = 0.05) mean(x < alpha, na.rm = TRUE)
by <- list(df$delta, df$n, df$scenario)
ind.beta <- grep(pattern = glob2rx("beta.est.data.*"), x = colnames(df))
df2 <- df[, -c(1:4, ind.beta)]
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
trajs.plot <- function(n, scenarios, alpha = 0.05, deltas = 0,
                       stat = "CvM", n.proj = 1:50, M = 1e3, main,
                       cex.main = 1) {

  # Colors
  col.05 <- "darkorchid4"
  col.10 <- "darkorange3"
  col.01 <- "firebrick3"

  if (identical(deltas, 0)) {

    par(mfrow = c(length(stat), length(n)))

  } else {

    par(mfrow = c(1, length(stat) * length(deltas)))

  }

  # Indices scenarios (aggregate returns them sorted)
  ind.sce <- match(x = scenarios, table = sort(unique(df$scenario)))

  show.main <- missing(main)
  for (tt in rev(stat)) {

    for (nn in n) {

      # Caption
      if (show.main) {

        main <- paste(stat, ", n = ", nn, sep = "")

      }

      for (delta in deltas) {

        if (delta == 0) {

          A1 <- trajs.power(n = nn, delta = delta, alpha = 0.05, stat = stat)
          matplot(n.proj, A1[n.proj, ind.sce], type = "o", pch = 16, cex = 0.5,
                  col = col.05, ylim = c(0, 0.15), lty = 1, xlim = c(1, max(n.proj)),
                  xlab = "Number of projections", ylab = "Empirical rejection rate",
                  main = main, cex.main = cex.main)
          A2 <- trajs.power(n = nn, delta = delta, alpha = 0.10, stat = stat)
          matlines(n.proj, A2[n.proj, ind.sce], type = "o", pch = 16, cex = 0.5,
                   lty = 1, col = col.10)
          A3 <- trajs.power(n = nn, delta = delta, alpha = 0.01, stat = stat)
          matlines(n.proj, A3[n.proj, ind.sce], type = "o", pch = 16, cex = 0.5,
                   lty = 1, col = col.01)
          abline(h = ci(p = 0.05, M = M), lwd = 3, col = col.05, lty = 2)
          abline(h = ci(p = 0.10, M = M), lwd = 3, col = col.10, lty = 2)
          abline(h = ci(p = 0.01, M = M), lwd = 3, col = col.01, lty = 2)

        } else {

          A <- trajs.power(n = nn, delta = delta, alpha = alpha, stat = stat)
          matplot(n.proj, A[n.proj, ind.sce], type = "o", pch = 16, cex = 0.5,
                  pch = 1, ylim = c(0, 1), lty = 1, xlim = c(1, max(n.proj)),
                  xlab = "Number of projections", ylab = "Empirical rejection rate",
                  main = main, cex.main = cex.main)

        }

      }

    }

  }

}

## Graphs FDR p-value vs projections for level

M <- 1e3
scenarios <- c(1, 7, 3, 5, 6, 4, 8, 9, 12)
trajs.plot(n = c(50, 100, 250), M = M, stat = "CvM", scenarios = scenarios)
trajs.plot(n = c(50, 100, 250), M = M, stat = "KS", scenarios = scenarios)

## Graphs FDR p-values vs projections for power

M <- 1e3
scenarios <- c(1, 7, 3, 5, 6, 4, 8, 9, 12)
trajs.plot(n = 50, M = M, deltas = 1:2, stat = "CvM", scenarios = scenarios)
trajs.plot(n = 100, M = M, deltas = 1:2, stat = "CvM", scenarios = scenarios)
trajs.plot(n = 250, M = M, deltas = 1:2, stat = "CvM", scenarios = scenarios)
trajs.plot(n = 50, M = M, deltas = 1:2, stat = "KS", scenarios = scenarios)
trajs.plot(n = 100, M = M, deltas = 1:2, stat = "KS", scenarios = scenarios)
trajs.plot(n = 250, M = M, deltas = 1:2, stat = "KS", scenarios = scenarios)

## Individual plots for sizes and powers

M <- 5e2
save <- FALSE
pow <- FALSE
for (n in c(50, 100, 250)) {

  ## H0

  # CvM
  if (save) pdf(paste("sizes_cvm_", n, ".pdf", sep = ""))
  trajs.plot(n = n, M = M, deltas = 0, stat = "CvM", scenarios = scenarios,
             main = substitute(list(H[0], CvM, n == nn), list(nn = n)),
             cex.main = 1)
  if (save) dev.off()

  # KS
  if (save) pdf(paste("sizes_ks_", n, ".pdf", sep = ""))
  trajs.plot(n = n, M = M, deltas = 0, stat = "KS", scenarios = scenarios,
             main = substitute(list(H[0], KS, n == nn), list(nn = n)),
             cex.main = 1)
  if (save) dev.off()

  ## H1, H2
  if (pow) {

    for (delta in 1:2) {

      # CvM
      if (save) pdf(paste("powers_cvm_", n, "_", delta, ".pdf", sep = ""))
      trajs.plot(n = n, M = M, deltas = delta, stat = "CvM", scenarios = scenarios,
                 main = substitute(list(H[0], CvM, n == nn), list(nn = n)),
                 cex.main = 1)
      if (save) dev.off()

      # KS
      if (save) pdf(paste("powers_ks_", n, "_", delta, ".pdf", sep = ""))
      trajs.plot(n = n, M = M, deltas = delta, stat = "KS", scenarios = scenarios,
                 main = substitute(list(H[0], KS, n == nn), list(nn = n)),
                 cex.main = 1)
      if (save) dev.off()

    }

  }

}

## Table average dn

ind.p <- match(x = "p", table = colnames(df))
summary.p <- aggregate.data.frame(x = df[, ind.p], by = by,
                                  FUN = function(x) {
                                    c(mean(x, na.rm = TRUE), sd(x, na.rm = TRUE),
                                      quantile(x, probs = c(0.00, 0.25, 0.5,
                                                            0.75, 1.00),
                                               na.rm = TRUE))
                                  })
summary.p <- as.matrix(summary.p)
colnames(summary.p) <- c("delta", "n", "scenario", "mean.p", "sd.p", "min.p",
                         "25%.p", "50%.p", "75%.p", "max.p")
summary.p

## Estimated vs real betas for n = 100

par(mfrow = c(3, 3), mar = c(5, 4, 4, 2) + 0.1)
t <- seq(0, 1, l = 201)
ind1 <- df$n == 100 & df$delta == 0
ind2 <- grep(pattern = glob2rx("beta.est.data.*"), x = colnames(df))
for (scenario in scenarios) {

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

  for (scenario in seq_along(scenarios)) {

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

