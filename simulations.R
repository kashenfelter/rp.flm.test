
# Install from GitHub - auth token should work for anyone
# devtools::install_github("egarpor/rp.flm.test", username = NULL, ref = "master", subdir = NULL, auth_token = "5d85d98bdd4ad74f3897192760c4bb20e311173e", host = "api.github.com")

# Load packages
library(fda.usc)
library(rp.flm.test)
library(simTool)
library(doSNOW)

# Run full simulation study
run.simulations <- TRUE

######################
## Simulation study ##
######################

# For general simulation
user.simulation.function <- function(scenario = 1, n = 100, delta = 0, est.method = "pc",
                                     p.fixed = FALSE, n.proj = 10, type.basis = "bspline",
                                     composite = TRUE, flm = TRUE, p.criterion = "SICc",
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
  projs <- rdir.pc(n = n.proj, X.fdata = samp$X.fdata, ncomp = 0.95, norm = FALSE,
                   sd = 1, zero.mean = TRUE)
  # projs <- r.ou(n = n.proj, t = seq(0, 1, l = 201), mu = 0, alpha = 0.5, sigma = 1)
  
  # Random projection test
  test.rp.flm <- tryCatch(rp.flm.test(X.fdata = samp$X.fdata, Y = samp$Y,
                                      beta0.fdata = beta0, B = B, pmax = 20,
                                      est.method = est.method, n.proj = n.proj,
                                      p.criterion = p.criterion, p = p.fixed,
                                      F.code = TRUE, verbose = FALSE,
                                      same.rwild = FALSE, type.basis = type.basis,
                                      projs = projs), error = function(e) e)

  # Flm test
  if (flm & composite) {

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
scenarios <- c(1, 7, 3, 5, 6, 4, 8, 9, 12)

# Sample sizes
n <- c(50, 100, 250)

# Deviations
delta <- 0:2

if (run.simulations) {

# Test
system.time(res <- user.simulation.function(scenario = 10, n = 50))

# Fake data generation
fg <- function(scenario) scenario

# Create a cluster
cl <- makeCluster(45, type = "SOCK")

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
  fallback = "simus.proj.full.comp.dep095.sd1.b1e3.m1e3.s",
  clusterSeed = rep(12345678, 6)

))

# Stop cluster
stopCluster(cl)

# Alert by email
system(paste("echo 'Sent from R. Time: ", time[3], " secs. Check output attached.' | mutt -s 'Simulations FLM finished!' -- edugarpor@gmail.com", sep = ""))

}

##############
## Analysis ##
##############

if (!run.simulations) {

## Read results

# Set wd
setwd("~/Dropbox/10 Linear model/New code")

# Load brute data composite hypothesis
load("simus.proj.full.comp.dep095.sd0.b1e4.m1e3.s.Rdata")
load("simus.proj.full.comp.dep095.sd0.b1e3.m1e3.s.Rdata")
#load("simus.proj.full.comp.dep095.sd1.b1e4.m1e3.s.Rdata")
#load("simus.proj.full.comp.dep095.sd1.b1e3.m1e3.s.Rdata")
#load("simus.proj.full.comp.indep.b1e4.m1e3.s.Rdata")
#load("simus.proj.full.comp.indep.b1e3.m1e3.s.Rdata")

# Load brute data simple hypothesis
load("simus.proj.full.simp.dep095.sd0.b1e4.m1e3.s.Rdata")
load("simus.proj.full.simp.dep095.sd0.b1e3.m1e3.s.Rdata")

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
ind.CvM <- grep(pattern = glob2rx("CvM*"), x = colnames(pow.05))
ind.KS <- grep(pattern = glob2rx("KS*"), x = colnames(pow.05))
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
                       cex.main = 1, cex.lab = 1) {

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
                  main = main, cex.main = cex.main, cex.lab = cex.lab)
          A2 <- trajs.power(n = nn, delta = delta, alpha = 0.10, stat = stat)
          matlines(n.proj, A2[n.proj, ind.sce], type = "o", pch = 16, cex = 0.5,
                   lty = 1, col = col.10)
          A3 <- trajs.power(n = nn, delta = delta, alpha = 0.01, stat = stat)
          matlines(n.proj, A3[n.proj, ind.sce], type = "o", pch = 16, cex = 0.5,
                   lty = 1, col = col.01)
          abline(h = ci(p = 0.05, M = M), lwd = 3, col = col.05, lty = 2)
          abline(h = ci(p = 0.10, M = M), lwd = 3, col = col.10, lty = 2)
          abline(h = ci(p = 0.01, M = M), lwd = 3, col = col.01, lty = 2)
          abline(h = 0.05, lwd = 1, col = col.05, lty = 3)
          abline(h = 0.10, lwd = 1, col = col.10, lty = 3)
          abline(h = 0.01, lwd = 1, col = col.01, lty = 3)

        } else {

          A <- trajs.power(n = nn, delta = delta, alpha = alpha, stat = stat)
          matplot(n.proj, A[n.proj, ind.sce], type = "o", pch = 16, cex = 0.5, 
                  ylim = c(0, 1), lty = 1, xlim = c(1, max(n.proj)),
                  xlab = "Number of projections", ylab = "Empirical rejection rate",
                  main = main, cex.main = cex.main, cex.lab = cex.lab, 
                  col = rainbow(length(scenarios)))
          for (k in seq_along(scenarios)) {
            
            text(x = n.proj[1], y = A[1, ind.sce[k]], labels = k, pos = 2, 
                 col = rainbow(length(scenarios))[k])
            
          }


          
        }

      }

    }

  }

}

## Graphs FDR p-value vs projections for level

M <- 1e3
scenarios <- c(1, 7, 3, 5, 6, 4, 8, 9, 12)
n.proj <- 1:50
trajs.plot(n = c(50, 100, 250), M = M, stat = "CvM", scenarios = scenarios, n.proj = n.proj)
trajs.plot(n = c(50, 100, 250), M = M, stat = "KS", scenarios = scenarios, n.proj = n.proj)

## Graphs FDR p-values vs projections for power

M <- 1e3
scenarios <- c(1, 7, 3, 5, 6, 4, 8, 9, 12)
n.proj <- 1:50
trajs.plot(n = 50, M = M, deltas = 1:2, stat = "CvM", scenarios = scenarios, n.proj = n.proj)
trajs.plot(n = 100, M = M, deltas = 1:2, stat = "CvM", scenarios = scenarios, n.proj = n.proj)
trajs.plot(n = 250, M = M, deltas = 1:2, stat = "CvM", scenarios = scenarios, n.proj = n.proj)
trajs.plot(n = 50, M = M, deltas = 1:2, stat = "KS", scenarios = scenarios, n.proj = n.proj)
trajs.plot(n = 100, M = M, deltas = 1:2, stat = "KS", scenarios = scenarios, n.proj = n.proj)
trajs.plot(n = 250, M = M, deltas = 1:2, stat = "KS", scenarios = scenarios, n.proj = n.proj)

## Individual plots for sizes and powers

M <- 1e3
n.proj <- 1:50
save <- TRUE
pow <- TRUE
suffix <- ""
for (n in c(50, 100, 250)) {

  ## H0

  # CvM
  if (save) pdf(paste("sizes_cvm_", n, suffix, ".pdf", sep = ""))
  trajs.plot(n = n, M = M, deltas = 0, stat = "CvM", scenarios = scenarios,
             main = substitute(list(H[0], CvM, n == nn), list(nn = n)),
             cex.main = 1.5, cex.lab = 1.5, n.proj = n.proj)
  if (save) dev.off()

  # KS
  if (save) pdf(paste("sizes_ks_", n, suffix, ".pdf", sep = ""))
  trajs.plot(n = n, M = M, deltas = 0, stat = "KS", scenarios = scenarios,
             main = substitute(list(H[0], KS, n == nn), list(nn = n)),
             cex.main = 1.5, cex.lab = 1.5, n.proj = n.proj)
  if (save) dev.off()

  ## H1, H2
  
  if (pow) {

    for (delta in 1:2) {

      # CvM
      if (save) pdf(paste("powers_cvm_", n, "_", delta, suffix, ".pdf", sep = ""))
      trajs.plot(n = n, M = M, deltas = delta, stat = "CvM", scenarios = scenarios,
                 main = substitute(list(H[1], CvM, n == nn, d == dd), 
                                   list(nn = n, dd = delta)),
                 cex.main = 2, cex.lab = 1.5, n.proj = n.proj)
      if (save) dev.off()

      # KS
      if (save) pdf(paste("powers_ks_", n, "_", delta, suffix, ".pdf", sep = ""))
      trajs.plot(n = n, M = M, deltas = delta, stat = "KS", scenarios = scenarios,
                 main = substitute(list(H[1], KS, n == nn, d == dd), 
                                   list(nn = n, dd = delta)),
                 cex.main = 2, cex.lab = 1.5, n.proj = n.proj)
      if (save) dev.off()

    }

  }

}

## Table average dn

ind.p <- match(x = "p", table = colnames(df))
summary.p <- aggregate.data.frame(x = df[, ind.p], by = by,
                                  FUN = function(x) {
                                    c(mean(x, na.rm = TRUE), sd(x, na.rm = TRUE),
                                      quantile(x, probs = c(0.00, 0.25, 0.50,
                                                            0.75, 1.00),
                                               na.rm = TRUE))
                                  })
summary.p <- as.matrix(summary.p)
colnames(summary.p) <- c("delta", "n", "scenario", "mean.p", "sd.p", "min.p",
                         "25%.p", "50%.p", "75%.p", "max.p")
summary.p <- as.data.frame(summary.p)
summary.p

# Table
for (k in seq_along(scenarios)) {
  
  numbers <- sprintf("%.2f", 
                     summary.p$mean.p[subset = summary.p$scenario == scenarios[k]])
  numbers <- paste(numbers, " (", sprintf("%.2f", summary.p$sd.p[subset = summary.p$scenario == scenarios[k]]), ")", sep = "")
  cat(paste("$H_{", k, ",\\delta}$ & ", 
            paste(numbers, collapse = " & "),  sep = ""), "\\\\\n")
  
}

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
ind.k <- c(1, 5, 10)
block.1 <- function(d = 0, a = 0.05) {

  t(rbind(
  trajs.power(n = 100, delta = d, alpha = a, stat = "CvM")[ind.k, ],
  trajs.power(n = 100, delta = d, alpha = a, stat = "KS")[ind.k, ],
  trajs.power(n = 100, delta = d, alpha = a, stat = "PCvM"),
  trajs.power(n = 250, delta = d, alpha = a, stat = "CvM")[ind.k, ],
  trajs.power(n = 250, delta = d, alpha = a, stat = "KS")[ind.k, ],
  trajs.power(n = 250, delta = d, alpha = a, stat = "PCvM")
  ))

}
block.2 <- function(d = 0, a = 0.05) {
  
  t(rbind(
    trajs.power(n = 50, delta = d, alpha = a, stat = "CvM")[ind.k, ],
    trajs.power(n = 50, delta = d, alpha = a, stat = "KS")[ind.k, ],
    trajs.power(n = 50, delta = d, alpha = a, stat = "PCvM")
  ))
  
}
block.3 <- function(d = 0, a = 0.05) {
  
  ind.k <- c(1, 5)
  t(rbind(
    trajs.power(n = 50, delta = d, alpha = a, stat = "CvM")[ind.k, ],
    trajs.power(n = 50, delta = d, alpha = a, stat = "KS")[ind.k, ],
    trajs.power(n = 100, delta = d, alpha = a, stat = "CvM")[ind.k, ],
    trajs.power(n = 100, delta = d, alpha = a, stat = "KS")[ind.k, ],
    trajs.power(n = 250, delta = d, alpha = a, stat = "CvM")[ind.k, ],
    trajs.power(n = 250, delta = d, alpha = a, stat = "KS")[ind.k, ]
  ))
  
}
tab.1 <- rbind(block.1(d = 0), block.1(d = 1), block.1(d = 2))
tab.2 <- rbind(block.2(d = 0), block.2(d = 1), block.2(d = 2))
tab.3 <- rbind(block.3(d = 0), block.3(d = 1), block.3(d = 2))

# Print table
tab <- tab.3
tab <- tab.2
tab <- tab.1
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

# Average power gain of PCvM wrt to CvM5
power.gain <- function(n, d, k = 5) {
  
  pcvm <- trajs.power(n = n, delta = d, stat = "PCvM")[1, ]
  cvm <- trajs.power(n = n, delta = d, stat = "CvM")[k, ]
  if (d == 0) {
    
    return(summary(abs(pcvm - 0.05) - abs(cvm - 0.05)))
    
  } else {
    
    return(summary((cvm - pcvm) / pcvm))
  
  }
  
}
power.gain(n = 50, d = 0)
power.gain(n = 100, d = 0)
power.gain(n = 250, d = 0)
power.gain(n = 50, d = 1)
power.gain(n = 50, d = 2)
power.gain(n = 100, d = 1)
power.gain(n = 100, d = 2)
power.gain(n = 250, d = 1)
power.gain(n = 250, d = 2)


## Time simulation

# For general simulation
user.simulation.function.time <- function(n = 100, scenario = 1, delta = 0,
                                          est.method = "pc", n.proj = 5, B = 1e3,
                                          seed = 2345678, p = 10, 
                                          projs = rproc2fdata(5, t = seq(0, 1, l = 201)), 
                                          ...) {

  # Fix different seeds
  set.seed(seed)

  # Sample
  samp <- tryCatch(r.mod(n = n, scenario = scenario, delta = delta, R2 = 0.95,
                         composite = TRUE), error = function(e) NA)

  # Get the p
  if (is.null(p)) {
    
    aux <- tryCatch(rp.flm.test(X.fdata = samp$X.fdata, Y = samp$Y,
                                beta0.fdata = NULL, B = 2, pmax = 10,
                                est.method = est.method, n.proj = n.proj,
                                p = NULL, F.code = TRUE, verbose = FALSE,
                                same.rwild = FALSE, projs = 0.95),
                    error = function(e) e)
    p <- tryCatch(aux$p, error = function(e) NULL)

  }
  
  # Flm test
  t.flm <- proc.time()
  test.flm <- tryCatch(flm.test(X.fdata = samp$X.fdata, Y = samp$Y, B = B,
                                est.method = est.method, p = p,
                                verbose = FALSE, plot.it = FALSE),
                       error = function(e) NA)
  t.flm <- proc.time() - t.flm

  # Random projection test
  t.rp <- proc.time()
  test.rp.flm <- tryCatch(rp.flm.test(X.fdata = samp$X.fdata, Y = samp$Y,
                                      beta0.fdata = NULL, B = B, pmax = 10,
                                      est.method = est.method, n.proj = n.proj,
                                      p = p, F.code = TRUE, verbose = FALSE,
                                      same.rwild = FALSE, projs = projs),
                          error = function(e) e)
  t.rp <- proc.time() - t.rp
  
  # beta0
  beta0 <- samp$beta.fdata
  
  # Flm test simple
  t.flm.s <- proc.time()
  test.flm <- tryCatch(flm.test(X.fdata = samp$X.fdata, Y = samp$Y, 
                                beta0.fdata = beta0, B = B,
                                est.method = est.method, p = p,
                                verbose = FALSE, plot.it = FALSE),
                       error = function(e) NA)
  t.flm.s <- proc.time() - t.flm.s
  
  # Random projection test simple
  t.rp.s <- proc.time()
  test.rp.flm <- tryCatch(rp.flm.test(X.fdata = samp$X.fdata, Y = samp$Y,
                                      beta0.fdata = beta0, B = B, pmax = 10,
                                      est.method = est.method, n.proj = n.proj,
                                      p = p, F.code = TRUE, verbose = FALSE,
                                      same.rwild = FALSE, projs = projs),
                          error = function(e) e)
  t.rp.s <- proc.time() - t.rp.s

  # Return results
  return(c("RP" = t.rp[3], "FLM" = t.flm[3], "RP.S" = t.rp.s[3], "FLM.S" = t.flm.s[3]))

}

# Simulation
M <- 100
nn <- 2^(4:11)

# Fake data generation
fg <- function(n) n

# Create a cluster
cl <- makeCluster(5, type = "SOCK")

# Simulation study
eg <- evalGrids(
  
  dataGrid = expandGrid("fg", n = nn),
  procGrid = expandGrid(proc = "user.simulation.function.time"),
  replications = M,
  progress = TRUE,
  discardGeneratedData = TRUE,
  cluster = cl,
  clusterLibraries = c("fda.usc", "rp.flm.test"),
  fallback = "times",
  clusterSeed = rep(12345678, 6)
  
)

# Stop cluster
stopCluster(cl)

# Alert by email
system("echo 'Sent from R. Check output attached.' | mutt -s 'Time simulations finished!' -- edugarpor@gmail.com")

# Set wd
#setwd("~/Dropbox/10 Linear model/New code")

# Load brute data
load("times.Rdata")
eg <- fallBackObj

# As data frame
times <- simTool:::as.data.frame.evalGrid(eg, summary.fun = mean)[, -c(1:3, 5:6)]

# Times
save <- TRUE
if (save) pdf("times.pdf")
plot(times[, c(1, 3)], type = "o", pch = 19, cex = 0.5, col = 3,
     xlab = "Sample size", ylab = "Time (seconds)", log = "xy")
lines(times[, 1:2], type = "o", pch = 19, cex = 0.5, col = 4)
lines(times[, c(1, 5)], type = "o", pch = 19, lty = 2, cex = 0.5, col = 3)
lines(times[, c(1, 4)], type = "o", pch = 19, lty = 2, cex = 0.5, col = 4)
legend("topleft", lwd = 2, col = rep(4:3, each = 2), lty = rep(1:2, 2), 
       legend = c("CvM and KS, composite", "CvM and KS, simple", 
                  "PCvM, composite", "PCvM, simple"))
if (save) dev.off()

}

