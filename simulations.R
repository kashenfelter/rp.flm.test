

# Install from GitHub - auth token should work for anyone
# pkgName <- "rp.flm.test"
# install_github(paste("egarpor/", pkgName, sep = ""), username = NULL, ref = "master", subdir = NULL, auth_token = "5d85d98bdd4ad74f3897192760c4bb20e311173e", host = "api.github.com")

# Load packages
library(fda.usc)
library(rp.flm.test)
library(simTool)
library(doSNOW)
library(sm)

# Skip chunks of code
donotrun <- TRUE

##########################
## New simulation study ##
##########################

# For general simulation
user.simulation.function <- function(model = 1, n = 100, delta = 0, est.method = "pc", 
                                     p.fixed = NA, n.proj, type.basis = "bspline",
                                     comp = TRUE, ...) {

  # Call of the test and seed
  call <- list("model" = model, "n" = n, "delta" = delta, "est.method" = est.method, 
               "p.fixed" = p.fixed, "n.proj" = n.proj, "comp" = comp, 
               "type.basis" = type.basis, "seed" = .Random.seed)

  # Set p.fixed to NULL if NA
  if (is.na(p.fixed)) {
    
    p.fixed <- NULL

  }
  
  # Sample
  samp <- tryCatch(r.mod(n = n, model = model, delta = delta, R2 = 0.95, 
                         comp = comp), error = function(e) NA)

  # Composite or simple hypothesis?
  if (comp) {
    
    beta0 <- NULL
    
  } else {
    
    argvals <- argvals(samp$X.fdata)
    beta0 <- fdata(mdata = rep(0, length(argvals)), argvals = argvals)
    
  }
  
  # Random projection test
  test.rp.flm <- tryCatch(rp.flm.test(X.fdata = samp$X.fdata, Y = samp$Y, 
                                      beta0.fdata = beta0, B = 1000, 
                                      est.method = est.method, n.proj = n.proj, 
                                      p = p.fixed, F.code = TRUE, verbose = FALSE, 
                                      same.rwild = FALSE, type.basis = type.basis,
                                      proj.gamma = 0.90), error = function(e) e)

  # Remove boot.statistics to save space (x 10 space if not removed)
  test.rp.flm$boot.statistics <- NULL

  # Flm test 
  # test.flm <- tryCatch(flm.test(X.fdata = samp$X.fdata, Y = samp$Y, B = 2000, 
  #                               est.method = c('pc', 'pls', 'basis')[est.method], 
  #                               p = p, verbose = FALSE, plot.it = FALSE), 
  #                      error = function(e) NA)

  # Attach call
  result <- test.rp.flm

  # Return results
  attr(result, "call") <- call
  return(result)

}

# Test
system.time(res <- user.simulation.function(model = 1, n = 100, n.proj = 50))
attr(res, "call")
plot(res$beta.est)

# Fake data generation
fg <- function(model) model

# Create a cluster
cl <- makeCluster(30, type = "SOCK")
 
# To save fallback.R in other directory
# setwd("~/Desktop/")

# Replicates
M <- 1000

# Models
models <- c(1:12)

# Simulation study
# Caution! If a NA is present in any of the parameters passed to expandGrid() 
# the default values will be used
time <- system.time(eg <- evalGrids(
  
  dataGrid = expandGrid("fg", model = models),
  procGrid = expandGrid(proc = "user.simulation.function", n = 100, delta = 0, 
                        est.method = "pc", type.basis = "bspline",
                        p.fixed = NA, n.proj = 50, comp = TRUE),
  replications = M,
  progress = TRUE,
  discardGeneratedData = TRUE,
  cluster = cl,
  clusterLibraries = c("fda.usc", "rp.flm.test"),
  fallback = "new-simus",
  clusterSeed = rep(12345678, 6)

))

# Stop cluster
stopCluster(cl)

# Alert by email
system(paste("echo 'Sent from R. Time: ", time[3], " secs. Check output attached.' | mutt -s 'Simulations FLM finished!' -a /home/egarcia/ScriptsFLM/new.simus.rp.flm.test.R -a /home/egarcia/ScriptsFLM/output.txt -- edugarpor@gmail.com", sep = ""))

# Alert by email (with potential large attachment)
system(paste("echo 'Sent from R. Time: ", time[3], " secs. Check RData and output attached.' | mutt -s 'Simulations FLM finished! - with attached results' -a /home/egarcia/ScriptsFLM/new.simus.rp.flm.test.R -a /home/egarcia/ScriptsFLM/output.txt -a /home/egarcia/ScriptsFLM/new-simus.Rdata -- edugarpor@gmail.com", sep = ""))

if (!donotrun) {

# Tried with p.fixed = 3 and only model = 1, which has a beta coefficient made of 3 PCs

# With M = 1000, n = 100, est.method = 1 and 10 clusters
#"Estimated replications per hour:  590"
# user   system  elapsed
# 95.706   21.176 6101.581

# Set wd
setwd("~/Dropbox/10 Linear model/New code")

# Load data
load("new-simus.Rdata")
eg <- fallBackObj

# Extract p-values
j <- 1
eg$simulation[[1]][[1]]$results[[j]]$method
pvalues.cvm <- pvalues.ks <- array(dim = c(length(models), 50, 3))
sapply(1:length(models), function(model) {
  
  pvalues.ks[model, , ] <<- sapply(c(0.1, 0.05, 0.01), function(a) {
    
    rowMeans(sapply(1:M, function(i) {
      eg$simulation[[model]][[i]]$results[[j]]$proj.p.values[, 2]
      }) <= a) # 1:CvM, 2:KS
    
  })
  
  pvalues.cvm[model, , ] <<- sapply(c(0.1, 0.05, 0.01), function(a) {
    
    rowMeans(sapply(1:M, function(i) {
      eg$simulation[[model]][[i]]$results[[j]]$proj.p.values[, 1]
      }) <= a) # 1:CvM, 2:KS
    
  })
  
})

# Half length of the 95% confidence interval for the proportion p from a sample of size M
ci <- function(p, M) {
  
  p + c(-1, 1) * qnorm(0.025) * sqrt(p * (1 - p) / M)
  
}

# Are p-values within CIs?
KS <- abs(pvalues.ks[, , 1] - 0.1) < ci(0.1, M = M)
KS
CvM <- abs(pvalues.cvm[, , 1] - 0.1) < ci(0.1, M = M)
CvM



# Colors
col.05 <- "darkorchid4"
col.10 <- "darkorange3"
col.01 <- "firebrick3"

matplot(t(pvalues.cvm[, , 1]), type = "l", lty = 1, col = col.10,
        xlab = "Number of projections", ylab = "Empirical rejection rate", 
        ylim = c(0, 0.15))
matlines(t(pvalues.cvm[, , 2]), type = "l", lty = 1, col = col.05)
matlines(t(pvalues.cvm[, , 3]), type = "l", lty = 1, col = col.01)
abline(h = ci(p = 0.10, M = 100), lwd = 3, col = col.10, lty = 2)
abline(h = ci(p = 0.05, M = 100), lwd = 3, col = col.05, lty = 2)
abline(h = ci(p = 0.01, M = 100), lwd = 3, col = col.01, lty = 2)



matplot(n.proj, t(A[ind.dev, ind.05[ind.test]]), type = "l", pch = 1, 
        cex = 0.5, col = col.05, ylim = ylim, lty = 1, 
        xlim = c(1, max(n.proj)), main = main, 
        )
matlines(n.proj, t(A[ind.dev, ind.10[ind.test]]), type = "l", pch = 1, 
         cex = 0.5, lty = 1, col = col.10)
matlines(n.proj, t(A[ind.dev, ind.01[ind.test]]), type = "l", pch = 1, 
         cex = 0.5, lty = 1, col = col.01)
abline(h = 0.05 + qnorm(0.975) * sqrt(0.05 * (1 - 0.05)/500) * c(-1, 1), 
       lwd = 3, col = col.05, lty = 2)
abline(h = 0.1 + qnorm(0.975) * sqrt(0.1 * (1 - 0.1)/500) * c(-1, 1), 
       lwd = 3, col = col.10, lty = 2)
abline(h = 0.01 + qnorm(0.975) * sqrt(0.01 * (1 - 0.01)/500) * c(-1, 1), 
       lwd = 3, col = col.01, lty = 2)






#############
## Results ##
#############


# Comparative plots
plots <- function(A, dev, test, n.proj, main) {
  
  # Colors
  col.05 <- "darkorchid4"
  col.10 <- "darkorange3"
  col.01 <- "firebrick3"
  
  # Indexes
  ind.05 <- seq(6, ncol(A) - 1, by = 3)
  ind.10 <- seq(5, ncol(A) - 1, by = 3)
  ind.01 <- seq(7, ncol(A) - 1, by = 3)
  
  # l.proj
  l.proj <- length(n.proj)
  
  # H0
  if (identical(dev, 1)) {
    
    par(mfrow = c(1, length(test)))
    
  } else {
    
    par(mfrow = c(1, length(test) * length(dev)))
    
  }
  for (tt in rev(test)) {
    
    ind.test <- (1 + (tt - 1) * l.proj):(l.proj * tt)
    
    for (d in dev) {
      
      ind.dev <- (1 + (d - 1) * 12):(12 * d)
      
      if (d == 1) {
        
        ylim <- c(0, 0.15)
        
      } else {
        
        ylim <- c(0, 1)
        
      }
      if (d == 1) {
        
        matplot(n.proj, t(A[ind.dev, ind.05[ind.test]]), type = "l", pch = 1, 
                cex = 0.5, col = col.05, ylim = ylim, lty = 1, 
                xlim = c(1, max(n.proj)), main = main, 
                xlab = "Number of projections", ylab = "Empirical rejection rate")
        matlines(n.proj, t(A[ind.dev, ind.10[ind.test]]), type = "l", pch = 1, 
                 cex = 0.5, lty = 1, col = col.10)
        matlines(n.proj, t(A[ind.dev, ind.01[ind.test]]), type = "l", pch = 1, 
                 cex = 0.5, lty = 1, col = col.01)
        abline(h = 0.05 + qnorm(0.975) * sqrt(0.05 * (1 - 0.05)/500) * c(-1, 1), 
               lwd = 3, col = col.05, lty = 2)
        abline(h = 0.1 + qnorm(0.975) * sqrt(0.1 * (1 - 0.1)/500) * c(-1, 1), 
               lwd = 3, col = col.10, lty = 2)
        abline(h = 0.01 + qnorm(0.975) * sqrt(0.01 * (1 - 0.01)/500) * c(-1, 1), 
               lwd = 3, col = col.01, lty = 2)
        
      } else {
        
        matplot(n.proj, t(A[ind.dev, ind.05[ind.test]]), type = "l", pch = 1, 
                ylim = ylim, lty = 1, main = main, xlim = c(1, max(n.proj)), 
                xlab = "Number of projections", ylab = "Empirical rejection rate")
        
      }
    }
  }
  
}


df <- as.data.frame(eg, summary.fun = mean)



# Join results p data-driven
join.results.pdata <- function(results) {
  
  results.pow <- aggregate(x = results[, c(9:(ncol(results) - 1))], by = list(mod = results$n.mod, 
                                                                              delta = results$n.delta, n = results$n.n, method = results$est.method), FUN = function(x) {
                                                                                res <- sapply(c(0.1, 0.05, 0.01), function(alpha) mean(x <= alpha, na.rm = TRUE))
                                                                                names(res) <- c(10, 5, 1)
                                                                                return(res)
                                                                              })
  mean.p <- aggregate(x = results[, ncol(results)], by = list(mod = results$n.mod, delta = results$n.delta, 
                                                              n = results$n.n, method = results$est.method), FUN = function(x) mean(x, na.rm = TRUE))[, 
                                                                                                                                                      5]
  results.pow <- as.data.frame(as.matrix(results.pow))
  results.pow <- cbind(results.pow, mean.p)
  colnames(results.pow)[ncol(results.pow)] <- "mean p"
  
  return(results.pow)
  
}



dontrun <- TRUE
if (!dontrun) {
  
  # Join results p data-driven
  join.results.pdata <- function(results) {
    
    results.pow <- aggregate(x = results[, c(9:(ncol(results) - 1))], by = list(mod = results$n.mod, 
                                                                                delta = results$n.delta, n = results$n.n, method = results$est.method), FUN = function(x) {
                                                                                  res <- sapply(c(0.1, 0.05, 0.01), function(alpha) mean(x <= alpha, na.rm = TRUE))
                                                                                  names(res) <- c(10, 5, 1)
                                                                                  return(res)
                                                                                })
    mean.p <- aggregate(x = results[, ncol(results)], by = list(mod = results$n.mod, delta = results$n.delta, 
                                                                n = results$n.n, method = results$est.method), FUN = function(x) mean(x, na.rm = TRUE))[, 
                                                                                                                                                        5]
    results.pow <- as.data.frame(as.matrix(results.pow))
    results.pow <- cbind(results.pow, mean.p)
    colnames(results.pow)[ncol(results.pow)] <- "mean p"
    
    return(results.pow)
    
  }
  
  # Join results p fixed
  join.results.pfixed <- function(results) {
    
    results.pow <- aggregate(x = results[, c(9:(ncol(results) - 1))], by = list(mod = results$n.mod, 
                                                                                delta = results$n.delta, n = results$n.n, method = results$est.method, p.fixed = results$p.fixed), 
                             FUN = function(x) {
                               res <- sapply(c(0.1, 0.05, 0.01), function(alpha) mean(x <= alpha, na.rm = TRUE))
                               names(res) <- c(10, 5, 1)
                               return(res)
                             })
    results.pow <- as.data.frame(as.matrix(results.pow))
    return(results.pow)
    
  }
  
  # Tables p data-driven
  A <- join.results.pdata(results)
  
  # Only 0.05
  A[, c(1:4, seq(6, ncol(A), by = 3), ncol(A))]
  
  # Only 0.10
  A[, c(1:4, seq(5, ncol(A), by = 3), ncol(A))]
  
  # Only 0.01
  A[, c(1:4, seq(7, ncol(A), by = 3), ncol(A))]
  
  # # Plot H0 pdf('Sizes_for_projs_diffs_Breplicates.pdf',width=14,height=7)
  # plots(A[A[,3]==100,],dev=2,test=1,n.proj=c(1:50),show.main=F) dev.off() # Plot H1
  # pdf('Powers_for_projs.pdf',width=24,height=7) plots(A,dev=2:3,show.main=F) dev.off()
  
  ## Individual plots ##
  {
    n.proj <- c(1:5, 10, 25, 50)
    ## H0 ##
    
    # H0, CvM, n=50
    pdf("sizes_cvm_50.pdf")
    plots(A[A[, 3] == 50, ], dev = 1, test = 1, n.proj = n.proj, main = expression(list(H[0], CvM, 
                                                                                        n == 50)))
    dev.off()
    # H0, CvM, n=100
    pdf("sizes_cvm_100.pdf")
    plots(A[A[, 3] == 100, ], dev = 1, test = 1, n.proj = n.proj, main = expression(list(H[0], CvM, 
                                                                                         n == 100)))
    dev.off()
    # H0, CvM, n=250
    pdf("sizes_cvm_250.pdf")
    plots(A[A[, 3] == 250, ], dev = 1, test = 1, n.proj = n.proj, main = expression(list(H[0], CvM, 
                                                                                         n == 250)))
    dev.off()
    
    # H0, KS, n=50
    pdf("sizes_ks_50.pdf")
    plots(A[A[, 3] == 50, ], dev = 1, test = 2, n.proj = n.proj, main = expression(list(H[0], KS, 
                                                                                        n == 50)))
    dev.off()
    # H0, KS, n=100
    pdf("sizes_ks_100.pdf")
    plots(A[A[, 3] == 100, ], dev = 1, test = 2, n.proj = n.proj, main = expression(list(H[0], KS, 
                                                                                         n == 100)))
    dev.off()
    # H0, KS, n=250
    pdf("sizes_ks_250.pdf")
    plots(A[A[, 3] == 250, ], dev = 1, test = 2, n.proj = n.proj, main = expression(list(H[0], KS, 
                                                                                         n == 250)))
    dev.off()
    
    ## H1 and H2 ##
    
    # H1, CvM, n=50
    pdf("powers_cvm_50_1.pdf")
    plots(A[A[, 3] == 50, ], dev = 2, test = 1, n.proj = n.proj, main = expression(list(H[1], CvM, 
                                                                                        n == 50)))
    dev.off()
    pdf("powers_cvm_50_2.pdf")
    plots(A[A[, 3] == 50, ], dev = 3, test = 1, n.proj = n.proj, main = expression(list(H[2], CvM, 
                                                                                        n == 50)))
    dev.off()
    # H1, CvM, n=100
    pdf("powers_cvm_100_1.pdf")
    plots(A[A[, 3] == 100, ], dev = 2, test = 1, n.proj = n.proj, main = expression(list(H[1], CvM, 
                                                                                         n == 100)))
    dev.off()
    pdf("powers_cvm_100_2.pdf")
    plots(A[A[, 3] == 100, ], dev = 3, test = 1, n.proj = n.proj, main = expression(list(H[2], CvM, 
                                                                                         n == 100)))
    dev.off()
    # H1, CvM, n=250
    pdf("powers_cvm_250_1.pdf")
    plots(A[A[, 3] == 250, ], dev = 2, test = 1, n.proj = n.proj, main = expression(list(H[1], CvM, 
                                                                                         n == 250)))
    dev.off()
    pdf("powers_cvm_250_2.pdf")
    plots(A[A[, 3] == 250, ], dev = 3, test = 1, n.proj = n.proj, main = expression(list(H[2], CvM, 
                                                                                         n == 250)))
    dev.off()
    
    # H1, KS, n=50
    pdf("powers_ks_50_1.pdf")
    plots(A[A[, 3] == 50, ], dev = 2, test = 2, n.proj = n.proj, main = expression(list(H[1], KS, 
                                                                                        n == 50)))
    dev.off()
    pdf("powers_ks_50_2.pdf")
    plots(A[A[, 3] == 50, ], dev = 3, test = 2, n.proj = n.proj, main = expression(list(H[2], KS, 
                                                                                        n == 50)))
    dev.off()
    # H1, KS, n=100
    pdf("powers_ks_100_1.pdf")
    plots(A[A[, 3] == 100, ], dev = 2, test = 2, n.proj = n.proj, main = expression(list(H[1], KS, 
                                                                                         n == 100)))
    dev.off()
    pdf("powers_ks_100_2.pdf")
    plots(A[A[, 3] == 100, ], dev = 3, test = 2, n.proj = n.proj, main = expression(list(H[2], KS, 
                                                                                         n == 100)))
    dev.off()
    # H1, KS, n=250
    pdf("powers_ks_250_1.pdf")
    plots(A[A[, 3] == 250, ], dev = 2, test = 2, n.proj = n.proj, main = expression(list(H[1], KS, 
                                                                                         n == 250)))
    dev.off()
    pdf("powers_ks_250_2.pdf")
    plots(A[A[, 3] == 250, ], dev = 3, test = 2, n.proj = n.proj, main = expression(list(H[2], KS, 
                                                                                         n == 250)))
    dev.off()
    
  }
  
  library(xtable)
  ind.05 <- seq(6, ncol(A) - 1, by = 3)
  ind.k <- c(51, 55, 1, 5)
  
  # First attemp
  mat <- cbind(A[A[, 3] == 50, c(ind.05[c(ind.k, length(ind.05))])], A[A[, 3] == 100, c(ind.05[c(ind.k, 
                                                                                                 length(ind.05))])], A[A[, 3] == 250, c(ind.05[c(ind.k, length(ind.05))])])
  xtable(mat, digits = c(0, rep(3, ncol(mat))))
  
  # Good one
  k <- 1
  for (dev in 1:3) {
    for (mod in 1:12) {
      H <- paste("$H_{", mod, ",", dev - 1, "}$", sep = "")
      cat(paste(H, paste("$", paste(sprintf(fmt = "%.3f", mat[k, ]), collapse = "$ & $"), "$", 
                         sep = ""), sep = " & "), "\\\\\n")
      k <- k + 1
    }
    if (dev == 3) {
      cat("\\bottomrule\\bottomrule\n")
    } else {
      cat("\\midrule\n")
    }
  }
  
}


##################### Time simulation ##

dontrun <- TRUE
if (!dontrun) {
  
  # Fix number of projections
  n.proj <- 5
  
  # For general simulation
  user.simulation.function <- function(model, n, delta, est.method, seed) {
    
    # Fix different seeds
    set.seed(seed)
    
    # Sample
    samp <- tryCatch(rmod(n = n, model = model, delta = delta, R2 = 0.95, comp = TRUE), error = function(e) NA)
    
    # Get the p
    aux <- tryCatch(rp.flm.test(X.fdata = samp$X.fdata, Y = samp$Y, B = 2, est.method = c("pc", "pls", 
                                                                                          "basis")[est.method], n.proj = n.proj, p = NULL, F.code = TRUE, verbose = FALSE, same.rwild = FALSE, 
                                sigma = "vexponential", par.list = list(theta = diff(range(samp$X.fdata$argvals))/20)), error = function(e) NA)
    p <- tryCatch(aux$p, error = function(e) NULL)
    
    # Random projection test
    t.rp <- proc.time()
    test.rp.flm <- tryCatch(rp.flm.test(X.fdata = samp$X.fdata, Y = samp$Y, B = 1000, est.method = c("pc", 
                                                                                                     "pls", "basis")[est.method], n.proj = n.proj, p = p, F.code = TRUE, verbose = FALSE, same.rwild = same.rwild, 
                                        sigma = "vexponential", par.list = list(theta = diff(range(samp$X.fdata$argvals))/20)), error = function(e) NA)
    t.rp <- proc.time() - t.rp
    
    # Flm test
    t.flm <- proc.time()
    test.flm <- tryCatch(flm.test(X.fdata = samp$X.fdata, Y = samp$Y, B = 1000, est.method = c("pc", 
                                                                                               "pls", "basis")[est.method], p = p, verbose = FALSE, plot.it = FALSE), error = function(e) NA)
    t.flm <- proc.time() - t.flm
    
    # Return results
    return(c(t.rp[3], t.flm[3]))
    
  }
  
  
  M <- 100
  nn <- 2^(3:10)
  res <- array(dim = c(length(nn), M, 2))
  pb <- txtProgressBar(style = 3)
  for (ni in seq_along(nn)) {
    cat(ni, "\n")
    for (i in 1:M) {
      m <- i%%12
      d <- i%%2
      if (m == 0) 
        m <- 12
      res[ni, i, ] <- user.simulation.function(model = m, n = nn[ni], delta = d, est.method = 1, 
                                               seed = round(log(i) * ni + 4335))
      setTxtProgressBar(pb = pb, value = i/M)
    }
    cat("Done\n")
  }
  
  # 1: RP 2: FLM
  load("times.RData")
  res <- times
  means <- apply(res, c(1, 3), mean, na.rm = T)
  means
  
  # Times
  pdf("times.pdf")
  plot(nn, means[, 2], type = "o", pch = 19, cex = 0.5, col = 3, xlab = "Sample size", ylab = "Time (seconds)")
  lines(nn, means[, 1], type = "o", pch = 19, cex = 0.5, col = 4)
  legend("topleft", lwd = 2, col = 3:4, legend = c("CvM and KS", "PCvM"))
  dev.off()
  
}






}















