rp.flm.test
===========

[![Travis-CI Build Status](https://travis-ci.org/egarpor/rp.flm.test.svg?branch=master)](https://travis-ci.org/egarpor/rp.flm.test)
[![License](https://img.shields.io/badge/license-MIT%20License-brightgreen.svg)](https://opensource.org/licenses/MIT)

## Overview

Package companion for the paper *Goodness-of-fit tests for the functional linear model based on randomly projected empirical processes* (Cuesta-Albertos *et al.*, 2017). It implements the proposed tests and allows to replicate the empirical results presented.

## Install

```
library(devtools)
install_github("egarpor/rp.flm.test")
```

Alternatively, see function `rp.flm.test` in the [`fda.usc`](http://cran.r-project.org/web/packages/fda.usc/) library (Febrero-Bande and Oviedo de la Fuente, 2012) for versions above 1.3.1.

## Usage

### Simulated data
```
# Generate data
set.seed(345678)
t <- seq(0, 1, l = 101)
n <- 100
X <- r.ou(n = n, t = t, alpha = 2, sigma = 0.5)
beta0 <- fdata(mdata = cos(2 * pi * t) - (t - 0.5)^2, argvals = t,
               rangeval = c(0,1))
Y1 <- inprod.fdata(X, beta0) + rnorm(n, sd = 0.1)
Y2 <- log(norm.fdata(X)) + rnorm(n, sd = 0.1)

# Do not reject FLM
rp.test1 <- rp.flm.test(X.fdata = X, Y = Y1, n.proj = c(1, 5))
rp.test1$p.values.fdr

# Reject FLM
rp.test2 <- rp.flm.test(X.fdata = X, Y = Y2, n.proj = c(1, 5))
rp.test2$p.values.fdr

# Estimations of beta
par(mfrow = c(1, 2))
plot(beta0, main = "beta")
lines(rp.test1$beta.est, col = 2)
plot(beta0, main = "beta")
lines(rp.test2$beta.est, col = 2)
```

### Tecator dataset
```
# Load data
data(tecator)
absorp <- tecator$absorp.fdata
ind <- 1:129 # or ind <- 1:215
x <- absorp[ind, ]
y <- tecator$y$Fat[ind]

# Composite hypothesis
rp.tecat <- rp.flm.test(X.fdata = x, Y = y)
rp.tecat$p.values.fdr[c(5, 10), ]

# Simple hypothesis
zero <- fdata(mdata = rep(0, length(x$argvals)), argvals = x$argvals,
              rangeval = x$rangeval)
rp.flm.test(X.fdata = x, Y = y, beta0.fdata = zero)$p.values.fdr
```

### AEMET dataset
```
# Load data
data(aemet)
wind.speed <- apply(aemet$wind.speed$data, 1, mean)
temp <- aemet$temp

# Remove the 5% of the curves with less depth (i.e. 4 curves)
par(mfrow = c(1, 1))
res.FM <- depth.FM(temp, draw = TRUE)
qu <- quantile(res.FM$dep, prob = 0.05)
l <- which(res.FM$dep <= qu)
lines(aemet$temp[l], col = 3)

# Data without outliers
wind.speed <- wind.speed[-l]
temp <- temp[-l]

# Composite hypothesis
rp.aemet <- rp.flm.test(X.fdata = temp, Y = wind.speed)
rp.aemet$p.values.fdr
apply(rp.aemet$p.values.fdr, 2, range)
```

## References

Cuesta-Albertos, J. A., García-Portugués, E., Febrero-Bande, M. and González-Manteiga, W. (2017). Goodness-of-fit tests for the functional linear model based on randomly projected empirical processes. *arXiv:1701.XXXXX*. <https://arxiv.org/abs/1701.XXXXX>

Febrero-Bande, M. and Oviedo de la Fuente, M. (2012). Statistical Computing in Functional Data Analysis: The R Package fda.usc. *Journal of Statistical Software*, 51(4), 1-28. <http://www.jstatsoft.org/v51/i04/>

