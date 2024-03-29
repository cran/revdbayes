## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(comment = "#>", collapse = TRUE)

## -----------------------------------------------------------------------------
library(revdbayes)
# Set the number of posterior samples required.
n <- 10000
set.seed(46)

## -----------------------------------------------------------------------------
thresh <- quantile(gom, probs = 0.65)
fp <- set_prior(prior = "flat", model = "gp", min_xi = -1)

## -----------------------------------------------------------------------------
t1 <- system.time(
gp1 <- rpost(n = n, model = "gp", prior = fp, thresh = thresh, data = gom,
             rotate = FALSE)
)[3]

## -----------------------------------------------------------------------------
t2 <- system.time(
gp2 <- rpost(n = n, model = "gp", prior = fp, thresh = thresh, data = gom)
)[3]

## -----------------------------------------------------------------------------
t3 <- system.time(
gp3 <- rpost(n = n, model = "gp", prior = fp, thresh = thresh, data = gom,
             rotate = FALSE, trans = "BC")
)[3]
t4 <- system.time(
gp4 <- rpost(n = n, model = "gp", prior = fp, thresh = thresh, data = gom,
             trans = "BC")
)[3]

## ----fig.show='hold', out.width="45%"-----------------------------------------
plot(gp1, ru_scale = FALSE, cex.main = 0.75, cex.lab = 0.75, 
  main = paste("no transformation \n pa = ", round(gp1$pa, 3), 
               ", time = ", round(t1, 2), "s"))
plot(gp2, ru_scale = TRUE, cex.main = 0.75, cex.lab = 0.75, 
  main = paste("rotation \n pa = ", round(gp2$pa, 3),
               ", time = ", round(t2, 2), "s"))
plot(gp3, ru_scale = TRUE, cex.main = 0.75, cex.lab = 0.75, 
  main = paste("Box-Cox \n pa = ", round(gp3$pa, 3),
               ", time = ", round(t3, 2), "s"))
plot(gp4, ru_scale = TRUE, cex.main = 0.75, cex.lab = 0.75,
  main = paste("Box-Cox and rotation \n pa = ", round(gp4$pa, 3),
               ", time = ", round(t4, 2), "s"))

## ----echo = TRUE--------------------------------------------------------------
thresh <- quantile(gom, probs = 0.95)
fp <- set_prior(prior = "flat", model = "gp", min_xi = -1)
t2 <- system.time(
gp2 <- rpost(n = n, model = "gp", prior = fp, thresh = thresh, data = gom)
)[3]
t4 <- system.time(
gp4 <- rpost(n = n, model = "gp", prior = fp, thresh = thresh, data = gom,
             trans = "BC")
)[3]

## ----fig.show='hold', echo = FALSE, out.width="45%"---------------------------
plot(gp2, ru_scale = TRUE, cex.main = 0.75, cex.lab = 0.75, 
  main = paste("rotation \n pa = ", round(gp2$pa, 3),
               ", time = ", round(t2, 2), "s"))
plot(gp4, ru_scale = TRUE, cex.main = 0.75, cex.lab = 0.75,
  main = paste("Box-Cox and rotation \n pa = ", round(gp4$pa, 3),
               ", time = ", round(t4, 2), "s"))

## -----------------------------------------------------------------------------
u_prior_fn <- function(x, ab) {
  #
  # Calculates the the log of the (improper) prior density for GP parameters 
  # (sigma_u, xi) in which log(sigma_u) is uniform on the real line and xi has 
  # a beta(alpha, beta)-type prior on the interval (-1, 1).
  #
  # Args:
  #   x  : A numeric vector.  GP parameter vector (sigma, xi).
  #   ab : A numeric vector.  Hyperparameter vector (alpha, beta), where
  #        alpha and beta must be positive.
  #
  # Returns : the value of the log-prior at x.
  #
  if (x[1] <= 0 | x[2] <= -1 | x[2] >= 1) {
    return(-Inf)
  }  
  return(-log(x[1]) + (ab[1] - 1) * log(1 + x[2]) +
           (ab[2] - 1) * log(1 - x[2]))
}
up <- set_prior(prior = u_prior_fn, ab = c(2, 2), model = "gp")
gp_u <- rpost(n = n, model = "gp", prior = up, thresh = thresh, data = gom)

## -----------------------------------------------------------------------------
data(portpirie)
mat <- diag(c(10000, 10000, 100))
pn <- set_prior(prior = "norm", model = "gev", mean = c(0,0,0), cov = mat)
p1 <- rpost(n = n, model = "gev", prior = pn, data = portpirie)
p1$pa

## ----fig.width = 7------------------------------------------------------------
plot(p1, cex.main = 0.75, cex.lab = 0.75, 
  main = "original scale")
plot(p1, ru_scale = TRUE, cex.main = 0.75, cex.lab = 0.75,
     main = "sampling scale")

## ----fig.width = 7------------------------------------------------------------
oldpar <- par(mfrow = c(1,3))
plot(density(p1$sim_vals[, 1], adj = 2), main = "", xlab = expression(mu))
plot(density(p1$sim_vals[, 2], adj = 2), main = "", xlab = expression(sigma))
plot(density(p1$sim_vals[, 3], adj = 2), main = "", xlab = expression(xi))
par(oldpar)

## -----------------------------------------------------------------------------
data(oxford)
prox <- set_prior(prior = "prob", model = "gev", quant = c(85, 88, 95), 
                  alpha = c(4, 2.5, 2.25, 0.25))
# revdbayes
ox_post <- rpost(n = 1000, model = "gev", prior = prox, data = oxford)

## ----fig.width = 7------------------------------------------------------------
oldpar <- par(mfrow = c(1,3))
plot(density(ox_post$sim_vals[, 1], adj = 2), main = "", xlab = expression(mu))
plot(density(ox_post$sim_vals[, 2], adj = 2), main = "", xlab = expression(sigma))
plot(density(ox_post$sim_vals[, 3], adj = 2), main = "", xlab = expression(xi))
par(oldpar)

## -----------------------------------------------------------------------------
data(rainfall)
rthresh <- 40
prrain <- set_prior(prior = "quant", model = "gev", prob = 10^-(1:3), 
                    shape = c(38.9, 7.1, 47), scale = c(1.5, 6.3, 2.6))

## ----fig.width = 7------------------------------------------------------------
# Annual maximum parameterisation (number of blocks = number of years of data)
rn0 <- rpost(n = n, model = "pp", prior = prrain, data = rainfall, 
  thresh = rthresh, noy = 54, rotate = FALSE)
plot(rn0)
rn0$pa
# Rotation about maximum a posteriori estimate (MAP)
system.time(
rn1 <- rpost(n = n, model = "pp", prior = prrain, data = rainfall, 
  thresh = rthresh, noy = 54)
)
plot(rn1, ru_scale = TRUE)
rn1$pa

## ----fig.width = 7------------------------------------------------------------
# Number of blocks = number of threshold excesses
n_exc <- sum(rainfall > rthresh, na.rm = TRUE)
rn2 <- rpost(n = n, model = "pp", prior = prrain, data = rainfall, 
  thresh = rthresh, noy = 54, use_noy = FALSE, rotate = FALSE)
rn2$pa
plot(rn2, ru_scale = TRUE)
# Number of blocks = number of threshold excesses, rotation about MAP
system.time(
rn3 <- rpost(n = n, model = "pp", prior = prrain, data = rainfall,
  thresh = rthresh, noy = 54, use_noy = FALSE)
)
plot(rn3, ru_scale = TRUE)
rn3$pa

## ----fig.width = 7------------------------------------------------------------
oldpar <- par(mfrow = c(1,3))
plot(density(rn1$sim_vals[, 1], adj = 2), main = "", xlab = expression(mu))
plot(density(rn1$sim_vals[, 2], adj = 2), main = "", xlab = expression(sigma))
plot(density(rn1$sim_vals[, 3], adj = 2), main = "", xlab = expression(xi))
par(oldpar)

## ----fig.width = 7------------------------------------------------------------
data(venice)
mat <- diag(c(10000, 10000, 100))
pv <- set_prior(prior = "norm", model = "gev", mean = c(0,0,0), cov = mat)
osv <- rpost(n = n, model = "os", prior = pv, data = venice)
plot(osv)

## ----fig.width = 7------------------------------------------------------------
plot(osv, ru_scale = TRUE)

