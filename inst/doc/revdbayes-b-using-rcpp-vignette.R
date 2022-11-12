## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(comment = "#>", collapse = TRUE)

## -----------------------------------------------------------------------------
library(revdbayes)
# Is the microbenchmark package available?
got_microbenchmark <- requireNamespace("microbenchmark", quietly = TRUE)
if (got_microbenchmark) {
  library(microbenchmark)
}
# Set the number of posterior samples required.
n <- 1000
set.seed(46)

## -----------------------------------------------------------------------------
u <- quantile(gom, probs = 0.65)
fp <- set_prior(prior = "flat", model = "gp", min_xi = -1)
if (got_microbenchmark) {
  res <- microbenchmark(
    rpost = rpost(n = n, model = "gp", prior = fp, thresh = u, data = gom),
    rpost_rcpp = rpost_rcpp(n = n, model = "gp", prior = fp, thresh = u, 
                            data =   gom)
  )
  print(res, signif = 3)
  options(microbenchmark.unit = "relative")
  print(res, signif = 2)
}  

## -----------------------------------------------------------------------------
mat <- diag(c(10000, 10000, 100))
pn <- set_prior(prior = "norm", model = "gev", mean = c(0,0,0), cov = mat)
if (got_microbenchmark) {
  res <- microbenchmark(
    rpost = rpost(n = n, model = "gev", prior = pn, data = portpirie),
    rpost_rcpp = rpost_rcpp(n = n, model = "gev", prior = pn, 
                            data = portpirie)
  )
}    
options(microbenchmark.unit = NULL)
print(res, signif = 3)
options(microbenchmark.unit = "relative")
print(res, signif = 2)

## -----------------------------------------------------------------------------
# Informative prior set using revdbayes
pr2 <- set_prior(prob = 10^-(1:3), shape = c(38.9, 7.1, 47),
                 scale = c(1.5, 6.3, 2.6), model = "gev", prior = "quant")
if (got_microbenchmark) {
  res <- microbenchmark(
    rpost = rpost(n = n, model = "pp", prior = pr2, data = rainfall,
                  thresh = 40, noy = 54),
    rpost_rcpp = rpost_rcpp(n = n, model = "pp", prior = pr2, 
                            data = rainfall, thresh = 40, noy = 54)
  )
}
options(microbenchmark.unit = NULL)
print(res, signif = 3)
options(microbenchmark.unit = "relative")
print(res, signif = 2)

## -----------------------------------------------------------------------------
# GP model, user-defined prior
ptr_gp_flat <- create_prior_xptr("gp_flat")
p_user <- set_prior(prior = ptr_gp_flat, model = "gp", min_xi = -1)
gpg <- rpost_rcpp(n = 1000, model = "gp", prior = p_user, thresh = u,
                  data = gom)

