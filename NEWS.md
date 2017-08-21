# revdbayes 1.2.1

## Bug fixes and minor improvements

* In some extreme cases (datasets with very small numbers of threshold excesses) calling `predict.evpost` with `type = "q"` and `x` close to 1 returns an imprecise value for the requested predictive quantiles.  This has been corrected by using `stats::uniroot` rather than `stats::nlminb`.

* A bug (missing `drop = FALSE` in subsetting a matrix) in `plot.evpred` produced an error message if `n_years` was scalar in the prior call to `predict.evpost`.  This bug has been corrected.

* The placing of ... in the function definitions of `rpost` and `rpost_rcpp` meant that it was not possible to supply the argument `r` to be passed to `rust::ru` or `rust::ru_rcpp` to change the ratio-of-uniforms tuning parameter `r`.  Furthermore, if `model = "os"` then trying to do this sets `ros` in error. This has been corrected.

* A bug meant that the values returned by `predict(evpost_object, type = "d")` being incorrect if `evpost_object` was returned from a call to `rpost` using `model = bingp`.  The values returned were too small: they differ from the correct values by a factor approximately equal to the proportion of observations that lie above the threshold.  This bug has been corrected.

# revdbayes 1.2.0

## New features

* Faster computation, owing to the use of packages Rcpp and RcppArmadillo in package rust (https://CRAN.R-project.org/package=rust).

* New function: `rpost_rcpp`.  

* New vignette. "Faster simulation using revdbayes".

* `set_prior` has been extended so that informative priors for GEV parameters can be specified using the arguments `prior = "prob"` or `prior = "quant"`.  It is no longer necessary to use the functions `prior.prob` and `prior.quant` from the evdbayes package to set these priors.

## Bug fixes and minor improvements

* The list returned from `set_prior` now contains default values for all the required arguments of a given in-built prior, if these haven't been specified by the user.  This simplifies the evaluation of prior densities using C++.

* The GEV functions `dgev`, `pgev`, `qgev`, `rgev` and the GP functions `dgp`, `pgp`, `qgp`, `rgp` have been rewritten to conform with the vectorised style of the standard functions for distributions, e.g. those found at `?Normal`.  This makes these functions more flexible, but also means that the user take care when calling them with vectors arguments or different lengths.

* The documentation for `rpost` has been corrected: previously it stated that the default for `use_noy` is `use_noy = FALSE`, when in fact it is `use_noy = TRUE`.

* Bug fixed in `plot.evpost` : previously, in the `d = 2` case, providing the graphical parameter `col` produced an error because `col = 8` was hard-coded in a call to `points`. Now the extra argument `points_par` enables the user to provide a list of arguments to `points`.

* All the (R, not C++) prior functions described in the documentation of `set_prior` are now exported.  This means that they can now be used in the function `posterior` in the `evdbayes` package.

* Unnecessary dependence on package `devtools` via Suggests is removed.

* Bugs fixed in the (R) prior functions `gp_norm`, `gev_norm` and `gev_loglognorm`.  The effect of the bug was negligible unless the prior variances are not chosen to be large.

* In a call to `rpost` or `rpost_rcpp` with `model = "os"` the user may provide `data` in the form of a vector of block maxima.  In this instance the output is equivalent to a call to these functions with `model = "gev"` with the same data. 

# revdbayes 1.1.0

## New features

* A new vignette (Posterior Predictive Extreme Value Inference using the revdbayes Package) provides an overview of most of the new features. Run browseVignettes("revdbayes") to access.

* S3 `predict()` method for class 'evpost' performs predictive inference about the largest observation observed in N years, returning an object of class `evpred`.
  
* S3 `plot()` for the `evpred` object returned by `predict.evpost`.

* S3 `pp_check()` method for class 'evpost' performs posterior predictive checks using the bayesplot package.

* Interface to the bayesplot package added in the S3 `plot.evpost` method.

* `model = bingp` can now be supplied to `rpost()` to add inferences about the probability of threshold exceedance to inferences about threshold excesses based on the Generalised Pareto (GP) model.  `set_bin_prior()` can be used to set a prior for this probability.

* `rprior_quant()`: to simulate from the prior distribution for GEV parameters proposed in Coles and Tawn (1996) [A Bayesian analysis of extreme rainfall data. Appl. Statist., 45, 463-478], based on independent gamma priors for differences between quantiles.  
   
* `prior_prob()`: to simulate from the prior distribution for GEV parameters based on Crowder (1992), in which independent beta priors are specified for ratios of probabilities (which is equivalent to a Dirichlet prior on differences between these probabilities).

## Bug fixes and minor improvements

* The spurious warning messages relating to checking that the model argument to `rpost()` is consistent with the prior set using `set-prior()` have been corrected.  These occurred when `model = "pp"` or `model = "os"`.
  
* The hyperparameter in the MDI prior was `a` in the documentation and `a_mdi` in the code.  Now it is `a` everywhere.
  
* In `set_prior` with `prior = "beta"` parameter vector `ab` has been corrected to `pq`.
  
* In the documentation of `rpost()` the description of the argument `noy` has been corrected.
  
* Package spatstat removed from the Imports field in description to avoid NOTE in CRAN checks.  
