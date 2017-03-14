
<!-- README.md is generated from README.Rmd. Please edit that file -->
revdbayes: Ratio-of-uniforms Sampling for Bayesian Extreme Value Analysis
-------------------------------------------------------------------------

### What does revdbayes do?

The `revdbayes` package uses the ratio-of-uniforms method to produce random samples from the posterior distributions that occur in some relatively simple Bayesian extreme value analyses. The functionality of revdbayes is similar to the evdbayes package <https://cran.r-project.org/package=evdbayes>, which uses Markov Chain Monte Carlo (MCMC) methods for posterior simulation.

### A simple example

The two main functions in `revdbayes` are `set_prior` and `rpost`. `set_prior` sets a prior for extreme value parameters. `rpost` samples from the posterior produced by updating this prior using the likelihood of observed data under an extreme value model. The following code sets a prior for Generalised Extreme Value (GEV) parameters based on a multivariate normal distribution and then simulates a random sample of size 1000 from the posterior distribution based on a dataset of annual maximum sea levels.

``` r
data(portpirie)
mat <- diag(c(10000, 10000, 100))
pn <- set_prior(prior = "norm", model = "gev", mean = c(0,0,0), cov = mat)
gevp  <- rpost(n = 1000, model = "gev", prior = pn, data = portpirie)
plot(gevp)
```

### Installation

To get the current released version from CRAN:

``` r
install.packages("revdbayes")
```

### Vignette

See `vignette("revdbayes-vignette", package = "revdbayes")` for an overview of the package and `vignette("revdbayes-predictive-vignette", package = "revdbayes")` for an outline of how to use revdbayes to perform posterior predictive extreme value inference.
