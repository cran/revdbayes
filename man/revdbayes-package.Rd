% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/revdbayes-package.R
\docType{package}
\name{revdbayes-package}
\alias{revdbayes}
\alias{revdbayes-package}
\title{revdbayes: Ratio-of-Uniforms Sampling for Bayesian Extreme Value Analysis}
\description{
Uses the multivariate generalized ratio-of-uniforms method to simulate
random samples from the posterior distributions commonly encountered in
Bayesian extreme value analyses.
}
\details{
The main functions in the revdbayes package are \code{\link{rpost}}
  and \code{\link{rpost_rcpp}}, which simulate random samples from the
  posterior distribution of extreme value model parameters using the
  functions \code{\link[rust]{ru}} and \code{\link[rust]{ru_rcpp}}
  from the rust package, respectively. The user chooses the extreme value
  model, the prior density for the parameters and provides the data.
  There are options to improve the probability of acceptance of the
  ratio-of-uniforms algorithm by working with transformation of the model
  parameters.

  The functions \code{\link{kgaps_post}} and \code{\link{dgaps_post}}
  simulate from the posterior distribution of the extremal index
  \eqn{\theta} based on the K-gaps model for threshold interexceedance
  times of Suveges and Davison (2010) and the similar D-gaps model of
  Holesovsky and Fusek (2020).  See also Attalides (2015).

  See \code{vignette("revdbayes-a-vignette", package = "revdbayes")} for an
  overview of the package and
  \code{vignette("revdbayes-b-using-rcpp-vignette", package = "revdbayes")}
   for an illustration of the improvements in efficiency produced using
   the Rcpp package.
   See
   \code{vignette("revdbayes-c-predictive-vignette", package = "revdbayes")}
   for an outline of how to use revdbayes to perform posterior predictive
   extreme value inference and
   \code{vignette("revdbayes-d-kgaps-vignette", package = "revdbayes")}
   considers Bayesian inference for the extremal index \eqn{\theta}
   using threshold inter-exceedance times.
}
\references{
Holesovsky, J. and Fusek, M. Estimation of the extremal index
  using censored distributions. Extremes 23, 197-213 (2020).
  \doi{10.1007/s10687-020-00374-3}

Northrop, P. J. (2016). rust: Ratio-of-Uniforms Simulation with
  Transformation. R package version 1.2.2.
  \url{https://cran.r-project.org/package=rust}.

Suveges, M. and Davison, A. C. (2010) Model
  misspecification in peaks over threshold analysis, \emph{The Annals of
  Applied Statistics}, \strong{4}(1), 203-221.
  \doi{10.1214/09-AOAS292}

Attalides, N. (2015) Threshold-based extreme value modelling,
  PhD thesis, University College London.
  \url{https://discovery.ucl.ac.uk/1471121/1/Nicolas_Attalides_Thesis.pdf}
}
\seealso{
\code{\link{set_prior}} to set a prior density for extreme value
  parameters.

\code{\link{rpost}} and \code{\link{rpost_rcpp}} to perform
  ratio-of-uniforms sampling from an extreme value posterior distribution.

\code{\link{kgaps_post}} and \code{\link{dgaps_post}} to sample
  from a posterior distribution for the extremal index based on
  inter-exceedance times.

The \code{\link[rust]{ru}} and \code{\link[rust]{ru_rcpp}}
  functions in the \code{\link[rust]{rust}} package for details of the arguments
  that can be passed to \code{ru} via \code{rpost} and for the form of the
  object (of class \code{"evpost"}) returned from \code{rpost}, which has the same
  structure as an object (of class \code{"ru"}) returned by \code{ru} and
  \code{ru_rcpp}.
}
\author{
\strong{Maintainer}: Paul J. Northrop \email{p.northrop@ucl.ac.uk} [copyright holder]

Other contributors:
\itemize{
  \item Scott D. Grimshaw [contributor]
}

}
