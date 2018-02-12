# ================================= ts_post ================================

#' Marginal Bayesian inference for time series extremes
#'
#' Uses the \code{\link[rust]{rust}} package to simulate from the posterior
#' distribution of the extremal index \eqn{\theta} based on the K-gaps model
#' for threshold interexceedance times of Suveges and Davison (2010).
#'
#' @param n A numeric scalar. The size of posterior sample required.
#' @param data A numeric vector of raw data.  No missing values are allowed.
#' @param thresh A numeric scalar.  Extreme value threshold applied to data.
#' @param decluster A logical scalar.  If \code{decluster = FALSE}
#'   (the default) then all threshold excesses are used, as in Fawcett and
#'   Walshaw (2012).  If \code{decluster = TRUE} then only cluster maxima
#'   are used.
#' @param k A numeric scalar.  Run parameter \eqn{K}, as defined in Suveges and
#'   Davison (2010).  Threshold inter-exceedances times that are not larger
#'   than \code{k} units are assigned to the same cluster, resulting in a
#'   \eqn{K}-gap equal to zero.  Specifically, the \eqn{K}-gap \eqn{S}
#'   corresponding to an inter-exceedance time of \eqn{T} is given by
#'   \eqn{S = max(T - K, 0)}.
#' @param inc_cens A logical scalar indicating whether or not to include
#'   contributions from censored inter-exceedance times relating to the
#'   first and last observation.  See Attalides (2015) for details.
#' @param alpha,beta Positive numeric scalars.  Parameters of a
#'   beta(\eqn{\alpha}, \eqn{\beta}) prior for \eqn{\theta}.
#' @param param A character scalar.  If \code{param = "logit"} (the default)
#'   then we simulate from the posterior distribution of
#'   \eqn{\phi = logit(\theta)} and then transform back to the
#'   \eqn{\theta}-scale.  If \code{param = "theta"} then we simulate
#'   directly from the posterior distribution of \eqn{\theta}, unless
#'   the sample K-gaps are all equal to zero or all positive, when we revert
#'   to \code{param = "logit"}.  This is to avoid sampling directly from a
#'   posterior with mode equal to 0 or 1.
#' @param use_rcpp A logical scalar.  If \code{TRUE} (the default) the
#'   rust function \code{\link[rust]{ru_rcpp}} is used for
#'   posterior simulation.  If \code{FALSE} the (slower) function
#'   \code{\link[rust]{ru}} is used.
#' @details Add details
#' @return Explain
#' @references Fawcett, L. and Walshaw, D. (2012) Estimating return levels
#'   from serially dependent extremes, \emph{Environmetrics}, \strong{23}(3),
#'   272-283. \url{http://dx.doi.org/10.1002/env.2133}
#' @references Suveges, M. and Davison, A. C. (2010) Model
#'   misspecification in peaks over threshold analysis, \emph{The Annals of
#'   Applied Statistics}, \strong{4}(1), 203-221.
#'   \url{http://dx.doi.org/10.1214/09-AOAS292}
#' @references Attalides, N. (2015) Threshold-based extreme value modelling,
#'   PhD thesis, University College London.
#'   \url{http://discovery.ucl.ac.uk/1471121/1/Nicolas_Attalides_Thesis.pdf}
#' @seealso \code{\link{kgaps_post}} for random sampling from K-gaps posterior
#'   distribution.
#' @seealso \code{\link{rpost}} and \code{\link{rpost_rcpp}} for sampling
#'   from an extreme value posterior distribution.
#' @examples
#' thresh <- quantile(newlyn, probs = 0.90)
#' x <- ts_post(newlyn, thresh)
ts_post <- function(data, thresh, decluster = FALSE, k = 1, n = 1000,
                    inc_cens = FALSE, alpha = 1, beta = 1,
                    param = c("logit", "theta"), use_rcpp = TRUE) {
  # Simulate a sample of size n from the posterior for the extremal index theta
  theta_sim <- kgaps_post(data, thresh, k, n, inc_cens, alpha, beta, param,
                          use_rcpp)
  theta <- theta_sim$sim_vals
  if (decluster) {
    # Decluster the data based on the value of k

  }
  return(theta_sim)
}

#' Declustering
#'
#' Describe
#' @param data A numeric vector of raw data.  No missing values are allowed.
#' @param thresh A numeric scalar.  Extreme value threshold applied to data.
#' @param k An integer scalar
#' @examples
#' thresh <- quantile(newlyn, probs = 0.99)
#' x <- find_clusters(newlyn, thresh)
find_clusters <- function(data, thresh, k = 1) {
  if (any(is.na(data))) {
    stop("No missing values are allowed in ''data''")
  }
  # Sample size, positions, number and proportion of exceedances
  nx <- length(data)
  exc_u <- (1:nx)[data > thresh]
  N_u <- length(exc_u)
  q_u <- N_u / nx
  # Inter-exceedances times and K-gaps
  T_u <- diff(exc_u)
  S_k <- pmax(T_u - k,0)
  # Threshold inter-exceedances times that are not larger than k units
  # (i.e. a K-gap equal to zero) are assigned to the same cluster
  print(exc_u)
  print(which(S_k == 0))
  print(which(S_k > 0))
  # A cluster starts with an exceedance that is preceded by a positive K-gaps
  # and ends with the next exceedance (possibly the exceedance itself) that
  # is followed by a positive K-gap
  clus <- numeric(nx)
  #
  exc_ind <- data > thresh
  print(exc_ind)
  if (k == 1) {
    cluster_durations <- rle(exc_ind)
  } else {
    print("HERE")
    print(zoo::rollsum(exc_ind, k = k))
    print(as.numeric(zoo::rollsum(exc_ind, k = k) > 0))
    cluster_durations <- rle(zoo::rollsum(exc_ind, k = k) > 0)
    which_to_reduce <- cluster_durations$values > 0
    cluster_durations$lengths[which_to_reduce] <-
      cluster_durations$lengths[which_to_reduce] - k + 1
  }
  return(cluster_durations)
}
