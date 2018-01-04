#' Multinomial cumulative distribution function (CDF)

#' Multinomial cumulative distribution function (CDF)
#'
#' i.e. If n balls are thrown into k bins, calculate
#' P(x1<=q1, x2<=q2, ..., xk<=qk) where xi are numbers of balls in each bin, and
#' qi are quantiles
#' Uses an approximation described in: \cr
#' Levin,  B.:  A  representation  for  multinomial  cumulative  distribution
#' functions. Ann. Stat. 9, 1123–1126 (1981)
#' @param q Vector of quantiles; If q has length 1 then the same quantile is
#'          used for each bin. If q is length k, each bin is assigned its own
#'          quantile
#' @param n Total number of balls
#' @param prob A vector of length k, describing the probability of a ball being
#'             placed into bin xk
#' @param log.p If TRUE, natural logarithm of probability is returned
#'        (Default: FALSE)
#' @param lower.tail logical; if TRUE (default), probabilities are P[X <= x],
#'                   otherwise, P[X > x].
#' @return Numeric; A probability value. Returns the string "<2.2e-16" for
#'         extremely small probabilities
#' @examples
#' # For 4 balls randomly thrown into 2 bins, probability that all bins end up
#' # with 2 or fewer balls
#' pmultinom(q = 2, n = 4, prob = c(0.5, 0.5))
#'
#' # For 6 balls randomly thrown into 3 bins, probability that at least one bin
#' # ends up with 3 or more balls
#' pmultinom(q = 2, n = 6, c(1/3, 1/3, 1/3))
#'
#' # Same as above, but with different probabilities for each bin
#' pmultinom(q = 2, n = 6, c(0.2, 0.3, 0.5))
#'
#' # For 4 balls randomly thrown into 2 bins, probability that the first bin has
#' # 0 balls and the second bin has 4
#' pmultinom(q = c(0, 4), n = 4, c(0.5, 0.5))
#' @seealso \code{\link{dmultinom}}, \code{\link{pbinom}}
#' @export
pmultinom <- function(q, n, prob, log.p = FALSE, lower.tail = TRUE) {
  if (length(q) == 1) {
    q <- rep(q, length(prob))
  }

  if (any(q > n)) {
    # If q >= n then xi <= q for all i
    p <- 0
  } else {
    p <- cdf_multinomial_lnP(K = length(prob), N = n, p = prob, n = q)

    if (is.nan(p)) {
      # If NaN is returned, p = 0
      p <- -Inf
    }

    if (p >= 0) {
      return(handle_overflow(lower.tail))
    }
  }

  if (!lower.tail) {
    # Calculate log(1 - exp(p)) without loss of precision
    p <- calculate_upper_tail(p)
  }

  if (!log.p) {
    p <- calculate_exp_p(p)
  }
  p
}

#' Convert the p-value for the lower tail to the upper tail p-value
#' @importFrom VGAM log1mexp
calculate_upper_tail <- function(plower) {
  log1mexp(-plower)
}

#' Calculate exp(p) and handle any underflow errors
calculate_exp_p <- function(logp) {
    p <- exp(logp)
    if (p < 0) {
      p <- "<2.2e-16"
    }
    p
}

#' Return the correct p-value in the case of an overflow error
handle_overflow <- function(lower.tail) {
  if (lower.tail) {
    p <- 1
  } else {
    p <- "<2.2e-16"
  }
  p
}
