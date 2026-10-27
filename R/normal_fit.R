#' Normal model for the distribution of true effects
#'
#' Fits the model in which the true effects are normally distributed,
#' \eqn{\theta_i \sim N(\mu, \tau^2)}, and each estimate is
#' \eqn{y_i \sim N(\theta_i, s_i^2)} with known standard error \eqn{s_i}.
#' The mean \eqn{\mu} and SD \eqn{\tau} of the true effects are estimated by
#' maximum likelihood, so that \eqn{y_i \sim N(\mu, \tau^2 + s_i^2)}. This is
#' the parametric counterpart of [em_smooth()] and a useful baseline.
#'
#' @param estimate Numeric vector of estimated effects.
#' @param se Numeric vector of standard errors (positive, same length as
#'   `estimate`), treated as known.
#' @return An object of class `"normal_fit"`, a list with the estimated `mu`
#'   and `tau`, their standard errors `se_mu` and `se_tau` (from the observed
#'   information, treating estimates as independent), `loglik`, and `data`.
#'   Use [summary.normal_fit()] for the mean, SD, quantiles and exceedance
#'   probabilities, computed exactly from the normal distribution.
#' @examples
#' set.seed(1)
#' se <- runif(300, 0.05, 0.4)
#' y <- rnorm(300, rnorm(300, 0.2, 0.2), se)
#' fit <- normal_fit(y, se)
#' fit
#' density(fit, at = c(0, 0.2, 0.4))
#' @export
normal_fit <- function(estimate, se) {
  check_inputs(estimate, se)
  nll <- function(p) {
    -sum(stats::dnorm(estimate, p[1], sqrt(exp(2 * p[2]) + se^2),
                      log = TRUE))
  }
  start <- c(mean(estimate), log(max(stats::sd(estimate), 1e-3)))
  o <- stats::optim(start, nll, method = "BFGS", hessian = TRUE)
  mu <- o$par[1]
  tau <- exp(o$par[2])
  vc <- tryCatch(solve(o$hessian), error = function(e) matrix(NA_real_, 2, 2))
  se_mu <- sqrt(vc[1, 1])
  se_tau <- tau * sqrt(vc[2, 2])           # delta method for exp()
  structure(list(mu = mu, tau = tau, se_mu = se_mu, se_tau = se_tau,
                 loglik = -o$value,
                 data = data.frame(estimate = estimate, se = se),
                 call = match.call()),
            class = "normal_fit")
}

#' Summary of a fitted normal model
#'
#' @param object An object returned by [normal_fit()].
#' @param probs Probabilities at which to report quantiles.
#' @param thresholds Values `t` at which to report the probability that a
#'   true effect exceeds `t`.
#' @param ... Unused.
#' @return A data frame with one row per statistic (mean, SD, quantiles,
#'   exceedance probabilities) and a column `estimate`, computed exactly from
#'   the fitted normal distribution. The format matches [summary.em_smooth()].
#' @export
summary.normal_fit <- function(object,
                               probs = c(0.10, 0.25, 0.50, 0.75, 0.90),
                               thresholds = 0, ...) {
  s <- c(mean = object$mu, sd = object$tau,
         stats::setNames(stats::qnorm(probs, object$mu, object$tau),
                         paste0("q", probs * 100)),
         stats::setNames(
           stats::pnorm(thresholds, object$mu, object$tau,
                        lower.tail = FALSE),
           paste0("p_gt_", thresholds)))
  data.frame(statistic = names(s), estimate = unname(s),
             row.names = NULL, stringsAsFactors = FALSE)
}

#' @export
print.normal_fit <- function(x, digits = 3, ...) {
  cat("Normal model for the distribution of true effects\n")
  cat("  estimates:", nrow(x$data), "\n")
  cat("  mu  =", round(x$mu, digits), "(SE", round(x$se_mu, digits), ")\n")
  cat("  tau =", round(x$tau, digits), "(SE", round(x$se_tau, digits),
      ")\n\n")
  s <- summary(x)
  print(stats::setNames(round(s$estimate, digits), s$statistic))
  invisible(x)
}

#' Density of the fitted normal distribution of true effects
#'
#' @param x An object returned by [normal_fit()].
#' @param at Numeric vector of values at which to evaluate the density.
#' @param ... Unused.
#' @return A numeric vector of densities, one per value of `at`.
#' @importFrom stats density
#' @export
density.normal_fit <- function(x, at, ...) {
  stats::dnorm(at, x$mu, x$tau)
}
