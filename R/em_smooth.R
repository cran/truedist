#' Smoothed EM estimate of the distribution of true effects
#'
#' Estimates the distribution \eqn{G} of latent true effects
#' \eqn{\theta_i} from estimates \eqn{y_i} that are normally distributed
#' around \eqn{\theta_i} with known standard errors \eqn{s_i}, which may
#' differ from one estimate to the next. \eqn{G} is a nonparametric mixing
#' distribution on a grid, estimated by the smoothed EM (EMS) algorithm of
#' Silverman, Jones, Wilson and Nychka (1990): each EM update of the grid
#' weights is followed by a Gaussian smoothing step.
#'
#' The smoothing step uses a Gaussian kernel with one constant standard
#' deviation, chosen from the data by a rule of thumb:
#' \deqn{\mathrm{SD} = c \times \mathrm{median}(s_i) \times
#' \sqrt{\log 500 / \log n}.}
#' The SD grows with the typical standard error and, slowly, as the sample
#' size \eqn{n} falls, following the logarithmic dependence on \eqn{n} that
#' theory gives for deconvolution (Fan, 1991). The constant \eqn{c} is 0.1 by
#' default, so that with \eqn{n = 500} and a typical standard error of 0.15
#' the SD is about 0.015. In simulations with \eqn{n = 500} and standard
#' errors between 0.05 and 0.4, this recovered the true CDF as well as or better
#' than 0.2, 0.33 or 0.5. This is a rule of thumb, not an optimal bandwidth: the
#' best amount of smoothing depends on how smooth the true distribution is,
#' which the data cannot reliably reveal. Features such as peak heights are
#' sensitive to it, while the mean, SD, quantiles and CDF are much less so;
#' consider rerunning with other values of `smoothing` as a sensitivity check.
#'
#' @param estimate Numeric vector of estimated effects.
#' @param se Numeric vector of standard errors (positive, same length as
#'   `estimate`), treated as known.
#' @param smoothing The constant \eqn{c} in the rule given in the Details
#'   section (at least 0.05).
#'   Larger values smooth more. With fewer than 3 estimates the function stops,
#'   and with fewer than 20 it warns that the result is very uncertain.
#'
#' @return An object of class `"em_smooth"`, a list with components
#'   `grid`, `weights` (probability mass on each grid point), `density`
#'   (`weights` divided by the grid step), `loglik`, `iterations`,
#'   `converged`, `smoothing` (the value of `smoothing` used), `smoothing_sd`
#'   (the standard deviation of the smoothing kernel), and `data` (the input
#'   estimates and standard errors).
#'   Use [summary.em_smooth()] for the mean, SD, quantiles and exceedance
#'   probabilities.
#'
#' @references
#' Silverman, B. W., Jones, M. C., Wilson, J. D., and Nychka, D. W. (1990).
#' A smoothed EM approach to indirect estimation problems, with particular
#' reference to stereology and emission tomography. *Journal of the Royal
#' Statistical Society, Series B*, 52(2), 271-303.
#' \doi{10.1111/j.2517-6161.1990.tb01788.x}
#'
#' Nychka, D. (1990). Some properties of adding a smoothing step to the EM
#' algorithm. *Statistics & Probability Letters*, 9(2), 187-193.
#'
#' Fan, J. (1991). On the optimal rates of convergence for nonparametric
#' deconvolution problems. *Annals of Statistics*, 19(3), 1257-1272.
#'
#' @examples
#' set.seed(1)
#' theta <- c(rnorm(300, 0, 0.1), rnorm(100, 0.8, 0.2))
#' se <- runif(400, 0.05, 0.4)
#' y <- rnorm(400, theta, se)
#' fit <- em_smooth(y, se)
#' fit
#' summary(fit, probs = c(0.1, 0.5, 0.9), thresholds = 0.2)
#' plot(fit)
#' @export
em_smooth <- function(estimate, se, smoothing = 0.1) {
  if (!is.numeric(smoothing) || length(smoothing) != 1 ||
      !is.finite(smoothing) || smoothing <= 0) {
    stop("`smoothing` must be a single positive number.", call. = FALSE)
  }
  if (smoothing < 0.05) {
    stop("`smoothing` must be at least 0.05; smaller values smooth too ",
         "little to be reliable.", call. = FALSE)
  }
  check_inputs(estimate, se)
  n <- length(estimate)
  if (n < 3) {
    stop("`em_smooth()` needs at least 3 estimates.", call. = FALSE)
  }
  if (n < 20) {
    warning("Only ", n, " estimates: the estimated distribution is very ",
            "uncertain and mostly reflects the amount of smoothing.",
            call. = FALSE)
  }

  # support grid: 200 points covering the estimates plus a margin
  pad <- 3 * max(se)
  grid <- seq(min(estimate) - pad, max(estimate) + pad, length.out = 200)
  step <- grid[2] - grid[1]

  smoothing_sd <- smoothing * stats::median(se) * sqrt(log(500) / log(n))
  L <- lik_matrix(estimate, se, grid)
  fit <- fit_weights(L, make_smoother(grid, smoothing_sd))
  w <- fit$weights

  structure(
    list(grid = grid, weights = w, density = w / step, loglik = fit$loglik,
         iterations = fit$iterations, converged = fit$converged,
         smoothing = smoothing, smoothing_sd = smoothing_sd,
         data = data.frame(estimate = estimate, se = se),
         call = match.call()),
    class = "em_smooth")
}

#' Posterior mean of each true effect
#'
#' Empirical Bayes posterior means of the true effects, using the estimated
#' distribution from [em_smooth()] as the prior.
#'
#' @param object An object returned by [em_smooth()].
#' @return A numeric vector with one posterior mean per estimate, in the
#'   order of the original data.
#' @examples
#' set.seed(1)
#' se <- runif(200, 0.05, 0.4)
#' y <- rnorm(200, rnorm(200, 0.2, 0.2), se)
#' fit <- em_smooth(y, se)
#' head(posterior_mean(fit))
#' @export
posterior_mean <- function(object) {
  stopifnot(inherits(object, "em_smooth"))
  L <- lik_matrix(object$data$estimate, object$data$se, object$grid)
  num <- sweep(L, 2, object$weights, `*`)
  as.vector(num %*% object$grid / rowSums(num))
}

# ---- internal helpers --------------------------------------------------

check_inputs <- function(estimate, se) {
  if (!is.numeric(estimate) || !is.numeric(se)) {
    stop("`estimate` and `se` must be numeric.", call. = FALSE)
  }
  if (length(estimate) != length(se)) {
    stop("`estimate` and `se` must have the same length.", call. = FALSE)
  }
  if (anyNA(estimate) || anyNA(se)) {
    stop("`estimate` and `se` must not contain missing values.",
         call. = FALSE)
  }
  if (any(!is.finite(estimate)) || any(!is.finite(se)) || any(se <= 0)) {
    stop("`estimate` must be finite and `se` must be positive.",
         call. = FALSE)
  }
  invisible(TRUE)
}

# n x K matrix of N(y_i; theta_k, se_i^2) densities
lik_matrix <- function(y, se, grid) {
  stats::dnorm(outer(y, grid, "-") / se) / se
}

# Gaussian smoothing matrix with standard deviation `sd`; each column sums
# to 1 so that mass is preserved.
make_smoother <- function(grid, sd) {
  S <- stats::dnorm(outer(grid, grid, "-"), 0, sd)
  sweep(S, 2, colSums(S), "/")
}

# Smoothed EM for the mixing weights on the grid
fit_weights <- function(L, smoother, max_iter = 1000, tol = 1e-8) {
  K <- ncol(L)
  w <- rep(1 / K, K)
  ll_old <- -Inf
  converged <- FALSE
  for (it in seq_len(max_iter)) {
    num <- sweep(L, 2, w, `*`)
    den <- pmax(rowSums(num), 1e-300)   # avoid 0/0 when an estimate has no support
    w_new <- as.vector(smoother %*% colMeans(num / den))
    w <- w_new / sum(w_new)
    ll <- sum(log(den))
    if (abs(ll - ll_old) < tol * (abs(ll) + 1)) {
      converged <- TRUE
      break
    }
    ll_old <- ll
  }
  list(weights = w, loglik = ll, iterations = it, converged = converged)
}

# Mean, SD, quantiles and exceedance probabilities of a grid distribution
dist_summary <- function(w, grid, probs, thresholds) {
  m <- sum(w * grid)
  s <- sqrt(sum(w * (grid - m)^2))
  cdf <- cumsum(w)
  q <- stats::approx(x = cdf, y = grid, xout = probs, ties = min,
                     rule = 2)$y
  p_gt <- vapply(thresholds, function(t) sum(w[grid > t]), numeric(1))
  c(mean = m, sd = s,
    stats::setNames(q, paste0("q", probs * 100)),
    stats::setNames(p_gt, paste0("p_gt_", thresholds)))
}
