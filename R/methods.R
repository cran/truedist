#' @export
print.em_smooth <- function(x, digits = 3, ...) {
  cat("Smoothed EM estimate of the distribution of true effects\n")
  cat("  estimates:", nrow(x$data), "| grid points:", length(x$grid), "\n")
  cat("  smoothing SD:", signif(x$smoothing_sd, 3),
      "| iterations:", x$iterations,
      if (x$converged) "(converged)" else "(NOT converged)", "\n\n")
  s <- summary(x)
  print(stats::setNames(round(s$estimate, digits), s$statistic))
  if (!x$converged) {
    warning("The EM algorithm did not converge.", call. = FALSE)
  }
  invisible(x)
}

#' Summary of an estimated true-effect distribution
#'
#' @param object An object returned by [em_smooth()].
#' @param probs Probabilities at which to report quantiles.
#' @param thresholds Values `t` at which to report the probability that a
#'   true effect exceeds `t`.
#' @param ... Unused.
#' @return A data frame with one row per statistic (mean, SD, quantiles,
#'   exceedance probabilities) and a column `estimate`.
#' @examples
#' set.seed(1)
#' se <- runif(200, 0.05, 0.4)
#' y <- rnorm(200, rnorm(200, 0.2, 0.2), se)
#' summary(em_smooth(y, se), probs = c(0.1, 0.5, 0.9), thresholds = 0.2)
#' @export
summary.em_smooth <- function(object,
                              probs = c(0.10, 0.25, 0.50, 0.75, 0.90),
                              thresholds = 0, ...) {
  s <- dist_summary(object$weights, object$grid, probs, thresholds)
  data.frame(statistic = names(s), estimate = unname(s),
             row.names = NULL, stringsAsFactors = FALSE)
}

#' Plot an estimated true-effect distribution
#'
#' @param x An object returned by [em_smooth()].
#' @param rug Draw a rug of the observed estimates.
#' @param xlab,ylab,main Axis labels and title.
#' @param ... Further arguments passed to [graphics::plot()].
#' @return `x`, invisibly.
#' @examples
#' set.seed(1)
#' se <- runif(200, 0.05, 0.4)
#' y <- rnorm(200, rnorm(200, 0.2, 0.2), se)
#' plot(em_smooth(y, se))
#' @export
plot.em_smooth <- function(x, rug = TRUE, xlab = "Effect", ylab = "Density",
                           main = "Estimated distribution of true effects",
                           ...) {
  graphics::plot(x$grid, x$density, type = "l", lwd = 2, col = "steelblue",
                 xlab = xlab, ylab = ylab, main = main,
                 ylim = c(0, max(x$density)), ...)
  if (rug) graphics::rug(x$data$estimate, col = "gray50")
  invisible(x)
}
