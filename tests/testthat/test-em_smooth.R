simulate_data <- function(n = 600, seed = 42) {
  set.seed(seed)
  theta <- c(rnorm(n * 0.6, 0.05, 0.05), rnorm(n * 0.4, 0.8, 0.15))
  se <- runif(n, 0.05, 0.35)
  list(theta = theta, se = se, y = rnorm(n, theta, se))
}

stat <- function(fit, name, ...) {
  s <- summary(fit, ...)
  s$estimate[s$statistic == name]
}

test_that("weights are a probability distribution on the grid", {
  d <- simulate_data()
  fit <- em_smooth(d$y, d$se)
  expect_s3_class(fit, "em_smooth")
  expect_equal(sum(fit$weights), 1)
  expect_true(all(fit$weights >= 0))
  expect_equal(length(fit$grid), length(fit$weights))
  expect_equal(sum(fit$density) * (fit$grid[2] - fit$grid[1]), 1)
})

test_that("recovers the mean and median of a bimodal true distribution", {
  d <- simulate_data()
  fit <- em_smooth(d$y, d$se)
  expect_lt(abs(stat(fit, "mean") - mean(d$theta)), 0.05)
  expect_lt(abs(stat(fit, "q50") - median(d$theta)), 0.1)
})

test_that("deconvolved distribution is less dispersed than the estimates", {
  d <- simulate_data()
  fit <- em_smooth(d$y, d$se)
  expect_lt(stat(fit, "sd"), sd(d$y))
})

test_that("the smoothing SD follows the rule in the documentation", {
  d <- simulate_data(150)
  fit <- em_smooth(d$y, d$se)
  expect_equal(fit$smoothing, 0.1)
  expect_equal(fit$smoothing_sd,
               0.1 * median(d$se) * sqrt(log(500) / log(length(d$se))))
  fit2 <- em_smooth(d$y, d$se, smoothing = 0.3)
  expect_equal(fit2$smoothing_sd, 3 * fit$smoothing_sd)
  expect_output(print(fit), "smoothing SD")
})

test_that("more smoothing gives a less peaked distribution", {
  d <- simulate_data()
  a <- em_smooth(d$y, d$se, smoothing = 0.05)
  b <- em_smooth(d$y, d$se, smoothing = 0.5)
  expect_gt(max(a$density), max(b$density))
})

test_that("tiny samples are rejected and small ones warn", {
  set.seed(6)
  expect_error(em_smooth(c(0.1, 0.2), c(0.1, 0.1)), "at least 3")
  se <- rep(0.1, 15)
  expect_warning(em_smooth(rnorm(15, 0.2, 0.2), se), "15 estimates")
})

test_that("smoothing is validated", {
  d <- simulate_data(100)
  expect_error(em_smooth(d$y, d$se, smoothing = 0.01), "at least 0.05")
  expect_error(em_smooth(d$y, d$se, smoothing = c(0.1, 0.2)), "single")
  expect_error(em_smooth(d$y, d$se, smoothing = -1), "positive")
})

test_that("summary accepts custom probs and thresholds", {
  d <- simulate_data(200)
  fit <- em_smooth(d$y, d$se)
  s <- summary(fit, probs = c(0.2, 0.8), thresholds = c(0, 0.5))
  expect_true(all(c("mean", "sd", "q20", "q80", "p_gt_0", "p_gt_0.5") %in%
                    s$statistic))
})

test_that("posterior means have one value per estimate and shrink", {
  d <- simulate_data(200)
  fit <- em_smooth(d$y, d$se)
  pm <- posterior_mean(fit)
  expect_length(pm, length(d$y))
  expect_lt(sd(pm), sd(d$y))
})

test_that("normal_fit recovers mu and tau and its summary is exact", {
  set.seed(3)
  n <- 2000
  se <- runif(n, 0.05, 0.4)
  y <- rnorm(n, rnorm(n, 0.3, 0.25), se)
  fit <- normal_fit(y, se)
  expect_lt(abs(fit$mu - 0.3), 0.03)
  expect_lt(abs(fit$tau - 0.25), 0.03)
  expect_equal(stat(fit, "q50"), fit$mu)
  expect_equal(stat(fit, "p_gt_0.3", thresholds = 0.3), 0.5, tolerance = 0.05)
  expect_equal(density(fit, at = fit$mu), dnorm(0, 0, fit$tau))
  expect_output(print(fit), "Normal model")
})

test_that("bad inputs are rejected with informative errors", {
  expect_error(em_smooth(c(1, 2), c(1, 2, 3)), "same length")
  expect_error(em_smooth(c(1, 2), c(1, -1)), "positive")
  expect_error(em_smooth(c(1, NA), c(1, 1)), "missing")
})

test_that("print and plot run without error", {
  d <- simulate_data(100)
  fit <- em_smooth(d$y, d$se)
  expect_output(print(fit), "Smoothed EM")
  pdf(NULL)
  on.exit(dev.off())
  expect_no_error(plot(fit))
})
