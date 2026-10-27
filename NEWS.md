# truedist 0.1.0

* First CRAN release.

* `em_smooth(estimate, se, smoothing = 0.1)`: smoothed EM estimate of the distribution of true effects. The smoothing uses one constant SD, equal to `smoothing` times the median standard error times a factor that grows slowly as the sample size falls.
* `summary()` computes the mean, SD, quantiles and exceedance probabilities on demand (`probs`, `thresholds`).
* `posterior_mean()`, and `print`, `summary` and `plot` methods.
* `normal_fit()` fits the normal model for the true effects.
