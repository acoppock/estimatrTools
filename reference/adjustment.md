# What adjustment did a fit actually use?

The fitting functions in this package fall back to a simpler
specification when covariate selection returns nothing or the adjusted
fit is degenerate. These accessors report what was actually run, so a
fallback is visible in the results rather than buried inside the
wrapper.

## Usage

``` r
adjustment(fit)

selected_covariates(fit)

fallback_reason(fit)
```

## Arguments

- fit:

  A fit returned by
  [`lm_lin_lasso`](https://alexandercoppock.com/estimatrTools/reference/lm_lin_lasso.md),
  [`lm_robust_lasso`](https://alexandercoppock.com/estimatrTools/reference/lm_robust_lasso.md),
  or
  [`lm_moderator_lasso`](https://alexandercoppock.com/estimatrTools/reference/lm_moderator_lasso.md).

## Value

`adjustment()` returns `"lin"`, `"robust"`, or `"none"`.
`selected_covariates()` returns the character vector of covariates
actually used, possibly empty. `fallback_reason()` returns
`NA_character_` when no fallback fired, and otherwise a description of
why it did.

## See also

Other fallback reporting:
[`fallback_log()`](https://alexandercoppock.com/estimatrTools/reference/fallback_log.md),
[`fallback_summary()`](https://alexandercoppock.com/estimatrTools/reference/fallback_summary.md)

## Examples

``` r
set.seed(1)
n <- 300
dat <- data.frame(Z = rep(0:1, n / 2), X1 = rnorm(n))
dat$Y <- 0.5 * dat$Z + rnorm(n)

# X1 is pure noise, so selection returns nothing and the fit falls back
fit <- lm_lin_lasso(Y ~ Z, ~ X1, data = dat)
adjustment(fit)
#> [1] "lin"
fallback_reason(fit)
#> [1] NA
```
