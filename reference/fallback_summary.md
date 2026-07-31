# Summarize what a collection of fits actually did

[`adjustment`](https://alexandercoppock.com/estimatrTools/reference/adjustment.md)
and
[`fallback_reason`](https://alexandercoppock.com/estimatrTools/reference/adjustment.md)
report on one fit at a time, which is not usable across a pipeline of
hundreds. This walks a list of fits and returns one row each, so a
fallback that quietly changed the estimator for a handful of
specifications is visible without inspecting every fit by hand.

## Usage

``` r
fallback_summary(fits, only_fallbacks = FALSE)
```

## Arguments

- fits:

  A list of fits from
  [`lm_lin_lasso`](https://alexandercoppock.com/estimatrTools/reference/lm_lin_lasso.md),
  [`lm_robust_lasso`](https://alexandercoppock.com/estimatrTools/reference/lm_robust_lasso.md),
  or
  [`lm_moderator_lasso`](https://alexandercoppock.com/estimatrTools/reference/lm_moderator_lasso.md).
  If the list is named, the names are used to label rows. Elements that
  were not produced by this package are reported with `adjustment = NA`.

- only_fallbacks:

  Logical. If `TRUE`, return only the rows where a fallback fired
  (default `FALSE`).

## Value

A data frame with columns `fit` (name or index), `adjustment`,
`n_selected`, `selected`, and `fallback_reason`.

## Details

A fallback is not an error: it means the requested adjustment could not
be produced and a simpler specification was used instead. That is
usually the right thing to do, but it changes what the estimate is, so
it belongs in the record rather than in the wrapper.

## See also

Other fallback reporting:
[`adjustment()`](https://alexandercoppock.com/estimatrTools/reference/adjustment.md),
[`fallback_log()`](https://alexandercoppock.com/estimatrTools/reference/fallback_log.md)

## Examples

``` r
set.seed(1)
n <- 400
dat <- data.frame(Z = rep(0:1, n / 2), X_sig = rnorm(n), X_noise = rnorm(n))
dat$Y <- 0.5 * dat$Z + 2 * dat$X_sig + rnorm(n)

fits <- list(
  informative = lm_lin_lasso(Y ~ Z, ~ X_sig + X_noise, data = dat),
  noise_only  = lm_lin_lasso(Y ~ Z, ~ X_noise, data = dat)
)
fallback_summary(fits)
#>           fit adjustment n_selected selected fallback_reason
#> 1 informative        lin          1    X_sig            <NA>
#> 2  noise_only        lin          1  X_noise            <NA>
fallback_summary(fits, only_fallbacks = TRUE)
#> [1] fit             adjustment      n_selected      selected       
#> [5] fallback_reason
#> <0 rows> (or 0-length row.names)
```
