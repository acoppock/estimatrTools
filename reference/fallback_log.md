# Record of what every fit did

[`adjustment`](https://alexandercoppock.com/estimatrTools/reference/adjustment.md)
and friends report on a fit you still hold.
[`fallback_summary`](https://alexandercoppock.com/estimatrTools/reference/fallback_summary.md)
reports on a list of fits you still hold. This reports on every fit made
in the session, including ones that were discarded, which is the usual
case in a pipeline that builds a fit inside `map()`, summarizes it, and
drops the model column.

## Usage

``` r
fallback_log(fallbacks_only = FALSE)

reset_fallback_log()
```

## Arguments

- fallbacks_only:

  Logical. Return only the calls where a fallback fired (default
  `FALSE`, which returns one row per call).

## Value

A data frame with one row per call: `call_index`, `fn`, `outcome`,
`treatment`, `adjustment`, `n_selected`, and `fallback_reason` (`NA`
when none fired).

`reset_fallback_log()` returns nothing and is called for its effect.

## Details

Falling back is not an error: it means the requested adjustment could
not be produced and a simpler specification was used. It does change
what the estimate is, though, so a pipeline that never looks at this log
has no way of knowing how many of its "adjusted" estimates are
unadjusted.

## Turning it off

Recording is on by default and costs a few scalars per call. Disable
with `options(estimatrTools.log = FALSE)`. `reset_fallback_log()` clears
what has accumulated, which is worth doing at the top of a pipeline so
the log describes that run rather than everything since the session
started.

## See also

[`adjustment`](https://alexandercoppock.com/estimatrTools/reference/adjustment.md)
for a single fit,
[`fallback_summary`](https://alexandercoppock.com/estimatrTools/reference/fallback_summary.md)
for a list of fits you have kept.

Other fallback reporting:
[`adjustment()`](https://alexandercoppock.com/estimatrTools/reference/adjustment.md),
[`fallback_summary()`](https://alexandercoppock.com/estimatrTools/reference/fallback_summary.md)

## Examples

``` r
reset_fallback_log()

set.seed(1)
n <- 300
dat <- data.frame(Z = rep(0:1, n / 2), X_sig = rnorm(n), X_noise = rnorm(n))
dat$Y <- 0.5 * dat$Z + 2 * dat$X_sig + rnorm(n)

invisible(lm_lin_lasso(Y ~ Z, ~ X_sig + X_noise, data = dat))
invisible(lm_lin_lasso(Y ~ Z, ~ X_noise, data = dat))

fallback_log()
#>   call_index           fn outcome treatment adjustment n_selected
#> 1          1 lm_lin_lasso       Y         Z        lin          1
#> 2          2 lm_lin_lasso       Y         Z        lin          1
#>   fallback_reason
#> 1            <NA>
#> 2            <NA>

# what fraction of fits were actually adjusted?
table(fallback_log()$adjustment)
#> 
#> lin 
#>   2 
```
