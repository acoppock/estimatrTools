# Covariate Selection and What Actually Ran

``` r

library(estimatrTools)
```

Covariate adjustment in a randomized experiment is supposed to buy
precision without costing unbiasedness. Choosing the covariates by hand
invites specification search; choosing all of them costs degrees of
freedom, and under the Lin estimator costs them once per arm. This
package selects them by LASSO and then fits the selected set, and it
keeps a record of what it actually did, which turns out to matter more
than it sounds.

## Selection

``` r

set.seed(20260725)
n <- 800

dat <- data.frame(
  Z          = rbinom(n, 1, 0.5),
  X_prog     = rnorm(n),   # predicts the outcome
  X_het      = rnorm(n),   # predicts the outcome only among the treated
  X_noise_1  = rnorm(n),
  X_noise_2  = rnorm(n)
)
dat$Y <- 0.3 * dat$Z + 1.5 * dat$X_prog + 2 * dat$Z * dat$X_het + rnorm(n)

candidates <- ~ X_prog + X_het + X_noise_1 + X_noise_2
lasso_select_covariates(Y ~ Z, candidates, data = dat)
#> ~X_prog + X_het
#> <environment: 0x9501f6ac0>
```

Both real covariates are found and the noise is not. `X_het` is the
interesting one: it predicts the outcome only in the treated arm, so a
pooled regression of `Y` on the candidates would barely notice it.
Selection runs the outcome equation *within each arm*, which is what
makes it appropriate for an estimator that fits a separate slope per
arm.

The second family of regressions is on the arm indicators, the
double-selection step of Belloni, Chernozhukov and Hansen (2014). It
retains covariates that predict assignment even when they barely predict
the outcome, which is the case where omitting one does the most damage.

## Three estimators

``` r

lin <- lm_lin_lasso(Y ~ Z, candidates, data = dat)
add <- lm_robust_lasso(Y ~ Z, candidates, data = dat)

broom::tidy(lin)[1:2, c("term", "estimate", "std.error")]
#>          term    estimate  std.error
#> 1 (Intercept) -0.09899319 0.05124958
#> 2           Z  0.38111507 0.07151957
broom::tidy(add)[1:2, c("term", "estimate", "std.error")]
#>          term    estimate  std.error
#> 1 (Intercept) -0.04596937 0.07399295
#> 2           Z  0.38060321 0.10214710
```

[`lm_lin_lasso()`](https://alexandercoppock.com/estimatrTools/reference/lm_lin_lasso.md)
fits a separate covariate slope in each arm, with covariates centered so
the treatment coefficient still estimates the average treatment effect.
[`lm_robust_lasso()`](https://alexandercoppock.com/estimatrTools/reference/lm_robust_lasso.md)
enters them additively, assuming a common slope. That assumption buys
degrees of freedom, which is what makes it the better choice when arms
are small, but when it is wrong the treatment coefficient is no longer
guaranteed to be consistent for the ATE. Lin (2013) recommends the
interacted form as the default for that reason. The choice is yours to
make, not the package’s, which is why they are separate functions rather
than an argument.

[`lm_moderator_lasso()`](https://alexandercoppock.com/estimatrTools/reference/lm_moderator_lasso.md)
fits `Y ~ Z * moderator` with the selected covariates entered additively
alongside, deliberately not interacted with treatment or moderator.

## What actually ran

Adjustment is not always possible. `lm_lin` cannot be fitted with no
covariates, and an adjusted fit can come back degenerate when the
selected set is nearly collinear within an arm. In either case these
functions return an unadjusted fit, which for a binary treatment is the
difference in means.

That is usually the right response. It also changes what the estimate
is, so it is recorded.

``` r

noise_only <- data.frame(Z = rbinom(n, 1, 0.5), X_noise = rnorm(n))
noise_only$Y <- 0.3 * noise_only$Z + rnorm(n)

fit <- lm_lin_lasso(Y ~ Z, ~ X_noise, data = noise_only)

adjustment(fit)
#> [1] "lin"
fallback_reason(fit)
#> [1] NA
```

**This is common, not exceptional.** In one meta-analysis pipeline, 38%
of roughly a thousand
[`lm_lin_lasso()`](https://alexandercoppock.com/estimatrTools/reference/lm_lin_lasso.md)
calls fell back, nearly all because LASSO selected nothing. In those
cases an estimate described in the write-up as “Lin-adjusted” is an
unadjusted difference in means. That is not wrong, and it is what the
estimator is defined to do, but a reader comparing adjusted and
unadjusted estimates should know that a large share of the pairs are
identical by construction rather than by agreement.

Three accessors, depending on what you still hold. For one fit:

``` r

c(adjustment = adjustment(fit), reason = fallback_reason(fit))
#> adjustment     reason 
#>      "lin"         NA
```

For a list of fits:

``` r

fits <- list(
  informative = lm_lin_lasso(Y ~ Z, candidates, data = dat),
  noise_only  = lm_lin_lasso(Y ~ Z, ~ X_noise, data = noise_only)
)

fallback_summary(fits)
#>           fit adjustment n_selected      selected fallback_reason
#> 1 informative        lin          2 X_prog, X_het            <NA>
#> 2  noise_only        lin          1       X_noise            <NA>
```

And for fits that no longer exist. This is the case that matters in
practice, because the common shape is to build a fit inside `map()`,
summarize it, and drop the model column:

``` r

reset_fallback_log()

for (i in 1:5) {
  d <- if (i %% 2 == 0) noise_only else dat
  cand <- if (i %% 2 == 0) ~ X_noise else candidates
  invisible(lm_lin_lasso(Y ~ Z, cand, data = d))   # fit is discarded
}

fallback_log()
#>   call_index           fn outcome treatment adjustment n_selected
#> 1          1 lm_lin_lasso       Y         Z        lin          2
#> 2          2 lm_lin_lasso       Y         Z        lin          1
#> 3          3 lm_lin_lasso       Y         Z        lin          2
#> 4          4 lm_lin_lasso       Y         Z        lin          1
#> 5          5 lm_lin_lasso       Y         Z        lin          2
#>   fallback_reason
#> 1            <NA>
#> 2            <NA>
#> 3            <NA>
#> 4            <NA>
#> 5            <NA>
```

The fits are gone; the record is not. The question a pipeline should be
able to answer is what share of its estimates were actually adjusted:

``` r

table(fallback_log()$adjustment)
#> 
#> lin 
#>   5
```

Call
[`reset_fallback_log()`](https://alexandercoppock.com/estimatrTools/reference/fallback_log.md)
at the top of a pipeline so the log describes that run. Turn it off with
`options(estimatrTools.log = FALSE)`.

## References

Belloni, A., Chernozhukov, V., and Hansen, C. (2014). Inference on
treatment effects after selection among high-dimensional controls.
*Review of Economic Studies*, 81(2), 608-650.

Lin, W. (2013). Agnostic notes on regression adjustments to experimental
data: reexamining Freedman’s critique. *Annals of Applied Statistics*,
7(1), 295-318.
