# Lin estimator on LASSO-selected covariates

Runs
[`lasso_select_covariates`](https://alexandercoppock.com/estimatrTools/reference/lasso_select_covariates.md)
and fits
[`lm_lin`](https://declaredesign.org/r/estimatr/reference/lm_lin.html)
on the selected set: a separate covariate slope in each treatment arm,
with covariates centered at their full-sample means so the treatment
coefficient remains an estimate of the average treatment effect.

## Usage

``` r
lm_lin_lasso(
  formula,
  covariates,
  data,
  weights = NULL,
  subset = NULL,
  clusters = NULL,
  se_type = NULL,
  ci = TRUE,
  alpha = 0.05,
  return_vcov = TRUE,
  try_cholesky = FALSE,
  lasso_args = list()
)
```

## Arguments

- formula:

  A two-sided formula `outcome ~ treatment`.

- covariates:

  Candidate covariates: a one-sided formula, a character vector of
  column names, or a bare expression.

- data:

  A data frame.

- weights, subset, clusters, se_type, ci, alpha, return_vcov,
  try_cholesky:

  Passed through to
  [`lm_lin`](https://declaredesign.org/r/estimatr/reference/lm_lin.html)
  or
  [`lm_robust`](https://declaredesign.org/r/estimatr/reference/lm_robust.html).

- lasso_args:

  A named list of further arguments for
  [`lasso_select_covariates`](https://alexandercoppock.com/estimatrTools/reference/lasso_select_covariates.md),
  such as `lambda_rule` or `seed`.

## Value

An `lm_robust` object, carrying the attributes described in
[`adjustment`](https://alexandercoppock.com/estimatrTools/reference/adjustment.md).

## Fallback

`lm_lin` cannot be fitted with no covariates, and an adjusted fit can
come back degenerate when the selected set is nearly collinear within an
arm. In either case this function returns an unadjusted
[`lm_robust`](https://declaredesign.org/r/estimatr/reference/lm_robust.html)
fit, which for a binary treatment is the difference in means. Because
that substitution changes the specification, it is recorded rather than
hidden: read it back with
[`adjustment`](https://alexandercoppock.com/estimatrTools/reference/adjustment.md)
and
[`selected_covariates`](https://alexandercoppock.com/estimatrTools/reference/adjustment.md).
The specific triggers are that selection returned `~1`, that `lm_lin`
threw, or that the treatment row of the adjusted fit has a non-finite
standard error.

## See also

[`lm_robust_lasso`](https://alexandercoppock.com/estimatrTools/reference/lm_robust_lasso.md)
for additive adjustment,
[`lm_moderator_lasso`](https://alexandercoppock.com/estimatrTools/reference/lm_moderator_lasso.md)
for a treatment-by-moderator model.

Other adjusted estimators:
[`lm_moderator_lasso()`](https://alexandercoppock.com/estimatrTools/reference/lm_moderator_lasso.md),
[`lm_robust_lasso()`](https://alexandercoppock.com/estimatrTools/reference/lm_robust_lasso.md)

## Examples

``` r
set.seed(1)
n <- 400
dat <- data.frame(Z = rep(0:1, n / 2), X1 = rnorm(n), X2 = rnorm(n))
dat$Y <- 0.5 * dat$Z + 1.5 * dat$X1 + rnorm(n)

fit <- lm_lin_lasso(Y ~ Z, ~ X1 + X2, data = dat)
adjustment(fit)
#> [1] "lin"
selected_covariates(fit)
#> [1] "X1"
```
