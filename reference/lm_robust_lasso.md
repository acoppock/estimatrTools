# Additive robust regression on LASSO-selected covariates

Runs
[`lasso_select_covariates`](https://alexandercoppock.com/estimatrTools/reference/lasso_select_covariates.md)
and fits
[`lm_robust`](https://declaredesign.org/r/estimatr/reference/lm_robust.html)
with the selected covariates entered additively, that is
`outcome ~ treatment + X1 + X2`.

## Usage

``` r
lm_robust_lasso(
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

## Details

Additive adjustment assumes a common covariate slope across arms. That
assumption buys degrees of freedom relative to
[`lm_lin_lasso`](https://alexandercoppock.com/estimatrTools/reference/lm_lin_lasso.md),
which fits a slope per arm, and it is what makes this the right choice
when arms are small enough that the interacted model is unstable. When
it is wrong, the treatment coefficient is no longer guaranteed to be
consistent for the average treatment effect, which is the reason Lin
(2013) recommends the interacted form as the default. Choose
deliberately.

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

## References

Lin, W. (2013). Agnostic notes on regression adjustments to experimental
data: reexamining Freedman's critique. *Annals of Applied Statistics*,
7(1), 295-318.
[doi:10.1214/12-AOAS583](https://doi.org/10.1214/12-AOAS583)

## See also

[`lm_lin_lasso`](https://alexandercoppock.com/estimatrTools/reference/lm_lin_lasso.md)
for arm-specific slopes.

Other adjusted estimators:
[`lm_lin_lasso()`](https://alexandercoppock.com/estimatrTools/reference/lm_lin_lasso.md),
[`lm_moderator_lasso()`](https://alexandercoppock.com/estimatrTools/reference/lm_moderator_lasso.md)

## Examples

``` r
set.seed(1)
n <- 400
dat <- data.frame(Z = rep(0:1, n / 2), X1 = rnorm(n), X2 = rnorm(n))
dat$Y <- 0.5 * dat$Z + 1.5 * dat$X1 + rnorm(n)

fit <- lm_robust_lasso(Y ~ Z, ~ X1 + X2, data = dat)
adjustment(fit)
#> [1] "robust"
```
