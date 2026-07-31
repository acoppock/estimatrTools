# Treatment-by-moderator model with LASSO-selected covariates

Fits `outcome ~ treatment * moderator + selected covariates`, where the
covariates are chosen by
[`lasso_select_covariates`](https://alexandercoppock.com/estimatrTools/reference/lasso_select_covariates.md)
from the focal `outcome ~ treatment` equation.

## Usage

``` r
lm_moderator_lasso(
  formula,
  moderator,
  data,
  covariates = NULL,
  clusters = NULL,
  ci = TRUE,
  alpha = 0.05,
  lasso_args = list()
)
```

## Arguments

- formula:

  A two-sided formula `outcome ~ treatment`, the focal estimand equation
  passed to selection. The moderator is supplied separately.

- moderator:

  Character scalar naming the moderator column.

- data:

  A data frame.

- covariates:

  Candidate covariates: a one-sided formula, a character vector of
  column names, or a bare expression.

- clusters, ci, alpha:

  Passed through to
  [`lm_robust`](https://declaredesign.org/r/estimatr/reference/lm_robust.html).

- lasso_args:

  A named list of further arguments for
  [`lasso_select_covariates`](https://alexandercoppock.com/estimatrTools/reference/lasso_select_covariates.md).

## Value

An `lm_robust` object, carrying the attributes described in
[`adjustment`](https://alexandercoppock.com/estimatrTools/reference/adjustment.md).
When the interacted fit with selected covariates cannot be produced, the
fallback is the same interaction model fitted on the full candidate
covariate set, not an unadjusted model.

## Details

The selected covariates enter additively and are deliberately not
interacted with treatment or with the moderator. Adding those
interactions would spend degrees of freedom on terms that are not the
estimand and, with a data-driven covariate set, would make the reported
interaction sensitive to which covariates happened to be selected.

The treatment-by-moderator coefficient describes how the treatment
effect varies across levels of an observed variable. The moderator is
not randomly assigned, so that variation is descriptive: it is not the
causal effect of the moderator, and it does not decompose the treatment
effect into a pathway.

## See also

Other adjusted estimators:
[`lm_lin_lasso()`](https://alexandercoppock.com/estimatrTools/reference/lm_lin_lasso.md),
[`lm_robust_lasso()`](https://alexandercoppock.com/estimatrTools/reference/lm_robust_lasso.md)

## Examples

``` r
set.seed(1)
n <- 600
dat <- data.frame(
  Z = rep(0:1, n / 2),
  X_pid = rep(c(0, 1), each = n / 2),
  X1 = rnorm(n), X2 = rnorm(n)
)
dat$Y <- 0.2 * dat$Z + 0.6 * dat$Z * dat$X_pid + 1.5 * dat$X1 + rnorm(n)

lm_moderator_lasso(Y ~ Z, moderator = "X_pid", data = dat, covariates = ~ X1 + X2)
#>                Estimate Std. Error    t value      Pr(>|t|)    CI Lower
#> (Intercept)  0.03729435 0.08620065  0.4326459  6.654288e-01 -0.13200019
#> Z            0.22121165 0.11711298  1.8888738  5.939482e-02 -0.00879345
#> X_pid       -0.21173121 0.12763984 -1.6588176  9.767956e-02 -0.46241063
#> X1           1.53575730 0.04090927 37.5405720 4.741440e-159  1.45541318
#> Z:X_pid      0.79325275 0.17013819  4.6624026  3.861417e-06  0.45910832
#>              CI Upper  DF
#> (Intercept) 0.2065889 595
#> Z           0.4512167 595
#> X_pid       0.0389482 595
#> X1          1.6161014 595
#> Z:X_pid     1.1273972 595
```
