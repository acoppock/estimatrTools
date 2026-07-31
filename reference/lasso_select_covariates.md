# Select covariates by post-double-selection LASSO

Chooses a covariate set for regression adjustment of an experiment by
running two families of LASSO regressions and taking their union:

## Usage

``` r
lasso_select_covariates(
  formula,
  covariates,
  data,
  lambda_rule = c("lambda.1se", "lambda.min"),
  seed = 999
)
```

## Arguments

- formula:

  A two-sided formula `outcome ~ treatment` giving the focal estimand.
  Only these two variables are read from it.

- covariates:

  A one-sided formula of candidate covariates (e.g. `~ X1 + X2`), a
  character vector of column names, or a bare expression such as
  `X1 + X2`.

- data:

  A data frame.

- lambda_rule:

  Either `"lambda.1se"` (default) or `"lambda.min"`. See the penalty
  section above.

- seed:

  Integer seed for the cross-validation fold assignment, so that
  selection is reproducible. The global random number state is restored
  on exit, so calling this function does not disturb the stream of the
  script that called it.

## Value

A one-sided formula naming the selected original covariates, or `~1`
when nothing was selected or selection could not be run. `~1` is a valid
input to
[`lm_robust`](https://declaredesign.org/r/estimatr/reference/lm_robust.html)
but not to
[`lm_lin`](https://declaredesign.org/r/estimatr/reference/lm_lin.html);
the fitting functions in this package handle that difference for you.

## Details

- \*\*Outcome equations\*\*: the outcome on all candidate covariates,
  separately within each treatment arm. Running these arm-by-arm rather
  than pooled is what makes the selection appropriate for the Lin (2013)
  estimator, whose adjustment is a separate slope per arm: a covariate
  that matters only for the treatment interaction is picked up here and
  would be missed by a pooled regression.

- \*\*Assignment equations\*\*: each arm indicator on all candidate
  covariates. This is the double-selection step of Belloni, Chernozhukov
  and Hansen (2014). It retains covariates that predict assignment even
  when they barely predict the outcome, which is what protects the
  estimate against omitted-variable bias from a covariate the outcome
  equations would have dropped.

Covariates are standardized before selection so the LASSO penalty
applies on a common scale, and selected columns of the model matrix are
mapped back to the original covariate names, so a factor is retained
whole rather than as a subset of its indicator columns.

## Choice of penalty

The default `lambda_rule = "lambda.1se"` is the one-standard-error rule:
the most parsimonious model within one standard error of the
cross-validation minimum. Parsimony matters more here than in plain OLS
because under
[`lm_lin_lasso`](https://alexandercoppock.com/estimatrTools/reference/lm_lin_lasso.md)
each retained covariate costs one parameter \*per arm\*, so
over-selection inflates standard errors faster than it would in an
additive model. Pass `lambda_rule = "lambda.min"` for the
cross-validation minimum instead.

## Rank deficiency in the interacted design

LASSO selects against the outcome and knows nothing about the design
that will consume its answer. Under
[`lm_lin_lasso`](https://alexandercoppock.com/estimatrTools/reference/lm_lin_lasso.md)
every selected covariate is interacted with every treatment arm, so a
set that is fine in a pooled regression can be rank deficient once
crossed with treatment.

When that happens this function **warns and returns the selection
unchanged**. It does not prune. Pruning would decide, silently and
inside the algorithm, something that belongs in the analysis script
where a reader can see it, and the right remedy is not the same in every
case:

- **Redundancy**: one covariate lies in the span of another, for
  instance a college indicator derived from an education factor.
  Dropping one of them at the call site is right, and the warning names
  which.

- **Empty cell**: a single factor level is unobserved in one arm, so one
  interaction column is all zero. Dropping the whole covariate to fix
  this is a poor trade, since a factor contributing fourteen columns
  loses thirteen good ones to remove one aliased one. Collapsing the
  sparse level in the cleaning script keeps the rest.

Left alone, `lm_lin` drops the aliased column itself and still returns
an estimate, so the warning is a prompt to make a decision rather than a
failure.

## Unusable candidates

Candidate columns that are absent from `data`, entirely missing, or
constant are dropped before the complete-case filter rather than after.
Dropping them afterwards is a trap: a single all-`NA` candidate empties
the complete-case set, and selection then silently returns nothing while
appearing to have run.

## References

Belloni, A., Chernozhukov, V., and Hansen, C. (2014). Inference on
treatment effects after selection among high-dimensional controls.
*Review of Economic Studies*, 81(2), 608-650.
[doi:10.1093/restud/rdt044](https://doi.org/10.1093/restud/rdt044)

Lin, W. (2013). Agnostic notes on regression adjustments to experimental
data: reexamining Freedman's critique. *Annals of Applied Statistics*,
7(1), 295-318.
[doi:10.1214/12-AOAS583](https://doi.org/10.1214/12-AOAS583)

## See also

Other covariate selection:
[`lasso_select_one()`](https://alexandercoppock.com/estimatrTools/reference/lasso_select_one.md)

## Examples

``` r
set.seed(1)
n <- 400
dat <- data.frame(
  Z  = rep(0:1, n / 2),
  X1 = rnorm(n),
  X2 = rnorm(n),
  X3 = rnorm(n)
)
dat$Y <- 0.5 * dat$Z + 1.5 * dat$X1 + rnorm(n)

# X1 predicts the outcome; X2 and X3 are noise
lasso_select_covariates(Y ~ Z, ~ X1 + X2 + X3, data = dat)
#> ~X1
#> <environment: 0x90ea70388>
```
