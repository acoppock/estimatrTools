# Select covariates that predict one variable

Runs a single LASSO regression of one variable on a pool of candidate
covariates and returns the covariates that survive. This is the narrow,
treatment-free counterpart to
[`lasso_select_covariates`](https://alexandercoppock.com/estimatrTools/reference/lasso_select_covariates.md):
it knows nothing about arms, assignment, or estimands, and answers only
"which of these candidates predict this column".

## Usage

``` r
lasso_select_one(
  response,
  covariates,
  data,
  subset = NULL,
  lambda_rule = c("lambda.1se", "lambda.min"),
  seed = 999
)
```

## Arguments

- response:

  Column to be predicted, as a character scalar or a bare column name.

- covariates:

  A one-sided formula of candidate covariates (e.g. `~ X1 + X2`), a
  character vector of column names, or a bare expression such as
  `X1 + X2`.

- data:

  A data frame.

- subset:

  Optional logical vector, the same length as `nrow(data)`, restricting
  the rows selection runs on. `NA` is treated as `FALSE`. Use it for
  equations estimated on part of the sample, such as an outcome equation
  among respondents who answered.

- lambda_rule:

  Either `"lambda.1se"` (default, the one-standard-error rule) or
  `"lambda.min"` (the cross-validation minimum).

- seed:

  Integer seed for the cross-validation fold assignment, so selection is
  reproducible. The global random number state is restored on exit, so
  calling this function does not disturb the stream of the script that
  called it.

## Value

A character vector of selected covariate names, in the order they were
supplied. Empty when nothing was selected or selection could not be run,
which are not distinguished: both mean "no covariates to use".

## Details

It exists because some procedures need selection equations that
`lasso_select_covariates` cannot express. That function takes
`outcome ~ treatment` and selects arm by arm and on each arm indicator,
which is the right thing for regression adjustment of an experiment. A
diagnostic for differential attrition, by contrast, needs two equations
with two different left-hand sides and no treatment in either: the
outcome among respondents who answered, and an indicator for having
dropped out among all those assigned. Calling this function twice
expresses that directly.

## Relation to the other selector

`lasso_select_covariates` returns a one-sided formula, because its
answer is consumed by
[`lm_lin`](https://declaredesign.org/r/estimatr/reference/lm_lin.html).
This function returns a character vector, because its answer is usually
combined with another selection before any model is fit. Wrap the result
in [`reformulate`](https://rdrr.io/r/stats/delete.response.html) if a
formula is what you need.

## Unusable candidates

Candidates that are absent from `data`, entirely missing, or constant
are dropped before the complete-case filter rather than after, and the
check is applied over the subset actually used. Dropping them afterwards
is a trap: a single all-`NA` candidate empties the complete-case set,
and selection then silently returns nothing while appearing to have run.
A covariate that is constant among the rows being used is equally
unusable even when it varies in the full data, which is why the check
follows `subset`.

Non-numeric candidates are expanded through `model.matrix` and mapped
back to their source covariate, so a factor is returned whole rather
than as a subset of its indicator columns.

## See also

[`lasso_select_covariates`](https://alexandercoppock.com/estimatrTools/reference/lasso_select_covariates.md)
for the post-double-selection procedure used to adjust an experiment.

Other covariate selection:
[`lasso_select_covariates()`](https://alexandercoppock.com/estimatrTools/reference/lasso_select_covariates.md)

## Examples

``` r
set.seed(1)
n <- 400
dat <- data.frame(X1 = rnorm(n), X2 = rnorm(n), X3 = rnorm(n))
dat$Y <- 1.5 * dat$X1 + rnorm(n)

# X1 predicts Y; X2 and X3 are noise
lasso_select_one("Y", ~ X1 + X2 + X3, data = dat)
#> [1] "X1"

# Selection restricted to half the sample
lasso_select_one("Y", ~ X1 + X2 + X3, data = dat,
                 subset = seq_len(n) <= n / 2)
#> [1] "X1"
```
