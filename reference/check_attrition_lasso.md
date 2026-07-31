# Test for differential attrition with LASSO-selected covariates

\*\*Experimental.\*\* Tests whether treatment predicts outcome
missingness, both unconditionally and allowing the pattern to differ
across covariates, using a double-LASSO step to keep the interacted test
from running out of degrees of freedom.

The fully interacted attrition model spends
`(n_arms - 1) x (1 + n_covariates)` degrees of freedom. With the whole
covariate pool that quickly drives events per variable below 10 for any
study with low missingness or many arms, and below that threshold the
model is unreliable in both directions: too sparse to detect real
differential dropout, and capable of spurious significance from
numerical instability. Selecting a parsimonious covariate set first is
what makes the interacted test usable.

Two selection equations, neither involving treatment, are run via
[`lasso_select_one`](https://alexandercoppock.com/estimatrTools/reference/lasso_select_one.md)
and their union enters the interacted model:

- \*\*Outcome equation\*\*: the outcome on all candidates, among
  respondents who answered. Retains the covariates where differential
  dropout would directly bias the treatment effect estimate.

- \*\*Dropout equation\*\*: an indicator for missing outcome on all
  candidates, among everyone assigned. Retains the covariates that
  characterize who drops out.

The covariate-free test is always computed and is the primary criterion
for arm-level differential dropout. The interacted test is secondary: it
detects subgroup-specific dropout that the simple test misses even when
overall rates are equal. Both are returned and the caller decides which
to act on.

## Usage

``` r
check_attrition_lasso(
  data,
  treatment,
  outcomes = NULL,
  covariates = NULL,
  epv_threshold = 10,
  lambda_rule = c("lambda.1se", "lambda.min"),
  alpha = 0.05,
  seed = 999,
  .method = estimatr::lm_robust,
  study_id = NULL,
  quiet = TRUE,
  ...
)
```

## Arguments

- data:

  A data frame or tibble. Rows with a missing treatment value are
  dropped before all calculations, so it may hold unassigned
  respondents.

- treatment:

  Unquoted name of the treatment variable.

- outcomes:

  Character vector of outcome variable names. If \`NULL\`, all columns
  beginning with \`"Y"\` are used.

- covariates:

  Candidate covariate pool: a one-sided formula, a character vector of
  column names, or a bare expression. If \`NULL\` or empty, only the
  covariate-free test is run.

- epv_threshold:

  Minimum events per variable required to run the interacted test.
  Default \`10\`. Below it, \`p_interacted\` is \`NA\` and
  \`epv_adequate\` is \`FALSE\`.

- lambda_rule:

  Either \`"lambda.1se"\` (default, conservative: fewer covariates,
  better events per variable) or \`"lambda.min"\`.

- alpha:

  Threshold at which the \`flag\_\*\` columns are set (default
  \`0.05\`).

- seed:

  Integer seed for the LASSO cross-validation folds, passed to
  \[lasso_select_one()\]. The global random number state is restored on
  exit.

- .method:

  Regression function for the tests (default \`estimatr::lm_robust\`).
  Must accept formula and data arguments.

- study_id:

  Optional character scalar. If provided, a \`study_id\` column holding
  this value is appended to the result, so results stack across studies.

- quiet:

  Logical. The default \`TRUE\` returns the result, which auto-prints at
  the console. \`FALSE\` prints it and returns invisibly.

- ...:

  Additional arguments passed to \`.method\`.

## Value

A tibble with one row per outcome and columns:

- outcome:

  Outcome variable name.

- n_assigned:

  Number of respondents with non-missing treatment.

- n_missing:

  Number of missing outcome observations.

- pct_missing:

  Proportion missing.

- p_simple:

  P-value from the covariate-free omnibus test, `I(missing) ~ Z`.

- n_outcome_eq, n_dropout_eq:

  Covariates selected by each equation.

- n_selected:

  Size of the union.

- selected_covariates:

  Comma-separated names of the selected set.

- df1:

  Numerator degrees of freedom the interacted test would spend,
  `(n_arms - 1) * (1 + n_selected)`.

- epv:

  Events per variable, `n_missing / df1`.

- epv_adequate:

  `epv >= epv_threshold` and `n_selected > 0`.

- p_interacted:

  P-value from the interacted test,
  `I(missing) ~ Z * (demeaned selected covariates)`, testing all
  treatment terms jointly. `NA` when `epv_adequate` is `FALSE`.

- estimable:

  `FALSE` when the missingness indicator does not vary, so there is no
  test to run at all and `p_simple` is `NA`.

- flag_simple, flag_interacted, flag:

  Whether each p-value, and either, falls at or below `alpha`.

## Experimental

This is an original procedure motivated by events-per-variable concerns,
not a reference implementation from a published paper. The double-LASSO
step is adapted from Belloni, Chernozhukov and Hansen (2014) but applied
to an attrition diagnostic rather than a treatment effect estimator, and
the second equation is not theirs (see below). The threshold of 10
events per variable follows the rule of thumb of Peduzzi et al. (1996)
for logistic regression; the same concern applies to a linear
probability model.

## Why selecting and testing on the same data is still valid here

Choosing covariates with the same data used to test them usually
invalidates the reference distribution. It does not here, because
selection never looks at treatment: the outcome equation uses the
outcome and the covariates, the dropout equation uses missingness and
the covariates. Under the null of no differential attrition, missingness
is independent of treatment, so conditioning on a selected set that is a
function of the covariates and missingness leaves treatment randomly
assigned and the test of the treatment terms approximately valid.

That argument depends on treatment being independent of the selection,
so it does not survive a blocked or clustered design analysed as if it
were simple random assignment.

Note also that the outcome equation is not the second equation of
Belloni, Chernozhukov and Hansen, which would be treatment on the
covariates. Since treatment is randomized, that equation selects nothing
asymptotically, which is why it is replaced here. The substitution is
deliberate, and it means the citation motivates the procedure rather
than licensing it.

## p_interacted is sensitive to the cross-validation fold draw

The covariate set is chosen by cross-validated LASSO, and which
covariates survive depends on how the folds happen to be drawn. When the
candidate pool is large relative to the number of events, that
dependence is strong enough to move `p_interacted` across any threshold
you might apply to it. On one real study with eight candidates, the same
outcome gave `n_selected` of 0 and `p_interacted` of `NA` under one
seed, 1 and 0.031 under another, and 1 and 0.423 under a third.

`seed` makes a given run reproducible; it does not make the answer
stable. Two consequences follow. Do not report `p_interacted` as though
it were a fixed property of the study without checking its sensitivity
across seeds. And prefer `p_simple` as the criterion for acting on a
study, since it involves no selection step and is invariant to all of
this.

## Clustered designs

The interacted test is a multi-degree-of-freedom Wald test, and such
tests over-reject badly under cluster randomization even with
cluster-robust standard errors: about 12 percent at a nominal 5 percent
with 30 clusters. The cause is not the denominator degrees of freedom,
so no correction repairs it. `p_simple` is unaffected. Treat
`p_interacted` as descriptive when assignment was clustered.

## When there is no p-value, and what that means

`status` records why, because the reasons mean opposite things.
`"tested"` is a computed p-value. `"no_attrition"` means nobody dropped
out, which is a **pass**: with no attrition there can be no differential
attrition. `"all_missing"` means everybody did, which is uninformative.
`"not_estimable"` means the indicator varied but no statistic came back.
`estimable` equals `status == "tested"`.

The distinction decides a denominator, and the two choices differ by a
lot. To report how much differential attrition a corpus has, count
`"no_attrition"` rows as passes and exclude only the uninformative ones.
To assess whether the computed p-values are uniform, as a valid design
implies, use `"tested"` rows only, because a row with no attrition
produces no draw from Uniform(0, 1) to compare against. On one real
corpus those two rates were 2.5 percent and 7.4 percent, so quoting the
second while describing the first overstates the failure rate roughly
threefold.

Reporting such a row as `p_value = 1` makes the mirror error: it puts a
spike at 1 into the uniformity diagnostic, which then reports a badly
non-uniform collection when nothing is wrong.

## References

Belloni, A., Chernozhukov, V., and Hansen, C. (2014). Inference on
treatment effects after selection among high-dimensional controls.
*Review of Economic Studies*, 81(2), 608-650.
[doi:10.1093/restud/rdt044](https://doi.org/10.1093/restud/rdt044)

Lin, W. (2013). Agnostic notes on regression adjustments to experimental
data: reexamining Freedman's critique. *Annals of Applied Statistics*,
7(1), 295-318.
[doi:10.1214/12-AOAS583](https://doi.org/10.1214/12-AOAS583)

Peduzzi, P., Concato, J., Kemper, E., Holford, T. R., and Feinstein, A.
R. (1996). A simulation study of the number of events per variable in
logistic regression analysis. *Journal of Clinical Epidemiology*,
49(12), 1373-1379.
[doi:10.1016/S0895-4356(96)00236-3](https://doi.org/10.1016/S0895-4356%2896%2900236-3)

## Examples

``` r
set.seed(42)
n <- 500
dat <- data.frame(
  Z = rep(c(0L, 1L), n / 2),
  X_age = rnorm(n, 50, 10),
  X_income = rnorm(n, 50000, 10000)
)
dat$Y_outcome <- 0.3 * dat$Z + 0.5 * scale(dat$X_age) + rnorm(n)
p_miss <- ifelse(dat$Z == 1, 0.05 + 0.01 * scale(dat$X_age), 0.05)
dat$Y_outcome[which(rbinom(n, 1, pmax(0, pmin(1, p_miss))) == 1)] <- NA

check_attrition_lasso(dat, Z,
  outcomes   = "Y_outcome",
  covariates = c("X_age", "X_income"))
#> # A tibble: 1 × 18
#>   outcome   n_assigned n_missing pct_missing p_simple n_outcome_eq n_dropout_eq
#>   <chr>          <int>     <int>       <dbl>    <dbl>        <int>        <int>
#> 1 Y_outcome        500        25        0.05    0.539            1            0
#> # ℹ 11 more variables: n_selected <int>, selected_covariates <chr>, df1 <int>,
#> #   epv <dbl>, epv_adequate <lgl>, p_interacted <dbl>, status <chr>,
#> #   estimable <lgl>, flag_simple <lgl>, flag_interacted <lgl>, flag <lgl>
```
