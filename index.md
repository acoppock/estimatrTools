# estimatrTools

Companion tools for [estimatr](https://declaredesign.org/r/estimatr/):
data-driven covariate selection for randomized experiments, and the
regression-adjusted estimators that consume it.

## Installation

``` r

# install.packages("remotes")
remotes::install_github("acoppock/estimatrTools")
```

## Selecting covariates

[`lasso_select_covariates()`](https://alexandercoppock.com/estimatrTools/reference/lasso_select_covariates.md)
runs two families of LASSO regressions and takes their union:

- **Outcome equations**, fit *within each treatment arm*. Running these
  arm-by-arm rather than pooled is what makes the selection appropriate
  for an estimator with a separate slope per arm: a covariate that
  matters only for the treatment interaction is picked up here and would
  be missed by a pooled fit.
- **Assignment equations**, one per arm indicator. This is the
  double-selection step of Belloni, Chernozhukov and Hansen (2014),
  which retains covariates that predict assignment even when they barely
  predict the outcome.

``` r

library(estimatrTools)

lasso_select_covariates(Y ~ Z, ~ X_age + X_educ + X_income, data = dat)
#> ~X_age
```

Selected model-matrix columns are mapped back to source covariates by
longest matching name, so a factor is kept whole rather than as a subset
of its indicator columns, and one covariate whose name is a prefix of
another does not drag the other in.

## Fitting

``` r

lm_lin_lasso(Y ~ Z, ~ X_age + X_educ, data = dat)     # Lin: slope per arm
lm_robust_lasso(Y ~ Z, ~ X_age + X_educ, data = dat)  # additive
lm_moderator_lasso(Y ~ Z, moderator = "X_pid", data = dat, covariates = ~ X_age)
```

[`lm_lin_lasso()`](https://alexandercoppock.com/estimatrTools/reference/lm_lin_lasso.md)
is the default choice.
[`lm_robust_lasso()`](https://alexandercoppock.com/estimatrTools/reference/lm_robust_lasso.md)
assumes a common covariate slope across arms, which buys degrees of
freedom but is not guaranteed consistent for the average treatment
effect when that assumption is wrong; Lin (2013) recommends the
interacted form for exactly that reason. Pick deliberately.

## Knowing what actually ran

Each fitting function falls back to a simpler specification when
selection returns nothing or the adjusted fit comes back degenerate.
That is usually the right response, but it changes what the estimate is,
so it is recorded rather than hidden.

**Selection returning nothing is common, not exceptional.** In one
meta-analysis pipeline, 38% of roughly a thousand
[`lm_lin_lasso()`](https://alexandercoppock.com/estimatrTools/reference/lm_lin_lasso.md)
calls fell back, nearly all because LASSO selected no covariates. In
those cases the “Lin-adjusted” estimate is an unadjusted difference in
means. A pipeline that never checks has no way of knowing what fraction
of its estimates that describes.

How to read it back depends on what you still have:

``` r

# one fit in hand
fit <- lm_lin_lasso(Y ~ Z, ~ X_age, data = dat)
adjustment(fit)           #> "lin", "robust", or "none"
selected_covariates(fit)
fallback_reason(fit)      #> NA when none fired

# a list of fits
fallback_summary(fits, only_fallbacks = TRUE)

# fits that were built and thrown away, the usual case inside map()
reset_fallback_log()
# ... run the pipeline ...
fallback_log()
table(fallback_log()$adjustment)
```

[`fallback_log()`](https://alexandercoppock.com/estimatrTools/reference/fallback_log.md)
exists because the common shape is `map(data, ~ lm_lin_lasso(...))`
piped into a summarizer and then `select(-model)`. By the time anyone
asks which fits were adjusted, the objects are gone. The log survives
them. Disable with `options(estimatrTools.log = FALSE)`.

## Vignette

[`vignette("covariate_selection")`](https://alexandercoppock.com/estimatrTools/articles/covariate_selection.md)
covers selection, the difference between the three estimators, and
reading the fallback record.

## Related packages

[excheckr](https://github.com/acoppock/excheckr) for experimental design
diagnostics, and [metaprep](https://github.com/acoppock/metaprep) for
carrying estimates and their covariance through a meta-analysis.

## References

Belloni, A., Chernozhukov, V., and Hansen, C. (2014). Inference on
treatment effects after selection among high-dimensional controls.
*Review of Economic Studies*, 81(2), 608-650.

Lin, W. (2013). Agnostic notes on regression adjustments to experimental
data: reexamining Freedman’s critique. *Annals of Applied Statistics*,
7(1), 295-318.

## License

MIT
