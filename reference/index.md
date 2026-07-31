# Package index

## Selecting covariates

Post-double-selection LASSO: the outcome within each treatment arm, and
each arm indicator on the candidates.

- [`lasso_select_covariates()`](https://alexandercoppock.com/estimatrTools/reference/lasso_select_covariates.md)
  : Select covariates by post-double-selection LASSO
- [`lasso_select_one()`](https://alexandercoppock.com/estimatrTools/reference/lasso_select_one.md)
  : Select covariates that predict one variable

## Adjusted estimators

Fit the selected set. Each falls back to a simpler specification when
selection returns nothing or the adjusted fit is degenerate.

- [`lm_lin_lasso()`](https://alexandercoppock.com/estimatrTools/reference/lm_lin_lasso.md)
  : Lin estimator on LASSO-selected covariates
- [`lm_robust_lasso()`](https://alexandercoppock.com/estimatrTools/reference/lm_robust_lasso.md)
  : Additive robust regression on LASSO-selected covariates
- [`lm_moderator_lasso()`](https://alexandercoppock.com/estimatrTools/reference/lm_moderator_lasso.md)
  : Treatment-by-moderator model with LASSO-selected covariates

## Attrition

Whether treatment predicts who is missing, allowing the pattern to
differ across covariates.

- [`check_attrition_lasso()`](https://alexandercoppock.com/estimatrTools/reference/check_attrition_lasso.md)
  : Test for differential attrition with LASSO-selected covariates

## What actually ran

A fallback changes what the estimate is, so it is recorded rather than
hidden. Which accessor you need depends on whether you still hold the
fit.

- [`adjustment()`](https://alexandercoppock.com/estimatrTools/reference/adjustment.md)
  [`selected_covariates()`](https://alexandercoppock.com/estimatrTools/reference/adjustment.md)
  [`fallback_reason()`](https://alexandercoppock.com/estimatrTools/reference/adjustment.md)
  : What adjustment did a fit actually use?
- [`fallback_summary()`](https://alexandercoppock.com/estimatrTools/reference/fallback_summary.md)
  : Summarize what a collection of fits actually did
- [`fallback_log()`](https://alexandercoppock.com/estimatrTools/reference/fallback_log.md)
  [`reset_fallback_log()`](https://alexandercoppock.com/estimatrTools/reference/fallback_log.md)
  : Record of what every fit did

## Package

- [`estimatrTools`](https://alexandercoppock.com/estimatrTools/reference/estimatrTools-package.md)
  [`estimatrTools-package`](https://alexandercoppock.com/estimatrTools/reference/estimatrTools-package.md)
  : estimatrTools: Data-Driven Covariate Adjustment for the 'estimatr'
  Estimators
