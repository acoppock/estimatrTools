# estimatrTools: Data-Driven Covariate Adjustment for the 'estimatr' Estimators

Companion tools for 'estimatr'. Selects covariates by LASSO within each
treatment arm and on each treatment-arm indicator (post-double
selection, following Belloni, Chernozhukov and Hansen (2014)
[doi:10.1093/restud/rdt044](https://doi.org/10.1093/restud/rdt044) ),
then fits the selected set with the Lin (2013) interacted estimator
[doi:10.1214/12-AOAS583](https://doi.org/10.1214/12-AOAS583) , an
additive robust regression, or a treatment-by-moderator interaction
model. Each fitting function falls back to an unadjusted specification
when selection returns nothing or the adjusted fit is degenerate, so a
fit is always returned and the fallback is reported rather than hidden.

## Selecting covariates

[`lasso_select_covariates`](https://alexandercoppock.com/estimatrTools/reference/lasso_select_covariates.md)
chooses a covariate set by running LASSO of the outcome on candidates
*within each treatment arm*, and of each arm indicator on candidates.
The first is what makes the selection appropriate for an estimator with
a separate slope per arm; the second is the double-selection step that
retains covariates predicting assignment even when they barely predict
the outcome.

## Fitting

[`lm_lin_lasso`](https://alexandercoppock.com/estimatrTools/reference/lm_lin_lasso.md)
fits the Lin (2013) interacted estimator on the selected set,
[`lm_robust_lasso`](https://alexandercoppock.com/estimatrTools/reference/lm_robust_lasso.md)
fits them additively, and
[`lm_moderator_lasso`](https://alexandercoppock.com/estimatrTools/reference/lm_moderator_lasso.md)
fits a treatment-by-moderator model with the selected covariates entered
additively alongside.

## Knowing what actually ran

Each fitting function falls back to a simpler specification when
selection returns nothing or the adjusted fit is degenerate. This is
usually the right response, but it changes what the estimate is, so it
is recorded rather than hidden. How to read it back depends on what you
still have:

- one fit in hand:
  [`adjustment`](https://alexandercoppock.com/estimatrTools/reference/adjustment.md),
  [`selected_covariates`](https://alexandercoppock.com/estimatrTools/reference/adjustment.md),
  [`fallback_reason`](https://alexandercoppock.com/estimatrTools/reference/adjustment.md)

- a list of fits:
  [`fallback_summary`](https://alexandercoppock.com/estimatrTools/reference/fallback_summary.md)

- fits that were built and discarded, which is the usual case inside a
  `map()` pipeline:
  [`fallback_log`](https://alexandercoppock.com/estimatrTools/reference/fallback_log.md)

The last matters more than it sounds. Selection returning nothing is
common, not exceptional, and when it happens the "adjusted" estimate is
an unadjusted difference in means. A pipeline that never reads the log
has no way of knowing what fraction of its estimates that describes.

## See also

Useful links:

- <https://alexandercoppock.com/estimatrTools/>

- <https://github.com/acoppock/estimatrTools>

- Report bugs at <https://github.com/acoppock/estimatrTools/issues>

## Author

**Maintainer**: Alexander Coppock <acoppock@gmail.com>

Authors:

- Alexander Coppock <acoppock@gmail.com>
