# Non-standard evaluation of weights, clusters, and subset.
#
# estimatr resolves these against `data` before the calling frame, which makes
# two failure modes possible and they pull in opposite directions. Forwarding an
# unsupplied argument as a symbol lets a same-named data column capture it; this
# is the bug that silently weighted every Lin fit in meta_propaganda on any study
# carrying a `weights` column, while the paired difference-in-means fits went
# unweighted. Evaluating the argument here instead breaks the ordinary usage of
# naming a column. Both directions are tested.

nse_data <- function(n = 400, seed = 1) {
  set.seed(seed)
  d <- data.frame(Z = rep(0:1, n / 2), X_sig = rnorm(n), X_noise = rnorm(n))
  d$Y <- 0.5 * d$Z + 2 * d$X_sig + rnorm(n)
  d
}

test_that("a data column named 'weights' does not silently weight the fit", {
  plain <- nse_data()
  trap <- plain
  trap$weights <- runif(nrow(trap), 0.5, 1.5)

  a <- lm_lin_lasso(Y ~ Z, ~ X_sig + X_noise, data = trap)
  b <- lm_lin_lasso(Y ~ Z, ~ X_sig + X_noise, data = plain)

  expect_false(isTRUE(a$weighted))
  expect_equal(broom::tidy(a)$estimate, broom::tidy(b)$estimate)
  expect_equal(broom::tidy(a)$std.error, broom::tidy(b)$std.error)
})

test_that("the same trap is avoided by lm_robust_lasso and lm_moderator_lasso", {
  plain <- nse_data()
  plain$X_mod <- rep(0:1, each = nrow(plain) / 2)
  trap <- plain
  trap$weights <- runif(nrow(trap), 0.5, 1.5)

  expect_false(isTRUE(lm_robust_lasso(Y ~ Z, ~ X_sig, data = trap)$weighted))
  expect_false(isTRUE(
    lm_moderator_lasso(Y ~ Z, moderator = "X_mod", data = trap, covariates = ~ X_sig)$weighted
  ))
})

test_that("weights supplied as a bare column name still resolve against data", {
  d <- nse_data()
  d$w <- runif(nrow(d), 0.5, 1.5)

  weighted <- lm_lin_lasso(Y ~ Z, ~ X_sig + X_noise, data = d, weights = w)
  unweighted <- lm_lin_lasso(Y ~ Z, ~ X_sig + X_noise, data = d)

  expect_true(isTRUE(weighted$weighted))
  expect_false(isTRUE(all.equal(
    broom::tidy(weighted)$estimate, broom::tidy(unweighted)$estimate
  )))
})

test_that("a data column named 'clusters' does not silently cluster the fit", {
  plain <- nse_data()
  trap <- plain
  trap$clusters <- rep(1:20, length.out = nrow(trap))

  a <- lm_lin_lasso(Y ~ Z, ~ X_sig + X_noise, data = trap)
  b <- lm_lin_lasso(Y ~ Z, ~ X_sig + X_noise, data = plain)

  expect_null(a$clusters)
  expect_equal(broom::tidy(a)$std.error, broom::tidy(b)$std.error)
})

test_that("subset resolves against data and actually restricts the sample", {
  d <- nse_data()
  d$keep <- d$X_noise > 0

  restricted <- lm_lin_lasso(Y ~ Z, ~ X_sig, data = d, subset = keep)
  full <- lm_lin_lasso(Y ~ Z, ~ X_sig, data = d)

  expect_lt(restricted$nobs, full$nobs)
  expect_equal(restricted$nobs, sum(d$keep))
})

test_that("selection ignores weights, by design", {
  # Selection is run unweighted even when the fit is weighted. Pinned here
  # because it is a modelling decision rather than an oversight: the LASSO
  # chooses which covariates are prognostic, which is a question about the
  # sample, not about the population the weights target.
  d <- nse_data()
  d$w <- runif(nrow(d), 0.5, 1.5)

  a <- lasso_select_covariates(Y ~ Z, ~ X_sig + X_noise, data = d)
  fit_w <- lm_lin_lasso(Y ~ Z, ~ X_sig + X_noise, data = d, weights = w)

  expect_equal(sort(all.vars(a)), sort(selected_covariates(fit_w)))
})
