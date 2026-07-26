fit_data <- function(n = 600, seed = 1) {
  set.seed(seed)
  dat <- data.frame(
    Z     = rep(0:1, n / 2),
    X_sig = rnorm(n),
    X_n1  = rnorm(n),
    X_n2  = rnorm(n)
  )
  dat$Y <- 0.5 * dat$Z + 2 * dat$X_sig + rnorm(n)
  dat
}

noise_data <- function(n = 300, seed = 11) {
  set.seed(seed)
  dat <- data.frame(Z = rep(0:1, n / 2), X1 = rnorm(n), X2 = rnorm(n))
  dat$Y <- 0.5 * dat$Z + rnorm(n)
  dat
}

test_that("lm_lin_lasso recovers the treatment effect and reports its adjustment", {
  fit <- lm_lin_lasso(Y ~ Z, ~ X_sig + X_n1 + X_n2, data = fit_data())

  expect_s3_class(fit, "lm_robust")
  expect_equal(adjustment(fit), "lin")
  expect_true("X_sig" %in% selected_covariates(fit))
  expect_true(is.na(fallback_reason(fit)))

  z_row <- broom::tidy(fit)[broom::tidy(fit)$term == "Z", ]
  expect_lt(z_row$conf.low, 0.5)
  expect_gt(z_row$conf.high, 0.5)
})

test_that("adjustment shrinks the standard error when a covariate predicts Y", {
  dat <- fit_data()
  adjusted <- lm_lin_lasso(Y ~ Z, ~ X_sig + X_n1 + X_n2, data = dat)
  unadjusted <- estimatr::lm_robust(Y ~ Z, data = dat)

  se_adj <- broom::tidy(adjusted)$std.error[2]
  se_un  <- broom::tidy(unadjusted)$std.error[2]
  expect_lt(se_adj, se_un)
})

test_that("lm_lin_lasso falls back and says why when nothing is selected", {
  fit <- lm_lin_lasso(Y ~ Z, ~ X1 + X2, data = noise_data())

  expect_equal(adjustment(fit), "none")
  expect_equal(selected_covariates(fit), character(0))
  expect_equal(fallback_reason(fit), "no covariates selected")
})

test_that("the fallback fit equals plain lm_robust", {
  dat <- noise_data()
  fit <- lm_lin_lasso(Y ~ Z, ~ X1 + X2, data = dat)
  plain <- estimatr::lm_robust(Y ~ Z, data = dat)

  expect_equal(broom::tidy(fit)$estimate, broom::tidy(plain)$estimate)
  expect_equal(broom::tidy(fit)$std.error, broom::tidy(plain)$std.error)
})

test_that("lm_robust_lasso adjusts additively", {
  fit <- lm_robust_lasso(Y ~ Z, ~ X_sig + X_n1 + X_n2, data = fit_data())

  expect_equal(adjustment(fit), "robust")
  expect_true("X_sig" %in% selected_covariates(fit))

  terms <- broom::tidy(fit)$term
  expect_true("X_sig" %in% terms)
  # additive, so no arm-specific interaction terms
  expect_false(any(grepl(":", terms)))

  z_row <- broom::tidy(fit)[broom::tidy(fit)$term == "Z", ]
  expect_lt(z_row$conf.low, 0.5)
  expect_gt(z_row$conf.high, 0.5)
})

test_that("lm_lin_lasso and lm_robust_lasso differ in their term structure", {
  dat <- fit_data()
  lin_terms <- broom::tidy(lm_lin_lasso(Y ~ Z, ~ X_sig + X_n1, data = dat))$term
  add_terms <- broom::tidy(lm_robust_lasso(Y ~ Z, ~ X_sig + X_n1, data = dat))$term

  expect_true(any(grepl(":", lin_terms)))
  expect_false(any(grepl(":", add_terms)))
})

test_that("lm_robust_lasso falls back when nothing is selected", {
  fit <- lm_robust_lasso(Y ~ Z, ~ X1 + X2, data = noise_data())
  expect_equal(adjustment(fit), "none")
  expect_equal(fallback_reason(fit), "no covariates selected")
})

test_that("a data column named 'clusters' does not silently cluster the fit", {
  dat <- fit_data()
  dat$clusters <- rep(1:30, length.out = nrow(dat))

  fit <- lm_lin_lasso(Y ~ Z, ~ X_sig + X_n1, data = dat)
  bare <- lm_lin_lasso(Y ~ Z, ~ X_sig + X_n1, data = fit_data())

  expect_equal(broom::tidy(fit)$std.error, broom::tidy(bare)$std.error)
  expect_null(fit$clusters)
})

test_that("clusters is honoured when supplied", {
  dat <- fit_data()
  dat$cl <- rep(1:30, length.out = nrow(dat))

  clustered <- lm_lin_lasso(Y ~ Z, ~ X_sig + X_n1, data = dat, clusters = cl)
  expect_equal(clustered$se_type, "CR2")
})

test_that("weights are passed through", {
  dat <- fit_data()
  dat$w <- runif(nrow(dat), 0.5, 1.5)

  weighted <- lm_lin_lasso(Y ~ Z, ~ X_sig + X_n1, data = dat, weights = w)
  unweighted <- lm_lin_lasso(Y ~ Z, ~ X_sig + X_n1, data = dat)
  expect_false(isTRUE(all.equal(
    broom::tidy(weighted)$estimate, broom::tidy(unweighted)$estimate
  )))
})

test_that("lasso_args reaches select_covariates", {
  dat <- fit_data()
  liberal <- lm_lin_lasso(Y ~ Z, ~ X_sig + X_n1 + X_n2, data = dat,
                          lasso_args = list(lambda_rule = "lambda.min"))
  parsimonious <- lm_lin_lasso(Y ~ Z, ~ X_sig + X_n1 + X_n2, data = dat)

  expect_gte(length(selected_covariates(liberal)),
             length(selected_covariates(parsimonious)))
})

test_that("lm_int_lasso fits the treatment-by-moderator interaction", {
  set.seed(2)
  n <- 800
  dat <- data.frame(
    Z = rep(0:1, n / 2),
    X_pid = rep(c(0, 1), each = n / 2),
    X_sig = rnorm(n),
    X_n1 = rnorm(n)
  )
  dat$Y <- 0.2 * dat$Z + 0.6 * dat$Z * dat$X_pid + 2 * dat$X_sig + rnorm(n)

  fit <- lm_int_lasso(Y ~ Z, moderator = "X_pid", data = dat,
                      covariates = ~ X_sig + X_n1)

  tidied <- broom::tidy(fit)
  int_row <- tidied[grepl(":", tidied$term), ]
  expect_equal(nrow(int_row), 1)
  expect_lt(int_row$conf.low, 0.6)
  expect_gt(int_row$conf.high, 0.6)
  expect_true("X_sig" %in% selected_covariates(fit))
})

test_that("lm_int_lasso works with no candidate covariates", {
  set.seed(2)
  n <- 400
  dat <- data.frame(Z = rep(0:1, n / 2), X_pid = rep(c(0, 1), each = n / 2))
  dat$Y <- 0.2 * dat$Z + 0.6 * dat$Z * dat$X_pid + rnorm(n)

  fit <- lm_int_lasso(Y ~ Z, moderator = "X_pid", data = dat)
  expect_equal(adjustment(fit), "none")
  expect_equal(selected_covariates(fit), character(0))
})

test_that("accessors return NULL for a fit this package did not produce", {
  plain <- estimatr::lm_robust(Y ~ Z, data = fit_data())
  expect_null(adjustment(plain))
  expect_null(selected_covariates(plain))
  expect_null(fallback_reason(plain))
})
