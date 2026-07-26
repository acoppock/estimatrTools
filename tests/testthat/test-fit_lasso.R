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

test_that("lasso_args reaches lasso_select_covariates", {
  dat <- fit_data()
  liberal <- lm_lin_lasso(Y ~ Z, ~ X_sig + X_n1 + X_n2, data = dat,
                          lasso_args = list(lambda_rule = "lambda.min"))
  parsimonious <- lm_lin_lasso(Y ~ Z, ~ X_sig + X_n1 + X_n2, data = dat)

  expect_gte(length(selected_covariates(liberal)),
             length(selected_covariates(parsimonious)))
})

test_that("lm_moderator_lasso fits the treatment-by-moderator interaction", {
  set.seed(2)
  n <- 800
  dat <- data.frame(
    Z = rep(0:1, n / 2),
    X_pid = rep(c(0, 1), each = n / 2),
    X_sig = rnorm(n),
    X_n1 = rnorm(n)
  )
  dat$Y <- 0.2 * dat$Z + 0.6 * dat$Z * dat$X_pid + 2 * dat$X_sig + rnorm(n)

  fit <- lm_moderator_lasso(Y ~ Z, moderator = "X_pid", data = dat,
                      covariates = ~ X_sig + X_n1)

  tidied <- broom::tidy(fit)
  int_row <- tidied[grepl(":", tidied$term), ]
  expect_equal(nrow(int_row), 1)
  expect_lt(int_row$conf.low, 0.6)
  expect_gt(int_row$conf.high, 0.6)
  expect_true("X_sig" %in% selected_covariates(fit))
})

test_that("lm_moderator_lasso works with no candidate covariates", {
  set.seed(2)
  n <- 400
  dat <- data.frame(Z = rep(0:1, n / 2), X_pid = rep(c(0, 1), each = n / 2))
  dat$Y <- 0.2 * dat$Z + 0.6 * dat$Z * dat$X_pid + rnorm(n)

  fit <- lm_moderator_lasso(Y ~ Z, moderator = "X_pid", data = dat)
  expect_equal(adjustment(fit), "none")
  expect_equal(selected_covariates(fit), character(0))
})

test_that("accessors return NULL for a fit this package did not produce", {
  plain <- estimatr::lm_robust(Y ~ Z, data = fit_data())
  expect_null(adjustment(plain))
  expect_null(selected_covariates(plain))
  expect_null(fallback_reason(plain))
})


# fallback_summary ----

test_that("fallback_summary reports one row per fit", {
  fits <- list(
    informative = lm_lin_lasso(Y ~ Z, ~ X_sig + X_n1, data = fit_data()),
    noise_only  = lm_lin_lasso(Y ~ Z, ~ X1 + X2, data = noise_data())
  )
  s <- fallback_summary(fits)

  expect_equal(nrow(s), 2)
  expect_equal(s$fit, c("informative", "noise_only"))
  expect_equal(s$adjustment, c("lin", "none"))
  expect_true(is.na(s$fallback_reason[1]))
  expect_equal(s$fallback_reason[2], "no covariates selected")
  expect_true(s$n_selected[1] > 0)
  expect_equal(s$n_selected[2], 0)
})

test_that("only_fallbacks filters to the substitutions", {
  fits <- list(
    a = lm_lin_lasso(Y ~ Z, ~ X_sig + X_n1, data = fit_data()),
    b = lm_lin_lasso(Y ~ Z, ~ X1 + X2, data = noise_data())
  )
  s <- fallback_summary(fits, only_fallbacks = TRUE)
  expect_equal(nrow(s), 1)
  expect_equal(s$fit, "b")
})

test_that("unnamed lists are labelled by index", {
  fits <- list(lm_lin_lasso(Y ~ Z, ~ X1 + X2, data = noise_data()))
  expect_equal(fallback_summary(fits)$fit, "1")
})

test_that("fits from outside the package are reported as NA, not dropped", {
  fits <- list(
    ours   = lm_lin_lasso(Y ~ Z, ~ X_sig, data = fit_data()),
    theirs = estimatr::lm_robust(Y ~ Z, data = fit_data())
  )
  s <- fallback_summary(fits)
  expect_equal(nrow(s), 2)
  expect_true(is.na(s$adjustment[2]))
  expect_true(is.na(s$n_selected[2]))
})

test_that("fallback_summary handles an empty list and rejects a non-list", {
  expect_equal(nrow(fallback_summary(list())), 0)
  expect_error(fallback_summary(1:3), "expected a list")
})

test_that("fallback_summary covers lm_robust_lasso and lm_moderator_lasso too", {
  set.seed(4)
  n <- 400
  dat <- data.frame(Z = rep(0:1, n / 2), X_pid = rep(c(0, 1), each = n / 2),
                    X_sig = rnorm(n))
  dat$Y <- 0.3 * dat$Z + 2 * dat$X_sig + rnorm(n)

  fits <- list(
    add = lm_robust_lasso(Y ~ Z, ~ X_sig, data = dat),
    int = lm_moderator_lasso(Y ~ Z, moderator = "X_pid", data = dat, covariates = ~ X_sig)
  )
  s <- fallback_summary(fits)
  expect_equal(nrow(s), 2)
  expect_true(all(s$adjustment %in% c("lin", "robust", "none")))
})


# fallback_log ----

test_that("the log records every call, including discarded fits", {
  reset_fallback_log()

  # Fits built anonymously and thrown away, which is the dominant calling shape
  invisible(lapply(1:2, function(i) lm_lin_lasso(Y ~ Z, ~ X_sig + X_n1, data = fit_data())))
  invisible(lm_lin_lasso(Y ~ Z, ~ X1 + X2, data = noise_data()))

  log <- fallback_log()
  expect_equal(nrow(log), 3)
  expect_equal(log$fn, rep("lm_lin_lasso", 3))
  expect_equal(log$outcome, rep("Y", 3))
  expect_equal(log$treatment, rep("Z", 3))
  expect_equal(log$adjustment, c("lin", "lin", "none"))
  expect_equal(sum(!is.na(log$fallback_reason)), 1)
  expect_equal(log$call_index, 1:3)
})

test_that("fallbacks_only narrows to the substitutions", {
  reset_fallback_log()
  invisible(lm_lin_lasso(Y ~ Z, ~ X_sig, data = fit_data()))
  invisible(lm_lin_lasso(Y ~ Z, ~ X1 + X2, data = noise_data()))

  expect_equal(nrow(fallback_log(fallbacks_only = TRUE)), 1)
  expect_equal(fallback_log(fallbacks_only = TRUE)$adjustment, "none")
})

test_that("the adjustment rate is recoverable from the log", {
  reset_fallback_log()
  invisible(lm_lin_lasso(Y ~ Z, ~ X_sig, data = fit_data()))
  invisible(lm_lin_lasso(Y ~ Z, ~ X1 + X2, data = noise_data()))
  invisible(lm_lin_lasso(Y ~ Z, ~ X1 + X2, data = noise_data()))

  rate <- mean(fallback_log()$adjustment != "none")
  expect_equal(rate, 1 / 3)
})

test_that("reset clears the log and an empty log has the right shape", {
  reset_fallback_log()
  expect_equal(nrow(fallback_log()), 0)
  expect_named(fallback_log(), c("call_index", "fn", "outcome", "treatment",
                                 "adjustment", "n_selected", "fallback_reason"))
})

test_that("logging can be switched off", {
  reset_fallback_log()
  withr_opt <- options(estimatrTools.log = FALSE)
  on.exit(options(withr_opt), add = TRUE)

  invisible(lm_lin_lasso(Y ~ Z, ~ X_sig, data = fit_data()))
  expect_equal(nrow(fallback_log()), 0)
})

test_that("lm_robust_lasso and lm_moderator_lasso are logged under their own names", {
  reset_fallback_log()
  set.seed(4)
  n <- 400
  dat <- data.frame(Z = rep(0:1, n / 2), X_pid = rep(c(0, 1), each = n / 2),
                    X_sig = rnorm(n))
  dat$Y <- 0.3 * dat$Z + 2 * dat$X_sig + rnorm(n)

  invisible(lm_robust_lasso(Y ~ Z, ~ X_sig, data = dat))
  invisible(lm_moderator_lasso(Y ~ Z, moderator = "X_pid", data = dat, covariates = ~ X_sig))

  expect_equal(fallback_log()$fn, c("lm_robust_lasso", "lm_moderator_lasso"))
})
