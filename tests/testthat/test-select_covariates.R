signal_data <- function(n = 500, seed = 1) {
  set.seed(seed)
  dat <- data.frame(
    Z     = rep(0:1, n / 2),
    X_sig = rnorm(n),
    X_n1  = rnorm(n),
    X_n2  = rnorm(n),
    X_n3  = rnorm(n)
  )
  dat$Y <- 0.5 * dat$Z + 2 * dat$X_sig + rnorm(n)
  dat
}

test_that("a strong outcome predictor is selected and noise is not", {
  sel <- lasso_select_covariates(Y ~ Z, ~ X_sig + X_n1 + X_n2 + X_n3, data = signal_data())
  expect_true("X_sig" %in% all.vars(sel))
  expect_false("X_n1" %in% all.vars(sel))
})

test_that("a covariate that only predicts assignment is selected", {
  set.seed(3)
  n <- 600
  X_assign <- rnorm(n)
  dat <- data.frame(
    X_assign = X_assign,
    X_noise  = rnorm(n)
  )
  # Assignment depends strongly on X_assign; the outcome does not
  dat$Z <- rbinom(n, 1, plogis(2 * X_assign))
  dat$Y <- 0.5 * dat$Z + rnorm(n)

  sel <- lasso_select_covariates(Y ~ Z, ~ X_assign + X_noise, data = dat)
  expect_true("X_assign" %in% all.vars(sel))
})

test_that("a covariate that only matters in one arm is selected", {
  set.seed(5)
  n <- 800
  dat <- data.frame(Z = rep(0:1, n / 2), X_het = rnorm(n), X_noise = rnorm(n))
  # X_het predicts Y only among the treated: a pooled regression would miss it
  dat$Y <- 0.3 * dat$Z + 3 * dat$Z * dat$X_het + rnorm(n)

  sel <- lasso_select_covariates(Y ~ Z, ~ X_het + X_noise, data = dat)
  expect_true("X_het" %in% all.vars(sel))
})

test_that("pure noise selects nothing and returns ~1", {
  set.seed(11)
  n <- 300
  dat <- data.frame(Z = rep(0:1, n / 2), X1 = rnorm(n), X2 = rnorm(n))
  dat$Y <- 0.5 * dat$Z + rnorm(n)

  sel <- lasso_select_covariates(Y ~ Z, ~ X1 + X2, data = dat)
  expect_equal(all.vars(sel), character(0))
})

test_that("an all-NA candidate does not disable selection", {
  dat <- signal_data()
  dat$X_unmeasured <- NA_real_

  sel <- lasso_select_covariates(Y ~ Z, ~ X_sig + X_unmeasured, data = dat)
  expect_true("X_sig" %in% all.vars(sel))
  expect_false("X_unmeasured" %in% all.vars(sel))
})

test_that("a constant candidate does not disable selection", {
  dat <- signal_data()
  dat$X_const <- 1

  sel <- lasso_select_covariates(Y ~ Z, ~ X_sig + X_const, data = dat)
  expect_true("X_sig" %in% all.vars(sel))
})

test_that("an absent candidate column is ignored", {
  sel <- lasso_select_covariates(Y ~ Z, c("X_sig", "X_not_in_data"), data = signal_data())
  expect_true("X_sig" %in% all.vars(sel))
})

test_that("covariates accept a formula, a character vector, and a bare expression", {
  dat <- signal_data()
  from_formula <- lasso_select_covariates(Y ~ Z, ~ X_sig + X_n1, data = dat)
  from_chr     <- lasso_select_covariates(Y ~ Z, c("X_sig", "X_n1"), data = dat)
  from_expr    <- lasso_select_covariates(Y ~ Z, X_sig + X_n1, data = dat)

  expect_equal(all.vars(from_formula), all.vars(from_chr))
  expect_equal(all.vars(from_formula), all.vars(from_expr))
})

test_that("a factor is retained under its own name, not its indicator columns", {
  set.seed(13)
  n <- 600
  dat <- data.frame(
    Z = rep(0:1, n / 2),
    X_region = factor(sample(c("N", "S", "E"), n, replace = TRUE)),
    X_noise = rnorm(n)
  )
  dat$Y <- 0.5 * dat$Z + 2 * (dat$X_region == "S") + rnorm(n)

  sel <- lasso_select_covariates(Y ~ Z, ~ X_region + X_noise, data = dat)
  expect_true("X_region" %in% all.vars(sel))
  expect_false(any(grepl("X_regionS", all.vars(sel))))
})

test_that("one covariate name prefixing another does not over-include", {
  # X_ed's indicator columns must not be attributed to X_e
  set.seed(17)
  n <- 700
  dat <- data.frame(
    Z    = rep(0:1, n / 2),
    X_e  = rnorm(n),
    X_ed = factor(sample(c("a", "b", "c"), n, replace = TRUE))
  )
  dat$Y <- 0.5 * dat$Z + 2.5 * (dat$X_ed == "c") + rnorm(n)

  sel <- lasso_select_covariates(Y ~ Z, ~ X_e + X_ed, data = dat)
  expect_true("X_ed" %in% all.vars(sel))
  expect_false("X_e" %in% all.vars(sel))
})

test_that("degenerate inputs return ~1 rather than erroring", {
  dat <- signal_data()

  one_arm <- dat
  one_arm$Z <- 0
  expect_equal(all.vars(lasso_select_covariates(Y ~ Z, ~ X_sig, data = one_arm)), character(0))

  no_covs <- lasso_select_covariates(Y ~ Z, character(0), data = dat)
  expect_equal(all.vars(no_covs), character(0))

  empty <- dat[0, ]
  expect_equal(all.vars(lasso_select_covariates(Y ~ Z, ~ X_sig, data = empty)), character(0))

  missing_outcome <- lasso_select_covariates(NotHere ~ Z, ~ X_sig, data = dat)
  expect_equal(all.vars(missing_outcome), character(0))
})

test_that("selection is reproducible and does not disturb the caller's RNG", {
  dat <- signal_data()

  a <- lasso_select_covariates(Y ~ Z, ~ X_sig + X_n1 + X_n2, data = dat)
  b <- lasso_select_covariates(Y ~ Z, ~ X_sig + X_n1 + X_n2, data = dat)
  expect_equal(all.vars(a), all.vars(b))

  set.seed(99)
  expected <- rnorm(3)
  set.seed(99)
  invisible(lasso_select_covariates(Y ~ Z, ~ X_sig + X_n1 + X_n2, data = dat))
  expect_equal(rnorm(3), expected)
})

test_that("lambda.min selects at least as much as lambda.1se", {
  dat <- signal_data()
  parsimonious <- lasso_select_covariates(Y ~ Z, ~ X_sig + X_n1 + X_n2 + X_n3, data = dat)
  liberal <- lasso_select_covariates(Y ~ Z, ~ X_sig + X_n1 + X_n2 + X_n3, data = dat,
                               lambda_rule = "lambda.min")

  expect_gte(length(all.vars(liberal)), length(all.vars(parsimonious)))
})

test_that("lambda_rule is validated", {
  expect_error(
    lasso_select_covariates(Y ~ Z, ~ X_sig, data = signal_data(), lambda_rule = "nope"),
    "should be one of"
  )
})


# rank-deficiency diagnosis ----

test_that("an empty level-by-arm cell warns and names the cell", {
  set.seed(21)
  n <- 400
  dat <- data.frame(
    Z = rep(0:1, each = n / 2),
    X_sig = rnorm(n),
    X_fac = factor(c(sample(c("a", "b", "c"), n / 2, replace = TRUE),
                     sample(c("a", "b"), n / 2, replace = TRUE)))
  )
  dat$Y <- 0.5 * dat$Z + 2 * dat$X_sig + 1.5 * (dat$X_fac == "c") + rnorm(n)

  expect_warning(
    sel <- lasso_select_covariates(Y ~ Z, ~ X_sig + X_fac, data = dat),
    "rank deficient"
  )
  # the selection is returned unchanged: nothing is pruned
  expect_true("X_fac" %in% all.vars(sel))
})

test_that("multi-arm treatments select per arm and fit", {
  # meta_propaganda has three- and four-arm designs; the per-arm outcome loop
  # runs once per arm, so this exercises a path two-arm tests never reach.
  set.seed(31)
  n <- 600
  d <- data.frame(
    Z = factor(rep(c("C", "T1", "T2"), each = n / 3)),
    X_sig = rnorm(n), X_noise = rnorm(n)
  )
  d$Y <- 0.3 * (d$Z == "T1") + 0.6 * (d$Z == "T2") + 2 * d$X_sig + rnorm(n)

  sel <- lasso_select_covariates(Y ~ Z, ~ X_sig + X_noise, data = d)
  expect_true("X_sig" %in% all.vars(sel))
  expect_false("X_noise" %in% all.vars(sel))
})

test_that("a covariate observed in only one arm is still handled end to end", {
  # The deryugina and allen shapes met in the wild, run through the whole
  # pipeline rather than the selector alone: a usable estimate must come back.
  set.seed(32)
  n <- 400
  d <- data.frame(
    Z = rep(0:1, each = n / 2),
    X_sig = rnorm(n),
    X_fac = factor(c(sample(c("a", "b"), n / 2, replace = TRUE),
                     sample(c("a", "b", "c"), n / 2, replace = TRUE)))
  )
  d$Y <- 0.5 * d$Z + 2 * d$X_sig + rnorm(n)

  fit <- suppressWarnings(lm_lin_lasso(Y ~ Z, ~ X_sig + X_fac, data = d))
  z_row <- broom::tidy(fit)[broom::tidy(fit)$term == "Z", ]
  expect_equal(nrow(z_row), 1)
  expect_true(is.finite(z_row$estimate))
  expect_true(is.finite(z_row$std.error))
})

test_that("a well-conditioned selection warns about nothing", {
  expect_silent(lasso_select_covariates(Y ~ Z, ~ X_sig + X_n1 + X_n2, data = signal_data()))
})

test_that("the selection is never pruned by the rank check", {
  set.seed(21)
  n <- 400
  dat <- data.frame(
    Z = rep(0:1, each = n / 2),
    X_sig = rnorm(n),
    X_fac = factor(c(sample(c("a", "b", "c"), n / 2, replace = TRUE),
                     sample(c("a", "b"), n / 2, replace = TRUE)))
  )
  dat$Y <- 0.5 * dat$Z + 2 * dat$X_sig + 1.5 * (dat$X_fac == "c") + rnorm(n)

  sel <- suppressWarnings(lasso_select_covariates(Y ~ Z, ~ X_sig + X_fac, data = dat))
  # lm_lin still fits: it aliases the offending column itself
  fit <- suppressWarnings(estimatr::lm_lin(Y ~ Z, covariates = sel, data = dat))
  z_row <- broom::tidy(fit)[2, ]
  expect_true(is.finite(z_row$std.error))
})


# small and degenerate samples ----
#
# Every problem this package met in the wild arrived through a subgroup: a
# politics stratum, a country subset, an arm split. Small and degenerate inputs
# are the normal case there, not the exotic one.

test_that("a sample too small to cross-validate returns ~1 rather than erroring", {
  set.seed(41)
  d <- data.frame(Z = rep(0:1, 4), X_a = rnorm(8), X_b = rnorm(8))
  d$Y <- 0.5 * d$Z + rnorm(8)
  expect_equal(all.vars(lasso_select_covariates(Y ~ Z, ~ X_a + X_b, data = d)), character(0))
})

test_that("no complete cases returns ~1", {
  set.seed(42)
  n <- 100
  d <- data.frame(Z = rep(0:1, n / 2), X_a = rnorm(n), X_b = rnorm(n))
  d$Y <- 0.5 * d$Z + rnorm(n)
  # X_a and X_b are missing on disjoint halves, so no row is complete on both
  d$X_a[1:(n / 2)] <- NA
  d$X_b[(n / 2 + 1):n] <- NA
  expect_equal(all.vars(lasso_select_covariates(Y ~ Z, ~ X_a + X_b, data = d)), character(0))
})

test_that("a single candidate is kept without cross-validating", {
  # cv.glmnet needs at least two columns to choose between. With exactly one
  # there is nothing to select, so it is included. Pinned because it means a
  # one-candidate call never really runs a LASSO.
  set.seed(43)
  d <- data.frame(Z = rep(0:1, 6), X_a = rnorm(12))
  d$Y <- 0.5 * d$Z + rnorm(12)

  fit <- lm_lin_lasso(Y ~ Z, ~ X_a, data = d)
  expect_equal(adjustment(fit), "lin")
  expect_equal(selected_covariates(fit), "X_a")
  expect_true(is.finite(broom::tidy(fit)$estimate[2]))
})

test_that("a tiny subgroup with several candidates falls back rather than failing", {
  # Fewer than ten rows per arm, so the per-arm outcome equations bail out.
  set.seed(46)
  d <- data.frame(Z = rep(0:1, 6), X_a = rnorm(12), X_b = rnorm(12), X_c = rnorm(12))
  d$Y <- 0.5 * d$Z + rnorm(12)

  fit <- lm_lin_lasso(Y ~ Z, ~ X_a + X_b + X_c, data = d)
  expect_true(adjustment(fit) %in% c("lin", "none"))
  expect_true(is.finite(broom::tidy(fit)$estimate[2]))
})

test_that("haven_labelled candidates are handled", {
  skip_if_not_installed("haven")
  set.seed(44)
  n <- 400
  d <- data.frame(Z = rep(0:1, n / 2), X_noise = rnorm(n))
  d$X_lab <- haven::labelled(sample(1:3, n, replace = TRUE),
                             c(low = 1, mid = 2, high = 3))
  d$Y <- 0.5 * d$Z + 1.5 * as.numeric(d$X_lab) + rnorm(n)

  sel <- lasso_select_covariates(Y ~ Z, ~ X_lab + X_noise, data = d)
  expect_true("X_lab" %in% all.vars(sel))
})
