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
  sel <- select_covariates(Y ~ Z, ~ X_sig + X_n1 + X_n2 + X_n3, data = signal_data())
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

  sel <- select_covariates(Y ~ Z, ~ X_assign + X_noise, data = dat)
  expect_true("X_assign" %in% all.vars(sel))
})

test_that("a covariate that only matters in one arm is selected", {
  set.seed(5)
  n <- 800
  dat <- data.frame(Z = rep(0:1, n / 2), X_het = rnorm(n), X_noise = rnorm(n))
  # X_het predicts Y only among the treated: a pooled regression would miss it
  dat$Y <- 0.3 * dat$Z + 3 * dat$Z * dat$X_het + rnorm(n)

  sel <- select_covariates(Y ~ Z, ~ X_het + X_noise, data = dat)
  expect_true("X_het" %in% all.vars(sel))
})

test_that("pure noise selects nothing and returns ~1", {
  set.seed(11)
  n <- 300
  dat <- data.frame(Z = rep(0:1, n / 2), X1 = rnorm(n), X2 = rnorm(n))
  dat$Y <- 0.5 * dat$Z + rnorm(n)

  sel <- select_covariates(Y ~ Z, ~ X1 + X2, data = dat)
  expect_equal(all.vars(sel), character(0))
})

test_that("an all-NA candidate does not disable selection", {
  dat <- signal_data()
  dat$X_unmeasured <- NA_real_

  sel <- select_covariates(Y ~ Z, ~ X_sig + X_unmeasured, data = dat)
  expect_true("X_sig" %in% all.vars(sel))
  expect_false("X_unmeasured" %in% all.vars(sel))
})

test_that("a constant candidate does not disable selection", {
  dat <- signal_data()
  dat$X_const <- 1

  sel <- select_covariates(Y ~ Z, ~ X_sig + X_const, data = dat)
  expect_true("X_sig" %in% all.vars(sel))
})

test_that("an absent candidate column is ignored", {
  sel <- select_covariates(Y ~ Z, c("X_sig", "X_not_in_data"), data = signal_data())
  expect_true("X_sig" %in% all.vars(sel))
})

test_that("covariates accept a formula, a character vector, and a bare expression", {
  dat <- signal_data()
  from_formula <- select_covariates(Y ~ Z, ~ X_sig + X_n1, data = dat)
  from_chr     <- select_covariates(Y ~ Z, c("X_sig", "X_n1"), data = dat)
  from_expr    <- select_covariates(Y ~ Z, X_sig + X_n1, data = dat)

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

  sel <- select_covariates(Y ~ Z, ~ X_region + X_noise, data = dat)
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

  sel <- select_covariates(Y ~ Z, ~ X_e + X_ed, data = dat)
  expect_true("X_ed" %in% all.vars(sel))
  expect_false("X_e" %in% all.vars(sel))
})

test_that("degenerate inputs return ~1 rather than erroring", {
  dat <- signal_data()

  one_arm <- dat
  one_arm$Z <- 0
  expect_equal(all.vars(select_covariates(Y ~ Z, ~ X_sig, data = one_arm)), character(0))

  no_covs <- select_covariates(Y ~ Z, character(0), data = dat)
  expect_equal(all.vars(no_covs), character(0))

  empty <- dat[0, ]
  expect_equal(all.vars(select_covariates(Y ~ Z, ~ X_sig, data = empty)), character(0))

  missing_outcome <- select_covariates(NotHere ~ Z, ~ X_sig, data = dat)
  expect_equal(all.vars(missing_outcome), character(0))
})

test_that("selection is reproducible and does not disturb the caller's RNG", {
  dat <- signal_data()

  a <- select_covariates(Y ~ Z, ~ X_sig + X_n1 + X_n2, data = dat)
  b <- select_covariates(Y ~ Z, ~ X_sig + X_n1 + X_n2, data = dat)
  expect_equal(all.vars(a), all.vars(b))

  set.seed(99)
  expected <- rnorm(3)
  set.seed(99)
  invisible(select_covariates(Y ~ Z, ~ X_sig + X_n1 + X_n2, data = dat))
  expect_equal(rnorm(3), expected)
})

test_that("lambda.min selects at least as much as lambda.1se", {
  dat <- signal_data()
  parsimonious <- select_covariates(Y ~ Z, ~ X_sig + X_n1 + X_n2 + X_n3, data = dat)
  liberal <- select_covariates(Y ~ Z, ~ X_sig + X_n1 + X_n2 + X_n3, data = dat,
                               lambda_rule = "lambda.min")

  expect_gte(length(all.vars(liberal)), length(all.vars(parsimonious)))
})

test_that("lambda_rule is validated", {
  expect_error(
    select_covariates(Y ~ Z, ~ X_sig, data = signal_data(), lambda_rule = "nope"),
    "should be one of"
  )
})
