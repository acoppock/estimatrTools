make_dat <- function(n = 400, seed = 1) {
  set.seed(seed)
  dat <- data.frame(X1 = rnorm(n), X2 = rnorm(n), X3 = rnorm(n))
  dat$Y <- 1.5 * dat$X1 + rnorm(n)
  dat
}

test_that("selects the predictive covariate and drops the noise", {
  dat <- make_dat()
  expect_equal(lasso_select_one("Y", ~ X1 + X2 + X3, data = dat), "X1")
})

test_that("accepts a character vector, a formula, and a bare expression", {
  dat <- make_dat()
  expect_equal(
    lasso_select_one("Y", c("X1", "X2", "X3"), data = dat),
    lasso_select_one("Y", ~ X1 + X2 + X3, data = dat)
  )
  expect_equal(
    lasso_select_one("Y", X1 + X2 + X3, data = dat),
    lasso_select_one("Y", ~ X1 + X2 + X3, data = dat)
  )
})

test_that("accepts a bare column name and a variable holding a name", {
  dat <- make_dat()
  target <- "Y"
  expect_equal(lasso_select_one(Y, ~ X1 + X2, data = dat), "X1")
  expect_equal(lasso_select_one(target, ~ X1 + X2, data = dat), "X1")
})

test_that("returns the covariates in the order supplied", {
  set.seed(3)
  n <- 400
  dat <- data.frame(A = rnorm(n), B = rnorm(n), C = rnorm(n))
  dat$Y <- dat$A + dat$C + rnorm(n)
  expect_equal(lasso_select_one("Y", ~ C + B + A, data = dat), c("C", "A"))
})

test_that("subset restricts the rows selection runs on", {
  set.seed(4)
  n <- 600
  dat <- data.frame(X1 = rnorm(n), X2 = rnorm(n))
  # X2 predicts Y only in the second half
  half <- seq_len(n) > n / 2
  dat$Y <- 1.5 * dat$X1 + ifelse(half, 2 * dat$X2, 0) + rnorm(n)

  expect_true("X2" %in% lasso_select_one("Y", ~ X1 + X2, data = dat, subset = half))
  expect_false("X2" %in% lasso_select_one("Y", ~ X1 + X2, data = dat, subset = !half))
})

test_that("subset treats NA as FALSE and validates its length", {
  dat <- make_dat()
  keep <- rep(c(TRUE, NA), length.out = nrow(dat))
  expect_silent(lasso_select_one("Y", ~ X1 + X2, data = dat, subset = keep))
  expect_error(
    lasso_select_one("Y", ~ X1 + X2, data = dat, subset = c(TRUE, FALSE)),
    "length nrow"
  )
})

test_that("a factor covariate is returned whole, not as indicator columns", {
  set.seed(5)
  n <- 600
  dat <- data.frame(X1 = rnorm(n), Xf = factor(sample(c("a", "b", "c"), n, TRUE)))
  dat$Y <- 2 * (dat$Xf == "c") + rnorm(n)
  out <- lasso_select_one("Y", ~ X1 + Xf, data = dat)
  expect_true("Xf" %in% out)
  expect_false(any(grepl("^Xfb$|^Xfc$", out)))
})

test_that("unusable candidates are dropped rather than emptying the complete-case set", {
  dat <- make_dat()
  dat$X_allna <- NA_real_
  dat$X_const <- 1
  # Without the pre-filter the all-NA column would empty complete.cases and
  # selection would silently return nothing.
  expect_equal(lasso_select_one("Y", ~ X1 + X2 + X_allna + X_const, data = dat), "X1")
})

test_that("usability is judged on the subset, not the full data", {
  set.seed(6)
  n <- 400
  dat <- data.frame(X1 = rnorm(n), X2 = rnorm(n))
  half <- seq_len(n) > n / 2
  # X2 is constant among the rows actually used
  dat$X2[half] <- 0
  dat$Y <- 1.5 * dat$X1 + rnorm(n)
  expect_false("X2" %in% lasso_select_one("Y", ~ X1 + X2, data = dat, subset = half))
})

test_that("absent columns and an absent response return empty rather than erroring", {
  dat <- make_dat()
  expect_equal(lasso_select_one("Y", ~ X1 + nope, data = dat), "X1")
  expect_equal(lasso_select_one("not_a_column", ~ X1 + X2, data = dat), character(0))
  expect_equal(lasso_select_one("Y", character(0), data = dat), character(0))
})

test_that("pure noise selects nothing", {
  set.seed(7)
  n <- 400
  dat <- data.frame(X1 = rnorm(n), X2 = rnorm(n), X3 = rnorm(n))
  dat$Y <- rnorm(n)
  expect_length(lasso_select_one("Y", ~ X1 + X2 + X3, data = dat), 0)
})

test_that("the global random number state is restored", {
  dat <- make_dat()
  set.seed(99)
  before <- rnorm(1)
  set.seed(99)
  invisible(lasso_select_one("Y", ~ X1 + X2 + X3, data = dat))
  after <- rnorm(1)
  expect_identical(before, after)
})

test_that("selection is reproducible across calls and responds to the seed", {
  dat <- make_dat()
  a <- lasso_select_one("Y", ~ X1 + X2 + X3, data = dat)
  b <- lasso_select_one("Y", ~ X1 + X2 + X3, data = dat)
  expect_identical(a, b)
  # a different seed is allowed to differ, but must itself be reproducible
  c1 <- lasso_select_one("Y", ~ X1 + X2 + X3, data = dat, seed = 7)
  c2 <- lasso_select_one("Y", ~ X1 + X2 + X3, data = dat, seed = 7)
  expect_identical(c1, c2)
})

test_that("lambda.min is at least as permissive as lambda.1se", {
  dat <- make_dat()
  n_1se <- length(lasso_select_one("Y", ~ X1 + X2 + X3, data = dat,
                                   lambda_rule = "lambda.1se"))
  n_min <- length(lasso_select_one("Y", ~ X1 + X2 + X3, data = dat,
                                   lambda_rule = "lambda.min"))
  expect_gte(n_min, n_1se)
})
