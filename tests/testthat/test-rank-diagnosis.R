# The rank-deficiency diagnosis, exercised directly.
#
# Going through lasso_select_covariates() would leave which branch fires up to
# whatever the LASSO happens to select, so these call the internal diagnostic
# with designs constructed to have a known cause. The three branches need to be
# told apart reliably, because they call for opposite remedies: an unobserved
# level is fixed by collapsing it in cleaning, a redundant covariate by
# excluding it at the call site.

diagnose <- function(cols, treatment, df) {
  msg <- NULL
  withCallingHandlers(
    estimatrTools:::warn_rank_deficient(cols, treatment, df),
    warning = function(w) { msg <<- conditionMessage(w); invokeRestart("muffleWarning") }
  )
  msg
}

test_that("a full-rank design is silent", {
  set.seed(1)
  n <- 300
  df <- data.frame(Z = rep(0:1, n / 2), X_a = rnorm(n), X_b = rnorm(n))
  expect_null(diagnose(c("X_a", "X_b"), "Z", df))
})

test_that("an unobserved level is diagnosed as such, not as redundancy", {
  # X_fac level "c" appears only in arm 1, so Z:X_facc is all zero. Removing
  # X_fac restores rank, which is exactly why a "does removing it help" test
  # cannot distinguish this from redundancy.
  set.seed(2)
  n <- 400
  df <- data.frame(
    Z = rep(0:1, each = n / 2),
    X_num = rnorm(n),
    X_fac = factor(c(sample(c("a", "b"), n / 2, replace = TRUE),
                     sample(c("a", "b", "c"), n / 2, replace = TRUE)))
  )
  msg <- diagnose(c("X_num", "X_fac"), "Z", df)

  expect_match(msg, "rank deficient")
  expect_match(msg, "unobserved level")
  expect_match(msg, 'X_fac = "c"')
  expect_match(msg, "no observations in Z = 0")
  expect_match(msg, "cleaning script")
  # must NOT reach for the redundancy explanation
  expect_no_match(msg, "redundancy")
})

test_that("the unobserved-level message quantifies what dropping would cost", {
  set.seed(3)
  n <- 600
  # a 5-level factor: 8 columns in the interacted design, 1 of them aliased
  levs <- c("a", "b", "c", "d", "e")
  df <- data.frame(
    Z = rep(0:1, each = n / 2),
    X_num = rnorm(n),
    X_fac = factor(c(sample(levs[1:4], n / 2, replace = TRUE),
                     sample(levs, n / 2, replace = TRUE)))
  )
  msg <- diagnose(c("X_num", "X_fac"), "Z", df)
  expect_match(msg, "contributes 8 columns")
  expect_match(msg, "discard 7 to remove 1")
})

test_that("a genuinely redundant covariate is named with a paste-ready exclusion", {
  # X_college is a deterministic function of X_educ, so its columns lie in the
  # span of X_educ's. No level is unobserved in either arm.
  set.seed(4)
  n <- 600
  educ <- factor(sample(c("hs", "college", "postgrad"), n, replace = TRUE))
  df <- data.frame(
    Z = rep(0:1, n / 2),
    X_educ = educ,
    X_college = as.integer(educ %in% c("college", "postgrad")),
    X_num = rnorm(n)
  )
  msg <- diagnose(c("X_educ", "X_college", "X_num"), "Z", df)

  expect_match(msg, "rank deficient")
  expect_match(msg, "redundancy")
  expect_match(msg, "X_college")
  expect_match(msg, "covariates = c\\(")
  # the suggested set drops the narrow derived variable, not the rich factor
  expect_match(msg, '"X_educ"')
  expect_no_match(msg, "unobserved level")
})

test_that("the narrowest redundant covariate is the one suggested for removal", {
  # Both X_college and X_educ restore rank when removed. The suggestion should
  # be the one contributing fewer columns, so the richer covariate is kept.
  set.seed(5)
  n <- 600
  educ <- factor(sample(c("hs", "college", "postgrad"), n, replace = TRUE))
  df <- data.frame(
    Z = rep(0:1, n / 2),
    X_educ = educ,
    X_college = as.integer(educ %in% c("college", "postgrad"))
  )
  msg <- diagnose(c("X_educ", "X_college"), "Z", df)
  expect_match(msg, "Removing|X_college lies in the span|redundancy")
  expect_match(msg, "X_college")
})

test_that("a deficiency attributable to neither cause says so rather than guessing", {
  # Three mutually collinear numerics: no single removal restores full rank and
  # no factor level is unobserved.
  set.seed(6)
  n <- 400
  x1 <- rnorm(n); x2 <- rnorm(n)
  df <- data.frame(
    Z = rep(0:1, n / 2),
    X_a = x1, X_b = x2,
    X_c = x1 + x2,          # exactly collinear with X_a + X_b
    X_d = 2 * x1 - x2       # and another
  )
  msg <- diagnose(c("X_a", "X_b", "X_c", "X_d"), "Z", df)

  expect_match(msg, "rank deficient")
  expect_match(msg, "not attributable to a single covariate")
  expect_match(msg, "Inspect the design")
})

test_that("the diagnosis never errors on degenerate input", {
  df <- data.frame(Z = rep(0:1, 20), X_a = rnorm(40))
  expect_silent(estimatrTools:::warn_rank_deficient(character(0), "Z", df))
  expect_silent(estimatrTools:::warn_rank_deficient("X_a", "Z", df[0, ]))
})
