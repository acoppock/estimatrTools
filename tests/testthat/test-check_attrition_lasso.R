make_attrition_dat <- function(n = 500, seed = 42, differential = TRUE) {
  set.seed(seed)
  dat <- data.frame(
    Z = rep(0:1, n / 2),
    X_age = rnorm(n, 50, 10),
    X_income = rnorm(n, 50000, 10000)
  )
  dat$Y_outcome <- 0.3 * dat$Z + 0.5 * as.vector(scale(dat$X_age)) + rnorm(n)
  p_miss <- if (differential) 0.10 + 0.20 * dat$Z else rep(0.15, n)
  dat$Y_outcome[rbinom(n, 1, p_miss) == 1] <- NA
  dat
}

test_that("returns one row per outcome with the documented columns", {
  dat <- make_attrition_dat()
  dat$Y_second <- rnorm(nrow(dat))
  dat$Y_second[1:40] <- NA
  out <- check_attrition_lasso(dat, Z, covariates = c("X_age", "X_income"))

  expect_equal(nrow(out), 2)
  expect_true(all(c("outcome", "n_assigned", "n_missing", "pct_missing",
                    "p_simple", "n_outcome_eq", "n_dropout_eq", "n_selected",
                    "selected_covariates", "df1", "epv", "epv_adequate",
                    "p_interacted", "estimable", "flag_simple",
                    "flag_interacted", "flag") %in% names(out)))
})

test_that("detects differential attrition and clears a clean study", {
  bad <- check_attrition_lasso(make_attrition_dat(differential = TRUE), Z,
                               outcomes = "Y_outcome",
                               covariates = c("X_age", "X_income"))
  good <- check_attrition_lasso(make_attrition_dat(differential = FALSE), Z,
                                outcomes = "Y_outcome",
                                covariates = c("X_age", "X_income"))
  expect_lt(bad$p_simple, 0.05)
  expect_true(bad$flag)
  expect_gt(good$p_simple, 0.05)
  expect_false(good$flag_simple)
})

test_that("a bare treatment name and a string are equivalent", {
  dat <- make_attrition_dat()
  a <- check_attrition_lasso(dat, Z, outcomes = "Y_outcome",
                             covariates = c("X_age", "X_income"))
  b <- check_attrition_lasso(dat, "Z", outcomes = "Y_outcome",
                             covariates = c("X_age", "X_income"))
  expect_equal(a, b)
})

test_that("an unvarying missingness indicator is reported as unestimable", {
  dat <- make_attrition_dat()
  dat$Y_none <- rnorm(nrow(dat))            # nobody drops out
  dat$Y_all <- NA_real_                     # everybody does

  none <- check_attrition_lasso(dat, Z, outcomes = "Y_none",
                                covariates = c("X_age", "X_income"))
  all_missing <- check_attrition_lasso(dat, Z, outcomes = "Y_all",
                                      covariates = c("X_age", "X_income"))

  for (out in list(none, all_missing)) {
    expect_false(out$estimable)
    expect_true(is.na(out$p_simple))
    expect_true(is.na(out$p_interacted))
    expect_false(out$flag)
  }
  # p = 1 would put a spike at 1 into any uniform-reference diagnostic
  expect_false(isTRUE(none$p_simple == 1))
})

test_that("a covariate not present in the data is an error, not a silent drop", {
  dat <- make_attrition_dat()
  expect_error(
    check_attrition_lasso(dat, Z, outcomes = "Y_outcome",
                          covariates = c("X_age", "X_typo")),
    "not found in data"
  )
})

test_that("an absent treatment or outcome column errors", {
  dat <- make_attrition_dat()
  expect_error(check_attrition_lasso(dat, Znope, outcomes = "Y_outcome"),
               "treatment column")
  expect_error(check_attrition_lasso(dat, Z, outcomes = "Y_nope"),
               "outcome column")
})

test_that("factor and NA-bearing candidates do not disable selection", {
  dat <- make_attrition_dat()
  dat$X_party <- factor(sample(c("D", "R", "I"), nrow(dat), TRUE))
  dat$X_gappy <- dat$X_age
  dat$X_gappy[1:50] <- NA

  out <- check_attrition_lasso(
    dat, Z, outcomes = "Y_outcome",
    covariates = c("X_age", "X_income", "X_party", "X_gappy")
  )
  # X_age drives the outcome, so selection must still find it
  expect_gt(out$n_selected, 0)
  expect_true(grepl("X_age", out$selected_covariates))
})

test_that("the EPV gate suppresses the interacted test when events are thin", {
  set.seed(11)
  n <- 300
  dat <- data.frame(Z = rep(0:1, n / 2))
  for (j in 1:10) dat[[paste0("X_", j)]] <- rnorm(n)
  dat$Y_outcome <- rowSums(dat[, paste0("X_", 1:10)]) + rnorm(n)
  dat$Y_outcome[sample(n, 8)] <- NA          # only 8 events

  out <- check_attrition_lasso(dat, Z, outcomes = "Y_outcome",
                               covariates = paste0("X_", 1:10))
  expect_false(out$epv_adequate)
  expect_true(is.na(out$p_interacted))
  expect_false(out$flag_interacted)
  # the covariate-free test is still reported
  expect_false(is.na(out$p_simple))
})

test_that("epv_threshold and alpha are honoured", {
  dat <- make_attrition_dat()
  strict <- check_attrition_lasso(dat, Z, outcomes = "Y_outcome",
                                  covariates = c("X_age", "X_income"),
                                  epv_threshold = 1e6)
  expect_false(strict$epv_adequate)
  expect_true(is.na(strict$p_interacted))

  # differential attrition is significant at .05 but not at .0001
  expect_true(check_attrition_lasso(dat, Z, outcomes = "Y_outcome",
                                    covariates = c("X_age", "X_income"),
                                    alpha = 0.05)$flag_simple)
  expect_false(check_attrition_lasso(dat, Z, outcomes = "Y_outcome",
                                     covariates = c("X_age", "X_income"),
                                     alpha = 1e-8)$flag_simple)
})

test_that("no covariates runs the covariate-free test only", {
  dat <- make_attrition_dat()
  out <- check_attrition_lasso(dat, Z, outcomes = "Y_outcome")
  expect_equal(out$n_selected, 0)
  expect_false(out$epv_adequate)
  expect_true(is.na(out$p_interacted))
  expect_false(is.na(out$p_simple))
})

test_that("rows with a missing treatment value are dropped before counting", {
  dat <- make_attrition_dat()
  dat_extra <- rbind(dat, transform(dat[1:20, ], Z = NA))
  out <- check_attrition_lasso(dat_extra, Z, outcomes = "Y_outcome",
                               covariates = c("X_age", "X_income"))
  expect_equal(out$n_assigned, nrow(dat))
})

test_that("multi-arm treatments spend more numerator degrees of freedom", {
  dat <- make_attrition_dat()
  dat$Zm <- sample(c("C", "T1", "T2"), nrow(dat), TRUE)
  binary <- check_attrition_lasso(dat, Z, outcomes = "Y_outcome",
                                  covariates = c("X_age", "X_income"))
  multi <- check_attrition_lasso(dat, Zm, outcomes = "Y_outcome",
                                 covariates = c("X_age", "X_income"))
  expect_equal(multi$df1, 2L * binary$df1)
})

test_that("a covariate whose name begins with the treatment name is not swept in", {
  set.seed(12)
  n <- 600
  dat <- data.frame(Z = rep(0:1, n / 2), Zeal = rnorm(n), X_a = rnorm(n))
  dat$Y_outcome <- dat$Zeal + dat$X_a + rnorm(n)
  dat$Y_outcome[rbinom(n, 1, 0.3) == 1] <- NA

  out <- check_attrition_lasso(dat, Z, outcomes = "Y_outcome",
                               covariates = c("Zeal", "X_a"),
                               epv_threshold = 0)
  # df1 must be (arms - 1) * (1 + n_selected), counting Zeal as a covariate and
  # not as a treatment term
  expect_equal(out$df1, as.integer(1 * (1 + out$n_selected)))
})

test_that("study_id is appended when supplied", {
  dat <- make_attrition_dat()
  out <- check_attrition_lasso(dat, Z, outcomes = "Y_outcome",
                               covariates = c("X_age", "X_income"),
                               study_id = "smith_2024_study_1")
  expect_equal(unique(out$study_id), "smith_2024_study_1")
})

test_that("quiet = FALSE prints and returns invisibly", {
  dat <- make_attrition_dat()
  expect_output(
    res <- withVisible(check_attrition_lasso(dat, Z, outcomes = "Y_outcome",
                                             covariates = c("X_age", "X_income"),
                                             quiet = FALSE))
  )
  expect_false(res$visible)
  expect_true(withVisible(check_attrition_lasso(
    dat, Z, outcomes = "Y_outcome", covariates = c("X_age", "X_income")
  ))$visible)
})

test_that("results are reproducible and the caller's RNG is untouched", {
  dat <- make_attrition_dat()
  a <- check_attrition_lasso(dat, Z, outcomes = "Y_outcome",
                             covariates = c("X_age", "X_income"))
  b <- check_attrition_lasso(dat, Z, outcomes = "Y_outcome",
                             covariates = c("X_age", "X_income"))
  expect_equal(a, b)

  set.seed(99)
  before <- rnorm(1)
  set.seed(99)
  invisible(check_attrition_lasso(dat, Z, outcomes = "Y_outcome",
                                  covariates = c("X_age", "X_income")))
  expect_identical(before, rnorm(1))
})

test_that("outcomes default to the Y-prefixed columns", {
  dat <- make_attrition_dat()
  dat$Y_extra <- rnorm(nrow(dat))
  dat$Y_extra[1:30] <- NA
  dat$not_an_outcome <- rnorm(nrow(dat))
  out <- check_attrition_lasso(dat, Z, covariates = c("X_age", "X_income"))
  expect_setequal(out$outcome, c("Y_outcome", "Y_extra"))
})
