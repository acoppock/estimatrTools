make_boot_data <- function(n = 200) {
  dat <- data.frame(Z = rep(0:1, n / 2), X = rnorm(n), w = runif(n, 0.5, 2))
  dat$Y <- 0.5 * dat$Z + dat$X + rnorm(n)
  dat
}

test_that("lm_bootstrap returns the original fit and one row per replicate per term", {
  set.seed(123)
  dat <- make_boot_data()
  boot_fit <- lm_bootstrap(Y ~ Z + X, data = dat, times = 50)

  expect_s3_class(boot_fit, "lm_bootstrap")
  expect_s3_class(boot_fit$original_mod, "lm_robust")
  expect_equal(nrow(boot_fit$boot_results), 50 * 3)
  expect_equal(
    levels(boot_fit$boot_results$term),
    c("(Intercept)", "Z", "X")
  )
  expect_equal(max(boot_fit$boot_results$replicate), 50)
  expect_false(any(is.na(boot_fit$boot_results$estimate)))
})

test_that("replicates differ from each other and from the original fit", {
  set.seed(123)
  dat <- make_boot_data()
  boot_fit <- lm_bootstrap(Y ~ Z + X, data = dat, times = 50)
  z_draws <- boot_fit$boot_results$estimate[boot_fit$boot_results$term == "Z"]

  expect_gt(length(unique(z_draws)), 1)
  expect_gt(stats::sd(z_draws), 0)
})

test_that("tidy reports the original estimate with bootstrap uncertainty", {
  set.seed(123)
  dat <- make_boot_data()
  boot_fit <- lm_bootstrap(Y ~ Z + X, data = dat, times = 200)
  tidied <- broom::tidy(boot_fit)

  expect_named(
    tidied, c("term", "estimate", "std.error", "conf.low", "conf.high")
  )
  expect_equal(tidied$term, c("(Intercept)", "Z", "X"))
  # The point estimate is the fit on the original data, not the bootstrap mean.
  expect_equal(tidied$estimate, unname(coef(boot_fit$original_mod)))
  expect_true(all(tidied$conf.low < tidied$conf.high))

  # A row bootstrap of an independent sample should land near the HC2 answer.
  hc2 <- estimatr::lm_robust(Y ~ Z + X, data = dat)
  expect_equal(tidied$std.error, hc2$std.error, tolerance = 0.2)
})

test_that("alpha widens the interval", {
  set.seed(123)
  dat <- make_boot_data()
  boot_fit <- lm_bootstrap(Y ~ Z + X, data = dat, times = 200)

  wide <- broom::tidy(boot_fit, alpha = 0.01)
  narrow <- broom::tidy(boot_fit, alpha = 0.10)
  expect_true(all(wide$conf.high - wide$conf.low > narrow$conf.high - narrow$conf.low))
})

test_that("arguments naming columns are evaluated inside each resample", {
  set.seed(123)
  dat <- make_boot_data()
  boot_fit <- lm_bootstrap(Y ~ Z + X, data = dat, times = 20, weights = w)

  expect_equal(
    coef(boot_fit$original_mod),
    coef(estimatr::lm_robust(Y ~ Z + X, data = dat, weights = w, se_type = "none"))
  )
  expect_false(any(is.na(boot_fit$boot_results$estimate)))

  # Weights travelling with the drawn rows is the point: a replicate paired with
  # the original ordering of `w` would give different draws.
  set.seed(123)
  unweighted <- lm_bootstrap(Y ~ Z + X, data = dat, times = 20)
  weighted_z <- boot_fit$boot_results$estimate[boot_fit$boot_results$term == "Z"]
  unweighted_z <- unweighted$boot_results$estimate[unweighted$boot_results$term == "Z"]
  expect_false(isTRUE(all.equal(weighted_z, unweighted_z)))
})

test_that("a replicate that cannot be fitted contributes NA and is reported", {
  set.seed(2)
  # One treated unit means many resamples draw none of it, so Z drops out.
  dat <- data.frame(Z = c(1, rep(0, 29)), X = rnorm(30))
  dat$Y <- rnorm(30)

  expect_warning(
    boot_fit <- lm_bootstrap(Y ~ Z + X, data = dat, times = 50),
    "replicates could not be fitted"
  )
  expect_true(any(is.na(boot_fit$boot_results$estimate)))

  # The summary still comes back, computed from the replicates that survived.
  tidied <- broom::tidy(boot_fit)
  expect_equal(nrow(tidied), 3)
  expect_true(all(is.finite(tidied$std.error)))
})
