#' Bootstrap an lm_robust fit
#'
#' Refits \code{\link[estimatr]{lm_robust}} on \code{times} resamples of the
#' data and returns the fit on the original data alongside the coefficients from
#' every resample. Summarize it with \code{\link[broom]{tidy}}, which reports the
#' bootstrap standard deviation as a standard error and percentile confidence
#' bounds.
#'
#' Each replicate is fitted with \code{se_type = "none"}, since the replicate's
#' own standard errors are never used: the spread across replicates is the
#' standard error.
#'
#' @section What gets resampled:
#' Rows, independently, \code{nrow(data)} of them with replacement. That is the
#' right unit when observations are independent and the wrong one when they are
#' not: with clustered data the row bootstrap breaks the within-cluster
#' dependence and returns intervals that are too narrow. This function does not
#' detect that case, so for clustered designs either use \code{lm_robust} with
#' \code{clusters =} directly, or resample clusters before calling it.
#'
#' @section What the point estimate is:
#' The coefficient from the fit on the original data, not the mean of the
#' bootstrap draws. Only the uncertainty comes from the resamples. The interval
#' is the plain percentile interval, with no bias correction or acceleration, so
#' it is not centered on the point estimate when the bootstrap distribution is
#' skewed.
#'
#' @section Replicates that cannot be fitted:
#' A resample can be rank deficient even when the original data is not, most
#' often when a factor level or a rare treatment arm is drawn zero times. Such a
#' replicate contributes \code{NA} rather than being silently redrawn, and
#' \code{lm_bootstrap} warns with the count so a summary computed from fewer
#' replicates than requested is visible rather than hidden.
#'
#' @param formula A formula passed to \code{\link[estimatr]{lm_robust}}.
#' @param data A data frame.
#' @param times Number of bootstrap replicates.
#' @param ... Further arguments for \code{\link[estimatr]{lm_robust}}, such as
#'   \code{weights} or \code{fixed_effects}. Arguments naming columns are
#'   evaluated inside each resample, so a weight travels with the rows that were
#'   drawn.
#'
#' @return An object of class \code{lm_bootstrap}: a list with
#'   \code{original_mod}, the \code{lm_robust} fit on the original data, and
#'   \code{boot_results}, a tibble with one row per replicate per term
#'   (\code{replicate}, \code{term}, \code{estimate}).
#'
#' @examples
#' set.seed(1)
#' n <- 200
#' dat <- data.frame(Z = rep(0:1, n / 2), X = rnorm(n))
#' dat$Y <- 0.5 * dat$Z + dat$X + rnorm(n)
#'
#' boot_fit <- lm_bootstrap(Y ~ Z + X, data = dat, times = 200)
#' broom::tidy(boot_fit)
#'
#' @seealso \code{\link{tidy.lm_bootstrap}} for the summary.
#' @importFrom stats coef complete.cases setNames
#' @export
lm_bootstrap <- function(formula, data, times = 1000, ...) {
  original_mod <- estimatr::lm_robust(
    formula = formula, data = data, se_type = "none", ...
  )
  term_levels <- broom::tidy(original_mod)$term

  n <- nrow(data)
  draws <- lapply(seq_len(times), function(b) {
    resample <- data[sample.int(n, n, replace = TRUE), , drop = FALSE]
    fit <- tryCatch(
      estimatr::lm_robust(
        formula = formula, data = resample, se_type = "none", ...
      ),
      error = function(e) NULL
    )
    if (is.null(fit)) {
      stats::setNames(rep(NA_real_, length(term_levels)), term_levels)
    } else {
      # Index by name: a rank-deficient resample can drop a term entirely, and
      # positional indexing would then shift every later coefficient onto the
      # wrong term.
      stats::setNames(unname(stats::coef(fit)[term_levels]), term_levels)
    }
  })

  boot_matrix <- do.call(rbind, draws)
  n_failed <- sum(!stats::complete.cases(boot_matrix))
  if (n_failed > 0) {
    warning(
      "lm_bootstrap: ", n_failed, " of ", times,
      " replicates could not be fitted for every term and contribute NA.",
      call. = FALSE
    )
  }

  boot_results <- tibble::tibble(
    replicate = rep(seq_len(times), each = length(term_levels)),
    term = factor(rep(term_levels, times = times), levels = term_levels),
    # Row-major, so each replicate's terms stay together in the order above.
    estimate = as.vector(t(boot_matrix))
  )

  structure(
    list(original_mod = original_mod, boot_results = boot_results),
    class = "lm_bootstrap"
  )
}

#' Summarize a bootstrapped lm_robust fit
#'
#' Reduces the replicates held by \code{\link{lm_bootstrap}} to one row per
#' term: the estimate from the original fit, the standard deviation across
#' replicates as its standard error, and percentile confidence bounds.
#'
#' Replicates that could not be fitted are excluded term by term, so a term that
#' survived every resample uses every replicate even when another term did not.
#' \code{lm_bootstrap} has already warned about the count.
#'
#' @param x An object of class \code{lm_bootstrap}.
#' @param alpha One minus the confidence level, so \code{0.05} gives a 95
#'   percent interval.
#' @param ... Ignored.
#'
#' @return A tibble with columns \code{term}, \code{estimate},
#'   \code{std.error}, \code{conf.low}, and \code{conf.high}, in the term order
#'   of the original fit.
#'
#' @examples
#' set.seed(1)
#' n <- 200
#' dat <- data.frame(Z = rep(0:1, n / 2), X = rnorm(n))
#' dat$Y <- 0.5 * dat$Z + dat$X + rnorm(n)
#'
#' boot_fit <- lm_bootstrap(Y ~ Z + X, data = dat, times = 200)
#' broom::tidy(boot_fit)
#' broom::tidy(boot_fit, alpha = 0.10)
#'
#' @importFrom broom tidy
#' @importFrom stats quantile sd
#' @export
tidy.lm_bootstrap <- function(x, alpha = 0.05, ...) {
  point <- broom::tidy(x$original_mod)
  by_term <- split(x$boot_results$estimate, x$boot_results$term)
  usable <- lapply(by_term, function(e) e[is.finite(e)])

  tibble::tibble(
    term = point$term,
    estimate = point$estimate,
    std.error = vapply(usable, stats::sd, numeric(1)),
    conf.low = vapply(
      usable, stats::quantile, numeric(1), probs = alpha / 2, names = FALSE
    ),
    conf.high = vapply(
      usable, stats::quantile, numeric(1), probs = 1 - alpha / 2, names = FALSE
    )
  )
}
