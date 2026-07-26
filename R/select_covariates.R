#' Select covariates by post-double-selection LASSO
#'
#' Chooses a covariate set for regression adjustment of an experiment by
#' running two families of LASSO regressions and taking their union:
#'
#' \itemize{
#'   \item **Outcome equations**: the outcome on all candidate covariates,
#'     separately within each treatment arm. Running these arm-by-arm rather
#'     than pooled is what makes the selection appropriate for the Lin (2013)
#'     estimator, whose adjustment is a separate slope per arm: a covariate
#'     that matters only for the treatment interaction is picked up here and
#'     would be missed by a pooled regression.
#'   \item **Assignment equations**: each arm indicator on all candidate
#'     covariates. This is the double-selection step of Belloni,
#'     Chernozhukov and Hansen (2014). It retains covariates that predict
#'     assignment even when they barely predict the outcome, which is what
#'     protects the estimate against omitted-variable bias from a covariate
#'     the outcome equations would have dropped.
#' }
#'
#' Covariates are standardized before selection so the LASSO penalty applies
#' on a common scale, and selected columns of the model matrix are mapped back
#' to the original covariate names, so a factor is retained whole rather than
#' as a subset of its indicator columns.
#'
#' @section Choice of penalty:
#' The default \code{lambda_rule = "lambda.1se"} is the one-standard-error
#' rule: the most parsimonious model within one standard error of the
#' cross-validation minimum. Parsimony matters more here than in plain OLS
#' because under \code{\link{lm_lin_lasso}} each retained covariate costs one
#' parameter *per arm*, so over-selection inflates standard errors faster than
#' it would in an additive model. Pass \code{lambda_rule = "lambda.min"} for
#' the cross-validation minimum instead.
#'
#' @section Rank deficiency in the interacted design:
#' LASSO selects against the outcome and knows nothing about the design that
#' will consume its answer. Under \code{\link{lm_lin_lasso}} every selected
#' covariate is interacted with every treatment arm, so a set that is fine in a
#' pooled regression can be rank deficient once crossed with treatment.
#'
#' When that happens this function \strong{warns and returns the selection
#' unchanged}. It does not prune. Pruning would decide, silently and inside the
#' algorithm, something that belongs in the analysis script where a reader can
#' see it, and the right remedy is not the same in every case:
#' \itemize{
#'   \item \strong{Redundancy}: one covariate lies in the span of another, for
#'     instance a college indicator derived from an education factor. Dropping
#'     one of them at the call site is right, and the warning names which.
#'   \item \strong{Empty cell}: a single factor level is unobserved in one arm,
#'     so one interaction column is all zero. Dropping the whole covariate to
#'     fix this is a poor trade, since a factor contributing fourteen columns
#'     loses thirteen good ones to remove one aliased one. Collapsing the sparse
#'     level in the cleaning script keeps the rest.
#' }
#' Left alone, \code{lm_lin} drops the aliased column itself and still returns
#' an estimate, so the warning is a prompt to make a decision rather than a
#' failure.
#'
#' @section Unusable candidates:
#' Candidate columns that are absent from \code{data}, entirely missing, or
#' constant are dropped before the complete-case filter rather than after.
#' Dropping them afterwards is a trap: a single all-\code{NA} candidate empties
#' the complete-case set, and selection then silently returns nothing while
#' appearing to have run.
#'
#' @param formula A two-sided formula \code{outcome ~ treatment} giving the
#'   focal estimand. Only these two variables are read from it.
#' @param covariates A one-sided formula of candidate covariates
#'   (e.g. \code{~ X1 + X2}), a character vector of column names, or a bare
#'   expression such as \code{X1 + X2}.
#' @param data A data frame.
#' @param lambda_rule Either \code{"lambda.1se"} (default) or
#'   \code{"lambda.min"}. See the penalty section above.
#' @param seed Integer seed for the cross-validation fold assignment, so that
#'   selection is reproducible. The global random number state is restored on
#'   exit, so calling this function does not disturb the stream of the script
#'   that called it.
#'
#' @return A one-sided formula naming the selected original covariates, or
#'   \code{~1} when nothing was selected or selection could not be run.
#'   \code{~1} is a valid input to \code{\link[estimatr]{lm_robust}} but not to
#'   \code{\link[estimatr]{lm_lin}}; the fitting functions in this package
#'   handle that difference for you.
#'
#' @references
#' Belloni, A., Chernozhukov, V., and Hansen, C. (2014). Inference on treatment
#' effects after selection among high-dimensional controls. \emph{Review of
#' Economic Studies}, 81(2), 608-650. \doi{10.1093/restud/rdt044}
#'
#' Lin, W. (2013). Agnostic notes on regression adjustments to experimental
#' data: reexamining Freedman's critique. \emph{Annals of Applied Statistics},
#' 7(1), 295-318. \doi{10.1214/12-AOAS583}
#'
#' @examples
#' set.seed(1)
#' n <- 400
#' dat <- data.frame(
#'   Z  = rep(0:1, n / 2),
#'   X1 = rnorm(n),
#'   X2 = rnorm(n),
#'   X3 = rnorm(n)
#' )
#' dat$Y <- 0.5 * dat$Z + 1.5 * dat$X1 + rnorm(n)
#'
#' # X1 predicts the outcome; X2 and X3 are noise
#' lasso_select_covariates(Y ~ Z, ~ X1 + X2 + X3, data = dat)
#'
#' @importFrom stats as.formula complete.cases model.matrix reformulate var coef terms setNames
#' @family covariate selection
#' @export
lasso_select_covariates <- function(formula, covariates, data,
                                    lambda_rule = c("lambda.1se", "lambda.min"),
                                    seed = 999) {
  lambda_rule <- match.arg(lambda_rule)

  covariate_cols <- resolve_covariate_names(substitute(covariates), covariates)

  outcome   <- as.character(formula[[2]])
  treatment <- as.character(formula[[3]])

  # Drop unusable candidates BEFORE complete.cases: an absent, all-NA, or
  # constant column would otherwise empty the complete-case set and silently
  # disable selection entirely.
  covariate_cols <- covariate_cols[covariate_cols %in% names(data)]
  usable <- vapply(covariate_cols, function(cl) {
    x <- data[[cl]]
    !all(is.na(x)) && length(unique(x[!is.na(x)])) > 1
  }, logical(1))
  covariate_cols <- covariate_cols[usable]

  if (length(covariate_cols) == 0) return(~1)
  if (!all(c(outcome, treatment) %in% names(data))) return(~1)

  df <- data[stats::complete.cases(data[, c(outcome, treatment, covariate_cols), drop = FALSE]), ,
             drop = FALSE]
  if (nrow(df) == 0) return(~1)

  y <- df[[outcome]]
  z <- df[[treatment]]
  if (length(unique(z)) < 2) return(~1)

  mm_formula <- stats::reformulate(covariate_cols)
  mm <- tryCatch(stats::model.matrix(mm_formula, data = df), error = function(e) NULL)
  if (is.null(mm) || ncol(mm) < 2) return(~1)

  X <- mm[, -1, drop = FALSE]

  # Map each model-matrix column to the covariate that produced it, using the
  # `assign` attribute rather than matching on names. R records the mapping
  # exactly: assign[j] is the index of the term that generated column j. Name
  # matching cannot recover this reliably, because a factor expands to columns
  # named <covariate><level> with no separator, so any prefix test guesses.
  mm_assign <- attr(mm, "assign")[-1]
  mm_terms <- attr(stats::terms(mm_formula, data = df), "term.labels")
  x_col_to_original <- stats::setNames(mm_terms[mm_assign], colnames(X))

  X <- scale(X)

  restore_seed <- capture_seed()
  on.exit(restore_seed(), add = TRUE)
  set.seed(seed)

  arms <- unique(z)

  selected_y <- unique(as.character(unlist(lapply(arms, function(arm) {
    lasso_nonzero(X[z == arm, , drop = FALSE], y[z == arm], lambda_rule)
  }))))

  selected_z <- unique(as.character(unlist(lapply(arms, function(arm) {
    lasso_nonzero(X, as.integer(z == arm), lambda_rule)
  }))))

  selected <- union(selected_y, selected_z)
  if (length(selected) == 0) return(~1)

  # Attribute selected columns to source covariates via the exact `assign`
  # mapping built above, preserving the original candidate order.
  selected_sources <- unique(unname(x_col_to_original[selected]))
  selected_sources <- selected_sources[!is.na(selected_sources)]
  original_cols <- covariate_cols[covariate_cols %in% selected_sources]

  if (length(original_cols) == 0) return(~1)

  warn_rank_deficient(original_cols, treatment, df)

  stats::reformulate(original_cols)
}


#' Warn when the interacted Lin design is rank deficient, and say why
#'
#' Distinguishes the two causes, because they call for different fixes. If
#' removing a single covariate restores full rank, that covariate is redundant
#' given the others and dropping it at the call site is the right response; the
#' narrowest such covariate is suggested, since a derived variable (a college
#' indicator) contributes fewer columns than the factor it came from. If no
#' single removal restores rank, the cause is a sparse cell rather than a
#' redundant covariate, and the empty level-by-arm combinations are named
#' instead.
#'
#' @param cols Selected covariate names.
#' @param treatment Name of the treatment variable.
#' @param df The complete-case data frame selection was run on.
#' @return Nothing; called for the warning.
#' @keywords internal
#' @noRd
warn_rank_deficient <- function(cols, treatment, df) {
  design <- function(v) {
    if (length(v) == 0) return(NULL)
    rhs <- paste(c(treatment, v, paste0(treatment, ":", v)), collapse = " + ")
    mm <- tryCatch(stats::model.matrix(stats::as.formula(paste("~", rhs)), data = df),
                   error = function(e) NULL)
    if (is.null(mm)) return(NULL)
    mm[stats::complete.cases(mm), , drop = FALSE]
  }

  mm <- design(cols)
  if (is.null(mm) || nrow(mm) == 0) return(invisible(NULL))
  full <- qr(mm)$rank
  if (full == ncol(mm)) return(invisible(NULL))

  deficiency <- ncol(mm) - full

  header <- sprintf(
    "lasso_select_covariates: the interacted Lin design is rank deficient (%d columns, rank %d).",
    ncol(mm), full
  )

  # Empty cells are checked FIRST. Removing the covariate that contains an
  # unobserved level always restores rank, so a "does removing it help" test
  # cannot tell the two causes apart and would report every empty cell as a
  # redundancy. Look for the empty cell directly instead.
  empties <- character(0)
  culprits <- character(0)
  for (v in cols) {
    x <- df[[v]]
    if (is.numeric(x)) next
    tb <- table(as.character(x), as.character(df[[treatment]]))
    z0 <- which(tb == 0, arr.ind = TRUE)
    if (nrow(z0) > 0) {
      culprits <- c(culprits, v)
      empties <- c(empties, sprintf('    %s = "%s" has no observations in %s = %s',
                                    v, rownames(tb)[z0[, 1]], treatment, colnames(tb)[z0[, 2]]))
    }
  }

  if (length(empties) > 0) {
    n_cols <- vapply(culprits, function(v) {
      m <- tryCatch(stats::model.matrix(stats::reformulate(v), data = df), error = function(e) NULL)
      if (is.null(m)) NA_integer_ else 2L * (ncol(m) - 1L)
    }, integer(1))
    body <- paste0(
      "  Cause: an unobserved level, not a redundant covariate.\n",
      paste(empties, collapse = "\n"), "\n",
      sprintf("  %s contributes %s columns to this design, so dropping it would discard %s to remove %d.\n",
              paste(culprits, collapse = " and "), paste(n_cols, collapse = "/"),
              paste(n_cols - deficiency, collapse = "/"), deficiency),
      "  Collapsing the sparse level in the cleaning script keeps the rest of the covariate."
    )
  } else {
    restorers <- Filter(function(v) {
      m <- design(setdiff(cols, v))
      !is.null(m) && nrow(m) > 0 && qr(m)$rank == ncol(m)
    }, cols)

    if (length(restorers) > 0) {
      width <- vapply(restorers, function(v) {
        m <- tryCatch(stats::model.matrix(stats::reformulate(v), data = df), error = function(e) NULL)
        if (is.null(m)) NA_integer_ else ncol(m) - 1L
      }, integer(1))
      suggest <- restorers[which.min(width)]
      keep <- setdiff(cols, suggest)
      body <- paste0(
        "  Cause: redundancy. ", suggest, " lies in the span of the others.\n",
        "  Exclude it at the call site so the choice is visible in the script:\n",
        "       covariates = c(", paste0('"', keep, '"', collapse = ", "), ")"
      )
    } else {
      body <- paste0(
        "  Cause: not attributable to a single covariate or an unobserved level; ",
        deficiency, " columns are collinear across the set.\n",
        "  Inspect the design before choosing what to exclude at the call site."
      )
    }
  }

  warning(paste0(header, "\n", body,
                 "\n  Left as is, lm_lin drops the aliased column and still returns an estimate."),
          call. = FALSE)
  invisible(NULL)
}


#' Names of the columns with non-zero LASSO coefficients
#'
#' @param x_mat Numeric design matrix.
#' @param response Numeric response.
#' @param lambda_rule Penalty selection rule.
#' @return Character vector of selected column names, possibly empty.
#' @keywords internal
#' @noRd
lasso_nonzero <- function(x_mat, response, lambda_rule) {
  if (nrow(x_mat) < 10 || ncol(x_mat) == 0) return(character(0))

  # scale() emits NaN for a zero-variance column, and var() then returns NA,
  # so coerce NA to FALSE rather than letting it propagate into the subscript.
  keep <- apply(x_mat, 2, function(x) stats::var(x, na.rm = TRUE) > 0)
  keep[is.na(keep)] <- FALSE
  x_mat <- x_mat[, keep, drop = FALSE]
  if (ncol(x_mat) == 0) return(character(0))

  ok <- stats::complete.cases(x_mat)
  x_mat <- x_mat[ok, , drop = FALSE]
  response <- response[ok]
  if (nrow(x_mat) < 10) return(character(0))
  if (isTRUE(stats::var(response, na.rm = TRUE) == 0)) return(character(0))

  # cv.glmnet needs at least two columns; with exactly one there is nothing to
  # select between, so keep it.
  if (ncol(x_mat) == 1) return(colnames(x_mat))

  fit <- tryCatch(
    glmnet::cv.glmnet(x_mat, response, alpha = 1, family = "gaussian"),
    error = function(e) NULL
  )
  if (is.null(fit)) return(character(0))

  coeffs <- as.matrix(stats::coef(fit, s = lambda_rule))
  rownames(coeffs)[coeffs[, 1] != 0 & rownames(coeffs) != "(Intercept)"]
}


#' Resolve a covariate specification to column names
#'
#' Accepts a one-sided formula, a character vector, or a bare expression.
#'
#' @param cov_expr The unevaluated \code{covariates} argument.
#' @param covariates The evaluated \code{covariates} argument.
#' @return Character vector of covariate names.
#' @keywords internal
#' @noRd
resolve_covariate_names <- function(cov_expr, covariates) {
  evaluated <- tryCatch(covariates, error = function(e) NULL)

  if (inherits(evaluated, "formula")) return(all.vars(evaluated))
  if (is.character(evaluated)) return(evaluated)
  if (inherits(cov_expr, "formula")) return(all.vars(cov_expr))
  all.vars(stats::as.formula(paste("~", deparse1(cov_expr))))
}


#' Capture and restore the global random number state
#'
#' Selection needs a fixed seed so that cross-validation fold assignment is
#' reproducible, but setting one should not silently reposition the caller's
#' random number stream.
#'
#' @return A function that restores the state captured at call time.
#' @keywords internal
#' @noRd
capture_seed <- function() {
  if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
    old <- get(".Random.seed", envir = globalenv(), inherits = FALSE)
    function() assign(".Random.seed", old, envir = globalenv())
  } else {
    function() {
      if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
        rm(".Random.seed", envir = globalenv())
      }
    }
  }
}
