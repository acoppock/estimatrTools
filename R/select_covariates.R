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
#' select_covariates(Y ~ Z, ~ X1 + X2 + X3, data = dat)
#'
#' @importFrom stats as.formula complete.cases model.matrix reformulate var coef
#' @family covariate selection
#' @export
select_covariates <- function(formula, covariates, data,
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

  X <- tryCatch(
    stats::model.matrix(stats::reformulate(covariate_cols), data = df)[, -1, drop = FALSE],
    error = function(e) NULL
  )
  if (is.null(X) || ncol(X) == 0) return(~1)
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

  # Map each selected model-matrix column back to its source covariate. Factors
  # expand to columns named <covariate><level> with no separator (e.g.
  # "X_region3", "X_supporterRegime supporter"), so a bare startsWith
  # over-includes whenever one covariate name is a string prefix of another:
  # both would match the longer covariate's indicator columns. Attribute each
  # selected name to its LONGEST matching covariate name instead. Identical to
  # the bare-prefix test whenever no covariate name prefixes another, and
  # correct when one does.
  source_of <- function(nm) {
    hits <- covariate_cols[startsWith(nm, covariate_cols)]
    if (length(hits) == 0) NA_character_ else hits[which.max(nchar(hits))]
  }
  selected_sources <- unique(vapply(selected, source_of, character(1)))
  original_cols <- covariate_cols[covariate_cols %in% selected_sources]

  if (length(original_cols) == 0) return(~1)
  stats::reformulate(original_cols)
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
