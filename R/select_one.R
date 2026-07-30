#' Select covariates that predict one variable
#'
#' Runs a single LASSO regression of one variable on a pool of candidate
#' covariates and returns the covariates that survive. This is the narrow,
#' treatment-free counterpart to \code{\link{lasso_select_covariates}}: it knows
#' nothing about arms, assignment, or estimands, and answers only "which of these
#' candidates predict this column".
#'
#' It exists because some procedures need selection equations that
#' \code{lasso_select_covariates} cannot express. That function takes
#' \code{outcome ~ treatment} and selects arm by arm and on each arm indicator,
#' which is the right thing for regression adjustment of an experiment. A
#' diagnostic for differential attrition, by contrast, needs two equations with
#' two different left-hand sides and no treatment in either: the outcome among
#' respondents who answered, and an indicator for having dropped out among all
#' those assigned. Calling this function twice expresses that directly.
#'
#' @section Relation to the other selector:
#' \code{lasso_select_covariates} returns a one-sided formula, because its answer
#' is consumed by \code{\link[estimatr]{lm_lin}}. This function returns a
#' character vector, because its answer is usually combined with another
#' selection before any model is fit. Wrap the result in
#' \code{\link[stats]{reformulate}} if a formula is what you need.
#'
#' @section Unusable candidates:
#' Candidates that are absent from \code{data}, entirely missing, or constant are
#' dropped before the complete-case filter rather than after, and the check is
#' applied over the subset actually used. Dropping them afterwards is a trap: a
#' single all-\code{NA} candidate empties the complete-case set, and selection
#' then silently returns nothing while appearing to have run. A covariate that is
#' constant among the rows being used is equally unusable even when it varies in
#' the full data, which is why the check follows \code{subset}.
#'
#' Non-numeric candidates are expanded through \code{model.matrix} and mapped
#' back to their source covariate, so a factor is returned whole rather than as a
#' subset of its indicator columns.
#'
#' @param response Column to be predicted, as a character scalar or a bare
#'   column name.
#' @param covariates A one-sided formula of candidate covariates
#'   (e.g. \code{~ X1 + X2}), a character vector of column names, or a bare
#'   expression such as \code{X1 + X2}.
#' @param data A data frame.
#' @param subset Optional logical vector, the same length as \code{nrow(data)},
#'   restricting the rows selection runs on. \code{NA} is treated as
#'   \code{FALSE}. Use it for equations estimated on part of the sample, such as
#'   an outcome equation among respondents who answered.
#' @param lambda_rule Either \code{"lambda.1se"} (default, the one-standard-error
#'   rule) or \code{"lambda.min"} (the cross-validation minimum).
#' @param seed Integer seed for the cross-validation fold assignment, so
#'   selection is reproducible. The global random number state is restored on
#'   exit, so calling this function does not disturb the stream of the script
#'   that called it.
#'
#' @return A character vector of selected covariate names, in the order they were
#'   supplied. Empty when nothing was selected or selection could not be run,
#'   which are not distinguished: both mean "no covariates to use".
#'
#' @seealso \code{\link{lasso_select_covariates}} for the post-double-selection
#'   procedure used to adjust an experiment.
#'
#' @examples
#' set.seed(1)
#' n <- 400
#' dat <- data.frame(X1 = rnorm(n), X2 = rnorm(n), X3 = rnorm(n))
#' dat$Y <- 1.5 * dat$X1 + rnorm(n)
#'
#' # X1 predicts Y; X2 and X3 are noise
#' lasso_select_one("Y", ~ X1 + X2 + X3, data = dat)
#'
#' # Selection restricted to half the sample
#' lasso_select_one("Y", ~ X1 + X2 + X3, data = dat,
#'                  subset = seq_len(n) <= n / 2)
#'
#' @importFrom stats complete.cases model.matrix reformulate terms setNames
#' @family covariate selection
#' @export
lasso_select_one <- function(response, covariates, data, subset = NULL,
                             lambda_rule = c("lambda.1se", "lambda.min"),
                             seed = 999) {
  lambda_rule <- match.arg(lambda_rule)

  # Accept a bare column name, a string, or a variable holding a string. Test the
  # unevaluated expression first: a bare name has no binding outside `data`, so
  # evaluating it up front would error.
  response_expr <- substitute(response)
  response_name <- if (is.character(response_expr)) {
    response_expr
  } else {
    evaluated <- tryCatch(response, error = function(e) NULL)
    if (is.character(evaluated) && length(evaluated) == 1L) {
      evaluated
    } else {
      deparse1(response_expr)
    }
  }

  covariate_cols <- resolve_covariate_names(substitute(covariates), covariates)

  if (!response_name %in% names(data)) return(character(0))

  # Apply the row restriction first, so every check below reflects the rows
  # selection will actually see.
  if (!is.null(subset)) {
    if (length(subset) != nrow(data)) {
      stop("lasso_select_one: subset must have length nrow(data).", call. = FALSE)
    }
    subset[is.na(subset)] <- FALSE
    data <- data[subset, , drop = FALSE]
  }
  if (nrow(data) == 0) return(character(0))

  # Drop unusable candidates BEFORE complete.cases, and judge usability on the
  # subset rather than the full data: a covariate constant among the rows in use
  # carries no information here even if it varies elsewhere.
  covariate_cols <- covariate_cols[covariate_cols %in% names(data)]
  usable <- vapply(covariate_cols, function(cl) {
    x <- data[[cl]]
    !all(is.na(x)) && length(unique(x[!is.na(x)])) > 1
  }, logical(1))
  covariate_cols <- covariate_cols[usable]
  if (length(covariate_cols) == 0) return(character(0))

  df <- data[stats::complete.cases(data[, c(response_name, covariate_cols), drop = FALSE]), ,
             drop = FALSE]
  if (nrow(df) == 0) return(character(0))

  mm_formula <- stats::reformulate(covariate_cols)
  mm <- tryCatch(stats::model.matrix(mm_formula, data = df), error = function(e) NULL)
  if (is.null(mm) || ncol(mm) < 2) return(character(0))

  x_mat <- mm[, -1, drop = FALSE]

  # Map each model-matrix column to the covariate that produced it via the
  # `assign` attribute, which records the mapping exactly. Name matching cannot
  # recover it: a factor expands to columns named <covariate><level> with no
  # separator, so any prefix test guesses.
  mm_assign <- attr(mm, "assign")[-1]
  mm_terms <- attr(stats::terms(mm_formula, data = df), "term.labels")
  x_col_to_original <- stats::setNames(mm_terms[mm_assign], colnames(x_mat))

  restore_seed <- capture_seed()
  on.exit(restore_seed(), add = TRUE)
  set.seed(seed)

  selected <- lasso_nonzero(scale(x_mat), df[[response_name]], lambda_rule)
  if (length(selected) == 0) return(character(0))

  selected_sources <- unique(unname(x_col_to_original[selected]))
  selected_sources <- selected_sources[!is.na(selected_sources)]

  covariate_cols[covariate_cols %in% selected_sources]
}
