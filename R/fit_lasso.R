#' Lin estimator on LASSO-selected covariates
#'
#' Runs \code{\link{lasso_select_covariates}} and fits
#' \code{\link[estimatr]{lm_lin}} on the selected set: a separate covariate
#' slope in each treatment arm, with covariates centered at their full-sample
#' means so the treatment coefficient remains an estimate of the average
#' treatment effect.
#'
#' @section Fallback:
#' \code{lm_lin} cannot be fitted with no covariates, and an adjusted fit can
#' come back degenerate when the selected set is nearly collinear within an
#' arm. In either case this function returns an unadjusted
#' \code{\link[estimatr]{lm_robust}} fit, which for a binary treatment is the
#' difference in means. Because that substitution changes the specification,
#' it is recorded rather than hidden: read it back with
#' \code{\link{adjustment}} and \code{\link{selected_covariates}}. The specific
#' triggers are that selection returned \code{~1}, that \code{lm_lin} threw, or
#' that the treatment row of the adjusted fit has a non-finite standard error.
#'
#' @param formula A two-sided formula \code{outcome ~ treatment}.
#' @param covariates Candidate covariates: a one-sided formula, a character
#'   vector of column names, or a bare expression.
#' @param data A data frame.
#' @param weights,subset,clusters,se_type,ci,alpha,return_vcov,try_cholesky
#'   Passed through to \code{\link[estimatr]{lm_lin}} or
#'   \code{\link[estimatr]{lm_robust}}.
#' @param lasso_args A named list of further arguments for
#'   \code{\link{lasso_select_covariates}}, such as \code{lambda_rule} or
#'   \code{seed}.
#'
#' @return An \code{lm_robust} object, carrying the attributes described in
#'   \code{\link{adjustment}}.
#'
#' @examples
#' set.seed(1)
#' n <- 400
#' dat <- data.frame(Z = rep(0:1, n / 2), X1 = rnorm(n), X2 = rnorm(n))
#' dat$Y <- 0.5 * dat$Z + 1.5 * dat$X1 + rnorm(n)
#'
#' fit <- lm_lin_lasso(Y ~ Z, ~ X1 + X2, data = dat)
#' adjustment(fit)
#' selected_covariates(fit)
#'
#' @seealso \code{\link{lm_robust_lasso}} for additive adjustment,
#'   \code{\link{lm_moderator_lasso}} for a treatment-by-moderator model.
#' @family adjusted estimators
#' @export
lm_lin_lasso <- function(formula, covariates, data,
                         weights = NULL, subset = NULL, clusters = NULL,
                         se_type = NULL, ci = TRUE, alpha = 0.05,
                         return_vcov = TRUE, try_cholesky = FALSE,
                         lasso_args = list()) {
  cov_names <- resolve_covariate_names(substitute(covariates), covariates)

  selected <- do.call(
    lasso_select_covariates,
    c(list(formula = formula, covariates = cov_names, data = data), lasso_args)
  )
  sel_vars <- all.vars(selected)

  nse_args <- list(weights = substitute(weights), subset = substitute(subset),
                   clusters = substitute(clusters))
  plain_args <- list(formula = formula, data = data, se_type = se_type,
                     ci = ci, alpha = alpha, return_vcov = return_vcov,
                     try_cholesky = try_cholesky)
  envir <- parent.frame()

  fallback <- function(reason) {
    fit <- call_estimatr(estimatr::lm_robust, nse_args, plain_args, envir)
    tag_fit(fit, adjustment = "none", selected = character(0), fallback = reason,
            fn = "lm_lin_lasso", formula = formula)
  }

  if (length(sel_vars) == 0) return(fallback("no covariates selected"))

  fit_adj <- tryCatch(
    call_estimatr(estimatr::lm_lin, nse_args,
                  c(plain_args, list(covariates = selected)), envir),
    error = function(e) NULL
  )
  if (is.null(fit_adj)) return(fallback("lm_lin failed to fit"))

  if (!finite_treatment_se(fit_adj, as.character(formula[[3]]))) {
    return(fallback("adjusted fit had a non-finite treatment standard error"))
  }

  tag_fit(fit_adj, adjustment = "lin", selected = sel_vars, fallback = NA_character_,
          fn = "lm_lin_lasso", formula = formula)
}


#' Additive robust regression on LASSO-selected covariates
#'
#' Runs \code{\link{lasso_select_covariates}} and fits
#' \code{\link[estimatr]{lm_robust}} with the selected covariates entered
#' additively, that is \code{outcome ~ treatment + X1 + X2}.
#'
#' Additive adjustment assumes a common covariate slope across arms. That
#' assumption buys degrees of freedom relative to \code{\link{lm_lin_lasso}},
#' which fits a slope per arm, and it is what makes this the right choice when
#' arms are small enough that the interacted model is unstable. When it is
#' wrong, the treatment coefficient is no longer guaranteed to be consistent
#' for the average treatment effect, which is the reason Lin (2013) recommends
#' the interacted form as the default. Choose deliberately.
#'
#' @inheritSection lm_lin_lasso Fallback
#'
#' @inheritParams lm_lin_lasso
#'
#' @return An \code{lm_robust} object, carrying the attributes described in
#'   \code{\link{adjustment}}.
#'
#' @examples
#' set.seed(1)
#' n <- 400
#' dat <- data.frame(Z = rep(0:1, n / 2), X1 = rnorm(n), X2 = rnorm(n))
#' dat$Y <- 0.5 * dat$Z + 1.5 * dat$X1 + rnorm(n)
#'
#' fit <- lm_robust_lasso(Y ~ Z, ~ X1 + X2, data = dat)
#' adjustment(fit)
#'
#' @references
#' Lin, W. (2013). Agnostic notes on regression adjustments to experimental
#' data: reexamining Freedman's critique. \emph{Annals of Applied Statistics},
#' 7(1), 295-318. \doi{10.1214/12-AOAS583}
#'
#' @seealso \code{\link{lm_lin_lasso}} for arm-specific slopes.
#' @importFrom stats as.formula
#' @family adjusted estimators
#' @export
lm_robust_lasso <- function(formula, covariates, data,
                            weights = NULL, subset = NULL, clusters = NULL,
                            se_type = NULL, ci = TRUE, alpha = 0.05,
                            return_vcov = TRUE, try_cholesky = FALSE,
                            lasso_args = list()) {
  cov_names <- resolve_covariate_names(substitute(covariates), covariates)

  selected <- do.call(
    lasso_select_covariates,
    c(list(formula = formula, covariates = cov_names, data = data), lasso_args)
  )
  sel_vars <- all.vars(selected)

  outcome   <- as.character(formula[[2]])
  treatment <- as.character(formula[[3]])

  nse_args <- list(weights = substitute(weights), subset = substitute(subset),
                   clusters = substitute(clusters))
  plain_args <- list(data = data, se_type = se_type, ci = ci, alpha = alpha,
                     return_vcov = return_vcov, try_cholesky = try_cholesky)
  envir <- parent.frame()

  fallback <- function(reason) {
    fit <- call_estimatr(estimatr::lm_robust, nse_args,
                         c(list(formula = formula), plain_args), envir)
    tag_fit(fit, adjustment = "none", selected = character(0), fallback = reason,
            fn = "lm_robust_lasso", formula = formula)
  }

  if (length(sel_vars) == 0) return(fallback("no covariates selected"))

  adj_formula <- stats::as.formula(paste0(
    "`", outcome, "` ~ `", treatment, "` + ",
    paste0("`", sel_vars, "`", collapse = " + ")
  ))

  fit_adj <- tryCatch(
    call_estimatr(estimatr::lm_robust, nse_args,
                  c(list(formula = adj_formula), plain_args), envir),
    error = function(e) NULL
  )
  if (is.null(fit_adj)) return(fallback("adjusted lm_robust failed to fit"))

  if (!finite_treatment_se(fit_adj, treatment)) {
    return(fallback("adjusted fit had a non-finite treatment standard error"))
  }

  tag_fit(fit_adj, adjustment = "robust", selected = sel_vars, fallback = NA_character_,
          fn = "lm_robust_lasso", formula = formula)
}


#' Treatment-by-moderator model with LASSO-selected covariates
#'
#' Fits \code{outcome ~ treatment * moderator + selected covariates}, where the
#' covariates are chosen by \code{\link{lasso_select_covariates}} from the focal
#' \code{outcome ~ treatment} equation.
#'
#' The selected covariates enter additively and are deliberately not
#' interacted with treatment or with the moderator. Adding those interactions
#' would spend degrees of freedom on terms that are not the estimand and, with
#' a data-driven covariate set, would make the reported interaction sensitive
#' to which covariates happened to be selected.
#'
#' The treatment-by-moderator coefficient describes how the treatment effect
#' varies across levels of an observed variable. The moderator is not
#' randomly assigned, so that variation is descriptive: it is not the causal
#' effect of the moderator, and it does not decompose the treatment effect into
#' a pathway.
#'
#' @param formula A two-sided formula \code{outcome ~ treatment}, the focal
#'   estimand equation passed to selection. The moderator is supplied
#'   separately.
#' @param moderator Character scalar naming the moderator column.
#' @param data A data frame.
#' @param covariates Candidate covariates: a one-sided formula, a character
#'   vector of column names, or a bare expression.
#' @param clusters,ci,alpha Passed through to
#'   \code{\link[estimatr]{lm_robust}}.
#' @param lasso_args A named list of further arguments for
#'   \code{\link{lasso_select_covariates}}.
#'
#' @return An \code{lm_robust} object, carrying the attributes described in
#'   \code{\link{adjustment}}. When the interacted fit with selected
#'   covariates cannot be produced, the fallback is the same interaction model
#'   fitted on the full candidate covariate set, not an unadjusted model.
#'
#' @examples
#' set.seed(1)
#' n <- 600
#' dat <- data.frame(
#'   Z = rep(0:1, n / 2),
#'   X_pid = rep(c(0, 1), each = n / 2),
#'   X1 = rnorm(n), X2 = rnorm(n)
#' )
#' dat$Y <- 0.2 * dat$Z + 0.6 * dat$Z * dat$X_pid + 1.5 * dat$X1 + rnorm(n)
#'
#' lm_moderator_lasso(Y ~ Z, moderator = "X_pid", data = dat, covariates = ~ X1 + X2)
#'
#' @importFrom stats as.formula coef reformulate
#' @family adjusted estimators
#' @export
lm_moderator_lasso <- function(formula, moderator, data, covariates = NULL,
                         clusters = NULL, ci = TRUE, alpha = 0.05,
                         lasso_args = list()) {
  cov_names <- if (is.null(covariates)) {
    character(0)
  } else {
    resolve_covariate_names(substitute(covariates), covariates)
  }

  selected <- if (length(cov_names) > 0) {
    do.call(
      lasso_select_covariates,
      c(list(formula = formula, covariates = cov_names, data = data), lasso_args)
    )
  } else {
    ~1
  }
  sel_vars <- all.vars(selected)

  outcome   <- as.character(formula[[2]])
  treatment <- as.character(formula[[3]])

  interaction_formula <- function(covs) {
    rhs <- paste0("`", treatment, "` * `", moderator, "`")
    if (length(covs) > 0) {
      rhs <- paste(rhs, "+", paste0("`", covs, "`", collapse = " + "))
    }
    stats::as.formula(paste0("`", outcome, "` ~ ", rhs))
  }

  nse_args <- list(clusters = substitute(clusters))
  envir <- parent.frame()

  call_lm_robust <- function(fml) {
    call_estimatr(estimatr::lm_robust, nse_args,
                  list(formula = fml, data = data, ci = ci, alpha = alpha),
                  envir)
  }

  fit <- tryCatch(call_lm_robust(interaction_formula(sel_vars)), error = function(e) NULL)

  if (!is.null(fit) && any(is.finite(stats::coef(fit)))) {
    return(tag_fit(fit,
                   adjustment = if (length(sel_vars) > 0) "robust" else "none",
                   selected = sel_vars,
                   fallback = NA_character_,
                   fn = "lm_moderator_lasso", formula = formula))
  }

  fit_full <- call_lm_robust(interaction_formula(cov_names))
  tag_fit(fit_full, adjustment = "robust", selected = cov_names,
          fallback = "selected-covariate interaction fit was degenerate; used the full candidate set",
          fn = "lm_moderator_lasso", formula = formula)
}


#' What adjustment did a fit actually use?
#'
#' The fitting functions in this package fall back to a simpler specification
#' when covariate selection returns nothing or the adjusted fit is degenerate.
#' These accessors report what was actually run, so a fallback is visible in
#' the results rather than buried inside the wrapper.
#'
#' @param fit A fit returned by \code{\link{lm_lin_lasso}},
#'   \code{\link{lm_robust_lasso}}, or \code{\link{lm_moderator_lasso}}.
#'
#' @return \code{adjustment()} returns \code{"lin"}, \code{"robust"}, or
#'   \code{"none"}. \code{selected_covariates()} returns the character vector
#'   of covariates actually used, possibly empty. \code{fallback_reason()}
#'   returns \code{NA_character_} when no fallback fired, and otherwise a
#'   description of why it did.
#'
#' @examples
#' set.seed(1)
#' n <- 300
#' dat <- data.frame(Z = rep(0:1, n / 2), X1 = rnorm(n))
#' dat$Y <- 0.5 * dat$Z + rnorm(n)
#'
#' # X1 is pure noise, so selection returns nothing and the fit falls back
#' fit <- lm_lin_lasso(Y ~ Z, ~ X1, data = dat)
#' adjustment(fit)
#' fallback_reason(fit)
#'
#' @family fallback reporting
#' @export
adjustment <- function(fit) attr(fit, "estimatrTools_adjustment")

#' @rdname adjustment
#' @export
selected_covariates <- function(fit) attr(fit, "estimatrTools_selected")

#' @rdname adjustment
#' @export
fallback_reason <- function(fit) attr(fit, "estimatrTools_fallback")


#' Summarize what a collection of fits actually did
#'
#' \code{\link{adjustment}} and \code{\link{fallback_reason}} report on one fit
#' at a time, which is not usable across a pipeline of hundreds. This walks a
#' list of fits and returns one row each, so a fallback that quietly changed the
#' estimator for a handful of specifications is visible without inspecting every
#' fit by hand.
#'
#' A fallback is not an error: it means the requested adjustment could not be
#' produced and a simpler specification was used instead. That is usually the
#' right thing to do, but it changes what the estimate is, so it belongs in the
#' record rather than in the wrapper.
#'
#' @param fits A list of fits from \code{\link{lm_lin_lasso}},
#'   \code{\link{lm_robust_lasso}}, or \code{\link{lm_moderator_lasso}}. If the list
#'   is named, the names are used to label rows. Elements that were not produced
#'   by this package are reported with \code{adjustment = NA}.
#' @param only_fallbacks Logical. If \code{TRUE}, return only the rows where a
#'   fallback fired (default \code{FALSE}).
#'
#' @return A data frame with columns \code{fit} (name or index),
#'   \code{adjustment}, \code{n_selected}, \code{selected}, and
#'   \code{fallback_reason}.
#'
#' @examples
#' set.seed(1)
#' n <- 400
#' dat <- data.frame(Z = rep(0:1, n / 2), X_sig = rnorm(n), X_noise = rnorm(n))
#' dat$Y <- 0.5 * dat$Z + 2 * dat$X_sig + rnorm(n)
#'
#' fits <- list(
#'   informative = lm_lin_lasso(Y ~ Z, ~ X_sig + X_noise, data = dat),
#'   noise_only  = lm_lin_lasso(Y ~ Z, ~ X_noise, data = dat)
#' )
#' fallback_summary(fits)
#' fallback_summary(fits, only_fallbacks = TRUE)
#'
#' @family fallback reporting
#' @export
fallback_summary <- function(fits, only_fallbacks = FALSE) {
  if (!is.list(fits)) {
    stop("fallback_summary: expected a list of fits, got ", class(fits)[1], ".")
  }

  empty <- data.frame(
    fit = character(0), adjustment = character(0), n_selected = integer(0),
    selected = character(0), fallback_reason = character(0),
    stringsAsFactors = FALSE
  )
  if (length(fits) == 0) return(empty)

  labels <- names(fits)
  if (is.null(labels)) labels <- as.character(seq_along(fits))
  labels[labels == ""] <- as.character(seq_along(fits))[labels == ""]

  rows <- lapply(seq_along(fits), function(i) {
    fit <- fits[[i]]
    sel <- selected_covariates(fit)
    data.frame(
      fit = labels[i],
      adjustment = if (is.null(adjustment(fit))) NA_character_ else adjustment(fit),
      n_selected = if (is.null(sel)) NA_integer_ else length(sel),
      selected = if (is.null(sel)) NA_character_ else paste(sel, collapse = ", "),
      fallback_reason = if (is.null(fallback_reason(fit))) NA_character_ else fallback_reason(fit),
      stringsAsFactors = FALSE
    )
  })

  out <- do.call(rbind, rows)
  if (only_fallbacks) out <- out[!is.na(out$fallback_reason), , drop = FALSE]
  rownames(out) <- NULL
  out
}


#' Record what a fit actually did
#'
#' @param fit A fitted model.
#' @param adjustment One of "lin", "robust", "none".
#' @param selected Character vector of covariates used.
#' @param fallback Reason string, or NA when no fallback fired.
#' @return The fit, with attributes attached.
#' @keywords internal
#' @noRd
tag_fit <- function(fit, adjustment, selected, fallback, fn = NA_character_,
                    formula = NULL) {
  attr(fit, "estimatrTools_adjustment") <- adjustment
  attr(fit, "estimatrTools_selected") <- selected
  attr(fit, "estimatrTools_fallback") <- fallback
  # Also record it where it will outlive the object: the dominant calling shape
  # discards the fit, so the attributes alone leave the fallback invisible.
  record_fit(fn, formula, adjustment, selected, fallback)
  fit
}


#' Call an estimatr fitting function, preserving non-standard evaluation
#'
#' \code{estimatr} captures \code{clusters}, \code{weights}, and \code{subset}
#' by non-standard evaluation and resolves them against \code{data}, so a user
#' can write \code{clusters = village_id} for a column of \code{data}. Two
#' things have to be true at once:
#'
#' \itemize{
#'   \item A supplied argument must reach \code{estimatr} unevaluated, so a
#'     bare column name still resolves. Evaluating it here would fail, because
#'     the column does not exist in this frame.
#'   \item An argument the caller never supplied must be omitted entirely
#'     rather than passed as \code{NULL}. Passing \code{NULL} is not neutral:
#'     \code{estimatr} looks the name up in \code{data} first, so a data frame
#'     with a column literally named \code{clusters} would silently cluster a
#'     fit that was never meant to be clustered.
#' }
#'
#' Building the call from the caller's own expressions and evaluating it in the
#' caller's frame satisfies both.
#'
#' @param fun The \code{estimatr} function to call.
#' @param nse_args Named list of unevaluated expressions captured with
#'   \code{substitute}; entries that are \code{NULL} are dropped.
#' @param plain_args Named list of ordinary values; \code{NULL} entries are
#'   dropped.
#' @param envir Frame in which to evaluate the constructed call.
#' @return The fitted model.
#' @keywords internal
#' @noRd
call_estimatr <- function(fun, nse_args, plain_args, envir) {
  nse_args <- nse_args[!vapply(nse_args, is.null, logical(1))]
  plain_args <- plain_args[!vapply(plain_args, is.null, logical(1))]

  the_call <- as.call(c(list(fun), plain_args, nse_args))
  eval(the_call, envir)
}


#' Does the treatment row of a fit have a usable standard error?
#'
#' @param fit A fitted model.
#' @param treatment Name of the treatment variable.
#' @return TRUE when a treatment term exists with a finite standard error.
#' @keywords internal
#' @noRd
finite_treatment_se <- function(fit, treatment) {
  tidied <- tryCatch(broom::tidy(fit), error = function(e) NULL)
  if (is.null(tidied) || !"term" %in% names(tidied)) return(FALSE)
  rows <- tidied[startsWith(tidied$term, treatment), , drop = FALSE]
  nrow(rows) > 0 && is.finite(rows$std.error[1])
}
