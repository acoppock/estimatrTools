# Package-level record of what each fit actually did.
#
# The attributes attached by tag_fit() answer the question for one fit, which is
# enough when a fit is assigned to a name and kept. The dominant calling shape in
# practice is not that: fits are built anonymously inside map(), passed straight
# to a summarizing function, and then dropped with select(-model). By the time
# anyone asks which fits were adjusted, the objects are gone.
#
# So the record is kept here as well, one row per call, and survives the fit.
# This is deliberate package state. It is bounded (a handful of scalars per
# call), it can be turned off with options(estimatrTools.log = FALSE), and it is
# cleared with reset_fallback_log().

.fallback_env <- new.env(parent = emptyenv())
.fallback_env$rows <- list()

#' Record of what every fit did
#'
#' \code{\link{adjustment}} and friends report on a fit you still hold.
#' \code{\link{fallback_summary}} reports on a list of fits you still hold. This
#' reports on every fit made in the session, including ones that were discarded,
#' which is the usual case in a pipeline that builds a fit inside \code{map()},
#' summarizes it, and drops the model column.
#'
#' Falling back is not an error: it means the requested adjustment could not be
#' produced and a simpler specification was used. It does change what the
#' estimate is, though, so a pipeline that never looks at this log has no way of
#' knowing how many of its "adjusted" estimates are unadjusted.
#'
#' @param fallbacks_only Logical. Return only the calls where a fallback fired
#'   (default \code{FALSE}, which returns one row per call).
#'
#' @return A data frame with one row per call: \code{call_index}, \code{fn},
#'   \code{outcome}, \code{treatment}, \code{adjustment}, \code{n_selected}, and
#'   \code{fallback_reason} (\code{NA} when none fired).
#'
#' @section Turning it off:
#' Recording is on by default and costs a few scalars per call. Disable with
#' \code{options(estimatrTools.log = FALSE)}. \code{reset_fallback_log()} clears
#' what has accumulated, which is worth doing at the top of a pipeline so the
#' log describes that run rather than everything since the session started.
#'
#' @examples
#' reset_fallback_log()
#'
#' set.seed(1)
#' n <- 300
#' dat <- data.frame(Z = rep(0:1, n / 2), X_sig = rnorm(n), X_noise = rnorm(n))
#' dat$Y <- 0.5 * dat$Z + 2 * dat$X_sig + rnorm(n)
#'
#' invisible(lm_lin_lasso(Y ~ Z, ~ X_sig + X_noise, data = dat))
#' invisible(lm_lin_lasso(Y ~ Z, ~ X_noise, data = dat))
#'
#' fallback_log()
#'
#' # what fraction of fits were actually adjusted?
#' table(fallback_log()$adjustment)
#'
#' @seealso \code{\link{adjustment}} for a single fit,
#'   \code{\link{fallback_summary}} for a list of fits you have kept.
#' @family fallback reporting
#' @export
fallback_log <- function(fallbacks_only = FALSE) {
  rows <- .fallback_env$rows

  empty <- data.frame(
    call_index = integer(0), fn = character(0), outcome = character(0),
    treatment = character(0), adjustment = character(0),
    n_selected = integer(0), fallback_reason = character(0),
    stringsAsFactors = FALSE
  )
  if (length(rows) == 0) return(empty)

  out <- do.call(rbind, lapply(rows, as.data.frame, stringsAsFactors = FALSE))
  rownames(out) <- NULL
  if (fallbacks_only) out <- out[!is.na(out$fallback_reason), , drop = FALSE]
  rownames(out) <- NULL
  out
}

#' @rdname fallback_log
#' @return \code{reset_fallback_log()} returns nothing and is called for its
#'   effect.
#' @export
reset_fallback_log <- function() {
  .fallback_env$rows <- list()
  invisible(NULL)
}

#' Append one call to the log
#'
#' @param fn Name of the fitting function.
#' @param formula The model formula, used to name the outcome and treatment.
#' @param adjustment One of "lin", "robust", "none".
#' @param selected Character vector of covariates used.
#' @param fallback Reason string, or NA.
#' @return Nothing; called for its effect.
#' @keywords internal
#' @noRd
record_fit <- function(fn, formula, adjustment, selected, fallback) {
  if (!isTRUE(getOption("estimatrTools.log", TRUE))) return(invisible(NULL))

  outcome <- tryCatch(as.character(formula[[2]])[1], error = function(e) NA_character_)
  treatment <- tryCatch(as.character(formula[[3]])[1], error = function(e) NA_character_)

  .fallback_env$rows[[length(.fallback_env$rows) + 1L]] <- list(
    call_index      = length(.fallback_env$rows) + 1L,
    fn              = fn,
    outcome         = outcome,
    treatment       = treatment,
    adjustment      = adjustment,
    n_selected      = length(selected),
    fallback_reason = if (is.null(fallback)) NA_character_ else fallback
  )
  invisible(NULL)
}
