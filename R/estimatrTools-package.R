#' @keywords internal
"_PACKAGE"

#' @section Selecting covariates:
#' \code{\link{lasso_select_covariates}} chooses a covariate set by running LASSO of
#' the outcome on candidates \emph{within each treatment arm}, and of each arm
#' indicator on candidates. The first is what makes the selection appropriate
#' for an estimator with a separate slope per arm; the second is the
#' double-selection step that retains covariates predicting assignment even when
#' they barely predict the outcome.
#'
#' @section Fitting:
#' \code{\link{lm_lin_lasso}} fits the Lin (2013) interacted estimator on the
#' selected set, \code{\link{lm_robust_lasso}} fits them additively, and
#' \code{\link{lm_moderator_lasso}} fits a treatment-by-moderator model with the
#' selected covariates entered additively alongside.
#'
#' @section Knowing what actually ran:
#' Each fitting function falls back to a simpler specification when selection
#' returns nothing or the adjusted fit is degenerate. This is usually the right
#' response, but it changes what the estimate is, so it is recorded rather than
#' hidden. How to read it back depends on what you still have:
#' \itemize{
#'   \item one fit in hand: \code{\link{adjustment}},
#'     \code{\link{selected_covariates}}, \code{\link{fallback_reason}}
#'   \item a list of fits: \code{\link{fallback_summary}}
#'   \item fits that were built and discarded, which is the usual case inside a
#'     \code{map()} pipeline: \code{\link{fallback_log}}
#' }
#' The last matters more than it sounds. Selection returning nothing is common,
#' not exceptional, and when it happens the "adjusted" estimate is an unadjusted
#' difference in means. A pipeline that never reads the log has no way of
#' knowing what fraction of its estimates that describes.
#'
#' @name estimatrTools-package
NULL
