#' Test for differential attrition with LASSO-selected covariates
#'
#' @description
#' **Experimental.** Tests whether treatment predicts outcome missingness, both
#' unconditionally and allowing the pattern to differ across covariates, using a
#' double-LASSO step to keep the interacted test from running out of degrees of
#' freedom.
#'
#' The fully interacted attrition model spends
#' \code{(n_arms - 1) x (1 + n_covariates)} degrees of freedom. With the whole
#' covariate pool that quickly drives events per variable below 10 for any study
#' with low missingness or many arms, and below that threshold the model is
#' unreliable in both directions: too sparse to detect real differential dropout,
#' and capable of spurious significance from numerical instability. Selecting a
#' parsimonious covariate set first is what makes the interacted test usable.
#'
#' Two selection equations, neither involving treatment, are run via
#' \code{\link{lasso_select_one}} and their union enters the interacted model:
#'
#' \itemize{
#'   \item **Outcome equation**: the outcome on all candidates, among respondents
#'     who answered. Retains the covariates where differential dropout would
#'     directly bias the treatment effect estimate.
#'   \item **Dropout equation**: an indicator for missing outcome on all
#'     candidates, among everyone assigned. Retains the covariates that
#'     characterize who drops out.
#' }
#'
#' The covariate-free test is always computed and is the primary criterion for
#' arm-level differential dropout. The interacted test is secondary: it detects
#' subgroup-specific dropout that the simple test misses even when overall rates
#' are equal. Both are returned and the caller decides which to act on.
#'
#' @section Experimental:
#' This is an original procedure motivated by events-per-variable concerns, not a
#' reference implementation from a published paper. The double-LASSO step is
#' adapted from Belloni, Chernozhukov and Hansen (2014) but applied to an
#' attrition diagnostic rather than a treatment effect estimator, and the second
#' equation is not theirs (see below). The threshold of 10 events per variable
#' follows the rule of thumb of Peduzzi et al. (1996) for logistic regression;
#' the same concern applies to a linear probability model.
#'
#' @section Why selecting and testing on the same data is still valid here:
#' Choosing covariates with the same data used to test them usually invalidates
#' the reference distribution. It does not here, because selection never looks at
#' treatment: the outcome equation uses the outcome and the covariates, the
#' dropout equation uses missingness and the covariates. Under the null of no
#' differential attrition, missingness is independent of treatment, so
#' conditioning on a selected set that is a function of the covariates and
#' missingness leaves treatment randomly assigned and the test of the treatment
#' terms approximately valid.
#'
#' That argument depends on treatment being independent of the selection, so it
#' does not survive a blocked or clustered design analysed as if it were simple
#' random assignment.
#'
#' Note also that the outcome equation is not the second equation of Belloni,
#' Chernozhukov and Hansen, which would be treatment on the covariates. Since
#' treatment is randomized, that equation selects nothing asymptotically, which is
#' why it is replaced here. The substitution is deliberate, and it means the
#' citation motivates the procedure rather than licensing it.
#'
#' @section p_interacted is sensitive to the cross-validation fold draw:
#' The covariate set is chosen by cross-validated LASSO, and which covariates
#' survive depends on how the folds happen to be drawn. When the candidate pool is
#' large relative to the number of events, that dependence is strong enough to move
#' \code{p_interacted} across any threshold you might apply to it. On one real
#' study with eight candidates, the same outcome gave \code{n_selected} of 0 and
#' \code{p_interacted} of \code{NA} under one seed, 1 and 0.031 under another, and
#' 1 and 0.423 under a third.
#'
#' \code{seed} makes a given run reproducible; it does not make the answer stable.
#' Two consequences follow. Do not report \code{p_interacted} as though it were a
#' fixed property of the study without checking its sensitivity across seeds. And
#' prefer \code{p_simple} as the criterion for acting on a study, since it involves
#' no selection step and is invariant to all of this.
#'
#' @section Clustered designs:
#' The interacted test is a multi-degree-of-freedom Wald test, and such tests
#' over-reject badly under cluster randomization even with cluster-robust standard
#' errors: about 12 percent at a nominal 5 percent with 30 clusters. The cause is
#' not the denominator degrees of freedom, so no correction repairs it.
#' \code{p_simple} is unaffected. Treat \code{p_interacted} as descriptive when
#' assignment was clustered.
#'
#' @param data A data frame or tibble. Rows with a missing treatment value are
#'   dropped before all calculations, so it may hold unassigned respondents.
#' @param treatment Unquoted name of the treatment variable.
#' @param outcomes Character vector of outcome variable names. If `NULL`, all
#'   columns beginning with `"Y"` are used.
#' @param covariates Candidate covariate pool: a one-sided formula, a character
#'   vector of column names, or a bare expression. If `NULL` or empty, only the
#'   covariate-free test is run.
#' @param epv_threshold Minimum events per variable required to run the
#'   interacted test. Default `10`. Below it, `p_interacted` is `NA` and
#'   `epv_adequate` is `FALSE`.
#' @param lambda_rule Either `"lambda.1se"` (default, conservative: fewer
#'   covariates, better events per variable) or `"lambda.min"`.
#' @param alpha Threshold at which the `flag_*` columns are set (default `0.05`).
#' @param seed Integer seed for the LASSO cross-validation folds, passed to
#'   [lasso_select_one()]. The global random number state is restored on exit.
#' @param .method Regression function for the tests (default
#'   `estimatr::lm_robust`). Must accept formula and data arguments.
#' @param study_id Optional character scalar. If provided, a `study_id` column
#'   holding this value is appended to the result, so results stack across
#'   studies.
#' @param quiet Logical. The default `TRUE` returns the result, which auto-prints
#'   at the console. `FALSE` prints it and returns invisibly.
#' @param ... Additional arguments passed to `.method`.
#'
#' @return A tibble with one row per outcome and columns:
#'   \describe{
#'     \item{outcome}{Outcome variable name.}
#'     \item{n_assigned}{Number of respondents with non-missing treatment.}
#'     \item{n_missing}{Number of missing outcome observations.}
#'     \item{pct_missing}{Proportion missing.}
#'     \item{p_simple}{P-value from the covariate-free omnibus test,
#'       \code{I(missing) ~ Z}.}
#'     \item{n_outcome_eq, n_dropout_eq}{Covariates selected by each equation.}
#'     \item{n_selected}{Size of the union.}
#'     \item{selected_covariates}{Comma-separated names of the selected set.}
#'     \item{df1}{Numerator degrees of freedom the interacted test would spend,
#'       \code{(n_arms - 1) * (1 + n_selected)}.}
#'     \item{epv}{Events per variable, \code{n_missing / df1}.}
#'     \item{epv_adequate}{\code{epv >= epv_threshold} and \code{n_selected > 0}.}
#'     \item{p_interacted}{P-value from the interacted test,
#'       \code{I(missing) ~ Z * (demeaned selected covariates)}, testing all
#'       treatment terms jointly. \code{NA} when \code{epv_adequate} is
#'       \code{FALSE}.}
#'     \item{estimable}{\code{FALSE} when the missingness indicator does not vary,
#'       so there is no test to run at all and \code{p_simple} is \code{NA}.}
#'     \item{flag_simple, flag_interacted, flag}{Whether each p-value, and either,
#'       falls at or below \code{alpha}.}
#'   }
#'
#' @section When there is no p-value, and what that means:
#' \code{status} records why, because the reasons mean opposite things.
#' \code{"tested"} is a computed p-value. \code{"no_attrition"} means nobody dropped
#' out, which is a \strong{pass}: with no attrition there can be no differential
#' attrition. \code{"all_missing"} means everybody did, which is uninformative.
#' \code{"not_estimable"} means the indicator varied but no statistic came back.
#' \code{estimable} equals \code{status == "tested"}.
#'
#' The distinction decides a denominator, and the two choices differ by a lot. To
#' report how much differential attrition a corpus has, count \code{"no_attrition"}
#' rows as passes and exclude only the uninformative ones. To assess whether the
#' computed p-values are uniform, as a valid design implies, use \code{"tested"}
#' rows only, because a row with no attrition produces no draw from Uniform(0, 1) to
#' compare against. On one real corpus those two rates were 2.5 percent and 7.4
#' percent, so quoting the second while describing the first overstates the failure
#' rate roughly threefold.
#'
#' Reporting such a row as \code{p_value = 1} makes the mirror error: it puts a spike
#' at 1 into the uniformity diagnostic, which then reports a badly non-uniform
#' collection when nothing is wrong.
#'
#' @references
#' Belloni, A., Chernozhukov, V., and Hansen, C. (2014). Inference on treatment
#' effects after selection among high-dimensional controls. \emph{Review of
#' Economic Studies}, 81(2), 608-650. \doi{10.1093/restud/rdt044}
#'
#' Lin, W. (2013). Agnostic notes on regression adjustments to experimental data:
#' reexamining Freedman's critique. \emph{Annals of Applied Statistics}, 7(1),
#' 295-318. \doi{10.1214/12-AOAS583}
#'
#' Peduzzi, P., Concato, J., Kemper, E., Holford, T. R., and Feinstein, A. R.
#' (1996). A simulation study of the number of events per variable in logistic
#' regression analysis. \emph{Journal of Clinical Epidemiology}, 49(12),
#' 1373-1379. \doi{10.1016/S0895-4356(96)00236-3}
#'
#' @examples
#' set.seed(42)
#' n <- 500
#' dat <- data.frame(
#'   Z = rep(c(0L, 1L), n / 2),
#'   X_age = rnorm(n, 50, 10),
#'   X_income = rnorm(n, 50000, 10000)
#' )
#' dat$Y_outcome <- 0.3 * dat$Z + 0.5 * scale(dat$X_age) + rnorm(n)
#' p_miss <- ifelse(dat$Z == 1, 0.05 + 0.01 * scale(dat$X_age), 0.05)
#' dat$Y_outcome[which(rbinom(n, 1, pmax(0, pmin(1, p_miss))) == 1)] <- NA
#'
#' check_attrition_lasso(dat, Z,
#'   outcomes   = "Y_outcome",
#'   covariates = c("X_age", "X_income"))
#'
#' @importFrom stats as.formula coef df.residual pf vcov
#' @family attrition diagnostics
#' @export
check_attrition_lasso <- function(data, treatment,
                                  outcomes      = NULL,
                                  covariates    = NULL,
                                  epv_threshold = 10,
                                  lambda_rule   = c("lambda.1se", "lambda.min"),
                                  alpha         = 0.05,
                                  seed          = 999,
                                  .method       = estimatr::lm_robust,
                                  study_id      = NULL,
                                  quiet         = TRUE,
                                  ...) {
  lambda_rule <- match.arg(lambda_rule)

  # Accept a bare column name or a string, without evaluating the promise: a
  # bare name has no binding outside `data`, so testing is.character(treatment)
  # first would error before the character branch could help.
  treatment_expr <- substitute(treatment)
  treatment_name <- if (is.character(treatment_expr)) {
    treatment_expr
  } else {
    deparse1(treatment_expr)
  }

  if (!treatment_name %in% names(data)) {
    stop("check_attrition_lasso: treatment column '", treatment_name,
         "' not found in data.", call. = FALSE)
  }

  y_cols <- if (is.null(outcomes)) {
    grep("^Y", names(data), value = TRUE)
  } else {
    outcomes
  }
  if (length(y_cols) == 0) {
    warning("check_attrition_lasso: no outcomes selected.")
    return(invisible(NULL))
  }
  absent_y <- setdiff(y_cols, names(data))
  if (length(absent_y) > 0) {
    stop("check_attrition_lasso: outcome column(s) not found in data: ",
         paste(absent_y, collapse = ", "), call. = FALSE)
  }

  # Drop unassigned respondents once, before anything is counted.
  data <- data[!is.na(data[[treatment_name]]), , drop = FALSE]

  cov_cols <- if (is.null(covariates)) {
    character(0)
  } else {
    resolve_covariate_names(substitute(covariates), covariates)
  }
  absent <- setdiff(cov_cols, names(data))
  if (length(absent) > 0) {
    stop("check_attrition_lasso: covariate column(s) not found in data: ",
         paste(absent, collapse = ", "),
         ". Silently testing a smaller set would misreport the diagnostic.",
         call. = FALSE)
  }
  has_covs <- length(cov_cols) > 0
  n_arms   <- length(unique(data[[treatment_name]]))

  results <- lapply(y_cols, function(yc) {
    miss       <- as.integer(is.na(data[[yc]]))
    n_assigned <- nrow(data)
    n_missing  <- sum(miss)

    if (length(unique(miss)) < 2L) {
      # Separate the two ways the indicator can fail to vary, because they mean
      # opposite things. Nobody missing is a pass: with no attrition there can be no
      # differential attrition, so the design question is answered in the
      # affirmative, and counting such a row as a missing test inflates any rate
      # computed over the rows that remain. Everybody missing is uninformative.
      return(tibble::tibble(
        outcome = yc, n_assigned = n_assigned, n_missing = as.integer(n_missing),
        pct_missing = mean(miss), p_simple = NA_real_,
        n_outcome_eq = 0L, n_dropout_eq = 0L, n_selected = 0L,
        selected_covariates = "", df1 = NA_integer_, epv = NA_real_,
        epv_adequate = FALSE, p_interacted = NA_real_,
        status = if (n_missing == 0L) "no_attrition" else "all_missing",
        estimable = FALSE,
        flag_simple = FALSE, flag_interacted = FALSE, flag = FALSE
      ))
    }

    # --- Covariate-free omnibus test ----
    fit_df     <- data
    fit_df$.miss <- miss
    simple_fit <- .method(stats::as.formula(paste(".miss ~", treatment_name)),
                          data = fit_df, ...)
    p_simple   <- broom::glance(simple_fit)$p.value

    # --- Two selection equations, neither involving treatment ----
    outcome_eq <- character(0)
    dropout_eq <- character(0)
    if (has_covs) {
      outcome_eq <- lasso_select_one(yc, cov_cols, data = data,
                                     subset = !is.na(data[[yc]]),
                                     lambda_rule = lambda_rule, seed = seed)
      dropout_eq <- lasso_select_one(".miss", cov_cols, data = fit_df,
                                     lambda_rule = lambda_rule, seed = seed)
    }
    selected <- union(outcome_eq, dropout_eq)
    n_sel    <- length(selected)

    # --- Events-per-variable gate ----
    df1    <- if (n_sel > 0) (n_arms - 1L) * (1L + n_sel) else NA_integer_
    epv    <- if (!is.na(df1) && df1 > 0) n_missing / df1 else NA_real_
    epv_ok <- !is.na(epv) && epv >= epv_threshold && n_sel > 0

    # --- Interacted test, only when the gate is passed ----
    p_interacted <- NA_real_
    if (epv_ok) {
      demeaned  <- demean_columns(data[, selected, drop = FALSE])
      int_df    <- cbind(fit_df[, c(".miss", treatment_name), drop = FALSE], demeaned)
      dm_names  <- colnames(demeaned)
      covar_rhs <- paste(paste0("`", dm_names, "`"), collapse = " + ")
      fmla      <- stats::as.formula(
        paste(".miss ~", treatment_name, "* (", covar_rhs, ")")
      )
      fit <- tryCatch(.method(fmla, data = int_df, ...), error = function(e) NULL)

      if (!is.null(fit)) {
        # In a fully interacted model every coefficient that is neither the
        # intercept nor a covariate main effect is a treatment term. Taking the
        # complement keeps a covariate whose name happens to begin with the
        # treatment name out of the test, which prefix matching would not.
        b_all      <- stats::coef(fit)
        bare       <- gsub("`", "", names(b_all), fixed = TRUE)
        test_terms <- names(b_all)[bare != "(Intercept)" & !bare %in% dm_names]
        q          <- length(test_terms)
        wald <- tryCatch({
          b <- b_all[test_terms]
          if (anyNA(b)) stop("rank deficient", call. = FALSE)
          v <- stats::vcov(fit)[test_terms, test_terms, drop = FALSE]
          as.numeric(t(b) %*% solve(v) %*% b)
        }, error = function(e) NULL)
        if (!is.null(wald)) {
          p_interacted <- stats::pf(wald / q, q, stats::df.residual(fit),
                                    lower.tail = FALSE)
        }
      }
    }

    flag_s <- !is.na(p_simple) && p_simple <= alpha
    flag_i <- epv_ok && !is.na(p_interacted) && p_interacted <= alpha

    tibble::tibble(
      outcome             = yc,
      n_assigned          = n_assigned,
      n_missing           = as.integer(n_missing),
      pct_missing         = mean(miss),
      p_simple            = p_simple,
      n_outcome_eq        = length(outcome_eq),
      n_dropout_eq        = length(dropout_eq),
      n_selected          = n_sel,
      selected_covariates = paste(selected, collapse = ", "),
      df1                 = df1,
      epv                 = epv,
      epv_adequate        = epv_ok,
      p_interacted        = p_interacted,
      status              = if (is.na(p_simple)) "not_estimable" else "tested",
      estimable           = !is.na(p_simple),
      flag_simple         = flag_s,
      flag_interacted     = flag_i,
      flag                = flag_s || flag_i
    )
  })

  out <- do.call(rbind, results)
  if (!is.null(study_id)) out$study_id <- study_id
  if (!quiet) {
    print(out)
    return(invisible(out))
  }
  out
}


#' Demean the columns of a covariate frame, following the Lin estimator
#'
#' Expands factor and character columns through \code{model.matrix}, then
#' subtracts each column's mean, so that the treatment coefficient in an
#' interacted model retains its average-effect interpretation.
#'
#' @param covariate_df A data frame of covariates.
#' @return A data frame of demeaned columns, names suffixed \code{"_dm"}.
#' @keywords internal
#' @noRd
demean_columns <- function(covariate_df) {
  mm_formula <- stats::reformulate(names(covariate_df))
  mm <- stats::model.matrix(
    mm_formula,
    data = stats::model.frame(mm_formula, covariate_df, na.action = stats::na.pass)
  )
  mm <- mm[, -1, drop = FALSE]
  centered <- sweep(mm, 2, colMeans(mm, na.rm = TRUE))
  colnames(centered) <- paste0(colnames(centered), "_dm")
  as.data.frame(centered)
}
