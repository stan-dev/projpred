#' Variable selection for generalized linear models
#'
#' Perform the projection predictive variable selection for a generalized
#' linear model fitted with rstanarm.
#' @param fit A \link[=stanreg-objects]{stanreg} object.
#' @param d_test A test dataset, which is used to evaluate model performance.
#' If not provided, training data is used.
#' @param method The method used in the variable selection. Possible options are
#' \code{'L1'} for L1-search and \code{'forward'} for forward selection.
#' Defaults to \code{'L1'}.
#' @param ... Optional arguments. Possible arguments and their defaults are:
#' \describe{
#'  \item{\code{ns = min(400, [number of draws])}}{
#'    Number of draws used in the variable selection.
#'    Cannot be larger than the number of draws in the full model.}
#'  \item{\code{nc = 0}}{
#'    If nonzero, a clustering with \code{nc} clusters is performed for
#'    the draws and the cluster centers are used in the variable selection
#'    instead of the actual draws.}
#'  \item{\code{nv = min(ncol(x) - 1, rankMatrix(x))}}{
#'    Maximum number of variables to be used in the projection (incl. intercept).
#'    Cannot be larger than \code{min(ncol(x) - 1, rankMatrix(x))}.}
#'  \item{\code{verbose = FALSE}}{
#'    If \code{verbose = TRUE}, prints information about the progress of the
#'    variable selection.}
#' }
#'
#' @return The original \link[=stanreg-objects]{stanreg} object augmented with an element 'varsel',
#' which is a list containing the following elements:
#' \describe{
#'  \item{\code{chosen}}{The order in which the variables were added to the submodel.}
#'  \item{\code{stats}}{An array with statistics of the submodel performance.}
#'  \item{\code{family}}{A \code{\link{family}}-object.}
#' }
#'
#' @examples
#' \dontrun{
#' ### Usage with stanreg objects
#' fit <- stan_glm(y~x, binomial())
#' fit_v <- varsel(fit)
#' plot_varsel(fit_v)
#' }
#'

#' @export
varsel <- function(fit, d_test = NA, method = 'L1', ns = 400L,
                   nv_max = NULL, intercept = NULL, verbose = F, ...) {
    UseMethod('varsel')
}

#' @export
varsel.stanreg <- function(fit, d_test = NA, method = 'L1', ns = 400L,
                           nv_max = NULL, intercept = NULL, verbose = F, ...) {
  .validate_for_varsel(fit)
  vars <- .extract_vars(fit)
  family_kl <- kl_helpers(family(fit))
  if(is.null(intercept)) intercept <- vars$intercept
  if(is.null(nv_max) || nv_max > NROW(vars$beta)) nv_max <- NROW(vars$beta)

  e <- get_data_and_parameters(vars, d_test, intercept, ns, family_kl)

  chosen <- select(method, e$p_full, e$d_train, family_kl, intercept, nv_max,
                   verbose)

  p_sub <- .get_submodels(chosen, c(0, seq_along(chosen)), family_kl, e$p_full,
                          e$d_train, intercept)
  sub_summaries <- .get_sub_summaries2(p_sub, chosen, e$d_eval, e$p_full, family_kl,
                                       intercept)
  full_summaries <- .get_full_summaries(e$d_eval, e$p_full, e$coef_full,
                                        family_kl, intercept)

  b_weights <- .get_bootstrap_ws(NROW(e$d_test$x))

  metrics <- .bootstrap_metrics(sub_summaries, full_summaries, e$d_eval,
                                family_kl, intercept, b_weights,
                                e$eval_is_test)
  kl <- .get_kl_array(p_sub)

  fit$varsel <- list(chosen = chosen, metrics = rbind(kl, metrics))

  fit
}

select <- function(method, p_full, d_train, family_kl, intercept, nv_max,
                   verbose) {
  #
  # Auxiliary function, performs variable selection with the given method,
  # and returns the variable ordering.
  #
  if (tolower(method) == 'l1') {
    chosen <- search_L1(p_full, d_train, family_kl, intercept, nv_max)
  } else if (tolower(method) == 'forward') {
    tryCatch(chosen <- search_forward(p_full, d_train, family_kl, intercept,
                                      nv_max, verbose),
             'error' = .varsel_errors)
  } else {
    stop(sprintf('Unknown search method: %s.', method))
  }
  return(chosen)
}

get_data_and_parameters <- function(vars, d_test, intercept, ns, family_kl) {
  # construct d_train/test, p_full and coef_full

  mu <- family_kl$mu_fun(vars$x, vars$alpha, vars$beta, vars$offset, intercept)

  d_train <- list(x = vars$x, weights = vars$weights, offset = vars$offset)

  # if test data doesn't exist, use training data to evaluate mse, r2, mlpd
  eval_is_test <- is.list(d_test)
  if(eval_is_test) {
    # check that d_test is of the correct form?
    d_eval <- d_test
    if(is.null(d_eval$weights)) d_eval$weights <- rep(1, nrow(d_eval$x))
    if(is.null(d_eval$offset)) d_eval$offset <- rep(0, nrow(d_eval$x))
  } else {
    d_eval <- vars[c('x', 'weights', 'offset', 'y')]
  }
  # indices of samples that are used in the projection
  s_ind <- round(seq(1, ncol(vars$beta), length.out  = ns))
  p_full <- list(mu = mu[, s_ind], dis = vars$dis[s_ind],
                 weights = rep(1/ns, ns))
  coef_full <- list(alpha = vars$alpha[s_ind], beta = vars$beta[, s_ind])

  # cluster the samples of the full model if clustering is wanted
  # for the variable selection
  do_clust <- F
  clust <- if(do_clust) get_p_clust(mu, vars$dis, ns) else NULL
  if(do_clust) p_full <- clust

  list(d_train = d_train, d_eval = d_eval, p_full = p_full,
       coef_full = coef_full, eval_is_test = eval_is_test)
}


