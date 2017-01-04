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
  if(ns > ncol(vars$beta)) {
    warning(paste0('Setting the number of samples to ', ncol(vars$beta),'.'))
    ns <- ncol(vars$beta)
  }

  if(is.null(intercept)) intercept <- vars$intercept
  if(is.null(nv_max) || nv_max > NROW(vars$beta)) nv_max <- NROW(vars$beta)

  e <- get_data_and_parameters(vars, d_test, intercept, ns, family_kl)

  chosen <- select(method, e$p_full, e$d_train, family_kl, intercept, nv_max,
                   verbose)

  p_sub <- .get_submodels(chosen, c(0, seq_along(chosen)), family_kl, e$p_full,
                          e$d_train, intercept)
  sub_summaries <- .get_sub_summaries2(chosen, e$p_full, e$data, p_sub, family_kl,
                                       intercept)
  full_summaries <- .get_full_summaries(e$p_full, e$data, e$coef_full,
                                        family_kl, intercept)

  b_weights <- .get_bootstrap_ws(NROW(e$d_test$x))

  metrics <- .bootstrap_metrics(sub_summaries, full_summaries, e$data,
                                family_kl, intercept, e$is_test, b_weights)
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
  # - Returns d_train, data, p_full, coef_full, is_test.
  # - If d_test is NA, data equals to d_train and is_test is FALSE.
  #   Otherwise data is set to d_train and is_test is TRUE.

  mu <- family_kl$mu_fun(vars$x, vars$alpha, vars$beta, vars$offset, intercept)

  d_train <- list(x = vars$x, weights = vars$weights, offset = vars$offset)

  # if test data doesn't exist, use training data to evaluate mse, r2, mlpd
  is_test <- is.list(d_test)
  if(is_test) {
    # check that d_test is of the correct form?
    data <- d_test
    if(is.null(data$weights)) data$weights <- rep(1, nrow(data$x))
    if(is.null(data$offset)) data$offset <- rep(0, nrow(data$x))
  } else {
    data <- vars[c('x', 'weights', 'offset', 'y')]
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

  list(d_train = d_train, data = data, p_full = p_full,
       coef_full = coef_full, is_test = is_test)
}


