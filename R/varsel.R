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

  e <- construct_objects(fit, vars, d_test, intercept, ns, family_kl)

  chosen <- select(method, e$p_full, e$d_train, family_kl, intercept, nv_max,
                   verbose)

  sub_list <- .get_sub_summaries2(chosen, e$d_train, e$d_test, e$p_full,
                                  family_kl, intercept)
  full_list <- .get_full_summaries(e$d_test, e$p_full, e$coef_full, family_kl)

  b_weights <- .gen_bootstrap_ws(length(e$d_test$y))

  stats <- stats_arr(sub_list, full_list, e$d_test, family_kl, intercept,
                     b_weights)
  fit$varsel <- list(chosen = chosen, stats = stats)

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

construct_objects <- function(fit, vars, d_test, intercept, ns, family_kl) {

  mu <- family_kl$mu_fun(vars$x, vars$alpha, vars$beta, vars$offset, intercept)

  d_train <- list(x = vars$x, weights = vars$weights, offset = vars$offset)

  # if test data doesn't exist, use training data to evaluate mse, r2, mlpd
  if(is.list(d_test)) {
    if(is.null(d_test$weights)) d_test$weights <- rep(1, nrow(d_test$x))
    if(is.null(d_test$offset)) d_test$offset <- rep(0, nrow(d_test$x))
    d_test$is_d_train <- F
  } else {
    d_test <- vars[c('x', 'weights', 'offset', 'y')]
    d_test$is_d_train <- T
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

  list(mu = mu, d_train = d_train, d_test = d_test,
       p_full = p_full, coef_full = coef_full)
}

stats_arr <- function(sub_list, full_list, d_test, family_kl,
                      intercept, b_weights) {

  kl_list <- sapply(sub_list, function(x) x$kl)
  stats_array <- data.frame(data = 'sel', size = seq_along(sub_list)-1,
                            delta = F, summary = 'kl', value = kl_list,
                            lq = NA, uq = NA)

  # evaluate performance on test data and
  # use bayesian bootstrap to get 95% credible intervals
  b_stats <- cbind(data = ifelse(d_test$is_d_train, 'train', 'test'),
                   .bootstrap_stats(sub_list, full_list, d_test, family_kl,
                                    b_weights, intercept))

  rbind(stats_array, b_stats, make.row.names = F)
}
