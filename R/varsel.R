#' Variable selection for generalized linear models
#'
#' Perform the projection predictive variable selection for a generalized
#' linear model fitted with rstanarm.
#' @param fit A \link[=stanreg-objects]{stanreg} object.
#' @param d_test A test dataset, which is used to evaluate model performance. If
#' not provided, training data is used.
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
varsel <- function(fit, d_test = NA, ...) {
  UseMethod('varsel')
}

#' @export
varsel.stanreg <- function(fit, d_test = NA, method='L1', ...) {

  .validate_for_varsel(fit)
  vars <- .extract_vars(fit)
  args <- .init_args(list(...), vars)
  family_kl <- kl_helpers(family(fit))
  # for models with no intercept, alpha = rep(0, nrow(beta))
  mu <- family_kl$mu_fun(vars$x, vars$alpha, vars$beta, vars$offset,
                         args$intercept)

  d_train <- list(x = vars$x, weights = vars$weights, offset = vars$offset)
  coef_init <- list(alpha = median(vars$alpha),
                    beta = matrix(apply(vars$beta, 1, median), ncol = 1))

  # indices of samples that are used in the projection
  s_ind <- round(seq(1, args$ns_total, length.out  = args$ns))
  p_full <- list(mu = mu[, s_ind], dis = vars$dis[s_ind],
                 cluster_w = rep(1/args$ns, args$ns))

  # cluster the samples of the full model if clustering is wanted
  # for the variable selection
  clust <- if(args$clust) .get_p_clust(mu, vars$dis, args$nc) else NULL
  p_sel <- if(args$clust) clust$p else p_full

  # the actual selection
  chosen <- select(method, p_sel, d_train, family_kl, args$intercept, args$nv,
                     args$regul, coef_init, args$verbose)

  # if test data doesn't exist, use training data to evaluate mse, r2, mlpd
  eval_data <- ifelse(is.list(d_test), 'test', 'train')
  if(eval_data == 'train') {
    d_test <- within(d_train, y <- vars$y)
  } else {
    if(is.null(d_test$weights)) d_test$weights <- rep(1, nrow(d_test$x))
    if(is.null(d_test$offset)) d_test$offset <- rep(0, nrow(d_test$x))
  }

  coef_full <- list(alpha = vars$alpha[s_ind], beta = vars$beta[, s_ind])

  # Perform the projection for the chosen variables and calculate lppds
  stats_list <- .summary_stats(chosen, d_train, d_test, p_full, family_kl,
                          args$intercept, args$regul, coef_init, coef_full)

  nv_list <- 1:length(stats_list$sub$mu) - args$intercept
  stats_array <- data.frame(data = 'sel', size = nv_list, delta = F,
                            summary = 'kl', value = stats_list$kl, lq = NA, uq = NA)

  # if function was called by cv_varsel, return also mu and lppd,
  if(args$cv) {
    res <- list(chosen = chosen, stats_list = stats_list, stats = stats_array)
    if(args$clust) res$cl <- clust$cl
    return(res)
  }

  # evaluate performance on test data and
  # use bayesian bootstrap to get 95% credible intervals
  b_weights <- .gen_bootstrap_ws(length(d_test$y), args$n_boot)
  b_stats <- .bootstrap_stats(stats_list, nv_list, d_test, family_kl, b_weights,
                              eval_data, args$intercept)

  stats_array <- rbind(stats_array, b_stats, make.row.names = F)
  res <- list(chosen = chosen, stats = stats_array)
  if(args$clust) res$cl <- clust$cl

  fit$varsel <- res
  fit
}



select <- function(method, p_sel, d_train, family_kl, intercept, pmax,
                   regul, coef_init, verbose) 
{
    #
    # Auxiliary function, performs variable selection with the given method,
    # and returns the variable ordering.
    #
    if (tolower(method) == 'l1') {
        chosen <- search_L1(p_sel, d_train, family_kl, intercept, pmax)
    } else if (tolower(method) == 'forward') {
        tryCatch(chosen <- fsel(p_sel, d_train, family_kl, intercept, pmax,
                                regul, coef_init, verbose),
                 'error' = .varsel_errors)
    } else {
        stop(sprintf('Unknown search method: %s', method))
    }
    return(chosen)
}


# function() {
#     # indices of samples that are used in the projection
#     s_ind <- round(seq(1, args$ns_total, length.out  = args$ns))
#     p_full <- list(mu = mu[, s_ind], dis = vars$dis[s_ind],
#                    cluster_w = rep(1/args$ns, args$ns))
#     
#     # cluster the samples of the full model if clustering is wanted
#     # for the variable selection
#     clust <- if(args$clust) .get_p_clust(mu, vars$dis, args$nc) else NULL
#     p_sel <- if(args$clust) clust$p else p_full
# }
