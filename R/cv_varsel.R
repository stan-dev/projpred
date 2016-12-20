#' Variable selection for generalized linear models with cross-validation
#'
#' Perform the projection predictive variable selection for a generalized
#' linear model fitted with rstanarm.
#' @param fit A \link[=stanreg-objects]{stanreg} object.
#' @param fits An array with cross-validated stanfits and the respective
#' test datasets returned by \link[=stanreg-objects]{cv_fit}(fit).
#' If not provided, \link[=stanreg-objects]{cv_fit}(fit) is called to
#' get the array.
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
#'  \item{\code{pctch}}{Percentage of cross-validation runs that included the given
#'    variable to a model of given size.}
#'  \item{\code{stats}}{An array with statistics of the submodel performance.}
#'  \item{\code{family}}{A \code{\link{family}}-object.}
#' }
#'
#' @examples
#' \dontrun{
#' ### Usage with stanreg objects
#' fit <- stan_glm(y~x, binomial())
#' fits <- kfold(fit)
#' fit_v <- cv_varsel(fit, fits)
#' plot_varsel(fit_v)
#' }
#'

#' @export
cv_varsel <- function(fit, fits = NULL, ...) {
  UseMethod('cv_varsel')
}

#' @export
cv_varsel.stanreg <- function(fit, k_fold = NULL, ...) {

  .validate_for_varsel(fit)
  if(is.null(k_fold)) {
    print(paste('k_fold not provided, performing 10-fold cross-validation',
                'for the stan model.'))
    k_fold <- glmproj::kfold(fit, save_fits = T)
  }

  if(!all(apply(k_fold$fits, 1, function(fits, fit) {
    .validate_for_varsel(fits$fit)
    is.vector(fits$omitted) && max(fits$omitted) <= nobs(fit) && all(fits$omitted > 0)
  }, fit))) stop('k_fold does not have the correct form.')

  k <- attr(k_fold, 'K')
  vars <- .extract_vars(fit)
  args <- .init_args(c(list(...), cv = T), vars, family(fit))

  d_test <- lapply(k_fold$fits[,'omitted'], function(omitted, d_full) {
    list(x = d_full$x[omitted,], y = d_full$y[omitted],
         weights = d_full$weights[omitted], offset = d_full$offset[omitted])
  }, vars)

  # max number of variables to be projected
  args$nv <- min(c(sapply(k_fold$fits[,'fit'], function(fit)
    rankMatrix(get_x(fit))), args$nv))

  msgs <- paste('Forward selection for the',
                c('full model.', paste0('fold number ', 1:k,'/',k,'.')))
  # perform the forward selection
  sel <- mapply(function(fit, d_test, msg, args) {
    print(msg)
    do.call(varsel, c(list(fit = fit, d_test = d_test), args))
  }, c(list(full = fit), k_fold$fits[,'fit']), c(list(full = NA), d_test),
  msgs, MoreArgs = list(args))

  # combine cross validated results
  d_cv <- as.list(data.frame(apply(
    simplify2array(d_test)[c('y', 'weights', 'offset'),],
    1, function(x) do.call(c, x))))

  # extract and combine mu and lppd from sel[stats_list,-1]
  stop('Does not work at the moment.')

  # evaluate performance on test data and
  # use bayesian bootstrap to get 95% credible intervals
  b_weights <- .gen_bootstrap_ws(length(d_cv$y), args$n_boot)
  nv_list <- 1:length(stats$sub$mu) - args$intercept
  b_stats <- rbind(sel[['stats',1]],
    .bootstrap_stats(sel[['stats_list',1]], nv_list, vars, args$family_kl,
                     b_weights, 'train', args$intercept),
    .bootstrap_stats(mu_cv, lppd_cv, nv_list, d_cv, args$family_kl, b_weights,
                     'test', args$intercept), make.row.names = F)

  # find out how many of cross-validated forward selection iterations select
  # the same variables as the forward selection with all the data.
  chosen_full <- sel[['chosen',1]]
  pctch <- mapply(function(var_ind, ind, arr, k) sum(arr[1:min(ind+0, nrow(arr)), ] == var_ind)/k,
                  chosen_full, seq_along(chosen_full),
                  MoreArgs = list(do.call(cbind, sel['chosen',-1]), k))

  res <- list(chosen = chosen_full, pctch = pctch, stats = b_stats)
  if(args$clust) res$cl <- sel[['cl',1]]

  fit$varsel <- res
  fit
}


loo_varsel <- function(fit, method='L1', ...) {
    
    # TODO, ADD COMMENTS
    vars <- .extract_vars(fit)
    args <- .init_args(list(...), vars)
    fam <- kl_helpers(family(fit))
    mu <- fam$mu_fun(vars$x, vars$alpha, vars$beta, vars$offset, args$intercept)
    dis <- vars$dis
    
    # training data and the fit of the full model
    d_train <- list(x = vars$x, weights = vars$weights, offset = vars$offset)
    p_full <- list(mu = mu, dis = dis)
    
    # perform the clustering for the full model
    # TODO DUMMY SOLUTION, USE ONE CLUSTER
    cl <- rep(1,n)
    
    # compute the log-likelihood for the full model to obtain the LOO weights
    loglik <- log_lik(fit)
    lw <- psislw(-loglik)$lw_smooth
    
    n <- dim(lw)[2]
    print(n)
    print(dim(mu))
    print(dim(lw))
    # res <- rep(0,n)
    
    nv <- c(0:args$nv) # TODO IMPLEMENT THIS PROPERLY
    pmax <- max(nv) ## TODO
    
    tic()
    chosen_mat <- matrix(rep(0, n*pmax), nrow=n)
    loo <- matrix(nrow=n, ncol=length(nv))
    for (i in 1:n) {
    	
    	# reweight the clusters according to the is-loo weights
    	p_sel <- .get_p_clust(mu, dis, cl=cl, wsample=exp(lw[,i]))$p
    	
    	# res[i] <- p_sel$mu[i]-y[i]
    	
    	# perform selection
    	chosen <- select(method, p_sel, d_train, fam, args$intercept, pmax, args$regul, NA, args$verbose)
    	chosen_mat[i,] <- chosen
    	
    	# project onto the selected models and compute the difference between
    	# training and loo density for the left-out point
    	#psub <- .get_submodels(chosen, nv, fam, p_sel, d_train, args$intercept)
    	d_test = list(x=matrix(x[i,],nrow=1), y=y[i], offset=d_train$offset[i], weights=1.0)
    	summaries <- .get_sub_summaries(chosen, nv, d_train, d_test, p_sel, fam, args$intercept)
    	
    	
    	for (k in 1:length(nv)) {
    	    loo[i,k] <- summaries[[k]]$lppd
    	}
    	
    	
    	#coef_full <- list(alpha = vars$alpha[s_ind], beta = vars$beta[, s_ind])
    	#.summary_stats(chosen, d_train, d_test, p_full, fam,
    	#               args$intercept, args$regul, NA, coef_full) 
    	
    	
        
        # mu_loo_i <- mu %*% exp(lw[,i])
        # p_full <- list(mu=mu_loo_i)
        
        # chosen <- search_L1(p_full, d_train, fam, args$intercept, args$nv)
        # print(dim(  )) 
    	print(sprintf('i = %d', i))
    }
    toc()
    
    # p_full <- list(b = b[, s_ind], mu = mu[, s_ind], dis = dis[s_ind],
    #                cluster_w = rep(1/args$ns, args$ns))
    
    # mse <- mean((apply(mu,1,'mean')-y)^2)
    # mse_loo <- mean(res^2)
    # 
    # print(mse)
    # print(mse_loo)
    
    return(loo)
    
}


    
    
    
    
    
    
