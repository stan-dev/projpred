#' Variable selection for generalized linear models
#'
#' Perform the projection predictive variable selection for a generalized
#' linear model fitted with rstanarm.
#' @param fit Either a \link[=stanreg-objects]{stanreg}-object or an object returned
#' by \link[=init_refmodel]{init_refmodel}.
#' @param d_test A test dataset, which is used to evaluate model performance.
#' If not provided, training data is used.
#' @param method The method used in the variable selection. Possible options are
#' \code{'L1'} for L1-search and \code{'forward'} for forward selection.
#' Default is 'forward' if the number of variables in the full data is at most 20, and
#' \code{'L1'} otherwise.
#' @param ns Number of posterior draws used in the variable selection.
#'    Cannot be larger than the number of draws in the full model.
#'    Ignored if nc is set.
#' @param nc Number of clusters to use in the clustered projection.
#'    Overrides the \code{ns} argument. Defaults to 1.
#' @param nspred Number of samples used for prediction (after selection). Ignored if ncpred is given.
#' @param ncpred Number of clusters used for prediction (after selection). Default is 1.
#' @param nv_max Maximum number of varibles until which the selection is continued.
#'    Defaults to min(D, floor(0.4*n)) where n is the number of observations and
#'    D the number of variables.
#' @param intercept Whether to use intercept in the submodels. Defaults to TRUE.
#' @param verbose If TRUE, may print out some information during the selection.
#'    Defaults to FALSE.
#' @param regul Amount of regularization in the projection. Usually there is no need for 
#' regularization, but sometimes for some models the projection can be ill-behaved and we
#' need to add some regularization to avoid numerical problems. Default is 1e-9.
#' @param ... Currently ignored.
#'
#'
#' @return The original fit-object object augmented with a field 'varsel',
#' which is a list containing the following elements:
#' \describe{
#'  \item{\code{vind}}{The order in which the variables were added to the submodel.}
#'  \item{\code{kl}}{KL-divergence for each submodel size.}
#'  \item{\code{summaries}}{Summary statistics computed during the selection.}
#'  \item{\code{d_test}}{The data used to evaluate the summaries.}
#'  \item{\code{family_kl}}{A modified \link{family}-object.}
#' }
#'
#' @examples
#' \donttest{
#' ### Usage with stanreg objects
#' fit <- stan_glm(y~x, binomial())
#' fit_v <- varsel(fit)
#' plot_varsel(fit_v)
#' }
#'

#' @export
varsel <- function(fit, d_test = NULL, method = NULL, ns = NULL, nc = NULL, 
                   nspred = NULL, ncpred = NULL, nv_max = NULL, 
                   intercept = NULL, verbose = F, regul=1e-9, ...) {


  .validate_for_varsel(fit)

  if ((is.null(ns) && is.null(nc)) || tolower(method)=='l1')
  	# use one cluster for selection by default, and always with L1-search
  	nc <- 1
  if (is.null(nspred) && is.null(ncpred))
    # use 1 clusters for prediction by default
    ncpred <- 1

  vars <- .extract_vars(fit)
  family_kl <- vars$fam
  
  if (is.null(method)) {
    if (dim(vars$x)[2] <= 20)
      method <- 'forward'
    else
      method <- 'L1'
  }

  if(is.null(intercept))
    intercept <- vars$intercept
  if(is.null(nv_max) || nv_max > NCOL(vars$x)) {
  	nv_max_default <- floor(0.4*length(vars$y)) # a somewhat sensible default limit for nv_max
  	nv_max <- min(NCOL(vars$x), nv_max_default)
  }

  # training and test data
  d_train <- .get_traindata(fit)
  if (is.null(d_test)) {
  	d_test <- d_train
  	d_type <- 'train'
  } else {
  	d_test <- .fill_offset_and_weights(d_test)
  	d_type <- 'test'
  }

  p_sel <- .get_refdist(fit, ns, nc)
  p_pred <- .get_refdist(fit, nspred, ncpred)

  # perform the selection
  vind <- select(method, p_sel, d_train, family_kl, intercept, nv_max, verbose, regul)

  # statistics for the selected submodels
  p_sub <- .get_submodels(vind, c(0, seq_along(vind)), family_kl, p_pred,
                          d_train, intercept, regul)
  sub <- .get_sub_summaries(p_sub, d_test, family_kl)

  #
  if (d_type == 'train') {
      full <- .weighted_summary_means(d_test, family_kl, vars$wsample, vars$mu, vars$dis)
  } else {
      # TODO, FIGURE OUT HOW TO HANDLE THE TEST PREDICTIONS FOR FULL MODEL WHEN COEF_FULL ARE NOT AVAILABLE
      warning('Test predictions for the full model not yet implemented.')
      full <- NULL
  }

  # store the relevant fields into fit
  fit$varsel <- list(vind = setNames(vind, vars$coefnames[vind]),
                     kl = sapply(p_sub, function(x) x$kl),
                     d_test = c(d_test[c('y','weights')], type = d_type),
                     summaries = list(sub = sub, full = full),
                     family_kl = family_kl)

  fit
}


select <- function(method, p_sel, d_train, family_kl, intercept, nv_max,
                   verbose, regul) {
  #
  # Auxiliary function, performs variable selection with the given method,
  # and returns the variable ordering.
  #
  if (NCOL(d_train$x) == 1)
    # special case, only one variable, so no need for selection
    return(1)
  if (tolower(method) == 'l1') {
    vind <- search_L1(p_sel, d_train, family_kl, intercept, nv_max)
  } else if (tolower(method) == 'forward') {
    if ( NCOL(p_sel$mu) == 1)
      # only one mu column (one cluster or one sample), so use the optimized version of the forward search
      vind <- search_forward1(p_sel, d_train, family_kl, intercept, nv_max, verbose, regul)
    else
      # routine that can be used with several clusters
      tryCatch(vind <- search_forward(p_sel, d_train, family_kl, intercept, nv_max, verbose, regul),
               'error' = .varsel_errors)
  } else {
    stop(sprintf('Unknown search method: %s.', method))
  }
  return(vind)
}

