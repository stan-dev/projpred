#' Extract draws of the linear predictor and draw from the predictive
#' distribution of the projected submodel
#'
#' \code{proj_linpred} extracts draws of the linear predictor and
#' \code{proj_predict} draws from the predictive distribution of the projected
#' submodel or submodels. If the projection has not been performed, the
#' functions also perform the projection.
#'
#' @name proj-pred
#'
#' @param object Either an object returned by \link[=varsel]{varsel}, \link[=cv_varsel]{cv_varsel}
#' or \link[=init_refmodel]{init_refmodel}, or alternatively any object that can be converted to a reference model.
#' @param xnew The predictor values used in the prediction. If \code{vind} is
#' specified, then \code{xnew} should either be a dataframe containing column names
#' that correspond to \code{vind} or a matrix with the number and order of columns
#' corresponding to \code{vind}. If \code{vind} is unspecified, then \code{xnew} must
#' either be a dataframe containing all the column names as in the original data or a matrix
#' with the same columns at the same positions as in the original data.
#' @param ynew New (test) target variables. If given, then the log predictive density
#' for the new observations is computed.
#' @param offsetnew Offsets for the new observations. By default a vector of
#' zeros.
#' @param weightsnew Weights for the new observations. For binomial model,
#' corresponds to the number trials per observation. For \code{proj_linpred},
#' this argument matters only if \code{ynew} is specified. By default a vector
#' of ones.
#' @param transform Should the linear predictor be transformed using the
#' inverse-link function? Default is \code{FALSE}. For \code{proj_linpred} only.
#' @param integrated If \code{TRUE}, the output is averaged over the
#' parameters. Default is \code{FALSE}. For \code{proj_linpred} only.
#' @param nv Number of variables in the submodel (the variable combination is
#' taken from the variable selection information). If a vector with several values,
#' then results for all specified model sizes are returned. Ignored if \code{vind} is specified.
#' By default use the automatically suggested model size.
#' @param draws Number of draws to return from the predictive distribution of
#' the projection. The default is 1000.
#' For \code{proj_predict} only.
#' @param seed An optional seed to use for drawing from the projection.
#' For \code{proj_predict} only.
#' @param ... Additional argument passed to \link{project} if \code{object}
#' is an object returned by \link{varsel} or \link{cv_varsel}.
#'
#' @return If the prediction is done for one submodel only (\code{nv} has length one
#' or \code{vind} is specified) and ynew is unspecified, a matrix or vector of
#' predictions (depending on the value of \code{integrated}). If \code{ynew} is specified,
#' returns a list with elements pred (predictions) and lpd (log predictive densities).
#' If the predictions are done for several submodel sizes, returns a list with one element
#' for each submodel.
#'
#' @examples
#' \donttest{
## Usage with stanreg objects
#' fit <- stan_glm(y~x, binomial())
#' vs <- varsel(fit)
#'
## compute predictions with 4 variables at the training points
#' pred <- proj_linpred(vs, xnew=x, nv = 4)
#' pred <- proj_predict(vs, xnew=x, nv = 4)
#'
#' }
#'
NULL

## The 'helper' for proj_linpred and proj_predict, ie. does all the
## functionality that is common to them. It essentially checks all the arguments
## and sets them to their respective defaults and then loops over the
## projections. For each projection, it evaluates the fun-function, which
## calculates the linear predictor if called from proj_linpred and samples from
## the predictive distribution if called from proj_predict.
proj_helper_poc <- function(object, xnew, offsetnew, weightsnew, nv, seed,
                            proj_predict, ...) {

  if (is.null(offsetnew)) offsetnew <- rep(0, nrow(xnew))
  if (is.null(weightsnew)) weightsnew <- rep(1, nrow(xnew))

  if ('projection' %in% class(object) ||
      (length(object)>0 && 'projection' %in% class(object[[1]]))) {
    proj <- object
  } else {
    ## reference model or varsel object obtained, so run the projection
    proj <- project_poc(object = object, nv = nv, ...)
  }

  if (!.is_proj_list(proj)) {
    proj <- list(proj)
  } else {
    ## proj is not a projection object
    if(any(sapply(proj, function(x) !('family_kl' %in% names(x)))))
      stop(paste('proj_linpred only works with objects returned by',
                 ' varsel, cv_varsel or project'))
  }

  ## FIXME: proj$sub_fit always exists, but does $coefficients?
  projected_sizes <- sapply(proj, function(x) NCOL(x$sub_fit$coefficients))
  nv <- list(...)$nv %ORifNULL% projected_sizes

  if (!all(nv %in% projected_sizes))
    stop(paste0('Linear prediction requested for nv = ',
                paste(nv, collapse = ', '),
                ', but projection performed only for nv = ',
                paste(projected_sizes, collapse = ', '), '.'))

  projs <- Filter(function(x) NCOL(x$sub_fit$coefficients) %in% nv, proj)
  names(projs) <- nv

  xnew_df <- is.data.frame(xnew)
  if (xnew_df) {
    terms <- unique(unlist(lapply(projs, function(x) unlist(unname(x$vind)))))
    xnew <- .df_to_model_mat(xnew, terms)
  }

  if (!is.matrix(xnew))
    stop('xnew not provided in the correct format. See ?proj-pred.')

  vind <- list(...)$vind
  if (!is.null(vind) && NCOL(xnew) != length(vind))
    stop(paste('The number of columns in xnew does not match with the given',
               'number of variable indices (vind).'))

  ## set random seed but ensure the old RNG state is restored on exit
  rng_state_old <- rngtools::RNGseed()
  on.exit(rngtools::RNGseed(rng_state_old))
  set.seed(seed)

  preds <- lapply(projs, function(proj) {
    if (xnew_df) {
      xtemp <- xnew[, min(1, length(proj$vind)):length(proj$vind), drop = F]
    } else if (!is.null(vind)) {
      ## columns of xnew are assumed to match to the given variable indices
      xtemp <- xnew
    } else {
      ## fetch the right columns from the feature matrix
      xtemp <- xnew[, unique(unname(unlist(proj$vind))), drop = F]
    }
    mu <- proj$family_kl$mu_fun(proj$sub_fit, xnew=xtemp)

    proj_predict(proj, mu, offsetnew, weightsnew)
  })

  .unlist_proj(preds)
}

#' @rdname proj-pred
#' @export
proj_linpred_poc <- function(object, xnew, ynew = NULL, offsetnew = NULL,
                             weightsnew = NULL, nv = NULL, transform = FALSE,
                             integrated = FALSE, ...) {

  ## function to perform to each projected submodel
  proj_predict <- function(proj, mu, offset, weights) {
    pred <- t(mu)
    if (!transform) pred <- proj$family_kl$linkfun(pred)
    if (integrated) {
      ## average over the parameters
      pred <- as.vector( proj$weights %*% pred )
    } else if (!is.null(dim(pred)) && nrow(pred) == 1) {
      ## return a vector if pred contains only one row
      pred <- as.vector(pred)
    }

    if (!is.null(ynew)) {
      ## compute also the log-density
      target <- .get_standard_y(ynew, weights, proj$family_kl)
      ynew <- target$y
      weights <- target$weights
      lpd <- proj$family_kl$ll_fun(mu, proj$dis, ynew, weights)
      if (integrated && !is.null(dim(lpd))) {
        lpd <- as.vector(apply(lpd, 1, log_weighted_mean_exp, proj$weights))
      } else if (!is.null(dim(lpd))) {
        lpd <- t(lpd)
      }
      return(nlist(pred, lpd))
    } else {
      return(pred)
    }
  }

  ## proj_helper lapplies fun to each projection in object
  proj_helper_poc(object = object, xnew = xnew, offsetnew = offsetnew,
                  weightsnew = weightsnew, nv = nv, seed = NULL, proj_predict =
                                                                   proj_predict, ...)
}

#' @rdname proj-pred
#' @export
proj_predict_poc <- function(object, xnew, offsetnew = NULL, weightsnew = NULL,
                             nv = NULL, draws = 1000, seed = NULL, ...) {

  ## function to perform to each projected submodel
  proj_predict <- function(proj, mu, offset, weights) {
    draw_inds <- sample(x = seq_along(proj$weights), size = draws,
                        replace = TRUE, prob = proj$weights)

    t(sapply(draw_inds, function(i) {
      proj$family_kl$ppd(mu[,i], proj$dis[i], weights)
    }))
  }

  ## proj_helper lapplies fun to each projection in object
  proj_helper_poc(object = object, xnew = xnew, offsetnew = offsetnew,
                  weightsnew = weightsnew, nv = nv, seed = seed,
                  proj_predict = proj_predict, ...)
}
