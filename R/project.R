##' Projection to submodels
##'
##' Perform projection onto submodels of selected sizes or a specified feature
##' combination.
##'
##' @param object Either a \code{refmodel}-type object created by \link[=get_refmodel]{get_refmodel}
##' or \link[=init_refmodel]{init_refmodel}, or an object which can be converted to a reference model
##' using \link[=get_refmodel]{get_refmodel}.
##' @param nv Number of variables in the submodel (the variable combination is taken from the
##' \code{varsel} information). If a list, then the projection is performed for each model size.
##' Default is the model size suggested by the variable selection (see function \code{suggest_size}).
##'  Ignored if \code{vind} is specified.
##' @param vind Variable indices onto which the projection is done. If specified, \code{nv} is ignored.
##' @param cv_search If TRUE, then the projected coefficients after L1-selection are computed
##' without any penalization (or using only the regularization determined by \code{regul}). If FALSE, then
##' the coefficients are the solution from the L1-penalized projection. This option is relevant only
##' if L1-search was used. Default is TRUE for genuine reference models and FALSE if \code{object} is
##' datafit (see \link[=init_refmodel]{init_refmodel}).
##' @param ns Number of samples to be projected. Ignored if \code{nc} is specified. Default is 400.
##' @param nc Number of clusters in the clustered projection.
##' @param intercept Whether to use intercept. Default is \code{TRUE}.
##' @param seed A seed used in the clustering (if \code{nc!=ns}). Can be used
##' to ensure same results every time.
##' @param regul Amount of regularization in the projection. Usually there is no need for
##' regularization, but sometimes for some models the projection can be ill-behaved and we
##' need to add some regularization to avoid numerical problems.
##' @param ... Currently ignored.
##'
##' @return A list of submodels (or a single submodel if projection was performed onto
##' a single variable combination), each of which contains the following elements:
##' \describe{
##'  \item{\code{kl}}{The KL divergence from the reference model to the submodel.}
##'  \item{\code{weights}}{Weights for each draw of the projected model.}
##'  \item{\code{dis}}{Draws from the projected dispersion parameter.}
##'  \item{\code{alpha}}{Draws from the projected intercept.}
##'  \item{\code{beta}}{Draws from the projected weight vector.}
##'  \item{\code{vind}}{The order in which the variables were added to the submodel.}
##'  \item{\code{intercept}}{Whether or not the model contains an intercept.}
##'  \item{\code{family}}{A modified \code{\link{family}}-object.}
##' }
##'
##'
##' @examples
##' \donttest{
## Usage with stanreg objects
##' fit <- stan_glm(y~x, binomial())
##' vs <- varsel(fit)
##'
## project onto the best model with 4 variables
##' proj4 <- project(vs, nv = 4)
##'
## project onto an arbitrary variable combination (variable indices 3,4 and 8)
##' proj <- project(fit, vind=c(3,4,8))
##' }
##'

##' @export
project <- function(object, nv = NULL, vind = NULL, cv_search = TRUE, ns = 400, nc = NULL,
                        intercept = NULL, seed = NULL, regul=1e-4, ...) {

  if ( !('vsel' %in% class(object) || 'cvsel' %in% class(object)) && is.null(vind) )
    stop(paste('The given object is not a variable selection -object.',
               'Run the variable selection first, or provide the variable indices (vind).'))

  refmodel <- get_refmodel(object)

  if (cv_search) {
    ## use non-cv_searched solution for datafits by default
    cv_search <- !inherits(refmodel, "datafit")
  }

  if (!is.null(vind) &&
      any(object$vind[1:length(vind)] != vind)) {
    ## search path not found, or the given variable combination
    ## not in the search path, then we need to project the
    ## required variables
    cv_search <- TRUE
  }

  if (!is.null(vind)) {
    ## if vind is given, nv is ignored (project only onto the given submodel)
    vind <- object$vind[vind]
    if (length(vind) > count_terms_chosen(vind))
      nv <- count_terms_chosen(vind)
    else
      nv <- length(vind)
  } else {
    ## by default take the variable ordering from the selection
    vind <- object$vind
    if (is.null(nv)) {
      if (!is.null(object$suggested_size) && !is.na(object$suggested_size)) {
        ## by default, project onto the suggested model size
        nv <- object$suggested_size
      } else
        stop('No suggested model size found, please specify nv or vind')
    } else {
      if (!is.numeric(nv) || any(nv < 0))
        stop('nv must contain non-negative values.')
      if (max(nv) > length(vind)) {
        stop(paste('Cannot perform the projection with', max(nv), 'variables,',
                   'because variable selection was run only up to', length(vind),
                   'variables.'))
      }
    }
  }

  if (is.null(nc))
    ns <- min(ns, NCOL(refmodel$mu))

  if (is.null(intercept))
    intercept <- refmodel$intercept

  family <- refmodel$family

  ## get the clustering or subsample
  p_ref <- .get_refdist(refmodel, ns = ns, nc = nc, seed = seed)

  ## project onto the submodels
  subm <- .get_submodels(object$spath, nv, family, p_ref,
                             refmodel, intercept, regul, cv_search=cv_search)

  ## add family
  proj <- lapply(subm, function(model) {
    model <- c(model, nlist(family), list(p_type = is.null(ns)))
    class(model) <- 'projection'
    return(model)
  })

  ## If only one model size, just return the proj instead of a list of projs
  .unlist_proj(proj)
}
