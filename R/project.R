#' @export
project <- function(object, nv, ...) {
  UseMethod('project')
}

#' @export
project.stanreg <- function(object, nv, ns = 400L, intercept = NULL, ...) {
  # Missing: Clustering!
  if(is.null(nv)) stop('nv not provided')
  .validate_for_varsel(object)
  if(!('varsel' %in% names(object)))
    stop(paste('The stanreg object doesn\'t contain information about the ',
               'variable selection. Run the variable selection first.'))

  vars <- .extract_vars(object)
  if(ns > ncol(vars$beta)) {
    warning(paste0('Setting the number of samples to ', ncol(vars$beta),'.'))
    ns <- ncol(vars$beta)
  }
  if(is.null(intercept)) intercept <- vars$intercept

  family_kl <- kl_helpers(family(object))

  if(max(nv) > length(object$varsel$chosen))
    stop(paste('Cannot perform the projection with', max(nv), 'variables,',
               'because the variable selection has been run only up to',
               length(object$varsel$chosen), 'variables.'))

  e <- .get_data_and_parameters(vars, NULL, intercept, ns, family_kl)

  object$proj <- list(
    p_sub = .get_submodels(object$varsel$chosen, nv, family_kl,
                           e$p_full, e$d_train, intercept),
    intercept = intercept)

  object
}
