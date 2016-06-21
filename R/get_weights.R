#' Observation weights
#'
#' Get the observation weights from a stanreg object (currently only for binomial model).
#' @param \code{fit} A \link[=stanreg-objects]{stanreg} object.

get_weights <- function(fit) {

  if('(weights)' %in% colnames(fit$model)) {
    fit$model[,'(weights)']
  } else if(fit$family$family == 'binomial' &&
            is.matrix(fit$model[[1]]) && ncol(fit$model[[1]])>1) {
    rowSums(fit$model[[1]])
  } else {
    rep(1, nrow(get_x(fit)))
  }
}
