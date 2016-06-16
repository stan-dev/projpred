#' Dispersion parameter
#'
#' Get the parameter related to dispersion for some models from a stanreg object.
#' @param fit A stanreg object.

get_dispersion <- function(fit) {

  dis <- NA

  if(fit$family$family %in% c('gaussian','Gamma')) {
    par <- switch(fit$family$family,
                    'gaussian' = 'sigma',
                    'Gamma' = 'shape')
    dis <- apply(as.array(fit$stanfit), 3, rbind)[, par]
  }

  dis
}
