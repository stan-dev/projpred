#' Predict
#'
#' A predict method for glmproj objects.
#'
#' @param object A glmproj object.
#' @param x A model matrix.
#' @param d The number of features to be used in the prediction.
#'  If not provided, uses all features.
#'
#' @export

predict.glmproj <- function(object, x, d) {
  if(missing(d)) d <- length(object$b)
  featureinds <- object$chosen[1:d]

  object$family$linkinv(x[,featureinds, drop = F]%*%object$b[[d]])
}
