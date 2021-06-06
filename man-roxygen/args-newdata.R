#' @param newdata Passed to argument \code{newdata} of the reference model's
#'   \code{extract_model_data} function (see \code{\link{init_refmodel}}).
#'   Provides the predictor (and possibly also the response) data for the new
#'   observations.
#' @param offsetnew Passed to argument \code{orhs} of the reference model's
#'   \code{extract_model_data} function (see \code{\link{init_refmodel}}).
#'   Used to get the offsets for the (new) observations.
#' @param weightsnew Passed to argument \code{wrhs} of the reference model's
#'   \code{extract_model_data} function (see \code{\link{init_refmodel}}).
#'   Used to get the weights for the (new) observations.
