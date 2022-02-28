#' @param newdata Passed to argument `newdata` of the reference model's
#'   `extract_model_data` function (see [init_refmodel()]). Provides the
#'   predictor (and possibly also the response) data for the new (or old)
#'   observations. May also be `NULL` (see argument `extract_model_data` of
#'   [init_refmodel()]).
#' @param offsetnew Passed to argument `orhs` of the reference model's
#'   `extract_model_data` function (see [init_refmodel()]). Used to get the
#'   offsets for the new (or old) observations.
#' @param weightsnew Passed to argument `wrhs` of the reference model's
#'   `extract_model_data` function (see [init_refmodel()]). Used to get the
#'   weights for the new (or old) observations.
