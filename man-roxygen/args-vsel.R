#' @param object Either a \code{refmodel}-type object created by
#'   \link[=init_refmodel]{init_refmodel}, an object which can be converted to a
#'   reference model using \link[=get_refmodel]{get_refmodel}, or a \code{vsel}
#'   object resulting from \code{varsel} or \code{cv_varsel}.
#' @param method The method used in the variable selection. Possible options are
#'   \code{'L1'} for L1-search and \code{'forward'} for forward selection.
#'   Default is 'forward' if the number of variables in the full data is at most
#'   20, and \code{'L1'} otherwise.
#' @param cv_search If \code{TRUE}, then the projected coefficients after
#'   L1-selection are computed without any penalization (or using only the
#'   regularization determined by \code{regul}). If \code{FALSE}, then the
#'   coefficients are the solution from the L1-penalized projection. This option
#'   is relevant only if \code{method = 'L1'}. Default is \code{TRUE} for
#'   genuine reference models and \code{FALSE} if \code{object} is datafit (see
#'   \link[=init_refmodel]{init_refmodel}).
#' @param ndraws Number of posterior draws used in the variable selection.
#'   Cannot be larger than the number of draws in the reference model.
#'   \strong{Caution:} For \code{ndraws <= 20}, the value of \code{ndraws} is
#'   passed to \code{nclusters} (so that clustering is used). Ignored if
#'   \code{nclusters} is not \code{NULL} or if \code{method = "L1"} (L1 search
#'   uses always one cluster). See also section "Details" below.
#' @param nclusters Number of clusters of posterior draws used in the variable
#'   selection. Ignored if \code{method = "L1"} (L1 search uses always one
#'   cluster). For the meaning of \code{NULL}, see argument \code{ndraws}. See
#'   also section "Details" below.
#' @param ndraws_pred Number of posterior draws used for prediction (after
#'   selection). Cannot be larger than the number of draws in the reference
#'   model. \strong{Caution:} For \code{ndraws_pred <= 20}, the value of
#'   \code{ndraws_pred} is passed to \code{nclusters_pred} (so that clustering
#'   is used). Ignored if \code{nclusters_pred} is not \code{NULL}. See also
#'   section "Details" below.
#' @param nclusters_pred Number of clusters of posterior draws used for
#'   prediction (after selection). For the meaning of \code{NULL}, see argument
#'   \code{ndraws_pred}. See also section "Details" below.
#' @param nterms_max Maximum number of variables until which the selection is
#'   continued. Defaults to \code{min(20, D, floor(0.4 * n))} where \code{n} is
#'   the number of observations and \code{D} the number of variables. Note that
#'   \code{nterms_max} does not count the intercept, so use \code{nterms_max =
#'   0} for the intercept-only model.
#' @param penalty Vector determining the relative penalties or costs for the
#'   variables. A value of \code{0} means that those variables have no cost and
#'   will therefore be selected first, whereas \code{Inf} means those variables
#'   will never be selected. Currently works only if \code{method = 'L1'}. By
#'   default \code{1} for each variable.
#' @param lambda_min_ratio Ratio between the smallest and largest lambda in the
#'   L1-penalized search. This parameter essentially determines how long the
#'   search is carried out, i.e., how large submodels are explored. No need to
#'   change the default value unless the program gives a warning about this.
#' @param nlambda Number of values in the lambda grid for L1-penalized search.
#'   No need to change unless the program gives a warning about this.
#' @param thresh Convergence threshold when computing L1-path. Usually no need
#'   to change this.
#' @param regul Amount of regularization in the projection. Usually there is no
#'   need for regularization, but sometimes for some models the projection can
#'   be ill-behaved and we need to add some regularization to avoid numerical
#'   problems.
#' @param search_terms A custom character vector of terms to consider for
#'   selection. The intercept (\code{"1"}) needs to be included explicitly. The
#'   default considers all the terms in the reference model's formula.
#' @param verbose A single logical value indicating whether to print out
#'   additional information while running (\code{TRUE}) or not (\code{FALSE}).
#' @param ... Additional arguments to be passed to the \code{get_refmodel}
#'   function.
