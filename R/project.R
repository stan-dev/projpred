#' Projection onto submodel(s)
#'
#' Project the reference model onto a single submodel consisting of a specific
#' combination of predictor terms or onto a single or multiple submodels of
#' specific sizes.
#'
#' @name project
#'
#' @param object Either a \code{refmodel}-type object created by
#'   \link[=get_refmodel]{get_refmodel} or \link[=init_refmodel]{init_refmodel},
#'   or an object which can be converted to a reference model using
#'   \link[=get_refmodel]{get_refmodel}.
#' @param nterms Number of terms in the submodel (the variable combination is
#'   taken from the \code{varsel} information). If a numeric vector, then the
#'   projection is performed for each model size. If \code{NULL}, the model size
#'   suggested by the variable selection (see function \code{suggest_size}).
#'   Ignored if \code{solution_terms} is specified. Note that \code{nterms} does
#'   not count the intercept, so use \code{nterms = 0} for the intercept-only
#'   model.
#' @param solution_terms Variable indices onto which the projection is done. If
#'   specified, \code{nterms} is ignored.
#' @param cv_search If TRUE, then the projected coefficients after L1-selection
#'   are computed without any penalization (or using only the regularization
#'   determined by \code{regul}). If FALSE, then the coefficients are the
#'   solution from the L1-penalized projection. This option is relevant only if
#'   L1-search was used. Default is TRUE for genuine reference models and FALSE
#'   if \code{object} is datafit (see \link[=init_refmodel]{init_refmodel}).
#' @param ndraws Number of posterior draws to be projected. Cannot be larger
#'   than the number of draws in the reference model. \strong{Caution:} For
#'   \code{ndraws <= 20}, the value of \code{ndraws} is passed to
#'   \code{nclusters} (so that clustering is used). Ignored if \code{nclusters}
#'   is not \code{NULL} or if the reference model is of class \code{"datafit"}
#'   (in which case one cluster is used). See also section "Details" below.
#' @param nclusters Number of clusters of posterior draws to be projected.
#'   Ignored if the reference model is of class \code{"datafit"} (in which case
#'   one cluster is used). For the meaning of \code{NULL}, see argument
#'   \code{ndraws}. See also section "Details" below.
#' @param seed A seed used for clustering the reference model's posterior draws
#'   (if \code{!is.null(nclusters)}). Can be used to ensure reproducible
#'   results. If \code{NULL}, no seed is set and therefore, the results are not
#'   reproducible. See \code{\link{set.seed}} for details.
#' @param regul Amount of ridge regularization when fitting the models in the
#'   projection. Usually there is no need for regularization, but sometimes for
#'   some models the projection can be ill-behaved and we need to add some
#'   regularization to avoid numerical problems.
#' @param ... Arguments passed to \link[=get_refmodel]{get_refmodel}.
#'
#' @details Using less draws or clusters in \code{ndraws} or \code{nclusters}
#'   than posterior draws in the reference model may result in slightly
#'   inaccurate projection performance. Increasing these arguments linearly
#'   affects the computation time.
#'
#' @return If the projection is performed onto a single submodel (i.e.,
#'   \code{nterms} has length one or \code{solution_terms} is specified), an
#'   object of class \code{"projection"} which is a \code{list} containing the
#'   following elements:
#'   \describe{
#'     \item{\code{dis}}{Projected draws for the dispersion parameter.}
#'     \item{\code{kl}}{The KL divergence from the submodel to the reference
#'     model.}
#'     \item{\code{weights}}{Weights for the projected draws.}
#'     \item{\code{solution_terms}}{A character vector of the submodel's
#'     predictor terms, ordered the way in which the terms were added to the
#'     submodel.}
#'     \item{\code{sub_fit}}{The submodel's fitted model object.}
#'     \item{\code{family}}{A modified \code{\link{family}}-object.}
#'     \item{\code{p_type}}{A single logical value indicating whether the
#'     reference model's posterior draws have been clustered for the projection
#'     (\code{TRUE}) or not (\code{FALSE}).}
#'     \item{\code{intercept}}{A single logical value indicating whether the
#'     reference model (as well as the submodel) contains an intercept
#'     (\code{TRUE}) or not (\code{FALSE}).}
#'     \item{\code{extract_model_data}}{The \code{extract_model_data()} function
#'     from the reference model (see \code{\link{init_refmodel}}).}
#'     \item{\code{refmodel}}{The reference model object (see
#'     \code{\link{init_refmodel}}).}
#'   }
#'   If the projection is performed onto more than one submodel, the output from
#'   above is returned for each submodel, giving a \code{list} with one element
#'   for each submodel.
#'
#' @examples
#' \donttest{
#' if (requireNamespace("rstanarm", quietly = TRUE)) {
#'   ### Usage with stanreg objects
#'   n <- 30
#'   d <- 5
#'   x <- matrix(rnorm(n * d), nrow = n)
#'   y <- x[, 1] + 0.5 * rnorm(n)
#'   data <- data.frame(x, y)
#'
#'   fit <- rstanarm::stan_glm(y ~ X1 + X2 + X3 + X4 + X5, gaussian(),
#'     data = data, chains = 2, iter = 500)
#'   vs <- varsel(fit)
#'
#'   # project onto the best model with 4 variables
#'   proj4 <- project(vs, nterms = 4)
#'
#'   # project onto an arbitrary variable combination
#'   proj <- project(fit, solution_terms = c("X1", "X3", "X5"))
#' }
#' }
#'
NULL

#' @rdname project
#' @export
project <- function(object, nterms = NULL, solution_terms = NULL,
                    cv_search = TRUE, ndraws = 400, nclusters = NULL,
                    seed = NULL, regul = 1e-4, ...) {
  if (!inherits(object, "vsel") && is.null(solution_terms)) {
    stop("The given object is not an object of class \"vsel\". Run the ",
         "variable selection first, or provide argument `solution_terms`.")
  }
  if (!inherits(object, "vsel") && !cv_search) {
    stop("The given object is not an object of class \"vsel\". Run the ",
         "variable selection first, or provide argument `cv_search = TRUE`.")
  }

  refmodel <- get_refmodel(object, ...)

  if (inherits(refmodel, "datafit")) {
    ## use non-cv_searched solution for datafits
    cv_search <- FALSE
  }

  if (!is.null(solution_terms) &&
      any(object$solution_terms[1:length(solution_terms)] != solution_terms)) {
    ## search path not found, or the given variable combination
    ## not in the search path, then we need to project the
    ## required variables
    cv_search <- TRUE
  }

  if (!is.null(solution_terms)) {
    ## if solution_terms is given, nterms is ignored
    ## (project only onto the given submodel)
    if (!is.null(object$solution_terms)) {
      vars <- object$solution_terms
    } else {
      ## project only the given model on a fit object
      vars <- setdiff(
        split_formula(refmodel$formula,
          data = refmodel$fetch_data(),
          add_main_effects = FALSE
        ),
        "1"
      )
    }
    if (length(solution_terms) > length(vars)) {
      stop(
        "Argument 'solution_terms' contains more terms than the number of ",
        "terms in the reference model."
      )
    }

    if (!all(solution_terms %in% vars)) {
      warning("At least one element of `solution_terms` could not be found ",
              "among the terms in the reference model. This element (or these ",
              "elements) is/are ignored.")
    }

    solution_terms <- intersect(solution_terms, vars)
    nterms <- length(solution_terms)
  } else {
    ## by default take the variable ordering from the selection
    solution_terms <- object$solution_terms
    if (is.null(nterms)) {
      if (!is.null(object$suggested_size) && !is.na(object$suggested_size)) {
        ## by default, project onto the suggested model size
        nterms <- min(object$suggested_size, length(solution_terms))
      } else {
        stop(
          "No suggested model size found, please specify nterms or solution",
          "terms"
        )
      }
    } else {
      if (!is.numeric(nterms) || any(nterms < 0)) {
        stop("nterms must contain non-negative values.")
      }
      if (max(nterms) > length(solution_terms)) {
        stop(paste(
          "Cannot perform the projection with", max(nterms), "variables,",
          "because variable selection was run only up to",
          length(solution_terms),
          "variables."
        ))
      }
    }
  }

  stopifnot(!is.null(ndraws))
  ndraws <- min(NCOL(refmodel$mu), ndraws)

  if (is.null(nclusters) && ndraws <= 20) {
    nclusters <- ndraws
  }
  if (!is.null(nclusters)) {
    nclusters <- min(NCOL(refmodel$mu), nclusters)
  }

  if (inherits(refmodel, "datafit")) {
    nclusters <- 1
  }

  intercept <- refmodel$intercept
  if (!intercept) {
    stop("Reference models without an intercept are currently not supported.")
  }
  family <- refmodel$family

  ## get the clustering or subsample
  p_ref <- .get_refdist(refmodel,
                        ndraws = ndraws, nclusters = nclusters, seed = seed)

  ## project onto the submodels
  subm <- .get_submodels(
    search_path = nlist(
      solution_terms,
      p_sel = object$search_path$p_sel,
      sub_fits = object$search_path$sub_fits
    ),
    nterms = nterms, family = family, p_ref = p_ref, refmodel = refmodel,
    intercept = intercept, regul = regul, cv_search = cv_search
  )
  ## add family
  proj <- lapply(subm, function(model) {
    model <- c(model, nlist(family))
    model$p_type <- !is.null(nclusters)
    model$intercept <- intercept
    model$extract_model_data <- refmodel$extract_model_data
    model$refmodel <- refmodel
    class(model) <- "projection"
    return(model)
  })
  ## If only one model size, just return the proj instead of a list of projs
  .unlist_proj(proj)
}
