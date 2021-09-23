#' Projection onto submodel(s)
#'
#' Project the posterior of the reference model onto the parameter space of a
#' single submodel consisting of a specific combination of predictor terms or
#' (after variable selection) onto the parameter space of a single or multiple
#' submodels of specific sizes.
#'
#' @param object An object which can be used as input to [get_refmodel()] (in
#'   particular, objects of class `refmodel`).
#' @param nterms Only relevant if `object` is of class `vsel` (returned by
#'   [varsel()] or [cv_varsel()]). Ignored if `!is.null(solution_terms)`.
#'   Number of terms for the submodel (the corresponding combination of
#'   predictor terms is taken from `object`). If a numeric vector, then the
#'   projection is performed for each element of this vector. If `NULL` (and
#'   `is.null(solution_terms)`), then the value suggested by the variable
#'   selection is taken (see function [suggest_size()]). Note that `nterms` does
#'   not count the intercept, so use `nterms = 0` for the intercept-only model.
#' @param solution_terms If not `NULL`, then this needs to be a character vector
#'   of predictor terms for the submodel onto which the projection will be
#'   performed. Argument `nterms` is ignored in that case. For an `object` which
#'   is not of class `vsel`, `solution_terms` must not be `NULL`.
#' @param cv_search A single logical value indicating whether to fit the
#'   submodels (again) (`TRUE`) or to retrieve the fitted submodels from
#'   `object` (`FALSE`). For an `object` which is not of class `vsel`,
#'   `cv_search` must be `TRUE`.
#' @param ndraws Only relevant if `cv_search` is `TRUE`. Number of posterior
#'   draws to be projected. **Caution:** For `ndraws <= 20`, the value of
#'   `ndraws` is passed to `nclusters` (so that clustering is used). Ignored if
#'   `nclusters` is not `NULL` or if the reference model is of class `datafit`
#'   (in which case one cluster is used). See also section "Details" below.
#' @param nclusters Only relevant if `cv_search` is `TRUE`. Number of clusters
#'   of posterior draws to be projected. Ignored if the reference model is of
#'   class `datafit` (in which case one cluster is used). For the meaning of
#'   `NULL`, see argument `ndraws`. See also section "Details" below.
#' @param seed Pseudorandom number generation (PRNG) seed by which the same
#'   results can be obtained again if needed. If `NULL`, no seed is set and
#'   therefore, the results are not reproducible. See [set.seed()] for details.
#'   Here, this seed is used for clustering the reference model's posterior
#'   draws (if `!is.null(nclusters)`).
#' @inheritParams varsel
#' @param ... Arguments passed to [get_refmodel()] (if [get_refmodel()] is
#'   actually used; see argument `object`).
#'
#' @details
#'
#' Notes:
#' * Arguments `ndraws` and `nclusters` are automatically truncated at the
#' number of posterior draws in the reference model (which is `1` for
#' `datafit`s).
#' * Using less draws or clusters in `ndraws` or `nclusters` than posterior
#' draws in the reference model may result in slightly inaccurate projection
#' performance. Increasing these arguments affects the computation time
#' linearly.
#'
#' @return If the projection is performed onto a single submodel (i.e.,
#'   `length(nterms) == 1 || !is.null(solution_terms)`), an object of class
#'   `projection` which is a `list` containing the following elements:
#'   \describe{
#'     \item{`dis`}{Projected draws for the dispersion parameter.}
#'     \item{`kl`}{The KL divergence from the submodel to the reference
#'     model.}
#'     \item{`weights`}{Weights for the projected draws.}
#'     \item{`solution_terms`}{A character vector of the submodel's
#'     predictor terms, ordered in the way in which the terms were added to the
#'     submodel.}
#'     \item{`sub_fit`}{The submodel's fitted model object.}
#'     \item{`family`}{A modified [`family`] object.}
#'     \item{`p_type`}{A single logical value indicating whether the
#'     reference model's posterior draws have been clustered for the projection
#'     (`TRUE`) or not (`FALSE`).}
#'     \item{`intercept`}{A single logical value indicating whether the
#'     reference model (as well as the submodel) contains an intercept
#'     (`TRUE`) or not (`FALSE`).}
#'     \item{`extract_model_data`}{The `extract_model_data` function
#'     from the reference model (see [init_refmodel()]).}
#'     \item{`refmodel`}{The reference model object.}
#'   }
#'   If the projection is performed onto more than one submodel, the output from
#'   above is returned for each submodel, giving a `list` with one element for
#'   each submodel.
#'
#' @examples
#' if (requireNamespace("rstanarm", quietly = TRUE)) {
#'   # Data:
#'   dat_gauss <- data.frame(y = df_gaussian$y, df_gaussian$x)
#'
#'   # The "stanreg" fit which will be used as the reference model (with small
#'   # values for `chains` and `iter`, but only for technical reasons in this
#'   # example; this is not recommended in general):
#'   fit <- rstanarm::stan_glm(
#'     y ~ X1 + X2 + X3 + X4 + X5, family = gaussian(), data = dat_gauss,
#'     QR = TRUE, chains = 2, iter = 500, refresh = 0, seed = 9876
#'   )
#'
#'   # Variable selection (here without cross-validation and with small values
#'   # for `nterms_max`, `nclusters`, and `nclusters_pred`, but only for the
#'   # sake of speed in this example; this is not recommended in general):
#'   vs <- varsel(fit, nterms_max = 3, nclusters = 5, nclusters_pred = 10,
#'                seed = 5555)
#'
#'   # Projection onto the best submodel with 2 predictor terms (with a small
#'   # value for `nclusters`, but only for the sake of speed in this example;
#'   # this is not recommended in general):
#'   prj_from_vs <- project(vs, nterms = 2, nclusters = 10, seed = 9182)
#'
#'   # Projection onto an arbitrary combination of predictor terms (with a small
#'   # value for `nclusters`, but only for the sake of speed in this example;
#'   # this is not recommended in general):
#'   prj <- project(fit, solution_terms = c("X1", "X3", "X5"), nclusters = 10,
#'                  seed = 9182)
#' }
#'
#' @export
project <- function(object, nterms = NULL, solution_terms = NULL,
                    cv_search = TRUE, ndraws = 400, nclusters = NULL,
                    seed = NULL, regul = 1e-4, ...) {
  if (inherits(object, "datafit")) {
    stop("project() does not support an `object` of class \"datafit\".")
  }
  if (!inherits(object, "vsel") && is.null(solution_terms)) {
    stop("Please provide an `object` of class \"vsel\" or use argument ",
         "`solution_terms`.")
  }
  if (!inherits(object, "vsel") && !cv_search) {
    stop("Please provide an `object` of class \"vsel\" or use ",
         "`cv_search = TRUE`.")
  }

  refmodel <- get_refmodel(object, ...)

  if (cv_search && inherits(refmodel, "datafit")) {
    warning("Automatically setting `cv_search` to `FALSE` since the reference ",
            "model is of class \"datafit\".")
    cv_search <- FALSE
  }

  if (!cv_search &&
      !is.null(solution_terms) &&
      any(
        solution_terms(object)[seq_along(solution_terms)] != solution_terms
      )) {
    warning("The given `solution_terms` are not part of the search path (from ",
            "`solution_terms(object)`), so `cv_search` is automatically set ",
            "to `TRUE`.")
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
                      add_main_effects = FALSE),
        "1"
      )
    }
    if (length(solution_terms) > length(vars)) {
      stop("Argument 'solution_terms' contains more terms than the number of ",
           "terms in the reference model.")
    }

    if (!all(solution_terms %in% vars)) {
      warning(
        "At least one element of `solution_terms` could not be found in the ",
        "table of solution terms (which is either `object$solution_terms` or ",
        "the vector of terms in the reference model, depending on whether ",
        "`object$solution_terms` is `NULL` or not). Elements which cannot be ",
        "found are ignored. The table of solution terms is here: `c(\"",
        paste(vars, collapse = "\", \""), "\")`."
      )
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
        stop("No suggested model size found, please specify nterms or solution",
             "terms")
      }
    } else {
      if (!is.numeric(nterms) || any(nterms < 0)) {
        stop("nterms must contain non-negative values.")
      }
      if (max(nterms) > length(solution_terms)) {
        stop(paste(
          "Cannot perform the projection with", max(nterms), "variables,",
          "because variable selection was run only up to",
          length(solution_terms), "variables."
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
