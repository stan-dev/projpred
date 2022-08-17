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
#' @param refit_prj A single logical value indicating whether to fit the
#'   submodels (again) (`TRUE`) or to retrieve the fitted submodels from
#'   `object` (`FALSE`). For an `object` which is not of class `vsel`,
#'   `refit_prj` must be `TRUE`. Note that currently, `refit_prj = FALSE`
#'   requires some caution, see GitHub issues #168 and #211.
#' @param ndraws Only relevant if `refit_prj` is `TRUE`. Number of posterior
#'   draws to be projected. Ignored if `nclusters` is not `NULL` or if the
#'   reference model is of class `datafit` (in which case one cluster is used).
#'   If both (`nclusters` and `ndraws`) are `NULL`, the number of posterior
#'   draws from the reference model is used for `ndraws`. See also section
#'   "Details" below.
#' @param nclusters Only relevant if `refit_prj` is `TRUE`. Number of clusters
#'   of posterior draws to be projected. Ignored if the reference model is of
#'   class `datafit` (in which case one cluster is used). For the meaning of
#'   `NULL`, see argument `ndraws`. See also section "Details" below.
#' @param seed Pseudorandom number generation (PRNG) seed by which the same
#'   results can be obtained again if needed. Passed to argument `seed` of
#'   [set.seed()], but can also be `NA` to not call [set.seed()] at all. Here,
#'   this seed is used for clustering the reference model's posterior draws (if
#'   `!is.null(nclusters)`) and for drawing new group-level effects when
#'   predicting from a multilevel submodel (however, not yet in case of a GAMM).
#' @inheritParams varsel
#' @param ... Arguments passed to [get_refmodel()] (if [get_refmodel()] is
#'   actually used; see argument `object`) as well as to the divergence
#'   minimizer (if `refit_prj` is `TRUE`).
#'
#' @details Arguments `ndraws` and `nclusters` are automatically truncated at
#'   the number of posterior draws in the reference model (which is `1` for
#'   `datafit`s). Using less draws or clusters in `ndraws` or `nclusters` than
#'   posterior draws in the reference model may result in slightly inaccurate
#'   projection performance. Increasing these arguments affects the computation
#'   time linearly.
#'
#'   Note that if [project()] is applied to output from [cv_varsel()], then
#'   `refit_prj = FALSE` will take the results from the *full-data* search.
#'
#' @return If the projection is performed onto a single submodel (i.e.,
#'   `length(nterms) == 1 || !is.null(solution_terms)`), an object of class
#'   `projection` which is a `list` containing the following elements:
#'   \describe{
#'     \item{`dis`}{Projected draws for the dispersion parameter.}
#'     \item{`kl`}{The Kullback-Leibler (KL) divergence from the submodel to the
#'     reference model. Note that in case of the Gaussian family, this is not
#'     the actual KL divergence but merely a proxy.}
#'     \item{`weights`}{Weights for the projected draws.}
#'     \item{`solution_terms`}{A character vector of the submodel's
#'     predictor terms, ordered in the way in which the terms were added to the
#'     submodel.}
#'     \item{`submodl`}{A `list` containing the submodel fits (one fit per
#'     projected draw).}
#'     \item{`cl_ref`}{A numeric vector of length equal to the number of
#'     posterior draws in the reference model, containing the cluster indices of
#'     these draws.}
#'     \item{`wdraws_ref`}{A numeric vector of length equal to the number of
#'     posterior draws in the reference model, giving the weights of these
#'     draws. These weights should be treated as not being normalized (i.e.,
#'     they don't necessarily sum to `1`).}
#'     \item{`p_type`}{A single logical value indicating whether the
#'     reference model's posterior draws have been clustered for the projection
#'     (`TRUE`) or not (`FALSE`).}
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
                    refit_prj = TRUE, ndraws = 400, nclusters = NULL,
                    seed = sample.int(.Machine$integer.max, 1), regul = 1e-4,
                    ...) {
  if (inherits(object, "datafit")) {
    stop("project() does not support an `object` of class \"datafit\".")
  }
  if (!inherits(object, "vsel") && is.null(solution_terms)) {
    stop("Please provide an `object` of class \"vsel\" or use argument ",
         "`solution_terms`.")
  }
  if (!inherits(object, "vsel") && !refit_prj) {
    stop("Please provide an `object` of class \"vsel\" or use ",
         "`refit_prj = TRUE`.")
  }

  refmodel <- get_refmodel(object, ...)

  # Set seed, but ensure the old RNG state is restored on exit:
  if (exists(".Random.seed", envir = .GlobalEnv)) {
    rng_state_old <- get(".Random.seed", envir = .GlobalEnv)
    on.exit(assign(".Random.seed", rng_state_old, envir = .GlobalEnv))
  }
  if (!is.na(seed)) set.seed(seed)

  if (refit_prj && inherits(refmodel, "datafit")) {
    warning("Automatically setting `refit_prj` to `FALSE` since the reference ",
            "model is of class \"datafit\".")
    refit_prj <- FALSE
  }

  stopifnot(is.null(solution_terms) || is.vector(solution_terms, "character"))
  if (!refit_prj &&
      !is.null(solution_terms) &&
      any(
        solution_terms(object)[seq_along(solution_terms)] != solution_terms
      )) {
    warning("The given `solution_terms` are not part of the solution path ",
            "(from `solution_terms(object)`), so `refit_prj` is automatically ",
            "set to `TRUE`.")
    refit_prj <- TRUE
  }

  if (!refit_prj) {
    warning("Currently, `refit_prj = FALSE` requires some caution, see GitHub ",
            "issues #168 and #211.")
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
        stop("No suggested model size found, please specify `nterms` or ",
             "`solution_terms`.")
      }
    } else {
      if (!is.numeric(nterms) || any(nterms < 0)) {
        stop("Argument `nterms` must contain non-negative values.")
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

  if (inherits(refmodel, "datafit")) {
    nclusters <- 1
  }

  ## get the clustering or subsample
  p_ref <- .get_refdist(refmodel, ndraws = ndraws, nclusters = nclusters)

  ## project onto the submodels
  submodels <- .get_submodels(
    search_path = nlist(
      solution_terms,
      p_sel = object$search_path$p_sel,
      submodls = object$search_path$submodls
    ),
    nterms = nterms, p_ref = p_ref, refmodel = refmodel, regul = regul,
    refit_prj = refit_prj, ...
  )

  # Output:
  projs <- lapply(submodels, function(initsubmodl) {
    proj_k <- initsubmodl
    proj_k$p_type <- !is.null(nclusters)
    proj_k$refmodel <- refmodel
    class(proj_k) <- "projection"
    return(proj_k)
  })
  ## If only one model size, just return the proj instead of a list of projs
  return(.unlist_proj(projs))
}
