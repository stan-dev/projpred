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
#'   [varsel()] or [cv_varsel()]). Ignored if `!is.null(solution_terms)`. Number
#'   of terms for the submodel (the corresponding combination of predictor terms
#'   is taken from `object`). If a numeric vector, then the projection is
#'   performed for each element of this vector. If `NULL` (and
#'   `is.null(solution_terms)`), then the value suggested by [suggest_size()] is
#'   taken (with default arguments for [suggest_size()], implying that this
#'   suggested size is based on the ELPD). Note that `nterms` does not count the
#'   intercept, so use `nterms = 0` for the intercept-only model.
#' @param solution_terms If not `NULL`, then this needs to be a character vector
#'   of predictor terms for the submodel onto which the projection will be
#'   performed. Argument `nterms` is ignored in that case. For an `object` which
#'   is not of class `vsel`, `solution_terms` must not be `NULL`.
#' @param refit_prj A single logical value indicating whether to fit the
#'   submodels (again) (`TRUE`) or---if `object` is of class `vsel`---to re-use
#'   the submodel fits from the full-data search that was run when creating
#'   `object` (`FALSE`). For an `object` which is not of class `vsel`,
#'   `refit_prj` must be `TRUE`. See also section "Details" below.
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
#'   [set.seed()], but can also be `NA` to not call [set.seed()] at all. If not
#'   `NA`, then the PRNG state is reset (to the state before calling
#'   [project()]) upon exiting [project()]. Here, `seed` is used for clustering
#'   the reference model's posterior draws (if `!is.null(nclusters)`) and for
#'   drawing new group-level effects when predicting from a multilevel submodel
#'   (however, not yet in case of a GAMM) and having global option
#'   `projpred.mlvl_pred_new` set to `TRUE`. (Such a prediction takes place when
#'   calculating output elements `dis` and `ce`.)
#' @param verbose A single logical value indicating whether to print out
#'   additional information during the computations. More precisely, this gets
#'   passed as `projpred_verbose` to the divergence minimizer function of the
#'   `refmodel` object. For the built-in divergence minimizers, this only has an
#'   effect in case of sequential computations (not in case of parallel
#'   projection, which is described in [projpred-package]).
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
#'   If `refit_prj = FALSE` (which is only possible if `object` is of class
#'   `vsel`), [project()] retrieves the submodel fits from the full-data search
#'   that was run when creating `object`. Usually, the search relies on a rather
#'   coarse clustering or thinning of the reference model's posterior draws (by
#'   default, [varsel()] and [cv_varsel()] use `nclusters = 20`). Consequently,
#'   [project()] with `refit_prj = FALSE` then inherits this coarse clustering
#'   or thinning.
#'
#' @return If the projection is performed onto a single submodel (i.e.,
#'   `length(nterms) == 1 || !is.null(solution_terms)`), an object of class
#'   `projection` which is a `list` containing the following elements:
#'   \describe{
#'     \item{`dis`}{Projected draws for the dispersion parameter.}
#'     \item{`ce`}{The cross-entropy part of the Kullback-Leibler (KL)
#'     divergence from the reference model to the submodel. For some families,
#'     this is not the actual cross-entropy, but a reduced one where terms which
#'     would cancel out when calculating the KL divergence have been dropped. In
#'     case of the Gaussian family, that reduced cross-entropy is further
#'     modified, yielding merely a proxy.}
#'     \item{`wdraws_prj`}{Weights for the projected draws.}
#'     \item{`solution_terms`}{A character vector of the submodel's predictor
#'     terms.}
#'     \item{`outdmin`}{A `list` containing the submodel fits (one fit per
#'     projected draw). This is the same as the return value of the
#'     `div_minimizer` function (see [init_refmodel()]), except if [project()]
#'     was used with an `object` of class `vsel` based on an L1 search as well
#'     as with `refit_prj = FALSE`, in which case this is the output from an
#'     internal *L1-penalized* divergence minimizer.}
#'     \item{`cl_ref`}{A numeric vector of length equal to the number of
#'     posterior draws in the reference model, containing the cluster indices of
#'     these draws.}
#'     \item{`wdraws_ref`}{A numeric vector of length equal to the number of
#'     posterior draws in the reference model, giving the weights of these
#'     draws. These weights should be treated as not being normalized (i.e.,
#'     they don't necessarily sum to `1`).}
#'     \item{`const_wdraws_prj`}{A single logical value indicating whether the
#'     projected draws have constant weights (`TRUE`) or not (`FALSE`).}
#'     \item{`refmodel`}{The reference model object.}
#'   }
#'   If the projection is performed onto more than one submodel, the output from
#'   above is returned for each submodel, giving a `list` with one element for
#'   each submodel.
#'
#'   The elements of an object of class `projection` are not meant to be
#'   accessed directly but instead via helper functions (see the main vignette
#'   and [projpred-package]). An exception is element `wdraws_prj` which is
#'   currently needed to weight quantities derived from the projected draws in
#'   case of clustered projection, e.g., after applying [as.matrix.projection()]
#'   (which throws a warning in case of clustered projection to make users aware
#'   of this problem).
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
#'   # Run varsel() (here without cross-validation and with small values for
#'   # `nterms_max`, `nclusters`, and `nclusters_pred`, but only for the sake of
#'   # speed in this example; this is not recommended in general):
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
                    refit_prj = TRUE, ndraws = 400, nclusters = NULL, seed = NA,
                    verbose = getOption("projpred.verbose_project", TRUE),
                    regul = 1e-4, ...) {
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

  if (exists(".Random.seed", envir = .GlobalEnv)) {
    rng_state_old <- get(".Random.seed", envir = .GlobalEnv)
  }
  if (!is.na(seed)) {
    # Set seed, but ensure the old RNG state is restored on exit:
    if (exists(".Random.seed", envir = .GlobalEnv)) {
      on.exit(assign(".Random.seed", rng_state_old, envir = .GlobalEnv))
    }
    set.seed(seed)
  }

  if (refit_prj && inherits(refmodel, "datafit")) {
    warning("Automatically setting `refit_prj` to `FALSE` since the reference ",
            "model is of class \"datafit\".")
    refit_prj <- FALSE
  }

  stopifnot(is.null(solution_terms) || is.vector(solution_terms, "character"))
  if (!refit_prj &&
      !is.null(solution_terms) &&
      any(object$solution_terms[seq_along(solution_terms)] != solution_terms)) {
    warning("The given `solution_terms` are not part of the solution path ",
            "(from `object`), so `refit_prj` is automatically set to `TRUE`.")
    refit_prj <- TRUE
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

    if (!all(solution_terms %in% vars)) {
      warning(
        "The following element(s) of `solution_terms` could not be found in ",
        "the table of possible solution terms: `c(\"",
        paste(setdiff(solution_terms, vars), collapse = "\", \""), "\")`. ",
        "These elements are ignored. (The table of solution terms is either ",
        "`object$solution_terms` or the vector of terms in the reference ",
        "model, depending on whether `object$solution_terms` is `NULL` or ",
        "not. Here, the table of solution terms is: `c(\"",
        paste(vars, collapse = "\", \""), "\")`.)"
      )
    }

    solution_terms <- intersect(solution_terms, vars)
    nterms <- length(solution_terms)
  } else {
    ## by default take the variable ordering from the selection
    solution_terms <- object$solution_terms
    if (is.null(nterms)) {
      sgg_size <- try(suggest_size(object, warnings = FALSE), silent = TRUE)
      if (!inherits(sgg_size, "try-error") && !is.null(sgg_size) &&
          !is.na(sgg_size)) {
        ## by default, project onto the suggested model size
        nterms <- min(sgg_size, length(solution_terms))
      } else {
        stop("Could not suggest a submodel size automatically; please specify ",
             "`nterms` or `solution_terms`.")
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

  nterms_max <- max(nterms)
  nterms_all <- count_terms_in_formula(refmodel$formula) - 1L
  if (nterms_max == nterms_all &&
      formula_contains_group_terms(refmodel$formula) &&
      (refmodel$family$family == "gaussian" || refmodel$family$for_latent)) {
    warning(
      "In case of the Gaussian family (also in case of the latent projection) ",
      "and multilevel terms, the projection onto the full model can be ",
      "instable and even lead to an error, see GitHub issue #323."
    )
  }

  ## get the clustering or thinning
  if (refit_prj) {
    p_ref <- get_refdist(refmodel, ndraws = ndraws, nclusters = nclusters)
  }

  ## project onto the submodels
  submodls <- get_submodls(
    search_path = nlist(
      solution_terms,
      p_sel = object$search_path$p_sel,
      outdmins = object$search_path$outdmins
    ),
    nterms = nterms, p_ref = p_ref, refmodel = refmodel, regul = regul,
    refit_prj = refit_prj, projpred_verbose = verbose, ...
  )

  # Output:
  if (refit_prj) {
    refdist_obj <- p_ref
  } else {
    refdist_obj <- object$search_path$p_sel
  }
  projs <- lapply(submodls, function(submodl) {
    proj_k <- submodl
    proj_k$const_wdraws_prj <- length(unique(refdist_obj$wdraws_prj)) == 1
    proj_k$refmodel <- refmodel
    class(proj_k) <- "projection"
    return(proj_k)
  })
  ## If only one model size, just return the proj instead of a list of projs
  return(unlist_proj(projs))
}
