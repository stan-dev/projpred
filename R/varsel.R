#' Variable selection (without cross-validation)
#'
#' Perform the projection predictive variable selection for GLMs, GLMMs, GAMs,
#' and GAMMs. This variable selection consists of a *search* part and an
#' *evaluation* part. The search part determines the solution path, i.e., the
#' best submodel for each number of predictor terms (model size). The evaluation
#' part determines the predictive performance of the submodels along the
#' solution path.
#'
#' @param object An object of class `refmodel` (returned by [get_refmodel()] or
#'   [init_refmodel()]) or an object that can be passed to argument `object` of
#'   [get_refmodel()].
#' @param d_test For internal use only. A `list` providing information about the
#'   test set which is used for evaluating the predictive performance of the
#'   reference model. If not provided, the training set is used.
#' @param method The method for the search part. Possible options are `"L1"` for
#'   L1 search and `"forward"` for forward search. If `NULL`, then `"forward"`
#'   is used if the reference model has multilevel or additive terms and `"L1"`
#'   otherwise. See also section "Details" below.
#' @param cv_search A single logical value indicating whether to fit the
#'   submodels along the solution path again (`TRUE`) or to retrieve their fits
#'   from the search part (`FALSE`) before using those (re-)fits in the
#'   evaluation part.
#' @param ndraws Number of posterior draws used in the search part. **Caution:**
#'   For `ndraws <= 20`, the value of `ndraws` is passed to `nclusters` (so that
#'   clustering is used). Ignored if `nclusters` is not `NULL` or if `method`
#'   turns out as `"L1"` (L1 search uses always one cluster). See also section
#'   "Details" below.
#' @param nclusters Number of clusters of posterior draws used in the search
#'   part. Ignored if `method` turns out as `"L1"` (L1 search uses always one
#'   cluster). For the meaning of `NULL`, see argument `ndraws`. See also
#'   section "Details" below.
#' @param ndraws_pred Only relevant if `cv_search` is `TRUE`. Number of
#'   posterior draws used in the evaluation part. **Caution:** For `ndraws_pred
#'   <= 20`, the value of `ndraws_pred` is passed to `nclusters_pred` (so that
#'   clustering is used). Ignored if `nclusters_pred` is not `NULL`. See also
#'   section "Details" below.
#' @param nclusters_pred Only relevant if `cv_search` is `TRUE`. Number of
#'   clusters of posterior draws used in the evaluation part. For the meaning of
#'   `NULL`, see argument `ndraws_pred`. See also section "Details" below.
#' @param nterms_max Maximum number of predictor terms until which the search is
#'   continued. If `NULL`, then `min(19, D)` is used where `D` is the number of
#'   terms in the reference model (or in `search_terms`, if supplied). Note that
#'   `nterms_max` does not count the intercept, so use `nterms_max = 0` for the
#'   intercept-only model.
#' @param penalty Only relevant if `method` turns out as `"L1"`. A numeric
#'   vector determining the relative penalties or costs for the predictors. A
#'   value of `0` means that those predictors have no cost and will therefore be
#'   selected first, whereas `Inf` means those predictors will never be
#'   selected. If `NULL`, then `1` is used for each predictor.
#' @param lambda_min_ratio Only relevant if `method` turns out as `"L1"`. Ratio
#'   between the smallest and largest lambda in the L1-penalized search. This
#'   parameter essentially determines how long the search is carried out, i.e.,
#'   how large submodels are explored. No need to change this unless the program
#'   gives a warning about this.
#' @param nlambda Only relevant if `method` turns out as `"L1"`. Number of
#'   values in the lambda grid for L1-penalized search. No need to change this
#'   unless the program gives a warning about this.
#' @param thresh Only relevant if `method` turns out as `"L1"`. Convergence
#'   threshold when computing the L1 path. Usually, there is no need to change
#'   this.
#' @param regul A number giving the amount of ridge regularization when
#'   projecting onto (i.e., fitting) submodels which are GLMs. Usually there is
#'   no need for regularization, but sometimes we need to add some
#'   regularization to avoid numerical problems.
#' @param search_terms A custom character vector of terms to consider for the
#'   search. The intercept (`"1"`) needs to be included explicitly. The default
#'   considers all the terms in the reference model's formula.
#' @param verbose A single logical value indicating whether to print out
#'   additional information during the computations.
#' @param seed Pseudorandom number generation (PRNG) seed by which the same
#'   results can be obtained again if needed. If `NULL`, no seed is set and
#'   therefore, the results are not reproducible. See [set.seed()] for details.
#'   Here, this seed is used for clustering the reference model's posterior
#'   draws (if `!is.null(nclusters)`).
#' @param ... Arguments passed to [get_refmodel()].
#'
#' @details Arguments `ndraws`, `nclusters`, `nclusters_pred`, and `ndraws_pred`
#'   are automatically truncated at the number of posterior draws in the
#'   reference model (which is `1` for `datafit`s). Using less draws or clusters
#'   in `ndraws`, `nclusters`, `nclusters_pred`, or `ndraws_pred` than posterior
#'   draws in the reference model may result in slightly inaccurate projection
#'   performance. Increasing these arguments affects the computation time
#'   linearly.
#'
#'   For argument `method`, there are some restrictions: For a reference model
#'   with multilevel or additive formula terms, only the forward search is
#'   available.
#'
#'   L1 search is faster than forward search, but forward search may be more
#'   accurate. Furthermore, forward search may find a sparser model with
#'   comparable performance to that found by L1 search, but it may also start
#'   overfitting when more predictors are added.
#'
#'   An L1 search may select interaction terms before the corresponding main
#'   terms are selected. If this is undesired, choose the forward search
#'   instead.
#'
#' @return An object of class `vsel`. The elements of this object are not meant
#'   to be accessed directly but instead via helper functions (see the vignette
#'   or type `?projpred`).
#'
#' @seealso [cv_varsel()]
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
#'   # Now see, for example, `?print.vsel`, `?plot.vsel`, `?suggest_size.vsel`,
#'   # and `?solution_terms.vsel` for possible post-processing functions.
#' }
#'
#' @export
varsel <- function(object, ...) {
  UseMethod("varsel")
}

#' @rdname varsel
#' @export
varsel.default <- function(object, ...) {
  refmodel <- get_refmodel(object, ...)
  return(varsel(refmodel, ...))
}

#' @rdname varsel
#' @export
varsel.refmodel <- function(object, d_test = NULL, method = NULL,
                            ndraws = 20, nclusters = NULL, ndraws_pred = 400,
                            nclusters_pred = NULL,
                            cv_search = !inherits(object, "datafit"),
                            nterms_max = NULL, verbose = TRUE,
                            lambda_min_ratio = 1e-5, nlambda = 150,
                            thresh = 1e-6, regul = 1e-4, penalty = NULL,
                            search_terms = NULL, seed = NULL, ...) {
  refmodel <- object

  ## fetch the default arguments or replace them by the user defined values
  args <- parse_args_varsel(
    refmodel = refmodel, method = method, cv_search = cv_search,
    nterms_max = nterms_max, nclusters = nclusters, ndraws = ndraws,
    nclusters_pred = nclusters_pred, ndraws_pred = ndraws_pred,
    search_terms = search_terms
  )
  method <- args$method
  cv_search <- args$cv_search
  nterms_max <- args$nterms_max
  nclusters <- args$nclusters
  ndraws <- args$ndraws
  nclusters_pred <- args$nclusters_pred
  ndraws_pred <- args$ndraws_pred
  search_terms <- args$search_terms

  if (is.null(d_test)) {
    d_type <- "train"
    test_points <- seq_len(NROW(refmodel$y))
    d_test <- nlist(
      y = refmodel$y, test_points, data = NULL, weights = refmodel$wobs,
      type = d_type, offset = refmodel$offset
    )
  } else {
    d_type <- d_test$type
  }

  ## reference distributions for selection and prediction after selection
  p_sel <- .get_refdist(refmodel, ndraws, nclusters, seed = seed)
  p_pred <- .get_refdist(refmodel, ndraws_pred, nclusters_pred, seed = seed)

  ## perform the selection
  opt <- nlist(lambda_min_ratio, nlambda, thresh, regul)
  search_path <- select(
    method = method, p_sel = p_sel, refmodel = refmodel,
    nterms_max = nterms_max, penalty = penalty, verbose = verbose, opt = opt,
    search_terms = search_terms
  )
  solution_terms <- search_path$solution_terms

  ## statistics for the selected submodels
  submodels <- .get_submodels(search_path, c(0, seq_along(solution_terms)),
                              p_ref = p_pred, refmodel = refmodel,
                              regul = regul, cv_search = cv_search)
  sub <- .get_sub_summaries(
    submodels = submodels, test_points = seq_along(refmodel$y),
    refmodel = refmodel
  )

  ## predictive statistics of the reference model on test data. if no test data
  ## are provided,
  ## simply fetch the statistics on the train data
  if (inherits(refmodel, "datafit")) {
    ## no actual reference model, so we don't know how to predict test
    ## observations
    ntest <- NROW(refmodel$y)
    ref <- list(mu = rep(NA, ntest), lppd = rep(NA, ntest))
  } else {
    if (d_type == "train") {
      mu_test <- refmodel$family$linkinv(
        refmodel$family$linkfun(refmodel$mu) + refmodel$offset
      )
    } else {
      mu_test <- refmodel$family$linkinv(
        refmodel$ref_predfun(refmodel$fit, newdata = d_test$data) +
          d_test$offset
      )
      mu_test <- unname(mu_test)
    }
    ref <- .weighted_summary_means(
      y_test = d_test, family = refmodel$family, wsample = refmodel$wsample,
      mu = mu_test, dis = refmodel$dis
    )
  }

  ## store the relevant fields into the object to be returned
  vs <- nlist(
    refmodel,
    search_path,
    d_test,
    summaries = nlist(sub, ref),
    solution_terms = search_path$solution_terms,
    kl = sapply(submodels, function(x) x$kl),
    nterms_max,
    nterms_all = count_terms_in_formula(refmodel$formula),
    method = method,
    cv_method = NULL,
    validate_search = NULL,
    ndraws,
    ndraws_pred,
    nclusters,
    nclusters_pred
  )
  ## suggest model size
  class(vs) <- "vsel"
  vs$suggested_size <- suggest_size(vs, warnings = FALSE)
  summary <- summary(vs)
  vs$summary <- summary$selection

  return(vs)
}

select <- function(method, p_sel, refmodel, nterms_max, penalty, verbose, opt,
                   search_terms = NULL) {
  ##
  ## Auxiliary function, performs variable selection with the given method,
  ## and returns the search_path, i.e., a list with the followint entries (the
  ## last three
  ## are returned only if one cluster projection is used for selection):
  ##   solution_terms: the variable ordering
  ##   beta: coefficients along the solution path
  ##   alpha: intercepts along the solution path
  ##   p_sel: the reference distribution used in the selection (the input
  ##   argument p_sel)
  ##
  ## routine that can be used with several clusters
  if (method == "l1") {
    search_path <- search_L1(p_sel, refmodel, nterms_max - refmodel$intercept,
                             penalty, opt)
    search_path$p_sel <- p_sel
    return(search_path)
  } else if (method == "forward") {
    search_path <- search_forward(p_sel, refmodel, nterms_max, verbose, opt,
                                  search_terms = search_terms)
    search_path$p_sel <- p_sel
    return(search_path)
  }
}

parse_args_varsel <- function(refmodel, method, cv_search, nterms_max,
                              nclusters, ndraws, nclusters_pred, ndraws_pred,
                              search_terms) {
  ##
  ## Auxiliary function for parsing the input arguments for varsel.
  ## The arguments specified by the user (or the function calling this function)
  ## are treated as they are, but if some are not given, then this function
  ## fills them in with the default values. The purpose of this function is to
  ## avoid repeating the same code both in varsel and cv_varsel.
  ##
  if (is.null(search_terms)) {
    search_terms <- split_formula(refmodel$formula,
                                  data = refmodel$fetch_data())
  }
  has_group_features <- formula_contains_group_terms(refmodel$formula)
  has_additive_features <- formula_contains_additive_terms(refmodel$formula)

  if (is.null(method)) {
    if (has_group_features || has_additive_features) {
      method <- "forward"
    } else {
      method <- "l1"
    }
  } else {
    method <- tolower(method)
    if (method == "l1" && (has_group_features || has_additive_features)) {
      stop("L1 search is only supported for GLMs.")
    }
  }

  if (!(method %in% c("l1", "forward"))) {
    stop("Unknown search method")
  }

  stopifnot(!is.null(cv_search))
  if (cv_search && inherits(refmodel, "datafit")) {
    warning("For an `object` of class \"datafit\", `cv_search` is ",
            "automatically set to `FALSE`.")
    cv_search <- FALSE
  }

  stopifnot(!is.null(ndraws))
  ndraws <- min(NCOL(refmodel$mu), ndraws)

  if (is.null(nclusters) && ndraws <= 20) {
    nclusters <- ndraws
  }
  if (!is.null(nclusters)) {
    nclusters <- min(NCOL(refmodel$mu), nclusters)
  }

  if (method == "l1") {
    nclusters <- 1
  }

  stopifnot(!is.null(ndraws_pred))
  ndraws_pred <- min(NCOL(refmodel$mu), ndraws_pred)

  if (is.null(nclusters_pred) && ndraws_pred <= 20) {
    nclusters_pred <- ndraws_pred
  }
  if (!is.null(nclusters_pred)) {
    nclusters_pred <- min(NCOL(refmodel$mu), nclusters_pred)
  }

  if (!is.null(search_terms)) {
    max_nv_possible <- count_terms_chosen(search_terms, duplicates = TRUE)
  } else {
    max_nv_possible <- count_terms_in_formula(refmodel$formula)
  }
  if (is.null(nterms_max)) {
    nterms_max <- 19
  }
  nterms_max <- min(max_nv_possible, nterms_max + 1)

  return(nlist(
    method, cv_search, nterms_max, nclusters, ndraws, nclusters_pred,
    ndraws_pred, search_terms
  ))
}
