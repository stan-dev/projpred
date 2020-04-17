#' Variable selection for generalized linear models
#'
#' Perform the projection predictive variable selection for generalized linear
#' models using generic reference models.
#'
#' @param object Either a \code{refmodel}-type object created by
#'   \link[=get_refmodel]{get_refmodel} or \link[=init_refmodel]{init_refmodel},
#'   or an object which can be converted to a reference model using
#'   \link[=get_refmodel]{get_refmodel}.
#' @param d_test A test dataset, which is used to evaluate model performance. If
#'   not provided, training data is used. Currently this argument is for
#'   internal use only.
#' @param method The method used in the variable selection. Possible options are
#'   \code{'L1'} for L1-search and \code{'forward'} for forward selection.
#'   Default is 'forward' if the number of variables in the full data is at most
#'   20,' and \code{'L1'} otherwise.
#' @param cv_search If TRUE, then the projected coefficients after L1-selection
#'   are computed without any penalization (or using only the regularization
#'   determined by \code{regul}). If FALSE, then the coefficients are the
#'   solution from the' L1-penalized projection. This option is relevant only if
#'   \code{method}='L1'. Default is TRUE for genuine reference models and FALSE
#'   if \code{object} is datafit (see \link[=init_refmodel]{init_refmodel}).
#' @param number_samples Number of posterior draws used in the variable
#'   selection. Cannot be larger than the number of draws in the reference
#'   model. Ignored if number_clusters is set.
#' @param number_clusters Number of clusters to use in the clustered projection.
#'   Overrides the \code{number_samples} argument. Defaults to 1.
#' @param number_samples_pred Number of samples used for prediction (after
#'   selection). Ignored if number_clusters_pred is given.
#' @param number_clusters_pred Number of clusters used for prediction (after
#'   selection). Default is 5.
#' @param nv_max Maximum number of varibles until which the selection is
#'   continued. Defaults to min(20, D, floor(0.4*n)) where n is the number of
#'   observations and D the number of variables.
#' @param intercept Whether to use intercept in the submodels. Defaults to TRUE.
#' @param penalty Vector determining the relative penalties or costs for the
#'   variables. Zero means that those variables have no cost and will therefore
#'   be selected first, whereas Inf means those variables will never be
#'   selected. Currently works only if method == 'L1'. By default 1 for each
#'   variable.
#' @param verbose If TRUE, may print out some information during the selection.
#'   Defaults to FALSE.
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
#' @param ... Additional arguments to be passed to the
#'   \code{get_refmodel}-function.
#'
#' @return An object of type \code{vsel} that contains information about the
#'   feature selection. The fields are not #' meant to be accessed directly by
#'   the user but instead via the helper #' functions (see the vignettes or type
#'   ?projpred #' to see the main functions in the package.)
#'
#' @examples
#' \donttest{
#' ## Usage with stanreg objects
#' fit <- stan_glm(y ~ x, binomial())
#' vs <- varsel(fit)
#' varsel_plot(vs)
#' }
#'
#' @export
varsel <- function(object, d_test = NULL, method = NULL, number_samples = NULL,
                   number_clusters = NULL, number_samples_pred = NULL,
                   number_clusters_pred = NULL, cv_search = FALSE,
                   nv_max = NULL, intercept = TRUE, verbose = TRUE,
                   lambda_min_ratio = 1e-5, nlambda = 150, thresh = 1e-6,
                   regul = 1e-4, penalty = NULL, search_terms = NULL, ...) {
  refmodel <- get_refmodel(object, ...)
  family <- refmodel$family

  ## fetch the default arguments or replace them by the user defined values
  args <- parse_args_varsel(
    refmodel, method, cv_search, intercept, nv_max,
    number_clusters, number_samples,
    number_clusters_pred, number_samples_pred, search_terms
  )
  method <- args$method
  cv_search <- args$cv_search
  intercept <- args$intercept
  nv_max <- args$nv_max
  number_clusters <- args$number_clusters
  number_samples <- args$number_samples
  number_clusters_pred <- args$number_clusters_pred
  number_samples_pred <- args$number_samples_pred
  search_terms <- args$search_terms
  has_group_features <- formula_contains_group_terms(refmodel$formula)

  if (method == "l1" && has_group_features) {
    stop(
      "l1 search is not supported for multilevel models",
    )
  }

  if (is.null(d_test)) {
    d_type <- "train"
    test_points <- seq_along(NROW(refmodel$y))
    d_test <- nlist(
      y = refmodel$y, test_points,
      data = NULL, weights = refmodel$wobs,
      type = d_type
    )
  }

  ## reference distributions for selection and prediction after selection
  p_sel <- .get_refdist(refmodel, number_samples, number_clusters)
  p_pred <- .get_refdist(refmodel, number_samples_pred, number_clusters_pred)

  ## perform the selection
  opt <- nlist(lambda_min_ratio, nlambda, thresh, regul)
  search_path <- select(method, p_sel, refmodel, family, intercept, nv_max,
    penalty, verbose, opt,
    search_terms = search_terms
  )
  solution_terms <- search_path$solution_terms

  ## statistics for the selected submodels
  p_sub <- .get_submodels(search_path, c(0, seq_along(solution_terms)), family, p_pred,
    refmodel, intercept, regul,
    cv_search = cv_search
  )
  sub <- .get_sub_summaries(p_sub, seq_along(refmodel$y), refmodel, family)

  ## predictive statistics of the reference model on test data. if no test data
  ## are provided,
  ## simply fetch the statistics on the train data
  if ("datafit" %in% class(refmodel)) {
    ## no actual reference model, so we don't know how to predict test
    ## observations
    ntest <- NROW(refmodel$y)
    ref <- list(mu = rep(NA, ntest), lppd = rep(NA, ntest))
  } else {
    if (d_type == "train") {
      ref <- .weighted_summary_means(
        d_test, family, refmodel$wsample,
        refmodel$mu, refmodel$dis
      )
    } else {
      mu_test <- refmodel$predfun(refmodel$fit, newdata = d_test$data)
      ref <- .weighted_summary_means(
        d_test, family, refmodel$wsample,
        mu_test, refmodel$dis
      )
    }
  }

  ## store the relevant fields into the object to be returned
  vs <- nlist(
    refmodel,
    search_path,
    d_test,
    summaries = nlist(sub, ref),
    family,
    solution_terms = search_path$solution_terms,
    kl = sapply(p_sub, function(x) x$kl),
    nv_max,
    nv_all = count_terms_in_subformula(refmodel$formula)
  )
  ## suggest model size
  class(vs) <- "vsel"
  vs$suggested_size <- suggest_size(vs,
    warnings = FALSE,
    has_group_features = has_group_features,
    search_terms = search_terms
  )

  return(vs)
}


select <- function(method, p_sel, refmodel, family, intercept, nv_max,
                   penalty, verbose, opt, search_terms = NULL) {
  ##
  ## Auxiliary function, performs variable selection with the given method,
  ## and returns the search_path, i.e., a list with the followint entries (the
  ## last three
  ## are returned only if one cluster projection is used for selection):
  ##   solution_terms: the variable ordering
  ##   beta: coefficients along the search path
  ##   alpha: intercepts along the search path
  ##   p_sel: the reference distribution used in the selection (the input
  ##   argument p_sel)
  ##
  ## routine that can be used with several clusters
  if (method == "l1") {
    search_path <- search_L1(
      p_sel, refmodel, family, intercept,
      nv_max - intercept, penalty, opt
    )
    search_path$p_sel <- p_sel
    return(search_path)
  } else if (method == "forward") {
    search_path <- search_forward(p_sel, refmodel, family,
      intercept, nv_max, verbose, opt, search_terms = search_terms
    )
    search_path$p_sel <- p_sel
    return(search_path)
  }
}


parse_args_varsel <- function(refmodel, method, cv_search, intercept, nv_max,
                              number_clusters, number_samples,
                              number_clusters_pred, number_samples_pred,
                              search_terms) {
  ##
  ## Auxiliary function for parsing the input arguments for varsel.
  ## The arguments specified by the user (or the function calling this function)
  ## are treated as they are, but if some are not given, then this function
  ## fills them in with the default values. The purpose of this function is to
  ## avoid repeating the same code both in varsel and cv_varsel.
  ##
  if (is.null(search_terms)) {
    search_terms <- split_formula(refmodel$formula)
  }
  has_group_features <- formula_contains_group_terms(refmodel$formula)

  if (is.null(method)) {
    if (has_group_features) {
      method <- "forward"
    } else {
      method <- "l1"
    }
  } else {
    method <- tolower(method)
  }

  if (!(method %in% c("l1", "forward"))) {
    stop("Unknown search method")
  }

  if (is.null(cv_search)) {
    cv_search <- !inherits(refmodel, "datafit")
  }

  if ((is.null(number_samples) && is.null(number_clusters)) || method == "l1") {
    ## use one cluster for selection by default, and always with L1-search
    number_clusters <- 1
  }
  if (is.null(number_samples_pred) && is.null(number_clusters_pred)) {
    ## use 5 clusters for prediction by default
    number_clusters_pred <- min(NCOL(refmodel$mu), 5)
  }

  max_nv_possible <- count_terms_in_subformula(refmodel$formula)
  if (is.null(intercept)) {
    intercept <- refmodel$intercept
  }
  if (is.null(nv_max)) {
    nv_max <- min(max_nv_possible, 20)
  } else {
    nv_max <- min(max_nv_possible, nv_max + 1)
  }

  return(nlist(method, cv_search, intercept, nv_max, number_clusters,
               number_samples, number_clusters_pred, number_samples_pred,
               search_terms))
}
