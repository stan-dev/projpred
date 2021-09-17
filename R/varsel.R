#' Variable selection for generalized linear models
#'
#' Perform the projection predictive variable selection for generalized linear
#' models, generalized linear and additive multilevel models using generic
#' reference models.
#'
#' @template args-vsel
#' @param d_test A test dataset which is used to evaluate model performance. If
#'   not provided, training data is used. Currently this argument is for
#'   internal use only.
#' @param seed Random seed used when clustering the posterior draws.
#'
#' @details Using less draws or clusters in \code{ndraws}, \code{nclusters},
#'   \code{nclusters_pred}, or \code{ndraws_pred} than posterior draws in the
#'   reference model may result in slightly inaccurate projection performance.
#'   Increasing these arguments linearly affects the computation time.
#'
#' @return An object of type \code{vsel} that contains information about the
#'   feature selection. The fields are not meant to be accessed directly by the
#'   user but instead via the helper functions (see the vignettes or type
#'   \code{?projpred} to see the main functions in the package).
#'
#' @examples
#' \donttest{
#' if (requireNamespace('rstanarm', quietly=TRUE)) {
#'   ### Usage with stanreg objects
#'   n <- 30
#'   d <- 5
#'   x <- matrix(rnorm(n*d), nrow=n)
#'   y <- x[,1] + 0.5*rnorm(n)
#'   data <- data.frame(x,y)
#'   fit <- rstanarm::stan_glm(y ~ X1 + X2 + X3 + X4 + X5, gaussian(),
#'                             data=data, chains=2, iter=500)
#'   vs <- varsel(fit)
#'   plot(vs)
#' }
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
  family <- refmodel$family

  ## fetch the default arguments or replace them by the user defined values
  ## use the intercept as indicated by the refmodel
  intercept <- NULL
  args <- parse_args_varsel(
    refmodel, method, cv_search, intercept, nterms_max,
    nclusters, ndraws, nclusters_pred, ndraws_pred, search_terms
  )
  method <- args$method
  cv_search <- args$cv_search
  intercept <- args$intercept
  nterms_max <- args$nterms_max
  nclusters <- args$nclusters
  ndraws <- args$ndraws
  nclusters_pred <- args$nclusters_pred
  ndraws_pred <- args$ndraws_pred
  search_terms <- args$search_terms
  has_group_features <- formula_contains_group_terms(refmodel$formula)
  has_additive_features <- formula_contains_additive_terms(refmodel$formula)

  if (method == "l1" && (has_group_features || has_additive_features)) {
    stop("L1 search is only supported for GLMs.")
  }

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
    family = family, intercept = intercept, nterms_max = nterms_max,
    penalty = penalty, verbose = verbose, opt = opt, search_terms = search_terms
  )
  solution_terms <- search_path$solution_terms

  ## statistics for the selected submodels
  submodels <- .get_submodels(search_path, c(0, seq_along(solution_terms)),
                              family = family, p_ref = p_pred,
                              refmodel = refmodel, intercept = intercept,
                              regul = regul, cv_search = cv_search)
  sub <- .get_sub_summaries(
    submodels = submodels, test_points = seq_along(refmodel$y),
    refmodel = refmodel, family = family
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
      mu_test <- refmodel$mu
    } else {
      mu_test <- family$linkinv(refmodel$ref_predfun(refmodel$fit,
                                                     newdata = d_test$data) +
                                  d_test$offset)
      mu_test <- unname(mu_test)
    }
    ref <- .weighted_summary_means(
      y_test = d_test, family = family, wsample = refmodel$wsample,
      mu = mu_test, dis = refmodel$dis
    )
  }

  ## store the relevant fields into the object to be returned
  vs <- nlist(
    refmodel,
    search_path,
    d_test,
    summaries = nlist(sub, ref),
    family,
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


select <- function(method, p_sel, refmodel, family, intercept, nterms_max,
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
      nterms_max - intercept, penalty, opt
    )
    search_path$p_sel <- p_sel
    return(search_path)
  } else if (method == "forward") {
    search_path <- search_forward(p_sel, refmodel, family,
                                  intercept, nterms_max, verbose, opt,
                                  search_terms = search_terms)
    search_path$p_sel <- p_sel
    return(search_path)
  }
}


parse_args_varsel <- function(refmodel, method, cv_search, intercept,
                              nterms_max, nclusters, ndraws, nclusters_pred,
                              ndraws_pred, search_terms) {
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

  if (is.null(intercept)) {
    intercept <- refmodel$intercept
  }
  if (!intercept) {
    stop("Reference models without an intercept are currently not supported.")
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
    method, cv_search, intercept, nterms_max, nclusters,
    ndraws, nclusters_pred, ndraws_pred,
    search_terms
  ))
}
