#' Variable selection for generalized linear models
#'
#' Perform the projection predictive variable selection for generalized linear
#' models, generalized linear and additive multilevel models using generic
#' reference models.
#'
#' @param object Either a \code{refmodel}-type object created by
#'   \link[=get_refmodel]{get_refmodel}, a \link[=init_refmodel]{init_refmodel},
#'   an object which can be converted to a reference model using
#'   \link[=get_refmodel]{get_refmodel} or a \code{vsel} object resulting from
#'   \code{varsel} or \code{cv_varsel}.
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
#' @param ndraws Number of posterior draws used in the variable selection.
#'   Cannot be larger than the number of draws in the reference model. Ignored
#'   if nclusters is set. Default is 10. In other words, we project a single
#'   draw from each cluster.
#' @param nclusters Number of clusters used for selection. Defaults to 10 and
#'   ignored if method='L1' (L1-search uses always one cluster). If nclusters is
#'   null we use as many clusters as draws to project.
#' @param ndraws_pred Number of projected draws used for prediction (after
#'   selection). Ignored if nclusters_pred is given. Note that setting less
#'   draws or clusters than posterior draws in the reference model may result in
#'   slightly inaccurate projection performance, although increasing this
#'   argument linearly affects the computation time. Default is 400.
#' @param nclusters_pred Number of clusters used for prediction (after
#'   selection). Default is 400. If nclusters_pred is null, we use as many
#'   clusters for prediction as ndraws_pred.
#' @param nterms_max Maximum number of varibles until which the selection is
#'   continued. Defaults to min(20, D, floor(0.4*n)) where n is the number of
#'   observations and D the number of variables.
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
#' @param search_terms A custom list of terms to evaluate for variable
#'   selection. By default considers all the terms in the reference model's
#'   formula.
#' @param seed Random seed used in the subsampling LOO. By default uses a fixed
#'   seed.
#' @param ... Additional arguments to be passed to the
#'   \code{get_refmodel}-function.
#'
#' @return An object of type \code{vsel} that contains information about the
#'   feature selection. The fields are not meant to be accessed directly by
#'   the user but instead via the helper functions (see the vignettes or type
#'   ?projpred to see the main functions in the package.)
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
#'   fit <- rstanarm::stan_glm(y ~ X1 + X2 + X3 + X4 + X5, gaussian(), data=data,
#'     chains=2, iter=500)
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
                            ndraws = NULL, nclusters = NULL, ndraws_pred = NULL,
                            nclusters_pred = NULL, cv_search = TRUE,
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

  if (method == "l1" && has_group_features) {
    stop(
      "l1 search is not supported for multilevel models",
    )
  }

  if (is.null(d_test)) {
    d_type <- "train"
    test_points <- seq_len(NROW(refmodel$y))
    d_test <- nlist(
      y = refmodel$y, test_points, data = NULL, weights = refmodel$wobs,
      type = d_type
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
  p_sub <- .get_submodels(search_path, c(0, seq_along(solution_terms)),
    family = family, p_ref = p_pred, refmodel = refmodel, intercept = intercept,
    regul = regul, cv_search = cv_search
  )
  sub <- .get_sub_summaries(
    submodels = p_sub, test_points = seq_along(refmodel$y), refmodel = refmodel,
    family = family
  )

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
      mu_test <- refmodel$mu
    } else {
      mu_test <- family$linkinv(refmodel$ref_predfun(refmodel$fit,
        newdata = d_test$data
      ))
    }
    ref <- .weighted_summary_means(
      y_test = d_test, family = family, wsample = refmodel$wsample,
      mu = mu_test, dis = refmodel$dis
    )
  }

  ## warn the user if the projection performance does not match the reference
  ## model's.
  ref_elpd <- get_stat(ref$mu, ref$lppd, d_test, family, "elpd",
    weights = ref$w
  )
  summ <- sub[[length(sub)]]
  proj_elpd <- get_stat(summ$mu, summ$lppd, d_test, family, "elpd",
    weights = summ$w
  )

  ## store the relevant fields into the object to be returned
  vs <- nlist(
    refmodel,
    search_path,
    d_test,
    summaries = nlist(sub, ref),
    family,
    solution_terms = search_path$solution_terms,
    kl = sapply(p_sub, function(x) x$kl),
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
  vs$suggested_size <- suggest_size(vs,
    warnings = FALSE
  )
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
      search_terms = search_terms
    )
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
      data = refmodel$fetch_data()
    )
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

  if (is.null(cv_search)) {
    cv_search <- !inherits(refmodel, "datafit")
  }

  if (is.null(ndraws) && is.null(nclusters)) {
    ndraws <- nclusters <- min(NCOL(refmodel$mu), 20)
  } else if (is.null(ndraws)) {
    ndraws <- nclusters <- min(NCOL(refmodel$mu), nclusters)
  } else if (is.null(nclusters)) {
    nclusters <- ndraws <- min(NCOL(refmodel$mu), ndraws)
  }

  if (method == "l1") {
    ndraws <- nclusters <- 1
  }

  if (is.null(ndraws_pred) && is.null(nclusters_pred)) {
    ## use 5 clusters for prediction by default
    ndraws_pred <- nclusters_pred <- min(NCOL(refmodel$mu), 400)
  } else if (is.null(nclusters_pred)) {
    nclusters_pred <- ndraws_pred <- min(NCOL(refmodel$mu), ndraws_pred)
  } else if (is.null(ndraws_pred)) {
    nclusters_pred <- ndraws_pred <- min(NCOL(refmodel$mu), nclusters_pred)
  }

  if (is.null(intercept)) {
    intercept <- refmodel$intercept
  }
  max_nv_possible <- count_terms_in_formula(refmodel$formula)
  if (!is.null(search_terms)) {
    max_nv_possible <- count_terms_chosen(
      search_terms, duplicates = TRUE, intercept = intercept
    )
  }
  if (is.null(nterms_max)) {
    nterms_max <- min(max_nv_possible, 20)
  } else {
    nterms_max <- min(max_nv_possible, nterms_max + as.numeric(intercept))
  }

  return(nlist(
    method, cv_search, intercept, nterms_max, nclusters,
    ndraws, nclusters_pred, ndraws_pred,
    search_terms
  ))
}
