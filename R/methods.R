#' Extract draws of the linear predictor and draw from the predictive
#' distribution of the projected submodel
#'
#' \code{proj_linpred} extracts draws of the linear predictor and
#' \code{proj_predict} draws from the predictive distribution of the projected
#' submodel or submodels. If the projection has not been performed, the
#' functions also perform the projection.
#'
#' @name proj-pred
#'
#' @param object Either an object returned by \link[=varsel]{varsel},
#'   \link[=cv_varsel]{cv_varsel} or \link[=init_refmodel]{init_refmodel}, or
#'   alternatively any object that can be converted to a reference model.
#' @param newdata The predictor values used in the prediction. If
#'   \code{solution_terms} is specified, then \code{newdata} should either be a
#'   dataframe containing column names that correspond to \code{solution_terms}
#'   or a matrix with the number and order of columns corresponding to
#'   \code{solution_terms}. If \code{solution_terms} is unspecified, then
#'   \code{newdata} must either be a dataframe containing all the column names
#'   as in the original data or a matrix with the same columns at the same
#'   positions as in the original data.
#' @param ynew New (test) target variables. If given, then the log predictive
#'   density for the new observations is computed.
#' @param offsetnew Offsets for the new observations. By default a vector of
#'   zeros. By default we take the weights from newdata as in the original
#'   model. Either NULL or right hand side formula.
#' @param weightsnew Weights for the new observations. For binomial model,
#'   corresponds to the number trials per observation. For \code{proj_linpred},
#'   this argument matters only if \code{ynew} is specified. By default we take
#'   the weights from newdata as in the original model. Either NULL or right
#'   hand side formula.
#' @param transform Should the linear predictor be transformed using the
#'   inverse-link function? Default is \code{FALSE}. For \code{proj_linpred}
#'   only.
#' @param integrated If \code{TRUE}, the output is averaged over the parameters.
#'   Default is \code{FALSE}. For \code{proj_linpred} only.
#' @param nterms Number of variables in the submodel (the variable combination
#'   is taken from the variable selection information). If a vector with several
#'   values, then results for all specified model sizes are returned. Ignored if
#'   \code{solution_terms} is specified. By default use the automatically
#'   suggested model size.
#' @param draws Number of draws to return from the predictive distribution of
#'   the projection. The default is 1000. For \code{proj_predict} only.
#' @param seed An optional seed to use for drawing from the projection. For
#'   \code{proj_predict} only.
#' @param ... Additional argument passed to \link{project} if \code{object} is
#'   an object returned by \link{varsel} or \link{cv_varsel}.
#'
#' @return If the prediction is done for one submodel only (\code{nterms} has
#'   length one or \code{solution_terms} is specified) and ynew is unspecified,
#'   a matrix or vector of predictions (depending on the value of
#'   \code{integrated}). If \code{ynew} is specified, returns a list with
#'   elements pred (predictions) and lpd (log predictive densities). If the
#'   predictions are done for several submodel sizes, returns a list with one
#'   element for each submodel.
#'
#' @examples
#' \donttest{
#' ## Usage with stanreg objects
#' fit <- stan_glm(y ~ x, binomial())
#' vs <- varsel(fit)
#'
#' ## compute predictions with 4 variables at the training points
#' pred <- proj_linpred(vs, newdata = x, nterms = 4)
#' pred <- proj_predict(vs, newdata = x, nterms = 4)
#' }
#'
NULL

## The 'helper' for proj_linpred and proj_predict, ie. does all the
## functionality that is common to them. It essentially checks all the arguments
## and sets them to their respective defaults and then loops over the
## projections. For each projection, it evaluates the fun-function, which
## calculates the linear predictor if called from proj_linpred and samples from
## the predictive distribution if called from proj_predict.
proj_helper <- function(object, newdata, offsetnew, weightsnew, nterms, seed,
                        proj_predict, ...) {
  if (is.null(newdata) ||
    !(inherits(newdata, "data.frame") ||
      inherits(newdata, "matrix"))) {
    stop("newdata must be a data.frame or a matrix")
  }

  if (inherits(object, "projection") ||
    (length(object) > 0 && inherits(object[[1]], "projection"))) {
    proj <- object
  } else {
    ## reference model or varsel object obtained, so run the projection
    proj <- project(object = object, nterms = nterms, ...)
  }

  if (!.is_proj_list(proj)) {
    proj <- list(proj)
  } else {
    ## proj is not a projection object
    if (any(sapply(proj, function(x) !("family" %in% names(x))))) {
      stop(paste(
        "proj_linpred only works with objects returned by",
        " varsel, cv_varsel or project"
      ))
    }
  }

  projected_sizes <- sapply(proj, function(x) {
    if (length(x$solution_terms) > 1) {
      count_terms_chosen(x$solution_terms)
    } else {
      1
    }
  })
  nterms <- list(...)$nterms %ORifNULL% projected_sizes

  if (!all(nterms %in% projected_sizes)) {
    stop(paste0(
      "Linear prediction requested for nterms = ",
      paste(nterms, collapse = ", "),
      ", but projection performed only for nterms = ",
      paste(projected_sizes, collapse = ", "), "."
    ))
  }

  projs <- Filter(function(x) length(x$solution_terms) + 1 %in% nterms, proj)
  names(projs) <- nterms

  solution_terms <- list(...)$solution_terms
  if (!is.null(solution_terms) &&
      length(solution_terms) > NCOL(newdata)) {
    stop(paste(
      "The number of columns in newdata does not match with the given",
      "number of terms indices (solution_terms)."
    ))
  }
  ## set random seed but ensure the old RNG state is restored on exit
  rng_state_old <- rngtools::RNGseed()
  on.exit(rngtools::RNGseed(rng_state_old))
  set.seed(seed)

  preds <- lapply(projs, function(proj) {
    extract_model_data <- proj$extract_model_data
    w_o <- extract_model_data(NULL,
      newdata = newdata, weightsnew,
      offsetnew, extract_y = FALSE
    )
    weightsnew <- w_o$weights
    offsetnew <- w_o$offset
    mu <- proj$family$mu_fun(proj$sub_fit,
      newdata = newdata, offset = offsetnew,
      weights = weightsnew
    )

    proj_predict(proj, mu, weightsnew)
  })

  return(.unlist_proj(preds))
}

#' @rdname proj-pred
#' @export
proj_linpred <- function(object, newdata, ynew = NULL, offsetnew = NULL,
                         weightsnew = NULL, nterms = NULL, transform = FALSE,
                         integrated = FALSE, seed = NULL, ...) {

  ## function to perform to each projected submodel
  proj_predict <- function(proj, mu, weights) {
    pred <- t(mu)
    if (!transform) pred <- proj$family$linkfun(pred)
    if (integrated) {
      ## average over the parameters
      pred <- as.vector(proj$weights %*% pred)
    } else if (!is.null(dim(pred)) && nrow(pred) == 1) {
      ## return a vector if pred contains only one row
      pred <- as.vector(pred)
    }

    return(nlist(pred, lpd = compute_lpd(
      ynew = ynew, pred = pred, proj = proj,
      weights = weights, integrated = integrated, transform = transform
    )))
  }

  ## proj_helper lapplies fun to each projection in object
  proj_helper(
    object = object, newdata = newdata, offsetnew = offsetnew,
    weightsnew = weightsnew, nterms = nterms, seed = seed,
    proj_predict = proj_predict, ...
  )
}

compute_lpd <- function(ynew, pred, proj, weights, integrated = FALSE,
                        transform = FALSE) {
  if (!is.null(ynew)) {
    ## compute also the log-density
    target <- .get_standard_y(ynew, weights, proj$family)
    ynew <- target$y
    weights <- target$weights
    ## if !transform then we are passing linkfun(mu)
    if (!transform) pred <- proj$family$linkinv(pred)
    lpd <- proj$family$ll_fun(pred, proj$dis, ynew, weights)
    if (integrated && !is.null(dim(lpd))) {
      lpd <- as.vector(apply(lpd, 1, log_weighted_mean_exp, proj$weights))
    } else if (!is.null(dim(lpd))) {
      lpd <- t(lpd)
    }
    return(lpd)
  } else {
    return(NULL)
  }
}

#' @rdname proj-pred
#' @export
proj_predict <- function(object, newdata, offsetnew = NULL, weightsnew = NULL,
                         nterms = NULL, ndraws = 1000, seed = NULL, ...) {

  ## function to perform to each projected submodel
  proj_predict <- function(proj, mu, weights) {
    draw_inds <- sample(
      x = seq_along(proj$weights), size = ndraws,
      replace = TRUE, prob = proj$weights
    )

    t(sapply(draw_inds, function(i) {
      proj$family$ppd(mu[, i], proj$dis[i], weights)
    }))
  }

  ## proj_helper lapplies fun to each projection in object
  proj_helper(
    object = object, newdata = newdata, offsetnew = offsetnew,
    weightsnew = weightsnew, nterms = nterms, seed = seed,
    proj_predict = proj_predict, ...
  )
}

#' Plot summary statistics related to variable selection
#'
#' @inheritParams summary.vsel 
#'
#' @examples
#' \donttest{
#' ### Usage with stanreg objects
#' fit <- stan_glm(y ~ x, binomial())
#' vs <- cv_varsel(fit)
#' plot(vs)
#' }
#'
#' @method plot vsel
#' @export
plot.vsel <- function(object, nterms_max = NULL, stats = "elpd",
                      deltas = FALSE, alpha = 0.32, baseline = NULL,
                      ...) {
  .validate_vsel_object_stats(object, stats)
  baseline <- .validate_baseline(object$refmodel, baseline, deltas)

  ## compute all the statistics and fetch only those that were asked
  nfeat_baseline <- .get_nfeat_baseline(object, baseline, stats[1])
  tab <- rbind(
    .tabulate_stats(object, stats,
      alpha = alpha,
      nfeat_baseline = nfeat_baseline
    ),
    .tabulate_stats(object, stats, alpha = alpha)
  )
  stats_table <- subset(tab, tab$delta == deltas)
  stats_ref <- subset(stats_table, stats_table$size == Inf)
  stats_sub <- subset(stats_table, stats_table$size != Inf)
  stats_bs <- subset(stats_table, stats_table$size == nfeat_baseline)


  if (NROW(stats_sub) == 0) {
    stop(paste0(
      ifelse(length(stats) == 1, "Statistics ", "Statistic "),
      paste0(unique(stats), collapse = ", "), " not available."
    ))
  }

  if (is.null(nterms_max)) {
    nterms_max <- max(stats_sub$size)
  } else {
    # don't exceed the maximum submodel size
    nterms_max <- min(nterms_max, max(stats_sub$size))
    if (nterms_max < 1) {
      stop("nterms_max must be at least 1")
    }
  }
  ylab <- if (deltas) "Difference to the baseline" else "Value"

  # make sure that breaks on the x-axis are integers
  n_opts <- c(4, 5, 6)
  n_possible <- Filter(function(x) nterms_max %% x == 0, n_opts)
  n_alt <- n_opts[which.min(n_opts - (nterms_max %% n_opts))]
  nb <- ifelse(length(n_possible) > 0, min(n_possible), n_alt)
  by <- ceiling(nterms_max / min(nterms_max, nb))
  breaks <- seq(0, by * min(nterms_max, nb), by)
  minor_breaks <- if (by %% 2 == 0) {
    seq(by / 2, by * min(nterms_max, nb), by)
  } else {
    NULL
  }

  # plot submodel results
  pp <- ggplot(
    data = subset(stats_sub, stats_sub$size <= nterms_max),
    mapping = aes_string(x = "size")
  ) +
    geom_linerange(aes_string(ymin = "lq", ymax = "uq", alpha = 0.1)) +
    geom_line(aes_string(y = "value")) +
    geom_point(aes_string(y = "value"))

  if (!all(is.na(stats_ref$se))) {
    # add reference model results if they exist
    pp <- pp + geom_hline(aes_string(yintercept = "value"),
      data = stats_ref,
      color = "darkred", linetype = 2
    )
  }
  if (baseline != "ref") {
    # add the baseline result (if different from the reference model)
    pp <- pp + geom_hline(aes_string(yintercept = "value"),
      data = stats_bs,
      color = "black", linetype = 3
    )
  }
  pp <- pp +
    scale_x_continuous(
      breaks = breaks, minor_breaks = minor_breaks,
      limits = c(min(breaks), max(breaks))
    ) +
    labs(x = "Number of variables in the submodel", y = ylab) +
    theme(legend.position = "none") +
    facet_grid(statistic ~ ., scales = "free_y")
  return(pp)
}

#' Summary statistics related to variable selection
#'
#' @param object The object returned by \link[=varsel]{varsel} or
#'   \link[=cv_varsel]{cv_varsel}.
#' @param nterms_max Maximum submodel size for which the statistics are
#'   calculated. For \code{plot.vsel} it must be at least 1.
#' @param stats One or several strings determining which statistics to
#'   calculate. Available statistics are:
#' \itemize{
#'  \item{elpd:} {(Expected) sum of log predictive densities}
#'  \item{mlpd:} {Mean log predictive density, that is, elpd divided by the
#'   number of datapoints.} \item{mse:} {Mean squared error (gaussian family
#'   only)}
#'  \item{rmse:} {Root mean squared error (gaussian family only)}
#'  \item{acc/pctcorr:} {Classification accuracy (binomial family only)}
#'  \item{auc:} {Area under the ROC curve (binomial family only)}
#' }
#' Default is \code{"elpd"}.
#' @param type One or more items from 'mean', 'se', 'lower' and 'upper'
#'   indicating which of these to compute (mean, standard error, and lower and
#'   upper credible bounds). The credible bounds are determined so that
#'   \code{1-alpha} percent of the mass falls between them.
#' @param deltas If \code{TRUE}, the submodel statistics are estimated relative
#'   to the baseline model (see argument \code{baseline}) instead of estimating
#'   the actual values of the statistics. Defaults to \code{FALSE}.
#' @param alpha A number indicating the desired coverage of the credible
#'   intervals. For example \code{alpha=0.32} corresponds to 68\% probability
#'   mass within the intervals, that is, one standard error intervals.
#' @param baseline Either 'ref' or 'best' indicating whether the baseline is the
#'   reference model or the best submodel found. Default is 'ref' when the
#'   reference model exists, and 'best' otherwise.
#' @param ... Currently ignored.
#'
#' @examples
#' \donttest{
#' ### Usage with stanreg objects
#' fit <- stan_glm(y ~ x, binomial())
#' vs <- cv_varsel(fit)
#' summary(vs, stats = c("acc"), type = c("mean", "se"))
#' }
#'
#' @method summary vsel
#' @export
summary.vsel <- function(object, nterms_max = NULL, stats = "elpd",
                         type = c("mean", "se"), deltas = FALSE,
                         alpha = 0.32, baseline = NULL, ...) {
  .validate_vsel_object_stats(object, stats)
  baseline <- .validate_baseline(object$refmodel, baseline, deltas)

  ## fetch statistics
  if (deltas) {
    nfeat_baseline <- .get_nfeat_baseline(object, baseline, stats[1])
    tab <- .tabulate_stats(object, stats,
      alpha = alpha, nfeat_baseline = nfeat_baseline
    )
  } else {
    tab <- .tabulate_stats(object, stats, alpha = alpha)
  }
  stats_table <- subset(tab, tab$size != Inf)


  ## these are the corresponding names for mean, se, upper and lower in the
  ## stats_table, and their suffices in the table to be returned
  qty <- unname(sapply(type, function(t) {
    switch(t, mean = "value", upper = "uq", lower = "lq", se = "se")
  }))
  suffix <- unname(sapply(type, function(t) {
    switch(t, mean = "", upper = ".upper", lower = ".lower", se = ".se")
  }))

  ## loop through all the required statistics
  arr <- data.frame(
    size = unique(stats_table$size),
    solution_terms = c(NA, object$solution_terms)
  )
  for (i in seq_along(stats)) {
    temp <- subset(stats_table, stats_table$statistic == stats[i], qty)
    newnames <- sapply(suffix, function(s) paste0(stats[i], s))
    colnames(temp) <- newnames
    arr <- cbind(arr, temp)
  }

  if (is.null(nterms_max)) {
    nterms_max <- max(stats_table$size)
  }

  if ("pct_solution_terms_cv" %in% names(object)) {
    arr$pct_solution_terms_cv <- c(NA, diag(object$pctch[, -1]))
  }

  return(subset(arr, arr$size <= nterms_max))
}

#' Print methods for vsel/vsel objects
#'
#' The \code{print} methods for vsel/vsel objects created by
#' \code{\link{varsel}} or \code{\link{cv_varsel}}) rely on
#' \code{\link{summary.vsel}} to display a summary of the results of the
#' projection predictive variable selection.
#'
#' @name print-vsel
#'
#' @param x An object of class vsel/vsel.
#' @param digits Number of decimal places to be reported (2 by default).
#' @param ... Further arguments passed to \code{\link{summary.vsel}}.
#'
#' @return Returns invisibly the data frame produced by
#'   \code{\link{summary.vsel}}.
#'
#' @export
#' @method print vsel
print.vsel <- function(x, digits = 2, ...) {
  stats <- summary.vsel(x, ...)
  v <- match("solution_terms", colnames(stats))
  stats[, -v] <- round(stats[, -v], digits)
  print(stats[, -v])
  return(invisible(stats))
}

##' Suggest model size
##'
##' This function can be used for suggesting an appropriate model size based on
##' a certain default rule. Notice that the decision rules are heuristic and
##' should be interpreted as guidelines. It is recommended that the user studies
##' the results via \code{plot.vsel} and/or \code{summary.vsel} and makes the
##' final decision based on what is most appropriate for the given problem.
##'
##' @param object The object returned by \link[=varsel]{varsel} or
##' \link[=cv_varsel]{cv_varsel}.
##' @param stat Statistic used for the decision. Default is 'elpd'. See
##'   \code{summary.vsel} for other possible choices.
##' @param alpha A number indicating the desired coverage of the credible
##'   intervals based on which the decision is made. E.g. \code{alpha=0.32}
##'   corresponds to 68\% probability mass within the intervals (one standard
##'   error intervals). See details for more information.
##' @param pct Number indicating the relative proportion between baseline model
##'   and null model utilities one is willing to sacrifice. See details for more
##'   information. @param type Either 'upper' (default) or 'lower' determining
##'   whether the decisions are based on the upper or lower credible bounds. See
##'   details for more information.
##' @param baseline Either 'ref' or 'best' indicating whether the baseline is
##'   the reference model or the best submodel found. Default is 'ref' when the
##'   reference model exists, and 'best' otherwise.
##' @param warnings Whether to give warnings if automatic suggestion fails,
##'   mainly for internal use. Default is TRUE, and usually no reason to set to
##'   FALSE.
##' @param ... Currently ignored.
##'
##' @details The suggested model size is the smallest model for which either the
##'   lower or upper (depending on argument \code{type}) credible bound of the
##'   submodel utility \eqn{u_k} with significance level \code{alpha} falls
##'   above \deqn{u_base - pct*(u_base - u_0)} Here \eqn{u_base} denotes the
##'   utility for the baseline model and \eqn{u_0} the null model utility. The
##'   baseline is either the reference model or the best submodel found (see
##'   argument \code{baseline}). The lower and upper bounds are defined to
##'   contain the submodel utility with probability 1-alpha (each tail has mass
##'   alpha/2).
##'
##' By default \code{ratio=0}, \code{alpha=0.32} and \code{type='upper'} which
##'   means that we select the smallest model for which the upper tail exceeds
##'   the baseline model level, that is, which is better than the baseline model
##'   with probability 0.16 (and consequently, worse with probability 0.84). In
##'   other words, the estimated difference between the baseline model and
##'   submodel utilities is at most one standard error away from zero, so the
##'   two utilities are considered to be close.
##'
##' NOTE: Loss statistics like RMSE and MSE are converted to utilities by
##'   multiplying them by -1, so call such as \code{suggest_size(object,
##'   stat='rmse', type='upper')} should be interpreted as finding the smallest
##'   model whose upper credible bound of the \emph{negative} RMSE exceeds the
##'   cutoff level (or equivalently has the lower credible bound of RMSE below
##'   the cutoff level). This is done to make the interpretation of the argument
##'   \code{type} the same regardless of argument \code{stat}.
##'
##' @examples
##' \donttest{
##' ### Usage with stanreg objects
##' fit <- stan_glm(y~x, binomial())
##' vs <- cv_varsel(fit)
##' suggest_size(vs)
##' }
##' 
##' @export
suggest_size <- function(object, ...) {
  UseMethod("suggest_size")
}

#' @rdname suggest_size
#' @export
suggest_size.vsel <- function(object, stat = "elpd", alpha = 0.32, pct = 0.0,
                              type = "upper", baseline = NULL, warnings = TRUE,
                              ...) {
  .validate_vsel_object_stats(object, stat)
  if (length(stat) > 1) {
    stop("Only one statistic can be specified to suggest_size")
  }

  if (.is_util(stat)) {
    sgn <- 1
  } else {
    sgn <- -1
    if (type == "upper") {
      type <- "lower"
    } else {
      type <- "upper"
    }
  }
  bound <- paste0(stat, ".", type)
  stats <- summary.vsel(object,
    stats = stat, alpha = alpha, type = c("mean", "upper", "lower"),
    baseline = baseline, deltas = TRUE
  )
  util_null <- sgn * unlist(unname(subset(stats, stats$size == 0, stat)))
  util_cutoff <- pct * util_null
  res <- subset(stats, sgn * stats[, bound] >= util_cutoff, "size")

  if (nrow(res) == 0) {
    ## no submodel satisfying the criterion found
    if (object$nterms_max == object$nterms_all) {
      suggested_size <- object$nterms_max
    } else {
      suggested_size <- NA
      if (warnings) {
        warning(paste(
          "Could not suggest model size. Investigate plot.vsel to identify",
          "if the search was terminated too early. If this is the case,",
          "run variable selection with larger value for nterms_max."
        ))
      }
    }
  } else {
    suggested_size <- max(min(res), 1) # always include intercept
  }

  return(suggested_size)
}

##' @method as.matrix ridgelm
as.matrix.ridgelm <- function(x, ...) {
  return(coef(x))
}

##' @method as.matrix lm
as.matrix.lm <- function(x, ...) {
  return(coef(x))
}

##' @method as.matrix glm
as.matrix.glm <- function(x, ...) {
  return(coef(x))
}

##' @method as.matrix noquote
as.matrix.noquote <- function(x, ...) {
  return(coef(x))
}

##' @method as.matrix list
as.matrix.list <- function(x, ...) {
  return(do.call(cbind, lapply(x, as.matrix.glm)))
}

##' @method t glm
t.glm <- function(x, ...) {
  return(t(as.matrix(x)))
}

##' @method t lm
t.lm <- function(x, ...) {
  return(t(as.matrix(x)))
}

##' @method t ridgelm
t.ridgelm <- function(x, ...) {
  return(t(as.matrix(x)))
}

##' @method t list
t.list <- function(x, ...) {
  return(t(as.matrix.list(x)))
}

##' @method as.matrix projection
##' @export
as.matrix.projection <- function(x, ...) {
  if (x$p_type) {
    warning(paste0(
      "Note, that projection was performed using",
      "clustering and the clusters might have different weights."
    ))
  }
  if (inherits(x$sub_fit, "list")) {
    res <- t(do.call(cbind, lapply(x$sub_fit, as.matrix.lm)))
  } else {
    res <- t(as.matrix.lm(x$sub_fit))
  }
  if (x$intercept) {
    if ("1" %in% x$solution_terms) {
      colnames(res) <- gsub("^1", "Intercept", x$solution_terms)
    } else {
      colnames(res) <- c("Intercept", x$solution_terms)
    }
  }
  if (x$family$family == "gaussian") res <- cbind(res, sigma = x$dis)
  return(res)
}

##' Create cross-validation indices
##'
##' Divide indices from 1 to \code{n} into subsets for \code{k}-fold cross
##' validation. These functions are potentially useful when creating the
##' \code{cvfits} and \code{cvfun} arguments for
##' \link[=init_refmodel]{init_refmodel}. The returned value is different for
##' these two methods, see below for details.
##'
##' @name cv-indices
##'
##' @param n Number of data points.
##' @param k Number of folds. Must be at least 2 and not exceed \code{n}.
##' @param out Format of the output, either 'foldwise' (default) or 'indices'.
##'   See below for details.
##' @param seed Random seed so that the same division could be obtained again if
##'   needed.
##'
##' @return \code{cvfolds} returns a vector of length \code{n} such that each
##'   element is an integer between 1 and \code{k} denoting which fold the
##'   corresponding data point belongs to. The returned value of \code{cv_ids}
##'   depends on the \code{out}-argument. If \code{out}='foldwise', the returned
##'   value is a list with \code{k} elements, each having fields \code{tr} and
##'   \code{ts} which give the training and test indices, respectively, for the
##'   corresponding fold. If \code{out}='indices', the returned value is a list
##'   with fields \code{tr} and \code{ts} each of which is a list with \code{k}
##'   elements giving the training and test indices for each fold.
##' @examples
##' \donttest{
##' ### compute sample means within each fold
##' n <- 100
##' y <- rnorm(n)
##' cv <- cv_ids(n, k=5)
##' cvmeans <- lapply(cv, function(fold) mean(y[fold$tr]))
##' }
##'
NULL

##' @rdname cv-indices
##' @export
cvfolds <- function(n, k, seed = NULL) {
  .validate_num_folds(k, n)

  ## set random seed but ensure the old RNG state is restored on exit
  if (exists(".Random.seed")) {
    rng_state_old <- .Random.seed
    on.exit(assign(".Random.seed", rng_state_old, envir = .GlobalEnv))
  }
  set.seed(seed)

  ## create and shuffle the indices
  folds <- rep_len(1:k, length.out = n)
  folds <- sample(folds, n, replace = FALSE)

  return(folds)
}

##' @rdname cv-indices
##' @export
cv_ids <- function(n, k, out = c("foldwise", "indices"), seed = NULL) {
  .validate_num_folds(k, n)
  out <- match.arg(out)

  # set random seed but ensure the old RNG state is restored on exit
  if (exists(".Random.seed")) {
    rng_state_old <- .Random.seed
    on.exit(assign(".Random.seed", rng_state_old, envir = .GlobalEnv))
  }
  set.seed(seed)

  # shuffle the indices
  ind <- sample(1:n, n, replace = FALSE)

  if (out == "foldwise") {
    cv <- lapply(1:k, function(i) {
      ts <- sort(ind[seq(i, n, k)]) # test set
      tr <- setdiff(1:n, ts) # training set
      list(tr = tr, ts = ts)
    })
  } else if (out == "indices") {
    cv <- list()
    cv$tr <- list()
    cv$ts <- list()
    for (i in 1:k) {
      ts <- sort(ind[seq(i, n, k)]) # test set
      tr <- setdiff(1:n, ts) # training set
      cv$tr[[i]] <- tr
      cv$ts[[i]] <- ts
    }
  }

  return(cv)
}

is.vsel <- function(object) {
  inherits(object, "vsel")
}
