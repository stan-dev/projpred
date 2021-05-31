.onAttach <- function(...) {
  ver <- utils::packageVersion("projpred")
  msg <- paste0("This is projpred version ", ver, ".")
  packageStartupMessage(msg)
}

weighted.sd <- function(x, w, na.rm = FALSE) {
  if (na.rm) {
    ind <- !is.na(w) & !is.na(x)
    n <- sum(ind)
  } else {
    n <- length(x)
    ind <- rep(TRUE, n)
  }
  w <- w / sum(w[ind])
  m <- sum(x[ind] * w[ind])
  sqrt(n / (n - 1) * sum(w[ind] * (x[ind] - m)^2))
}

log_weighted_mean_exp <- function(x, w) {
  x <- x + log(w)
  max_x <- max(x)
  max_x + log(sum(exp(x - max_x)))
}

log_sum_exp <- function(x) {
  max_x <- max(x)
  max_x + log(sum(exp(x - max_x)))
}

auc <- function(x) {
  resp <- x[, 1]
  pred <- x[, 2]
  weights <- x[, 3]
  n <- nrow(x)
  ord <- order(pred, decreasing = TRUE)
  resp <- resp[ord]
  pred <- pred[ord]
  weights <- weights[ord]
  w0 <- w1 <- weights
  w0[resp == 1] <- 0 # true negative weights
  w1[resp == 0] <- 0 # true positive weights
  cum_w0 <- cumsum(w0)
  cum_w1 <- cumsum(w1)

  ## ignore tied predicted probabilities, keeping only the rightmost one
  rightmost.prob <- c(diff(pred) != 0, TRUE)
  fpr <- c(0, cum_w0[rightmost.prob]) / cum_w0[n]
  tpr <- c(0, cum_w1[rightmost.prob]) / cum_w1[n]
  delta_fpr <- c(diff(fpr), 0)
  delta_tpr <- c(diff(tpr), 0)

  ## sum the area of the rectangles that fall completely below the ROC curve
  ## plus half the area of the rectangles that are cut in two by the curve
  return(sum(delta_fpr * tpr) + sum(delta_fpr * delta_tpr) / 2)
}

bootstrap <- function(x, fun = mean, b = 1000, oobfun = NULL, seed = NULL,
                      ...) {
  #
  # bootstrap an arbitrary quantity fun that takes the sample x
  # as the first input. other parameters to fun can be passed in as ...
  # example: boostrap(x,mean)
  #

  # set random seed but ensure the old RNG state is restored on exit
  if (exists(".Random.seed")) {
    rng_state_old <- .Random.seed
    on.exit(assign(".Random.seed", rng_state_old, envir = .GlobalEnv))
  }
  set.seed(seed)

  seq_x <- seq.int(NROW(x))
  is_vector <- NCOL(x) == 1
  bsstat <- rep(NA, b)
  oobstat <- rep(NA, b)
  for (i in 1:b) {
    bsind <- sample(seq_x, replace = TRUE)
    bsstat[i] <- fun(if (is_vector) x[bsind] else x[bsind, ], ...)
    if (!is.null(oobfun)) {
      oobind <- setdiff(seq_x, unique(bsind))
      oobstat[i] <- oobfun(if (is_vector) x[oobind] else x[oobind, ], ...)
    }
  }
  if (!is.null(oobfun)) {
    return(list(bs = bsstat, oob = oobstat))
  } else {
    return(bsstat)
  }
}

# from rstanarm
`%ORifNULL%` <- function(a, b) if (is.null(a)) b else a

.is.wholenumber <- function(x) abs(x - round(x)) < .Machine$double.eps^0.5

.validate_num_folds <- function(k, n) {
  if (!is.numeric(k) || length(k) != 1 || !.is.wholenumber(k)) {
    stop("Number of folds must be a single integer value.")
  }
  if (k < 2) {
    stop("Number of folds must be at least 2.")
  }
  if (k > n) {
    stop("Number of folds cannot exceed n.")
  }
}

.validate_vsel_object_stats <- function(object, stats) {
  if (!inherits(object, c("vsel"))) {
    stop(
      "The object is not a variable selection object. ",
      "Run variable selection first"
    )
  }

  recognized_stats <- c(
    "elpd", "mlpd", "mse", "rmse", "acc",
    "pctcorr", "auc"
  )
  binomial_only_stats <- c("acc", "pctcorr", "auc")
  family <- object$family$family

  if (is.null(stats)) {
    stop("Statistic specified as NULL.")
  }
  for (stat in stats) {
    if (!(stat %in% recognized_stats)) {
      stop(sprintf("Statistic '%s' not recognized.", stat))
    }
    if (stat %in% binomial_only_stats && family != "binomial") {
      stop("Statistic '", stat, "' available only for the binomial family.")
    }
  }
}

.validate_baseline <- function(refmodel, baseline, deltas) {
  if (is.null(baseline)) {
    if (inherits(refmodel, "datafit")) {
      baseline <- "best"
    } else {
      baseline <- "ref"
    }
  } else {
    if (!(baseline %in% c("ref", "best"))) {
      stop("Argument 'baseline' must be either 'ref' or 'best'.")
    }
    if (baseline == "ref" && deltas == TRUE && inherits(refmodel, "datafit")) {
      # no reference model (or the results missing for some other reason),
      # so cannot compute differences between the reference model and submodels
      stop(
        "Cannot use deltas = TRUE and baseline = 'ref' when there is no ",
        "reference model."
      )
    }
  }
  return(baseline)
}

.get_standard_y <- function(y, weights, fam) {
  # return y and the corresponding observation weights into the 'standard' form:
  # for binomial family, y is transformed into a vector with values between 0
  # and 1, and weights give the number of observations at each x. for all other
  # families, y and weights are kept as they are (unless weights is a vector
  # with length zero in which case it is replaced by a vector of ones).
  if (NCOL(y) == 1) {
    if (length(weights) > 0) {
      weights <- unname(weights)
    } else {
      weights <- rep(1, length(y))
    }
    if (fam$family == "binomial") {
      if (is.factor(y)) {
        if (nlevels(y) > 2) {
          stop("y cannot contain more than two classes if specified as factor.")
        }
        y <- as.vector(y, mode = "integer") - 1 # zero-one vector
      }
    } else {
      if (is.factor(y)) {
        stop("y cannot be a factor for models other than the binomial model.")
      }
    }
  } else if (NCOL(y) == 2) {
    weights <- y[, 2]
    y <- y[, 1]
  } else {
    stop("y cannot have more than two columns.")
  }
  return(nlist(y, weights))
}

.get_refdist <- function(refmodel, ndraws = NULL, nclusters = NULL, seed = NULL,
                         thinning = TRUE) {
  #
  # Creates the reference distribution based on the refmodel-object, and the
  # desired number of clusters (nclusters) or number of subsamples (ndraws). If
  # nclusters is specified, then clustering is used and ndraws is ignored.
  # Returns a list with fields:
  #
  #   mu: n-by-s matrix, vector of expected values for y for each draw/cluster.
  #       here s means either the number of draws ndraws or clusters nclusters
  #       used, depending on which one is used.
  #   var: n-by-s matrix, vector of predictive variances for y for each
  #         draw/cluster which which are needed for projecting the dispersion
  #         parameter (note that this can be unintuitively zero for those
  #         families that do not have dispersion) weights: s-element vector of
  #   weights for the draws/clusters
  #   cl: cluster assignment for each posterior draw, that is, a vector that has
  #       length equal to the number of posterior draws and each value is an
  #       integer between 1 and s
  # set random seed but ensure the old RNG state is restored on exit
  if (exists(".Random.seed")) {
    rng_state_old <- .Random.seed
    on.exit(assign(".Random.seed", rng_state_old, envir = .GlobalEnv))
  }
  set.seed(seed)

  family <- refmodel$family
  S <- NCOL(refmodel$mu) # number of draws in the reference model

  if (is.null(ndraws)) {
    ndraws <- S
  }

  if (!is.null(nclusters)) {
    # use clustering (ignore ndraws argument)
    if (nclusters == 1) {
      # special case, only one cluster
      cl <- rep(1, S)
      p_ref <- .get_p_clust(family, refmodel$mu, refmodel$dis,
                            wobs = refmodel$wobs, cl = cl
      )
    } else if (nclusters == NCOL(refmodel$mu)) {
      # number of clusters equal to the number of samples, so return the samples
      return(.get_refdist(refmodel, ndraws = nclusters, seed = seed))
    } else {
      # several clusters
      if (nclusters > NCOL(refmodel$mu)) {
        stop(
          "The number of clusters nclusters cannot exceed the number of ",
          "columns in mu."
        )
      }
      p_ref <- .get_p_clust(family, refmodel$mu, refmodel$dis,
                            wobs = refmodel$wobs, nclusters = nclusters
      )
    }
  } else {
    # Perform thinning or subsampling (depending on argument `thinning`) from
    # the reference model.
    if (ndraws > NCOL(refmodel$mu)) {
      stop(
        "The number of draws ndraws cannot exceed the number of ",
        "columns in mu."
      )
    }
    if (thinning) {
      s_ind <- round(seq(from = 1, to = S, length.out = ndraws))
    } else {
      s_ind <- sample(seq_len(S), size = ndraws)
    }
    cl <- rep(NA, S)
    cl[s_ind] <- c(1:ndraws)
    predvar <- sapply(s_ind, function(j) {
      family$predvar(refmodel$mu[, j, drop = FALSE], refmodel$dis[j])
    })
    p_ref <- list(
      mu = refmodel$mu[, s_ind, drop = FALSE], var = predvar,
      dis = refmodel$dis[s_ind], weights = rep(1 / ndraws, ndraws),
      cl = cl
    )
  }

  return(p_ref)
}

.get_p_clust <- function(family, mu, dis, nclusters = 10,
                         wobs = rep(1, dim(mu)[1]),
                         wsample = rep(1, dim(mu)[2]), cl = NULL) {
  # Function for perfoming the clustering over the samples.
  #
  # cluster the samples in the latent space if no clustering provided
  if (is.null(cl)) {
    f <- family$linkfun(mu)
    out <- kmeans(t(f), nclusters, iter.max = 50)
    cl <- out$cluster # cluster indices for each sample
  } else if (typeof(cl) == "list") {
    # old clustering solution provided, so fetch the cluster indices
    if (is.null(cl$cluster)) {
      stop(
        "argument cl must be a vector of cluster indices or a clustering ",
        "object returned by k-means."
      )
    }
    cl <- cl$cluster
  }

  # (re)compute the cluster centers, because they may be different from the ones
  # returned by kmeans if the samples have differing weights
  # number of clusters (assumes labeling 1,...,nclusters)
  nclusters <- max(cl, na.rm = TRUE)
  centers <- matrix(0, nrow = nclusters, ncol = dim(mu)[1])
  wcluster <- rep(0, nclusters) # cluster weights
  eps <- 1e-10
  for (j in 1:nclusters) {
    # compute normalized weights within the cluster, 1-eps is for numerical
    # stability
    ind <- which(cl == j)
    ws <- wsample[ind] / sum(wsample[ind]) * (1 - eps)

    # cluster centers and their weights
    centers[j, ] <- mu[, ind, drop = FALSE] %*% ws
    wcluster[j] <- sum(wsample[ind]) # unnormalized weight for the jth cluster
  }
  wcluster <- wcluster / sum(wcluster)

  # predictive variances
  predvar <- sapply(1:nclusters, function(j) {
    # compute normalized weights within the cluster, 1-eps is for numerical
    # stability
    ind <- which(cl == j)
    ws <- wsample[ind] / sum(wsample[ind]) * (1 - eps)
    family$predvar(mu[, ind, drop = FALSE], dis[ind], ws)
  })

  # combine the results
  p <- list(
    mu = unname(t(centers)),
    var = predvar,
    weights = wcluster,
    cl = cl
  )
  return(p)
}

.is_proj_list <- function(proj) {
  !("family" %in% names(proj))
}

.unlist_proj <- function(p) if (length(p) == 1) p[[1]] else p

## create a named list using object names
nlist <- function(...) {
  m <- match.call()
  dots <- list(...)
  no_names <- is.null(names(dots))
  has_name <- if (no_names) FALSE else nzchar(names(dots))
  if (all(has_name)) {
    return(dots)
  }
  nms <- as.character(m)[-1]
  if (no_names) {
    names(dots) <- nms
  } else {
    names(dots)[!has_name] <- nms[!has_name]
  }
  dots
}

## ifelse operator
"%||%" <- function(x, y) {
  if (is.null(x)) x <- y
  x
}

#' Execute a Function Call
#'
#' Execute a function call similar to \code{\link{do.call}}, but without
#' deparsing function arguments.
#'
#' @param what Either a function or a non-empty character string naming the
#'   function to be called.
#' @param args A list of arguments to the function call. The names attribute of
#'   \code{args} gives the argument names.
#' @param pkg Optional name of the package in which to search for the
#'   function if \code{what} is a character string.
#'
#' @return The result of the (evaluated) function call.
#'
#' @keywords internal
#' @export
do_call <- function(what, args, pkg = NULL) {
  call <- ""
  if (length(args)) {
    if (!is.list(args)) {
      stop2("'args' must be a list.")
    }
    fun_args <- names(args)
    if (is.null(fun_args)) {
      fun_args <- rep("", length(args))
    } else {
      nzc <- nzchar(fun_args)
      fun_args[nzc] <- paste0("`", fun_args[nzc], "` = ")
    }
    names(args) <- paste0(".x", seq_along(args))
    call <- paste0(fun_args, names(args), collapse = ",")
  } else {
    args <- list()
  }
  if (is.function(what)) {
    args$.fun <- what
    what <- ".fun"
  } else {
    what <- paste0("`", as_one_character(what), "`")
    if (!is.null(pkg)) {
      what <- paste0(as_one_character(pkg), "::", what)
    }
  }
  call <- paste0(what, "(", call, ")")
  eval2(call, envir = args, enclos = parent.frame())
}

# like 'eval' but parses characters before evaluation
eval2 <- function(expr, envir = parent.frame(), ...) {
  if (is.character(expr)) {
    expr <- parse(text = expr)
  }
  eval(expr, envir, ...)
}

# coerce 'x' to a single character string
as_one_character <- function(x, allow_na = FALSE) {
  s <- substitute(x)
  x <- as.character(x)
  if (length(x) != 1L || anyNA(x) && !allow_na) {
    s <- deparse_combine(s, max_char = 100L)
    stop2("Cannot coerce '", s, "' to a single character value.")
  }
  x
}

stop2 <- function(...) {
  stop(..., call. = FALSE)
}

# combine deparse lines into one string
deparse_combine <- function(x, max_char = NULL) {
  out <- paste0(deparse(x), collapse = "")
  if (isTRUE(max_char > 0)) {
    out <- substr(out, 1L, max_char)
  }
  out
}

#' @importFrom magrittr %>%
#' @export
magrittr::`%>%`

`%:::%` <- function(pkg, fun) {
  # Note: `utils::getFromNamespace(fun, pkg)` could probably be used, too (but
  # its documentation is unclear about the inheritance from parent
  # environments).
  get(fun, envir = asNamespace(pkg), inherits = FALSE)
}

# Function where() is not exported by package tidyselect, so redefine it here to
# avoid a note in R CMD check which would occur for usage of
# tidyselect:::where():
where <- `%:::%`("tidyselect", "where")

get_as.matrix_cls_projpred <- function() {
  ### Only works when projpred is loaded via devtools::load_all():
  # as.matrix_meths_projpred <- methods("as.matrix")
  # as.matrix_meths_projpred <- as.matrix_meths_projpred[
  #   attr(as.matrix_meths_projpred, "info")$from == "projpred"
  # ]
  ###
  as.matrix_meths_projpred <- grep(
    "^as\\.matrix\\.",
    ls(envir = asNamespace("projpred")),
    value = TRUE
  )
  as.matrix_cls_projpred <- sub("^as\\.matrix\\.", "", as.matrix_meths_projpred)
  return(as.matrix_cls_projpred)
}

