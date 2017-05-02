#' Extract draws of the linear predictor and draw from the predictive
#' distribution of the projected submodel
#'
#' \code{proj_linpred} extracts draws of the linear predictor and
#' \code{proj_predict} draws from the predictive distribution of the projected
#' submodel or submodels. If the projection has not been performed, the
#' functions also perform the projection.
#'
#' @name proj_pred
#'
#' @param object The object returned by \link[=varsel]{varsel},
#' \link[=cv_varsel]{cv_varsel} or \link[=project]{project}.
#' @param xnew The predictor values used in the prediction. The number and 
#' order of the columns should be the same as in the original full data if
#' argument \code{nv} is specified (see below). However, if argument 
#' \code{vind} is specified, then the number and order of columns should
#' correspond to \code{vind}.
#' @param ynew New (test) target variables.
#' @param draws Number of draws to return for \code{proj_predict}. Default is
#' \code{nc}.
#' @param offsetnew Offsets for the new observations. Defaults to a vector of
#' 0s.
#' @param weightsnew Weights for the new observations. For binomial model,
#' corresponds to the number trials per observation. For \code{proj_linpred},
#' this argument matters only if \code{ynew} is specified. Defaults to a vector
#' of 1s.
#' @param transform Should the linear predictor be transformed using the
#' inverse-link function? Default is \code{FALSE}.
#' @param integrated If \code{TRUE}, the output is averaged over the
#' parameters. Default is \code{FALSE}.
#' @param nv Number of variables in the submodel (the variable combination is
#' taken from the \code{varsel} information). If a list, then results for all
#' specified model sizes are returned. Ignored if \code{vind} is specified.
#' @param vind Variable indices with which to predict. If specified, \code{nv}
#' is ignored.
#' @param ns Number of draws to be projected. Ignored if \code{nc} is specified.
#' @param nc Number of clusters in the clustered projection. Default is 50.
#' @param intercept Whether to use intercept. Default is \code{TRUE}.
NULL

proj_helper <- function(object, xnew, nv, vind, ns, nc = NULL) {
  if (is.null(xnew))
    stop('Please provide xnew.')
  if (!is.null(vind) && NCOL(xnew) != length(vind))
    stop(paste('The number of columns in xnew does not match with the given',
               'number of variable indices (vind).'))

  if (!is.null(vind))
    nv <- NULL # ensure nv is ignored if vind is set

  if( 'projection' %in% class(object) || (length(object)>0 && 'projection' %in% class(object[[1]])) ) {
      proj <- object
      if(any(sapply(list(nv, vind, ns, nc), Negate(is.null))))
          warning('nv, vind, ns and nc are ignored when object is a projection.')
  } else {
      # reference model obtained, so run the projection
      proj <- project(object, nv = nv, ns = ns, nc = nc, vind = vind)
  }

  if(!.is_proj_list(proj)) {
    proj <- list(proj)
  } else {
    # proj is not a projection object
    if(any(sapply(proj, function(x) !('family_kl' %in% names(x)))))
      stop(paste('proj_linpred only works with objects returned by',
                 ' varsel, cv_varsel or project'))
  }

  projected_sizes <- sapply(proj, function(x) NROW(x$beta))
  if(is.null(nv))
    nv <- projected_sizes

  if(!all(nv %in% projected_sizes))
    stop(paste0('Linear prediction requested for nv = ',
                paste(nv, collapse = ', '),
                ', but projection performed only for nv = ',
                paste(projected_sizes, collapse = ', '), '.'))

  projs <- Filter(function(x) NROW(x$beta) %in% nv, proj)
  names(projs) <- nv

  projs
}

#' @rdname proj_pred
#' @export
proj_linpred <- function(object, xnew, ynew = NULL, offsetnew = NULL,
                         weightsnew = NULL, transform = FALSE,
                         integrated = FALSE, nv = NULL, vind = NULL,
                         ns = NULL, nc = NULL) {
    if (is.null(offsetnew))
        offsetnew <- rep(0, nrow(xnew))
    if (is.null(weightsnew))
        weightsnew <- rep(1, nrow(xnew))
    projs <- proj_helper(object = object, xnew = xnew, nv = nv, vind = vind,
                         ns = ns, nc = nc)

    preds <- lapply(projs, function(proj) {
        if (!is.null(vind))
            # columns of xnew are assumed to match to the given variable indices
            xtemp <- xnew
        else
            xtemp <- xnew[, proj$ind, drop = F]
        mu <- proj$family_kl$mu_fun(xtemp, proj$alpha, proj$beta, offsetnew)
        if(transform)
            pred <- t(mu)
        else
            pred <- t(proj$family_kl$linkfun(mu))
        if (integrated)
            # average over the parameters
            pred <- as.vector( proj$weights %*% pred )
        else if (!is.null(dim(pred)) && dim(pred)[1]==1)
            # return a vector if pred contains only one row
            pred <- as.vector(pred)
        if (!is.null(ynew)) {
            # compute also the log-density
            temp <- .get_standard_y(ynew,weightsnew,proj$family_kl)
            ynew <- temp$y
            weightsnew <- temp$weights
            lpd <- proj$family_kl$ll_fun(mu, proj$dis, ynew, weightsnew)
            if (integrated && !is.null(dim(lpd)))
                lpd <- as.vector(apply(lpd, 1, log_weighted_mean_exp, proj$weights))
            else if (!is.null(dim(lpd)))
                lpd <- t(lpd)
            return(list(pred=pred, lpd=lpd))
        } else
            return(list(pred=pred))
    })

    .unlist_proj(preds)
}

#' @rdname proj_pred
#' @export
proj_predict <- function(object, xnew, offsetnew = NULL, weightsnew = NULL, 
                         nv = NULL, vind = NULL, ns = NULL, nc = NULL, 
                         draws = NULL, seed = NULL) {
  if (!is.null(seed))
    set.seed(seed)
  if (is.null(offsetnew))
    offsetnew <- rep(0, nrow(xnew))
  if (is.null(weightsnew))
    weightsnew <- rep(1, nrow(xnew))

  projs <- proj_helper(object = object, xnew = xnew, nv = nv, vind = vind,
                       ns = ns, nc = nc)
  if(is.null(draws))
    draws <- length(projs[[1]]$weights)
  
  preds <- lapply(projs, function(proj) {
    if (!is.null(vind)) {
      # columns of xnew are assumed to match to the given variable indices
      xtemp <- xnew
    } else {
      xtemp <- xnew[, proj$ind, drop = F]
    }
    mu <- proj$family_kl$mu_fun(xtemp, proj$alpha, proj$beta, offsetnew)
    dis <- proj$dis
    draw_inds <- sample(x = seq_along(proj$weights), size = draws, 
                        replace = TRUE, prob = proj$weights)
    t(sapply(draw_inds, function(i) {
      proj$family_kl$ppd_fun(mu[,i], dis[i], weightsnew)
    }))
  })
  
  .unlist_proj(preds)
}

#' Plotting or printing summary statistics related to variable selection
#'
#' \code{varsel_statistics} can be used to obtain summary statistics related to
#' variable selection. The same statistics can be plotted with
#' \code{varsel_plot}.
#'
#' @name varsel-statistics
#'
#' @param object The object returned by \link[=varsel]{varsel} or
#' \link[=cv_varsel]{cv_varsel}.
#' @param nv_max Maximum submodel size for which the statistics are calculated.
#' @param statistics A list of strings of statistics to calculate. Available
#' options are: kl, mse, mlpd, kl, (gaussian only), pctcorr (binomial only).
#' If \code{NULL}, set to varsel_plot plots only mlpd, but varsel_statistics
#' return all the statistics.
#' @param deltas If \code{TRUE}, the difference between the full model and the
#' submodel is returned instead of the actual value of the statistic.
#' Defaults to \code{FALSE}.
#' @param n_boot Number of bootstrap samples for calculating the credible
#' intervals of the statistics.
#' @param alpha A number indicating the desired coverage of the credible
#' intervals. Eg. \code{alpha=0.1} corresponds to 90\% probability mass
#' within the intervals. Defaults to \code{0.1}.
NULL

#' @rdname varsel-statistics
#' @export
varsel_plot <- function(object, ..., nv_max = NULL, statistics = NULL, deltas = T,
                        n_boot = 1000, alpha = 0.1) {
    if(!('varsel' %in% names(object)))
        stop(paste('The provided object doesn\'t contain information about the',
                   'variable selection. Run the variable selection first.'))

    boot_stats <- .bootstrap_stats(object$varsel, n_boot, alpha)
    
    #
    if(is.null(statistics)) statistics <- 'mlpd' #as.character(unique(stats$statistic))
    if(deltas) {
      full_stats <- data.frame(statistic = statistics, value = 0)
    } else {
      boot_vals <- subset(boot_stats, size == 0 & delta == F &
                            statistic %in% statistics, 'value', drop = T)
      boot_deltas <- subset(boot_stats, size == 0 & delta == T &
                            statistic %in% statistics, 'value', drop = T)
      if('kl' %in% statistics) {
        boot_deltas <- c(0, boot_deltas)
        boot_vals[1] <- 0
      }
      full_stats <- data.frame(
        statistic = subset(boot_stats, size == 0 & delta == F &
                             statistic %in% statistics, 'statistic'),
        value = boot_vals - boot_deltas)
    }
     
    stats <- subset(boot_stats, delta == deltas | statistic == 'kl')
    arr <- subset(stats, statistic %in% statistics)
    

    if(NROW(arr) == 0) {
        stop(paste0(ifelse(length(statistics)==1, 'Statistics ', 'Statistic '),
                    paste0(unique(statistics), collapse=', '), ' not available.'))
    }

    if(is.null(nv_max)) nv_max <- max(arr$size)
    ylab <- if(deltas) 'Difference to the full model' else 'value'

    # make sure that breaks on the x-axis are integers
    n_opts <- c(4,5,6)
    n_possible <- Filter(function(x) nv_max %% x == 0, n_opts)
    n_alt <- n_opts[which.min(n_opts - (nv_max %% n_opts))]
    nb <- ifelse(length(n_possible) > 0, min(n_possible), n_alt)
    by <- ceiling(nv_max/min(nv_max, nb))
    breaks <- seq(0, by*min(nv_max, nb), by)
    minor_breaks <- if(by%%2 == 0)
        seq(by/2, by*min(nv_max, nb), by)
    else
      NULL

    ggplot(data = subset(arr, size <= nv_max), mapping = aes(x = size)) +
        geom_linerange(aes(ymin = lq, ymax = uq, alpha=0.1)) +
        geom_line(aes(y = value)) +
        geom_point(aes(y = value)) +
        geom_hline(aes(yintercept = value), data = full_stats,
                   color = 'darkred', linetype=2) +
        scale_x_continuous(breaks = breaks, minor_breaks = minor_breaks,
                           limits = c(min(breaks), max(breaks))) +
        labs(x = 'Number of variables in the submodel', y = ylab) +
        theme(legend.position = 'none') +
        facet_grid(statistic ~ ., scales = 'free_y')
}

#' @rdname varsel-statistics
#' @export
varsel_statistics <- function(object, ..., nv_max = NULL, deltas = F) {
    if(!('varsel' %in% names(object)))
        stop(paste('The provided object doesn\'t contain information about the',
                   'variable selection. Run the variable selection first.'))

    stats <- subset(.bootstrap_stats(object$varsel, NULL, 0.5),
                    delta == deltas | statistic == 'kl')
    statistics <- as.character(unique(stats$statistic))

    arr <- data.frame(sapply(statistics, function(sname) {
        unname(subset(stats, statistic == sname, 'value'))
    }))
    arr <- cbind(size = unique(stats$size), arr)

    if(is.null(nv_max)) nv_max <- max(stats$size)

    arr$chosen <- c(NA, object$varsel$chosen)
    if('pctch' %in% names(object$varsel))
      arr$pctch <- c(NA, diag(object$varsel$pctch[,-1]))

    subset(arr, size <= nv_max)
}

#' Generic reference model initialization
#' 
#' Initializes a structure that can be used as a reference fit for the 
#' projective variable selection.
#' This function is provided to allow construction of the reference fit
#' using also other tools than \code{rstanarm}, because only a limited amount
#' of information is needed for the actual projection and the variable selection.
#' 
#' 
#' @name refmodel
#' 
#' @param x Predictor matrix of dimension \code{n}-by-\code{D} containing the candidate
#'  variables for selection
#' (i.e. variables from which to select the submodel). Rows denote the observations
#' and columns the different variables. Note that this matrix can be different from
#' the one used to construct the reference model.
#' @param y Vector of length \code{n} giving the target variable values.
#' @param family \link[=family]{family} object giving the model family
#' @param mu \code{n}-by-{S} matrix of expected values for \code{y}, each column corresponding
#' to one posterior draw for the reference model.
#' @param dis Vector of length \code{S} giving the posterior draws for the dispersion parameter
#' in the reference model if there is such a parameter in the model family. For Gaussian
#' observation model this is the noise std \code{sigma}.
#' @param offset Offset to be added to the linear predictor in the projection. (Same as in 
#' function \code{glm})
#' @param wobs Observation weights. The weights should sum to \code{n}. 
#' If omitted, equal weights are assumed. 
#' @param wsample vector of length \code{S} giving the weights for the posterior draws.
#' The weights should sum to one. If omitted, equal weights are assumed.
#' @param intercept Whether to use intercept. Default is \code{TRUE}.
#' @param loglik \code{S}-by-\code{n} matrix giving the log-likelihood values
#' for the reference model for each pair of \code{S} posterior draws and \code{n} observations.
#' Can be omitted but is mandatory for performing the LOO validation.
#' 
#' @return An object that can be passed to all the functions that
#' take the reference fit as the first argument, such as \code{varsel}, \code{cv_varsel},
#' \code{proj_predict} and \code{proj_linpred}.
NULL

#' @rdname refmodel
#' @export
init_refmodel <- function(x, y, family, mu=NULL, dis=NULL, offset=NULL, wobs=NULL, wsample=NULL,
                          intercept=TRUE, loglik=NULL) {
    
    # fill in the missing values with their defaults
    if (is.null(mu))
        mu <- y
    mu <- unname(as.matrix(mu))
    
    S <- NCOL(mu) # number of samples in the reference model
    n <- length(y)
    if (is.null(dis))
        dis <- rep(1, S)
    if (is.null(offset))
        offset <- rep(0, n)
    if (is.null(wobs))
        wobs <- rep(1, n)
    if (is.null(wsample))
        wsample <- rep(1/S, S)
    if (is.null(intercept))
        intercept <- TRUE
    
    # figure out column names for the variables
    if (!is.null(colnames(x)))
        coefnames <- colnames(x)
    else
        coefnames <- paste0('x',1:ncol(x))

    # y and the observation weights in a standard form
    temp <- .get_standard_y(y, wobs, family)
    
    fit <- list(x=x, y=temp$y, fam=kl_helpers(family), mu=mu, dis=dis, coefnames=coefnames,
                offset=offset, wobs=temp$weights, wsample=wsample, intercept=intercept, loglik=loglik)
    
    # define the class of the retuned object to be 'refmodel'
    class(fit) <- 'refmodel'
    return(fit)
}






