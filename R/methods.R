#' Linear predictor of a projected submodel
#'
#' Extract draws of the linear predictor from the projected submodel or 
#' submodels. If the projection has not been performed, the function
#' also performs the projection. 
#'
#' @param object The object returned by \link[=varsel]{varsel}, 
#' \link[=cv_varsel]{cv_varsel} or \link[=project]{project}.
#' @param xnew The predictor values used in the prediction. The number and order of the columns
#'  should be the same as in the original full data if argument \code{nv} is specified (see below).
#'  However, if argument \code{vind} is specified, then the number and order of columns should
#'  correspond to \code{vind}. 
#' @param ynew New (test) target variables.
#' @param offsetnew Offsets for the new observations.
#' @param weightsnew Weights for the new observations. This argument matters only if \code{ynew} is specified.
#' @param transform Should the linear predictor be transformed using the inverse-link function? 
#' Default is \code{FALSE}.
#' @param integrated If \code{TRUE}, the output is averaged over the
#' parameters. Default is \code{FALSE}.
#' @param nv Number of variables in the submodel (the variable combination is taken from the
#' \code{varsel} information). If a list, then results for all specified
#' model sizes are returned. Ignored if \code{vind} is specified.
#' @param vind Variable indices with which to predict. If specified, \code{nv} is ignored.
#' @param ns Number of samples to be projected. Ignored if \code{nc} is specified.
#' @param nc Number of clusters in the clustered projection. Default is 50.
#' @param intercept Whether to use intercept. Default is \code{TRUE}.

#' @export
proj_linpred <- function(object, xnew, ynew = NULL, offsetnew = NULL, 
                         weightsnew = NULL, transform = FALSE, 
                         integrated = FALSE, nv = NULL, vind = NULL,
                         ns = NULL, nc = NULL) {

    if (is.null(xnew))
        stop('Please provide xnew.')
    if (!is.null(vind) && NCOL(xnew) != length(vind))
        stop('The number of columns in xnew does not match with the given number of variable indices (vind).')

    if (!is.null(vind))
        nv <- NULL # ensure nv is ignored if vind is set

    if('stanreg' %in% class(object)) {
      if( !('proj' %in% names(object)) || !is.null(nc) || !is.null(ns) || !is.null(vind) ) {
        proj <- project(object, nv=nv, ns=ns, nc=nc, vind=vind, return_fit=FALSE, ...)
      } else {
        proj <- object$proj
      }
    }

    if(!.is_proj_list(proj))
        proj <- list(proj)

    if (is.null(offsetnew))
        offsetnew <- rep(0, nrow(xnew))

    # project only onto the model sizes specified in nv
    projected_sizes <- sapply(proj, function(x) NROW(x$beta))
    if(is.null(nv)) nv <- projected_sizes

    if(!all(nv %in% projected_sizes))
        stop(paste0('Linear prediction requested for nv = ',
                    paste(nv, collapse = ', '),
                    ', but projection performed only for nv = ',
                    paste(projected_sizes, collapse = ', '), '.'))

    projs <- Filter(function(x) NROW(x$beta) %in% nv, proj)
    names(projs) <- nv

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
            if (is.null(weightsnew))
                weightsnew <- rep(1, NROW(ynew))
            temp <- .get_standard_y(ynew,weightsnew)
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

    if (length(preds)==1)
        return(preds[[1]])
    else
        return(preds)
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

    stats <- subset(.bootstrap_stats(object$varsel, n_boot, alpha),
                    delta == deltas | statistic == 'kl')
    if(is.null(statistics)) statistics <- 'mlpd' #as.character(unique(stats$statistic))
    arr <- subset(stats, statistic %in% statistics)

    if(NROW(arr) == 0) {
        stop(paste0(ifelse(length(statistics)==1, 'Statistics ', 'Statistic '),
                    paste0(unique(statistics), collapse=', '), ' not available.'))
    }

    if(is.null(nv_max)) nv_max <- max(arr$size)
    ylab <- if(deltas) 'Difference to the full model' else 'value'

    ggplot(data = subset(arr, size <= nv_max), mapping = aes(x = size)) +
        # geom_ribbon(aes(ymin = lq, ymax = uq), alpha = 0.3) +
        geom_errorbar(aes(ymin = lq, ymax = uq, width=0.2, alpha=0.1)) +
        geom_line(aes(y = value)) +
        geom_point(aes(y = value)) +
        geom_hline(aes(yintercept = value), subset(arr, size == max(size)),
                   color = 'darkred', linetype=2) +
        coord_cartesian(xlim = c(0, nv_max)) +
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

#' Initialize the reference model
#' @export
init_refmodel <- function(x, y, family, mu=NULL, dis=NULL, offset=NULL, wobs=NULL, wsample=NULL,
                          intercept=TRUE, loglik=NULL) {

    # fill in the missing values with their defaults
    if (is.null(mu))
        mu <- y
    mu <- as.matrix(mu)
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

    fit <- list(x=x, y=y, fam=kl_helpers(family), mu=mu, dis=dis, offset=offset,
                wobs=wobs, wsample=wsample, intercept=intercept, loglik=loglik)
    return(fit)
}