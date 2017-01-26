# All functions that users will use to extract model parameters,
# plot variable selection statistics etc.

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

#' @export
proj_linpred <- function(object, transform = FALSE, xnew = NULL, ynew = NULL, offsetnew = NULL, 
						 newdata = NULL, nv = NULL, integrated = FALSE, ...) {
  
  # TODO, IMPLEMENT THE PROJECTION/PREDICTION WITH AN ARBITRARY VARIABLE COMBINATION 
	
  if(!('proj' %in% names(object)))
  	object <- project(object, nv=nv, ...)

  vars <- .extract_vars(object)
  family_kl <- vars$fam

  if (is.null(xnew)) {
  	xnew <- vars$x
  	if (is.null(ynew))
  	  ynew <- vars$y
  }
  nt <- nrow(xnew)
  if (is.null(offsetnew))
  	offsetnew <- rep(0,nt)
  
  # project only model the sizes of which are specified in nv
  projected_sizes <- sapply(object$proj$p_sub, function(psub) NROW(psub$beta))
  if(is.null(nv)) nv <- projected_sizes

  if(!all(nv %in% projected_sizes))
    stop(paste0('Linear prediction requested for nv = ',
                paste(nv, collapse = ', '),
                ', but projection performed only for nv = ',
                paste(projected_sizes, collapse = ', '), '.'))

  projs <- Filter(function(psub) NROW(psub$beta) %in% nv, object$proj$p_sub)
  names(projs) <- nv
  
  preds <- mapply(function(proj, nv) {
    ch <- object$varsel$chosen[min(nv,1):nv]
    mu <- family_kl$mu_fun(xnew[,ch,drop=F], proj$alpha, proj$beta, offsetnew)
    if(transform)
    	pred <- t(mu)
    else
    	pred <- t(family_kl$linkfun(mu))
    if (integrated)
    	# average over the parameters
    	pred <- as.vector( proj$weights %*% pred )
    else if (!is.null(dim(pred)) && dim(pred)[1]==1)
    	# return a vector if pred contains only one row
    	pred <- as.vector(pred)
    if (!is.null(ynew)) {
    	# compute also the log-density
    	lpd <- family_kl$ll_fun(mu, proj$dis, ynew)
    	if (integrated && !is.null(dim(lpd)))
    		lpd <- as.vector(apply(lpd, 1, log_weighted_mean_exp, proj$weights))
    	else if (!is.null(dim(lpd)))
    		lpd <- t(lpd)
    	return(list(pred=pred, lpd=lpd))
    } else
    	return(list(pred=pred))
  }, projs, nv, SIMPLIFY = F)

  if (length(preds)==1)
  	return(preds[[1]])
  else
  	return(preds)
}

#' @export
proj_coef <- function(object, ...) {
  
  if(!('proj' %in% names(object)))
    stop(paste('The provided object doesn\'t contain information about the projection.',
               'Run the projection first.'))

  fun <- function(b, w) drop(b%*%w)
  proj_coef_helper(object, fun)
}

#' @export
proj_se <- function(object, ...) {
  
  if(!('proj' %in% names(object)))
    stop(paste('The provided object doesn\'t contain information about the projection.',
               'Run the projection first.'))
  # weighted standard deviation (using cluster weights)
  fun <- function(b, w) {
    drop(sqrt(((b - rowMeans(b))^2)%*%w / ((sum(w>0)-1)/sum(w>0)*sum(w))))
  }
  proj_coef_helper(object, fun)
}

proj_coef_helper <- function(object, fun) {
  # calculates 'fun' for each projected weight vector b and sample weights w
  coefnames <- names(coef(object))
  if(object$proj$intercept) coefnames <- coefnames[-1]

  lapply(object$proj$p_sub, function(proj) {
    b <- proj$beta
    if(NROW(b) == 0) return(0)
    rownames(b) <- coefnames[object$varsel$chosen[1:NROW(b)]]
    if(object$proj$intercept) {
      b <- rbind(proj$alpha, b)
      rownames(b)[1] <- names(coef(object))[1]
    }

    fun(b, proj$weights)
  })
}

#' @export
proj_sigma <- function(object, ...) {
  # only gaussian family supported currently
  
  if(!('proj' %in% names(object)))
    stop(paste('The provided object doesn\'t contain information about the projection.',
               'Run the projection first.'))
  vars <- .extract_vars(object)
  if(!(vars$fam$family %in% c('gaussian')))
    stop('Sigma available only for the gaussian family.')

  lapply(object$proj$p_sub, function(proj) {
    if(family(object)$family == 'gaussian') {
      drop(sqrt(proj$dis^2%*%proj$weights))
    }
  })
}

#' @export
varsel_plot <- function(x, ..., nv_max = NULL, statistics = NULL, deltas = T,
                        n_boot = 1000, alpha = 0.05) {
  if(!('varsel' %in% names(x)))
    stop(paste('The provided object doesn\'t contain information about the',
               'variable selection. Run the variable selection first.'))

  stats <- subset(.bootstrap_stats(x$varsel, n_boot, alpha),
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
    facet_grid(statistic ~ ., scales = 'free_y')
}

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
  if('pctch' %in% names(object$varsel)) arr$pctch <- c(NA, object$varsel$pctch)

  subset(arr, size <= nv_max)
}

