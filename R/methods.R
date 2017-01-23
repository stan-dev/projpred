# All functions that users will use to extract model parameters,
# plot variable selection statistics etc.

#' @export
init_refmodel <- function(x, y, family, mu=NULL, dis=NULL, offset=NULL, wobs=NULL, wsample=NULL, intercept=TRUE) {
    
    # fill in the missing values with their defaults
    if (is.null(mu))
        mu <- y
    mu <- as.matrix(mu)
    S <- NCOL(mu) # number of samples in the reference model
    n <- length(y)
    if (is.null(dis))
        dis <- rep(NA, S)
    if (is.null(offset))
        offset <- rep(0, n)
    if (is.null(wobs))
        wobs <- rep(1, n)
    if (is.null(wsample))
        wsample <- rep(1/S, S)
    if (is.null(intercept))
        intercept <- TRUE
    
    fit <- list(x=x, y=y, fam=kl_helpers(family), mu=mu, dis=dis, offset=offset,
                wobs=wobs, wsample=wsample, intercept=intercept)
    return(fit)
}

#' @export
proj_linpred <- function(object, transform = FALSE, newdata = NULL, offset = NULL, nv = NULL, integrated = FALSE, ...) {
  .validate_for_varsel(object)
  if(!('proj' %in% names(object)))
    stop(paste('The provided object doesn\'t contain information about the',
               'projection. Run the projection first.'))

  family_kl <- kl_helpers(family(object))

  data <- rstanarm:::pp_data(object, newdata, offset = offset)
  obj_intercept <- attr(object$terms,'intercept') %ORifNULL% 0
  if(obj_intercept) data$x <- data$x[,-1]


  # project only model the sizes of which are specified in nv
  projected_sizes <- sapply(object$proj$p_sub, function(x) NROW(x$beta))
  if(is.null(nv)) nv <- projected_sizes

  if(!all(nv %in% projected_sizes))
    stop(paste0('Linear prediction requested for nv = ',
                paste(nv, collapse = ', '),
                ', but projection performed only for nv = ',
                paste(projected_sizes, collapse = ', '), '.'))

  projs <- Filter(function(x) NROW(x$beta) %in% nv, object$proj$p_sub)
  names(projs) <- nv
  
  mapply(function(proj, nv) {
    ch <- object$varsel$chosen[min(nv,1):nv]
    mu <- t(family_kl$mu_fun(data$x[, ch, drop = F],
                             proj$alpha,
                             proj$beta, data$offset, object$proj$intercept))
    if(transform)
    	qty <- mu
    else
    	qty <- family_kl$linkfun(mu)
    if (integrated)
    	return(as.vector( proj$weights %*% qty ))
    else
    	return( qty )
  }, projs, nv, SIMPLIFY = F)

}

#' @export
proj_coef <- function(object, ...) {
  .validate_for_varsel(object)
  if(!('proj' %in% names(object)))
    stop(paste('The provided object doesn\'t contain information about the projection.',
               'Run the projection first.'))

  fun <- function(b, w) drop(b%*%w)
  proj_coef_helper(object, fun)
}

#' @export
proj_se <- function(object, ...) {
  .validate_for_varsel(object)
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
  # only gaussian family supported
  .validate_for_varsel(object)
  if(!('proj' %in% names(object)))
    stop(paste('The provided object doesn\'t contain information about the projection.',
               'Run the projection first.'))
  if(!(family(object)$family %in% 'gaussian'))
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
               'variable selection. Run the variable selection first!'))

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
               'variable selection. Run the variable selection first!'))

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

