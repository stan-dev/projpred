#' @export
project <- function(object, nv, ...) {
  UseMethod('project')
}

#' @export
#' @export
project.stanreg <- function(object, nv, ns = 400L, intercept = NULL, ...) {
  # Missing: Clustering!
  if(is.null(nv)) stop('nv not provided')
  .validate_for_varsel(object)
  if(!('varsel' %in% names(object)))
    stop(paste('The stanreg object doesn\'t contain information about the ',
               'variable selection. Run the variable selection first.'))

  vars <- .extract_vars(object)
  if(ns > ncol(vars$beta)) {
    warning(paste0('Setting the number of samples to ', ncol(vars$beta),'.'))
    ns <- ncol(vars$beta)
  }
  if(is.null(intercept)) intercept <- vars$intercept

  family_kl <- kl_helpers(family(object))

  if(max(nv) > length(object$varsel$chosen))
    stop(paste('Cannot perform the projection with', max(nv), 'variables,',
               'because the variable selection has been run only up to',
               length(object$varsel$chosen), 'variables.'))

  e <- get_data_and_parameters(vars, NA, intercept, ns, family_kl)
  # p_sub <- .get_submodels(object$varsel$chosen, nv, family_kl, e$p_full,
  #                         e$d_train, intercept)
  #
  # names <- names(coef(object))[object$varsel$chosen]
  # if(intercept) names <- names[-1]
  #
  # object$proj <- mapply(function(p_sub, nv) {
  #   p_sub$kl <- NULL
  #   p_sub$b <- p_sub$beta
  #   if(nv>0) rownames(p_sub$b) <- names[1:nv]
  #   if(intercept) {
  #     p_sub$b <- rbind(p_sub$alpha, p_sub$b)
  #     rownames(p_sub$b)[1] <- names(coef(object))[1]
  #   }
  #   p_sub$kl <- NULL
  #   p_sub$beta <- NULL
  #   p_sub$alpha <- NULL
  #   if(!(family_kl$family %in% c('gaussian', 'Gamma'))) p_sub$dis <- NULL
  #   p_sub$nv <- nv
  # }, p_sub, nv, SIMPLIFY = F)

  object$proj <- list(
    p_sub = .get_submodels(object$varsel$chosen, nv, family_kl,
                           e$p_full, e$d_train, intercept),
    intercept = intercept)

  object
}

#' @export
proj_linpred <- function(object, transform = FALSE, newdata = NULL, offset = NULL, nv = NULL, ...) {
  .validate_for_varsel(object)
  if(!('proj' %in% names(object)))
    stop(paste('The stanreg object doesn\'t contain information about the',
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
    mu <- family_kl$mu_fun(data$x[, ch, drop = F],
                           proj$alpha,
                           proj$beta, data$offset, object$proj$intercept)
    if(transform) mu else family_kl$linkfun(mu)
  }, projs, nv, SIMPLIFY = F)

}

#' @export
proj_coef <- function(object, ...) {
  .validate_for_varsel(object)
  if(!('proj' %in% names(object)))
    stop(paste('The stanreg object doesn\'t contain information about the projection.',
               'Run the projection first.'))

  fun <- function(b, w) drop(b%*%w)
  proj_coef_helper(object, fun)
}

#' @export
proj_se <- function(object, ...) {
  .validate_for_varsel(object)
  if(!('proj' %in% names(object)))
    stop(paste('The stanreg object doesn\'t contain information about the projection.',
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
    stop(paste('The stanreg object doesn\'t contain information about the projection.',
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
varsel_plot <- function(x, ..., nv_max = NULL, metrics = NULL, deltas = T,
                        data = 'test') {
  if(!('varsel' %in% names(x)))
    stop(paste('The stanreg object doesn\'t contain information about the',
               'variable selection. Run the variable selection first!'))

  data_remove <- if(data=='train') 'test' else 'train'
  if(is.null(metrics)) metrics <- as.character(unique(x$varsel$metrics$metric))
  arr <- subset(x$varsel$metrics, data != data_remove &
                  (delta == deltas | metric == 'kl') & metric %in% metrics)

  if(nrow(arr) == 0)
    stop(paste0(ifelse(length(metrics)==1, 'Summaries ', 'Summary '),
                paste0(unique(metrics), collapse = ', '),
                ' not evaluated on ', data, ' data.'))

  if(is.null(nv_max)) nv_max <- max(arr$size)
  ylab <- if(deltas) expression(Delta) else 'value'

  ggplot(data = subset(arr, size <= nv_max), mapping = aes(x = size)) +
    geom_ribbon(aes(ymin = lq, ymax = uq), alpha = 0.3) +
    geom_line(aes(y = value)) +
    geom_hline(aes(yintercept = value), subset(arr, size == max(size)),
               color = 'darkred') +
    coord_cartesian(xlim = c(0, nv_max)) +
    labs(x = 'Number of variables in the submodel', y = ylab) +
    facet_grid(metric ~ ., scales = 'free_y')
}

#' @export
varsel_summary <- function(object, ..., nv_max = NULL, deltas = F, data = 'test') {
  if(!('varsel' %in% names(object)))
    stop(paste('The stanreg object doesn\'t contain information about the',
               'variable selection. Run the variable selection first!'))

  if(is.null(nv_max)) nv_max <- max(object$varsel$metrics$size)
  data_remove <- if(data=='train') 'test' else 'train'

  metrics <- as.character(unique(object$varsel$metrics$metric))
  arr_list <- lapply(metrics, function(sname, metrics_arr, dr) {
    res <- subset(metrics_arr, metric == sname & (delta == deltas | metric == 'kl')
                  & data != dr, c('size', 'value'))
    setNames(res, c('size', sname))
  }, object$varsel$metrics, data_remove)
  # combine the arrays
  arr <- Reduce(merge, arr_list)
  # If no test data and results from training data are not wanted, return only KL
  if(nrow(arr) == 0)  arr <- setNames(subset(
    object$varsel$metrics, metric == 'kl', c('size','value')), c('size', 'kl'))

  # make sure that the list is ordered
  arr <- arr[order(arr$size),]

  arr$chosen <- c(NA, object$varsel$chosen)
  if('pctch' %in% names(object$varsel)) arr$pctch <- c(NA, object$varsel$pctch)

  subset(arr, size <= nv_max, c('size', intersect(colnames(arr),
                                              c(metrics, 'chosen', 'pctch'))))
}
