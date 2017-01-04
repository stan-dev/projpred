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

  data <- rstanarm:::pp_data(object, newdata, offset = offset)


  # project only model the sizes of which are specified in nv
  projected_sizes <- sapply(object$proj$p_sub, function(x) NROW(x$beta))
  if(is.null(nv)) nv <- projected_sizes

  if(!all(nv %in% projected_sizes))
    stop(paste0('Linear prediction requested for nv = ',
                paste(nv, collapse = ', '),
                ', but projection performed only for nv = ',
                paste(projected_sizes, collapse = ', '), '.'))

  projs <- Filter(function(x) NROW(x$beta) %in% nv, object$proj$p_sub)
  chosen <- object$varsel$chosen
  if(object$proj$intercept) chosen <- c(1, chosen + 1)

  lapply(projs, function(proj) {
    family_kl$mu_fun(x[,chosen])
    res <- t(dat$x[, chosen[1:nrow(proj$b)], drop = F]%*%proj$beta +  + dat$offset)
    if(transform) family(object)$linkinv(res) else res
  })

}

#' @export
proj_coef <- function(object, ...) {
  .validate_for_varsel(object)
  if(!('proj' %in% names(object)))
    stop(paste('The stanreg object doesn\'t contain information about the projection.',
               'Run the projection first.'))
  lapply(object$proj, function(x) drop(x$b%*%x$weights))
}

#' @export
proj_se <- function(object, ...) {
  .validate_for_varsel(object)
  if(!('proj' %in% names(object)))
    stop(paste('The stanreg object doesn\'t contain information about the projection.',
               'Run the projection first.'))
  # weighted standard deviation (using cluster weights)
  lapply(object$proj, function(x) {
    n0 <- sum(x$weights>0)
    drop(sqrt(((x$b - rowMeans(x$b))^2)%*%x$weights /
                ((n0-1)/n0*sum(x$weights))))
  })
}

#' @export
proj_sigma <- function(object, ...) {
  .validate_for_varsel(object)
  if(!('proj' %in% names(object)))
    stop(paste('The stanreg object doesn\'t contain information about the projection.',
               'Run the projection first.'))
  lapply(object$proj, function(x) if('dis' %in% names(x)) sqrt(mean(x$dis^2)) else 1)
}


#' @export
varsel_plot <- function(x, ..., nv = NULL, summaries = NULL, deltas = T,
                        data = 'test') {
  if(!('varsel' %in% names(x)))
    stop(paste('The stanreg object doesn\'t contain information about the variable',
               'selection. Run the variable selection first!'))

  data_remove <- if(data=='train') 'test' else 'train'
  if(is.null(summaries)) summaries <- as.character(unique(x$varsel$stats$summary))
  arr <- subset(x$varsel$stats, data != data_remove & (delta == deltas | summary == 'kl')
                & summary %in% summaries)

  if(nrow(arr) == 0)
    stop(paste0(ifelse(length(summaries)==1, 'Summaries ', 'Summary '),
                paste0(unique(summaries), collapse = ', '),
                ' not evaluated on ', data, ' data.'))

  if(is.null(nv)) nv <- max(arr$size)
  ylab <- if(deltas) expression(Delta) else 'value'

  ggplot(data = subset(arr, size <= nv), mapping = aes(x = size)) +
    geom_ribbon(aes(ymin = lq, ymax = uq), alpha = 0.3) +
    geom_line(aes(y = value)) +
    geom_hline(aes(yintercept = value), subset(arr, size == max(size)),
               color = 'darkred') +
    coord_cartesian(xlim = c(0, nv)) +
    labs(x = 'Number of variables in the submodel', y = ylab) +
    facet_grid(summary ~ ., scales = 'free_y')
}

#' @export
varsel_summary <- function(object, ..., nv = NULL, deltas = F, train = F) {
  if(!('varsel' %in% names(object)))
     stop(paste('Stanreg object doesn\'t contain information about the variable',
                'selection.\nRun the variable selection first!'))
  if(is.null(nv)) nv <- max(object$varsel$stats$size)
  data_remove <- if(train) 'test' else 'train'

  summaries <- as.character(unique(object$varsel$stats$summary))
  arr_list <- lapply(summaries, function(sname, stats, dr) {
    res <- subset(stats, summary == sname & (delta == deltas | summary == 'kl')
                  & data != dr, c('size', 'value'))
    setNames(res, c('size', sname))
  }, object$varsel$stats, data_remove)
  # combine the arrays
  arr <- Reduce(merge, arr_list)
  # If no test data and results from training data are not wanted, return only KL
  if(nrow(arr) == 0)  arr <- setNames(subset(
    object$varsel$stats, summary == 'kl', c('size','value')), c('size', 'kl'))

  # make sure that the list is ordered
  arr <- arr[order(arr$size),]

  arr$chosen <- c(object$varsel$chosen, NA)
  if('pctch' %in% names(object$varsel)) arr$pctch <- c(object$varsel$pctch, NA)

  subset(arr, size <= nv, c('size', intersect(colnames(arr),
                                              c(summaries, 'chosen', 'pctch'))))
}
