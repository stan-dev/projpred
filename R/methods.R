#' @export
project <- function(object, nv, ...) {
  UseMethod('project')
}

#' @export
#' @export
project.stanreg <- function(object, nv, ...) {
  if(is.null(nv)) stop('nv not provided')
  .validate_for_varsel(object)
  if(!('varsel' %in% names(object)))
    stop(paste('The stanreg object doesn\'t contain information about the ',
               'variable selection. Run the variable selection first.'))

  vars <- .extract_vars(object)
  args <- .init_args(list(...), vars)
  family_kl <- kl_helpers(family(object))

  if(max(nv) > length(object$varsel$chosen))
    stop(paste('Cannot perform the projection with', max(nv), 'variables,',
               'because the variable selection has been run only up to',
               length(object$varsel$chosen), 'variables.'))

  v_inds_max <- object$varsel$chosen[1:max(nv)]
  # the line above fails with nv=0 (ie. projection with only intercept)
  if(length(nv) == 1 && nv == 0) v_inds_max <- 0

  d_train <- list(x = vars$x[,v_inds_max],
                  weights = vars$weights,
                  offset = vars$offset)

  mu <- family_kl$mu_fun(vars$x, vars$alpha, vars$beta, vars$offset,
                         vars$intercept)
  dis <- vars$dis
  coef_init <- list(alpha = median(vars$alpha),
                    beta = matrix(apply(vars$beta, 1, median), ncol = 1))

  s_ind <- round(seq(1, args$ns_total, length.out  = args$ns))
  p_full <- list(mu = mu[, s_ind], dis = dis[s_ind],
                 weights = rep(1/args$ns, args$ns))

  projfun <- .get_proj_handle(family_kl)
  names <- names(coef(object))
  if(vars$intercept) names <- names[-1]

  object$proj <- lapply(nv, function(nv, names) {
    # if no intercept and 0 variables, return 'trivial'  result
    if(nv == 0 & vars$intercept == F) {
      mu_null <- family_kl$linkinv(matrix(0, NROW(p_full$mu), NCOL(p_full$mu)))
      return(list(weights = p_full$weights,
                  dis = family_kl$dis_fun(p_full, d_train, list(mu = mu_null)),
                  b = matrix(0, 1, NCOL(p_full$mu)),
                  intercept = 0,
                  nv = nv))
    }

    seq <- if(nv>0) 1:nv else 0
    proj <- projfun(seq, p_full, d_train, vars$intercept, args$regul,
                    coef_init)
    rownames(proj$beta) <- names[seq]
    proj$b <- proj$beta[seq, , drop = F]
    proj$intercept <- vars$intercept
    if(proj$intercept) {
      proj$b <- rbind(proj$alpha, proj$b)
      rownames(proj$b)[1] <- names(coef(object))[1]
    }
    proj$kl <- NULL
    proj$beta <- NULL
    proj$alpha <- NULL
    if(!(family_kl$family %in% c('gaussian', 'Gamma'))) proj$dis <- NULL
    proj$nv <- nv
    proj
  }, names[v_inds_max])

  object
}

#' @export
proj_linpred <- function(object, transform = FALSE, newdata = NULL, offset = NULL, nv = NULL, ...) {
  .validate_for_varsel(object)
  if(!('proj' %in% names(object)))
    stop(paste('The stanreg object doesn\'t contain information about the',
               'projection. Run the projection first.'))

  dat <- rstanarm:::pp_data(object, newdata, offset = offset)

  # project only model the sizes of which are specified in nv
  if(is.null(nv)) nv <- sapply(object$proj, function(x) x$nv)
  projected_sizes <- sapply(object$proj, function(x) x$nv)
  if(!all(nv %in% projected_sizes))
    stop(paste0('Linear prediction requested for nv = ',
                paste(nv, collapse = ', '),
                ', but projection performed only for nv = ',
                paste(projected_sizes, collapse = ', '), '.'))

  projs <- Filter(function(x) x$nv %in% nv, object$proj)
  chosen <- object$varsel$chosen
  if(projs[[1]]$intercept) chosen <- c(1, chosen + 1)

  lapply(projs, function(proj) {
    res <- t(dat$x[, chosen[1:nrow(proj$b)], drop = F]%*%proj$b + dat$offset)
    if(transform) family(object)$linkinv(res) else res
  })

}

#' @export
proj_coef <- function(object, ...) {
  .validate_for_varsel(object)
  if(!('proj' %in% names(object)))
    stop(paste('The stanreg object doesn\'t contain information about the projection.',
               'Run the projection first.'))
  lapply(object$proj, function(x) rowMeans(x$b))
}

#' @export
proj_se <- function(object, ...) {
  .validate_for_varsel(object)
  if(!('proj' %in% names(object)))
    stop(paste('The stanreg object doesn\'t contain information about the projection.',
               'Run the projection first.'))
  lapply(object$proj, function(x) apply(x$b, 1, sd))
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
