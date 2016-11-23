#' @export
project <- function(object, nv, ...) {
  UseMethod('project')
}

#' @export
project.stanreg <- function(object, nv, ...) {
  if(is.null(nv)) stop('nv not provided')
  .validate_for_varsel(object)
  if(!('varsel' %in% names(object)))
    stop(paste('The stanreg object doesn\'t contain information about the variable',
               'selection. Run the variable selection first.'))
  if(max(nv) >= length(object$varsel$chosen))
    stop(paste('Cannot perform the projection with', max(nv), 'variables, because the',
               'variable selection has been run only up to',
               length(object$varsel$chosen), 'variables.'))

  vars <- .extract_vars(object)
  ns_total <- ncol(vars$b)
  args <- .init_args(list(...), vars, family(object))
  if(args$intercept) nv <- nv + 1
  v_inds_max <- object$varsel$chosen[1:max(nv)]
  if(args$intercept) vars$x <- cbind(1, vars$x)

  d_train <- list(x = vars$x[,v_inds_max],
                  weights = vars$weights,
                  offset = vars$offset)

  mu <- args$family_kl$linkinv(vars$x%*%vars$b + d_train$offset)
  dis <- vars$dis
  b0 <- matrix(coef(object)[v_inds_max], ncol = 1)

  s_ind <- round(seq(1, ns_total, length.out  = args$ns))
  p_full <- list(mu = mu[, s_ind], dis = dis[s_ind], cluster_w = rep(1/args$ns, args$ns))

  projfun <- .get_proj_handle(args$family_kl)

  object$proj <- lapply(nv, function(nv, p_full, d_train, b0, args, names) {
    vars <- projfun(1:nv, p_full, d_train, unname(b0), args)
    rownames(vars$b) <- names[1:nrow(vars$b)]
    vars$kl <- NULL
    vars
  }, p_full, d_train, b0, args, names(coef(object))[v_inds_max])

  object
}

#' @export
proj_linpred <- function(object, transform = FALSE, newdata = NULL, offset = NULL, ...) {
  .validate_for_varsel(object)
  if(!('proj' %in% names(object)))
    stop(paste('The stanreg object doesn\'t contain information about the projection.',
               'Run the projection first.'))

  dat <- rstanarm:::pp_data(object, newdata, offset = offset)

  lapply(object$proj, function(proj, dat, chosen) {
    res <- t(dat$x[, chosen[1:nrow(proj$b)], drop = F]%*%proj$b + dat$offset)
    if(transform) family(object)$linkinv(res) else res
  }, dat, object$varsel$chosen)

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
varsel_plot <- function(x, ..., nv = NULL, summaries = NULL, deltas = T, train = F) {
  if(!('varsel' %in% names(x)))
    stop(paste('Stanreg object doesn\'t contain information about the variable',
               'selection. Run the variable selection first!'))

  data_remove <- if(train) 'test' else 'train'
  if(is.null(summaries)) summaries <- as.character(unique(x$varsel$stats$summary))
  arr <- subset(x$varsel$stats, data != data_remove & (delta == deltas | summary == 'kl')
                & summary %in% summaries)

  if(nrow(arr) == 0)
    stop(paste0('summaries must contain at least one of the following values: ',
                paste0(unique(x$varsel$stats$summary), collapse = ', '), '.'))

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
