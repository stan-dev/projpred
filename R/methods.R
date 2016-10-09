#' @export
project <- function(object, fit, nv, ...) {
  UseMethod('project')
}

#' @export
project.varsel <- function(object, fit, nv, ...) {
  if(is.null(nv)) stop('nv not provided')
  if(is.null(fit)) stop('fit not provided')
  vars <- .extract_vars(fit)
  ns_total <- ncol(vars$b)
  args <- .init_args(list(...), vars, family(fit))
  if(args$intercept) nv <- nv + 1
  nc_sel <- length(object$cl$size)
  v_inds_max <- object$chosen[1:max(nv)]
  if(args$intercept) vars$x <- cbind(1, vars$x)

  d_train <- list(x = vars$x[,v_inds_max],
                  weights = vars$weights, offset = vars$offset)

  mu <- args$family_kl$linkinv(vars$x%*%vars$b + d_train$offset)
  dis <- vars$dis
  b0 <- matrix(unname(coef(fit))[v_inds_max], ncol = 1)

  if(args$clust) {
    # if clustering has been already performed, use that
    cl <- if(args$nc == nc_sel) object$cl else NULL
    clust <- .get_p_clust(mu, dis, args, cl)
    p_full <- clust$p
  } else {
    s_ind <- round(seq(1, ns_total, length.out  = args$ns))
    p_full <- list(mu = mu[, s_ind], dis = dis[s_ind], cluster_w = rep(1/args$ns, args$ns))
  }

  proj <- .get_proj_handle(args$family_kl)

  p_sub <- lapply(nv, function(s, p_full, d_train, b0, args) {
    vars <- proj(NULL, 1:s, p_full, d_train, b0, args)
    vars$kl <- NULL
    vars
  }, p_full, d_train, b0, args)

  res <- list(proj_params = p_sub, chosen = v_inds_max, family = family(fit),
              cluster_w = p_full$cluster_w, intercept = vars$intercept)

  structure(res, class = 'proj')
}

#' @export
predict.varsel <- function(object, fit, nv, newdata, ...) {
  # does the projection first, if a projection object is not provided
  predict(project(object, fit, nv, ...), newdata, ...)
}

#' @export
predict.proj <- function(object, newdata, ...) {
  if(!is.list(newdata)) newdata <- list(x = newdata)
  if(object$intercept) newdata$x <- cbind(1, newdata$x)
  if(is.null(newdata$offset)) newdata$offset <- rep(0, nrow(newdata$x))

  lapply(object$proj_params, function(proj, newdata, inds, w, family) {
    drop(newdata$x[, object$chosen[1:nrow(proj$b)], drop = F]%*%(proj$b%*%w)) +
      newdata$offset
  }, newdata, object$chosen, object$cluster_w, object$family)

}

#' @export
plot.varsel <- function(x, summaries = NULL, deltas = T, train = F, nv = NULL, ...) {

  data_remove <- if(train) 'test' else 'train'
  if(is.null(summaries)) {
    arr <- subset(x$stats, data != data_remove & (delta == deltas | summary == 'kl'))
  } else {
    arr <- subset(x$stats, data != data_remove & (delta == deltas | summary == 'kl')
                   & summary %in% summaries)
    if(nrow(arr) == 0)
      stop(paste0('summaries must contain at least one of the following values: ',
                  paste0(unique(x$stats$summary), collapse = ','), '.'))
  }
  if(is.null(nv)) nv <- max(arr$size)
  ylab <- if(deltas) 'delta' else 'value'
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
summary.varsel <- function(object, nv = NULL, ..., digits = 3) {
  if(is.null(nv)) nv <- max(object$stats$size)
  if('test' %in% object$stats$data) {
    summaries <- setdiff(unique(object$stats$summary), c('kl'))
    # suffixes are to ensure unique column names when doing merge.
    arr <- Reduce(function(x, y) merge(x, y, by = 'size', suffixes = c(-ncol(x),ncol(x))),
           c(lapply(summaries, function(sum, stats) {
             subset(stats, data == 'test' & summary == sum & !delta, c('size', 'value'))
           }, object$stats),
           list(kl = subset(object$stats, summary == 'kl', c('size', 'value')))))
    arr <- cbind(arr,
                 chosen = c(object$chosen,NA),
                 pctch = c(object$pctch,NA))
    arr <- setNames(arr, c('size', summaries, 'kl', 'chosen', 'pctch'))
  } else {
    arr <- data.frame(subset(x$stats, summary == 'kl' & value > 0,
                             c('size','value')), chosen = x$chosen)
    if(!is.null(x$pctch)) arr$pctch <- x$pctch
    names(arr)[2] <- 'kl'
  }
  subset(arr, size <= nv)
}

#' @export
print.varsel <- function(x, digits = 3, nv = NULL, ...) {
  # switch from kl & value > 0 to size < max(size)
  cat('Table of the model size, the index',
      'of the variable added last')
  if(is.null(nv)) nv <- max(x$stats$size)
  if(!is.null(x$pctch)) {
    cat(',\nfraction of cv-runs that included the selected variable to a',
        '\nmodel of same size ')
  } else {cat('\n')}
  cat('and the respective KL divergence.\n\n')
  arr <- data.frame(
    subset(x$stats, summary == 'kl' & value > 0, c('size','value')))
  arr$chosen <- x$chosen
  if(!is.null(x$pctch)) arr$pctch <- x$pctch
  names(arr)[2] <- 'kl'
  print(subset(arr, size <= nv), digits = digits, width = 12, right = T, row.names = F)
}

#' @export
print.proj <- function(x, digits = 5, ...) {
  cat('Projected coefficients for the submodels.\n')
  lapply(x$proj_params, function(pars, weights, chosen, digits) {
    coefs <- round(drop(pars$b%*%weights), digits = digits)
    cat(paste0('\nModel size : ', NROW(pars$b) - x$intercept, '.\nChosen variables: ',
               paste0(chosen[1:NROW(pars$b)], collapse = ' '), '.\nCoefficients: ',
               paste0(coefs, collapse = ' '), '.\n'))
    if('dis' %in% names(pars))
      cat(paste0('Dispersion parameter: ', round(sqrt(pars$dis^2%*%weights), digits), '.\n'))
  }, x$cluster_w, x$chosen, digits)
  invisible(x)

}
