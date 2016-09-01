#' @export
project <- function(object, fit, size, ...) {
  UseMethod('project')
}

#' @export
project.varsel <- function(object, fit, size, ...) {
  vars <- .extract_vars(fit)
  ns_total <- ncol(vars$b)
  args <- .init_args(list(...), vars, family(fit))
  args$rank_x <- Inf # dont check if matrix is invertible anymore
  nc_sel <- length(object$cl$size)
  v_inds_max <- object$chosen[1:max(size)]

  d_train <- list(x = vars$x[,v_inds_max], w = vars$w, offset = vars$offset)

  b <- vars$b[v_inds_max, ]
  mu <- args$family_kl$linkinv(d_train$x%*%b + d_train$offset)
  dis <- vars$dis
  b0 <- matrix(rowMeans(b), ncol = 1)

  if(args$clust) {
    # if clustering has been already performed, use that
    cl <- if(args$nc == nc_sel) object$cl else NULL
    clust <- .get_p_clust(mu, dis, args, cl)
    p_full <- clust$p
  } else {
    s_ind <- round(seq(1, ns_total, length.out  = args$ns))
    p_full <- list(mu = mu[, s_ind], dis = dis[s_ind], cluster_w = rep(1, args$ns))
  }

  proj <- .get_proj_handle(args$family_kl$family)

  p_sub <- lapply(size, function(s, p_full, d_train, b0, args) {
    vars <- proj(NULL, 1:s, p_full, d_train, b0, args)
    vars$kl <- NULL
    vars
  }, p_full, d_train, b0, args)

  res <- list(proj_params = p_sub,
              chosen = v_inds_max,
              weights = p_full$cluster_w,
              intercept = vars$intercept,
              family = family(fit))

  structure(res, class = 'proj')
}

#' @export
predict.varsel <- function(object, fit, size, newdata, ...) {
  # does the projection first, if a projection object is not provided
  predict(project(object, fit, size, ...), newdata, ...)
}

#' @export
predict.proj <- function(object, newdata, ...) {
  if(!is.list(newdata)) newdata <- list(x = newdata)
  if(object$intercept) newdata$x <- cbind(1, newdata$x)
  if(!is.null(newdata$w)) newdata$x <- newdata$x*newdata$w

  lapply(object$proj_params, function(proj, newdata, inds, w, family) {
    eta <- drop(newdata$x[, object$chosen[1:nrow(proj$b)], drop = F]%*%(proj$b%*%w))
    if(!is.null(newdata$offset)) eta <- eta + newdata$offset
    family$linkinv(eta)
  }, newdata, object$chosen, object$weights, object$family)

}

#' @export
plot.varsel <- function(x, summaries = NULL, deltas = F, ...) {

  if(is.null(summaries)) {
    data <- subset(x$stats, data != 'train' & (delta == deltas | summary == 'kl'))
  } else {
    data <- subset(x$stats, data != 'train' & (delta == deltas | summary == 'kl')
                   & summary %in% summaries)
    if(nrow(data) == 0)
      stop(paste0('summaries must contain at least one of the following values: ',
                  paste0(unique(x$stats$summary), collapse = ','), '.'))
  }
  ylab <- if(deltas) 'value' else 'delta'
  ggplot(data = data, mapping = aes(x = nvar)) +
    geom_ribbon(aes(ymin = lq, ymax = uq), alpha = 0.3) +
    geom_line(aes(y = value)) +
    geom_hline(aes(yintercept = value), subset(data, nvar == max(nvar)),
               color = 'darkred') +
    labs(x = 'Number of variables in the submodel', y = ylab) +
    facet_grid(summary ~ ., scales = 'free_y')
}

#' @export
summary.varsel <- function(object, ..., digits = 3) {
  if('test' %in% object$stats$data) {
    summaries <- setdiff(unique(object$stats$summary), c('kl'))
    # suffixes are to ensure unique column names when doing merge.
    arr <- Reduce(function(x, y) merge(x, y, by = 'nvar', suffixes = c(-ncol(x),ncol(x))),
           c(lapply(summaries, function(sum, stats) {
             subset(stats, data == 'test' & summary == sum & delta, c('nvar', 'value'))
           }, object$stats),
           list(kl = subset(object$stats, summary == 'kl', c('nvar', 'value')))))
    arr <- cbind(arr,
                 chosen = c(object$chosen,NA),
                 pctch = c(object$pctch,NA))
    arr <- setNames(arr, c('size', summaries, 'kl', 'chosen', 'pctch'))
    print(arr, digits = digits, width = 12, right = T, row.names = F)
  } else {
    print(object, ..., digits = digits)
  }
}

#' @export
print.varsel <- function(x, digits = 3, ...) {
  cat('\nTable of the model size, the index',
      'of the variable added last')
  cols <- c('size','kl','chosen')
  if(!is.null(x$pctch)) {
    cols <- c(cols,'pctch')
    cat(',\nfraction of cv-runs that included the selected variable to a',
        '\nmodel of same size ')
  } else {cat('\n')}
  cat('and the respective KL divergence.\n\n')
  arr <- data.frame(
    subset(x$stats, summary == 'kl' & value > 0, c('nvar','value')),
    chosen = x$chosen)
  if(!is.null(x$pctch)) arr <- cbind(arr, pctch = x$pctch)
  names(arr)[1:2] <- c('size','kl')
  print(arr[,cols], digits = digits, width = 12, right = T, row.names = F)
}

#' @export
print.proj <- function(x, digits = 3, ...) {
  cat('\nProjected coefficients for the submodels.\n')
  lapply(x$proj_params, function(pars, weights, chosen, digits) {
    coefs <- round(drop(pars$b%*%weights)/length(weights), digits = digits)
    cat(paste0('\nModel size : ', NROW(pars$b), '.\nChosen variables: ',
               paste0(chosen[1:NROW(pars$b)], collapse = ' '), '.\nCoefs: ',
               paste0(coefs, collapse = ' '), '.\n'))
    # if model has dispersion, print dispersion...
  }, x$weights, x$chosen, digits)
  invisible(x)

}
