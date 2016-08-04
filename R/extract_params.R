
extract_params <- function(fit) {
  res <- list()
  res$x <- unname(get_x(fit))
  res$intercept <- attr(fit$terms,'intercept')

  e <- extract(fit$stanfit)
  res$b <- e$beta
  if(res$intercept) res$b <- unname(cbind(e$alpha, res$b))
  res$b <- t(res$b)
  res$dis <- unname(e$dispersion)
  if(is.null(res$dis)) res$dis <- NA

  # get weights from the fit
  if('(weights)' %in% colnames(fit$model)) {
    res$w <- fit$model[,'(weights)']
  } else if(length(weights(fit)) > 0) {
    res$w <- weights(fit)
  } else {
    res$w <- rep(1, nrow(res$x))
  }

  # get y
  y <- unname(get_y(fit))
  if(is.matrix(y)) y <- y[,1]
  res$y <- y

  res
}
