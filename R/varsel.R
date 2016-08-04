#' @export varsel

varsel <- function(fit, d_test = NA, ...) {
  UseMethod('varsel')
}

varsel.stanreg <- function(fit, d_test = NA, ...) {

  args <- list(...)

  params <- extract_params(fit)
  args$family <- kl_helpers(family(fit))
  if(is.list(d_test) && is.null(d_test$w)) d_test$w <- rep(1, nrow(d_test$x))

  if(ncol(params$x) < 2)
    stop('data must have at least 2 features.')
  if(!(family(fit)$family %in% c('gaussian','binomial','poisson')))
    stop(paste0(family(fit)$family, 'family not yet supported'))

  args$intercept <- params$intercept

  ns_total <- ncol(params$b)

  # set additional argumets to defaults if they are missing etc.
  args$n_s <- min(ifelse(is.null(args$n_s), 400, args$n_s), ns_total)
  args$d <- min(ncol(params$x)-1, args$d, rankMatrix(params$x))
  if(is.null(args$avg)) args$avg <- FALSE
  if(is.null(args$verbose)) args$verbose <- FALSE
  if(is.list(d_test) && is.null(d_test$w)) d_test$w <- 1

  d <- list(x = params$x, w = params$w)
  # Sample indices to be used with forward selection and final projection
  sind <- round(seq(1, ns_total, length.out  = args$n_s))
  b0 <- matrix(rowMeans(params$b), ncol = 1)
  p <- with(params, list(mu = args$family$linkinv(x%*%b[,sind]), dis = dis[sind]))

  if(args$avg) {
    # calculate the means of the variables
    p_means <- with(params, list(mu = matrix(rowMeans(args$family$linkinv(x%*%b)), ncol = 1),
                                 dis = sqrt(mean(dis^2))))
  } else {
    p_means <- NULL
  }


  sel <- fsel(p, d, d_test, p_means, b0, args)

  res <- list(submodel = as.data.frame(sel), full = list())

  if(is.list(d_test)) {
    test_stats <- test_stat(list(b = params$b[,sind], dis = params$dis[sind]),
                            1:ncol(d_test$x), d_test, args$family)
    res$full <- c(res$full, test_stats)
    res$full <- res$full[names(res$full) %in% names(res$submodel)]
  }

  res$family <- family(fit)

  structure(res, class = 'varsel')
}
