#' @export varsel varsel.stanreg

varsel <- function(fit, d_test = NA, ...) {
  UseMethod('varsel')
}

varsel.stanreg <- function(fit, d_test = NA, ...) {

  args <- list(...)

  params <- extract_params(fit)
  args$family_kl <- kl_helpers(family(fit))

  if(ncol(params$x) < 2)
    stop('data must have at least 2 features.')
  if(!(args$family_kl$family %in% c('gaussian','binomial','poisson')))
    stop(paste0(args$family_kl$family, 'family not yet supported'))

  args$intercept <- params$intercept

  ns_total <- ncol(params$b)

  # set additional argumets to defaults if they are missing etc.
  args$ns <- min(ifelse(is.null(args$ns), 400, args$ns), ns_total)
  args$rank_x <- rankMatrix(params$x)
  args$nv <- min(ncol(params$x) - 1, args$nv, args$cv_nv, args$rank_x - 1)
  if(is.null(args$avg)) args$avg <- FALSE
  if(is.null(args$verbose)) args$verbose <- FALSE
  args$cv <- ifelse(is.null(args$cv), F, args$cv)

  d <- list(x = params$x, w = params$w, offset = params$offset)
  # if test data doesn't exist, use training data to evaluate mse, r2, mlpd
  eval_data <- ifelse(is.list(d_test), 'test', 'train')
  if(eval_data == 'train') {
    d_test <- within(d, y <- params$y)
  } else {
    if(args$intercept && ncol(d$x) == ncol(d_test$x) + 1) d_test$x <- cbind(1, d_test$x)
    if(is.null(d_test$w)) d_test$w <- rep(1, nrow(d_test$x))
    if(is.null(d_test$offset)) d_test$offset <- rep(0, nrow(d_test$x))
  }

  # Sample indices to be used with forward selection and final projection
  s_ind <- round(seq(1, ns_total, length.out  = args$ns))
  b <- params$b[, s_ind]
  dis <- params$dis[s_ind]
  b0 <- matrix(rowMeans(b), ncol = 1)
  p <- list(mu = args$family_kl$linkinv(d$x%*%b + d$offset), dis = dis)

  if(args$avg) {
    # calculate the means of the variables
    p_means <- with(params, list(mu = matrix(rowMeans(args$family_kl$linkinv(x%*%b + offset)), ncol = 1),
                                 dis = sqrt(mean(dis^2))))
  } else {
    p_means <- NULL
  }

  sel <- fsel(p, d, d_test, p_means, b0, args)

  mu_full_samp <- args$family_kl$linkinv(d_test$x%*%b)
  mu_full <- list(rowMeans(mu_full_samp))
  lppd_full <- list(apply(args$family_kl$ll_fun(mu_full_samp, dis, d_test$y, d_test$w), 1, log_mean_exp))
  mu <- c(sel$mu, mu_full)
  lppd <- c(sel$lppd, lppd_full)

  nvs <- c(1:length(sel$mu), ncol(d$x)) - args$intercept
  n <- length(d_test$y)
  n_boot <- 1000
  b_weights <- matrix(rexp(n * n_boot, 1), ncol = n)
  b_weights <- b_weights/rowSums(b_weights)

  res <- list(chosen = sel$chosen,
              stats = rbind(summary_stats(mu, lppd, d_test, nvs, args$family_kl, b_weights, eval_data),
                            data.frame(data = 'train', nvar = nvs, summary = 'kl', value = c(sel$kl, 0), lq = NA, uq = NA)))

  if(args$cv) {
    res$mu <- mu
    res$lppd <- lppd
  } else {
    res$family <- family(fit)
  }

  structure(res, class = 'varsel')
}
