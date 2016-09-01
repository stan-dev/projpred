#' Variable selection for generalized linear models
#'
#' Perform the projection predictive variable selection for a generalized
#' linear model fitted with rstanarm.
#' @param fit A \link[=stanreg-objects]{stanreg} object.
#' @param d_test A test dataset, which is used to evaluate model performance. If
#' not provided, training data is used.
#' @param ... Optional arguments. Possible arguments and their defaults are:
#' \describe{
#'  \item{\code{ns = min(400, [number of samples])}}{
#'    Number of samples used in the variable selection.
#'    Cannot be larger than the number of samples.}
#'  \item{\code{nc = 0}}{
#'    If nonzero, samples are clustered and the cluster centers are
#'    used in the variable selection instead of the actual samples.
#'  }
#'  \item{\code{nv = min(ncol(x) - 1, rankMatrix(x))}}{
#'    Maximum number of features to be used in the projection (incl. intercept).
#'    Cannot be larger than \code{min(ncol(x) - 1, rankMatrix(x))}.}
#'  \item{\code{verbose = FALSE}}{
#'    If \code{verbose = TRUE}, prints information about the progress of the
#'    variable selection.}
#' }
#'
#' @return A list with class \code{'varsel'} containing the following elements:
#' \describe{
#'  \item{\code{chosen}}{The order in which the features were added to the submodel.}
#'  \item{\code{stats}}{An array with statistics of the submodel performance.}
#'  \item{\code{family}}{A \code{\link{family}}-object.}
#' }
#'
#' @examples
#' \dontrun{
#' ### Usage with stanreg objects
#' fit <- stan_glm(y~x, binomial())
#' vars <- varsel(fit)
#' plot(vars)
#' }
#'


#' @export
varsel <- function(fit, d_test = NA, ...) {
  UseMethod('varsel')
}

#' @export
varsel.stanreg <- function(fit, d_test = NA, ...) {

  vars <- .extract_vars(fit)
  args <- .init_args(list(...), vars, family(fit))

  if(ncol(vars$x) < 2)
    stop('data must have at least 2 features.')
  if(!(args$family_kl$family %in% c('gaussian','binomial','poisson')))
    stop(paste0(args$family_kl$family, 'family currently not supported'))

  d_train <- list(x = vars$x, w = vars$w, offset = vars$offset)

  # if test data doesn't exist, use training data to evaluate mse, r2, mlpd
  eval_data <- ifelse(is.list(d_test), 'test', 'train')
  if(eval_data == 'train') {
    d_test <- within(d_train, y <- vars$y)
  } else {
    if(args$intercept && any(d_test$x[,1] != 1)) d_test$x <- cbind(1, d_test$x)
    d_test$w <- d_test$w %ORifNULL% rep(1, nrow(d_test$x))
    d_test$offset <- d_test$offset %ORifNULL% rep(0, nrow(d_test$x))
  }

  b <- vars$b
  mu <- args$family_kl$linkinv(d_train$x%*%b + d_train$offset)
  dis <- vars$dis
  b0 <- matrix(rowMeans(b), ncol = 1)

  # Sample indices to be used with forward selection and final projection
  s_ind <- round(seq(1, args$ns_total, length.out  = args$ns))

  p_full <- list(mu = mu[, s_ind], dis = dis[s_ind], cluster_w = rep(1, args$ns))

  clust <- if(args$clust) .get_p_clust(mu, dis, args) else NULL

  sel <- fsel(p_full, d_train, d_test, clust$p, b0, args)

  # get parameters of the full and the submodel and combine them to
  # one list for bootstrapping the 95% intervals
  mu_full_samp <- args$family_kl$linkinv(d_test$x%*%b[, s_ind] + d_test$offset)
  mu_full <- list(rowMeans(mu_full_samp))
  lppd_full <- list(apply(
    args$family_kl$ll_fun(mu_full_samp, dis[s_ind], d_test$y, d_test$w),
    1, log_mean_exp))
  mu_sub_full <- c(sel$mu, mu_full)
  lppd_sub_full <- c(sel$lppd, lppd_full)

  # perform bootstrapping to get 95% intervals for the summaries
  nvs <- c(1:length(sel$mu), ncol(d_train$x)) - args$intercept
  n <- length(d_test$y)
  n_boot <- 1000
  b_weights <- matrix(rexp(n * n_boot, 1), ncol = n)
  b_weights <- b_weights/rowSums(b_weights)
  stats <- rbind(
    .bootstrap_stats(mu_sub_full, lppd_sub_full, d_test, nvs, args$family_kl,
                     b_weights, eval_data),
    data.frame(data = 'sel', nvar = nvs, delta = F, summary = 'kl',
               value = c(sel$kl, 0), lq = NA, uq = NA))

  res <- list(chosen = sel$chosen, stats = stats)

  # if clustering was performed, return it
  if(args$clust) res$cl <- clust$cl
  # if function was called by cv_varsel, return also mu and lppd,
  # otherwise return the family object
  if(args$cv) {
    res$mu <- mu_sub_full
    res$lppd <- lppd_sub_full
  } else {
    res$family <- family(fit)
  }

  structure(res, class = 'varsel')
}
