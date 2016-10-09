#' Variable selection for generalized linear models
#'
#' Perform the projection predictive variable selection for a generalized
#' linear model fitted with rstanarm.
#' @param fit A \link[=stanreg-objects]{stanreg} object.
#' @param d_test A test dataset, which is used to evaluate model performance. If
#' not provided, training data is used.
#' @param ... Optional arguments. Possible arguments and their defaults are:
#' \describe{
#'  \item{\code{ns = min(400, [number of draws])}}{
#'    Number of draws used in the variable selection.
#'    Cannot be larger than the number of draws in the full model.}
#'  \item{\code{nc = 0}}{
#'    If nonzero, a clustering with \code{nc} clusters is performed for
#'    the draws and the cluster centers are used in the variable selection
#'    instead of the actual draws.}
#'  \item{\code{nv = min(ncol(x) - 1, rankMatrix(x))}}{
#'    Maximum number of variables to be used in the projection (incl. intercept).
#'    Cannot be larger than \code{min(ncol(x) - 1, rankMatrix(x))}.}
#'  \item{\code{verbose = FALSE}}{
#'    If \code{verbose = TRUE}, prints information about the progress of the
#'    variable selection.}
#' }
#'
#' @return A list with class \code{'varsel'} containing the following elements:
#' \describe{
#'  \item{\code{chosen}}{The order in which the variables were added to the submodel.}
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

  .validate_for_varsel(fit)
  vars <- .extract_vars(fit)
  args <- .init_args(list(...), vars, family(fit))

  d_train <- list(x = if(args$intercept) cbind(1, vars$x) else vars$x,
                  weights = vars$weights, offset = vars$offset)

  # if test data doesn't exist, use training data to evaluate mse, r2, mlpd
  eval_data <- ifelse(is.list(d_test), 'test', 'train')
  if(eval_data == 'train') {
    d_test <- within(d_train, y <- vars$y)
  } else {
    if(args$intercept) d_test$x <- cbind(1, d_test$x)
    if(is.null(d_test$weights)) d_test$weights <- rep(1, nrow(d_test$x))
    if(is.null(d_test$offset)) d_test$offset <- rep(0, nrow(d_test$x))
  }

  b <- vars$b
  mu <- args$family_kl$linkinv(d_train$x%*%b + d_train$offset)
  dis <- vars$dis
  b0 <- matrix(unname(coef(fit)), ncol = 1)

  s_ind <- round(seq(1, args$ns_total, length.out  = args$ns))
  p_full <- list(mu = mu[, s_ind], dis = dis[s_ind], cluster_w = rep(1/args$ns, args$ns))
  clust <- if(args$clust) .get_p_clust(mu, dis, args) else NULL

  sel <- fsel(p_full, d_train, d_test, clust$p, b0, args)

  full_stats_temp <- .summary_stats(
    d_test, 1:nrow(vars$b), list(b = vars$b[, s_ind], dis = vars$dis[s_ind]), args)
  mu_all <- c(sel$mu, list(full_stats_temp$mu))
  lppd_all <- c(sel$lppd, list(full_stats_temp$lppd))
  nv <- c(1:length(sel$mu), nrow(b)) - args$intercept
  kl <- data.frame(data = 'sel', size = nv, delta = F, summary = 'kl',
                   value = c(sel$kl, 0), lq = NA, uq = NA)

  # if function was called by cv_varsel, return also mu and lppd,
  if(args$cv) {
    res <- list(chosen = sel$chosen, mu = mu_all, lppd = lppd_all, stats = kl)
    if(args$clust) res$cl <- clust$cl
    return(res)
  }

  # evaluate performance on test data and
  # use bayesian bootstrap to get 95% credible intervals
  b_weights <- .gen_bootstrap_ws(length(d_test$y), args$n_boot)
  stats <- rbind(kl, .bootstrap_stats(mu_all, lppd_all, nv, d_test, args$family_kl,
                                      b_weights, eval_data), make.row.names = F)

  res <- list(chosen = sel$chosen, stats = stats, family = family(fit))
  if(args$clust) res$cl <- clust$cl

  structure(res, class = 'varsel')
}
