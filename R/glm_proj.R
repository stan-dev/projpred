#' GLM projection to submodel.
#'
#' @param fit A \link[=stanreg-objects]{stanreg} object.
#' @param ... Optional arguments. Possible arguments and their defaults are:
#' \describe{
#'  \item{\code{n_out = min(400, [number of samples])}}{
#'    Number of samples used in the projection.
#'    Cannot be larger than the number of samples.}
#'  \item{\code{n_sel = min(800, [number of samples])}}{
#'    Number of samples used in the variable selection.
#'    Cannot be larger than the number of samples.}
#'  \item{\code{avg = FALSE}}{
#'    A logical value indicating whether selection is done using KL divergence
#'    between average of the samples instead of average KL divergence for each
#'    sample. If \code{avg = TRUE}, \code{n_sel} is ignored.}
#'  \item{\code{d = min(ncol(x) - 1, rankMatrix(x))}}{
#'    Maximum number of features to be used in the projection (incl. intercept).
#'    Cannot be larger than \code{min(ncol(x) - 1, rankMatrix(x))}.}
#'  \item{\code{glmproj.cores = getOption('glmproj.cores', parallel::detectCores())}}{
#'    Number of cores used when averaging KL divergence over samples. This can be
#'    set for an entire R session by \code{options('glmproj.cores' = NUMBER)}.}
#'  \item{\code{verbose = FALSE}}{
#'    If \code{verbose = TRUE}, prints information about the progress of the
#'    variable selection (and about the progress of the projection if only a single core is used).}
#' }
#'
#' @return A list with class \code{'glmproj'} containing the following elements:
#' \describe{
#'  \item{\code{chosen}}{The order in which the features were added.}
#'  \item{\code{kl}}{KL divergences with corresponding features added to the model.}
#'  \item{\code{b}}{Projected weight vectors.}
#'  \item{\code{dis}}{Projected dispersion parameters if the model has a dispersion parameter.}
#'  \item{\code{family}}{A \code{\link{family}}-object.}
#' }
#'
#' @examples
#' \dontrun{
#' ### Usage with stanreg objects
#' # use 400 samples to select the model, project results
#' # using 800 samples:
#' gproj <- glm_proj(fit, n_sel = 400, n_out = 800)
#'
#' # or use average of the samples in the selection
#' # (significantly faster, but potentially less accurate)
#' gproj <- glm_proj(fit, n_out = 800, avg = TRUE)
#' }
#'
#' @importFrom rstan extract
#' @importFrom rstanarm get_x
#'
#' @export glm_proj

glm_proj <- function(fit, ...) {
  UseMethod('glm_proj')
}

#' @export

glm_proj.stanreg <- function(fit, ...) {

  args <- list(...)
  x <- get_x(fit)
  family <- fit$family

  e <- extract(fit$stanfit)
  b_p <- t(e$beta)
  if('alpha' %in% names(e)) b_p <- rbind(drop(e$alpha), b_p)
  dis_p <- e$dispersion

  w <- get_weights(fit)

  res <- glm_proj_helper(x, b_p, w, dis_p, family, args)
  structure(res, class = 'glmproj')
}
