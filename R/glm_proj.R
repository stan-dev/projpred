#' GLM projection to submodel.
#'
#' @param fit Either a stanreg object or a list with the model matrix \code{x}, the observation
#'  weights \code{w}, a \code{\link{family}}-object and the dispersion parameter \code{dis}, if applicable.
#' @param ... Optional arguments. Possible arguments and their defaults are:
#' \describe{
#'  \item{\code{n_out = min(400, [number of samples])}}{
#'    Number of samples used in the projection.}
#'  \item{\code{n_sel = min(800, [number of samples])}}{
#'    Number of samples used in the variable selection.}
#'  \item{\code{avg = FALSE}}{
#'    A logical value indicating whether selection is done using KL divergence
#'    between average of the samples instead of average KL divergence for each
#'    sample. If \code{TRUE}, \code{n_sel} is ignored.}
#'  \item{\code{d = ncol(x) - 1}}{
#'    Maximum number of features to be used in the projection (incl. intercept).}
#'  \item{\code{glmproj.cores = getOption('glmproj.cores', parallel::detectCores())}}{
#'    Number of cores used when averaging KL divergence over samples. This can be
#'    set for an entire R session by \code{options('glmproj.cores' = NUMBER)}.}
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
#' @importFrom rstanarm get_x
#'
#' @export glm_proj

glm_proj <- function(fit, ...) {
  UseMethod('glm_proj')
}

#' @export

glm_proj.stanreg <- function(fit, ...) {

  args <- list(...)
  x <- unname(get_x(fit))
  b_p <- unname(t(apply(as.array(fit$stanfit), 3, rbind)[, colnames(get_x(fit)), drop = FALSE]))
  w <- get_weights(fit)
  dis_p <- get_dispersion(fit)
  family <- fit$family

  res <- glm_proj_fsel(x, b_p, w, dis_p, family, args)
  structure(res, class = 'glmproj')
}

#' @export

glm_proj.list <- function(fit, ...) {

  if(any(!(c('x','b','family') %in% names(fit))))
    stop('fit must contain elements x, b and family')

  args <- list(...)
  x <- fit$x
  b_p <- fit$b
  w <- fit$w
  dis_p <- fit$dis
  family <- fit$family

  if(is.null(dis_p)) dis_p <- NA
  if(is.null(w)) w <- rep(1, nrow(x))
  if(family$family %in% c('gaussian','Gamma') && is.na(dis_p))
    stop('dispersion parameter not provided')

  res <- glm_proj_fsel(x, b_p, w, dis_p, family, args)
  structure(res, class = 'glmproj')
}
