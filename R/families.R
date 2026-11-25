#' Extra family objects
#'
#' Family objects not in the set of default [`family`] objects.
#'
#' @name extra-families
#'
#' @param link Name of the link function. In contrast to the default [`family`]
#'   objects, this has to be a character string here.
#' @param nu Degrees of freedom for the Student-\eqn{t} distribution.
#'
#' @details [Student_t()] is deprecated and will be removed in a future
#'   release. [project()], [varsel()], and [cv_varsel()] do not support it.
#'   For a Student-\eqn{t} reference model, use the latent projection (see
#'   [extend_family()] and section ["Example: Student-t
#'   distribution"](https://mc-stan.org/projpred/articles/latent.html#studentex)
#'   of the latent-projection vignette).
#'
#' @return A family object analogous to those described in [`family`].
#'
NULL

# Define a student-t family object. Dispersion is defined to be the scale
# parameter of the distribution.

#' @rdname extra-families
#' @export
Student_t <- function(link = "identity", nu = 3) {
  warning("Student_t() is deprecated and will be removed in a future ",
          "release. For a Student-t reference model, use the latent ",
          "projection (see the latent-projection vignette).")
  if (!(link %in% c("identity", "log", "inverse"))) {
    stop("Non-supported link: ", link)
  }
  if (!is.character(link)) {
    stop("Link must be a string.")
  }

  ## fetch the link statistics
  stats <- make.link(link)

  ## variance function
  varfun <- function(mu) {
    if (nu > 2) {
      rep(nu / (nu - 2), length(mu))
    } else {
      rep(Inf, length(mu))
    }
  }

  ## create the object and append the relevant fields
  fam <- nlist(
    family = "Student_t",
    nu,
    link,
    linkfun = stats$linkfun,
    linkinv = stats$linkinv,
    variance = varfun,
    dev.resids = function(y, mu, wt, dis = 1) {
      wt * (nu + 1) * log(1 + 1 / nu * ((y - mu) / dis)^2)
    },
    aic = function(y, n, mu, wt, dev) {
      stop("aic not implemented for Student-t.")
    },
    mu.eta = stats$mu.eta,
    initialize = expression({
      stop("initialization for Student-t not implemented.")
    }),
    validmu = function(mu) {
      TRUE
    },
    valideta = stats$valideta
  )

  return(structure(fam, class = "family"))
}
