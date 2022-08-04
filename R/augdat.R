#__________________________________________________________________________
# Helper functions for the augmented-data projection
#__________________________________________________________________________

#' Augmented-data projection: Internals
#'
#' The augmented-data projection makes extensive use of *augmented-rows
#' matrices* and *augmented-length vectors*. In the following, \eqn{N},
#' \eqn{C_{\mathrm{cat}}}{C_cat}, \eqn{C_{\mathrm{lat}}}{C_lat},
#' \eqn{S_{\mathrm{ref}}}{S_ref}, and \eqn{S_{\mathrm{prj}}}{S_prj} from help
#' topic [refmodel-init-get] are used. Furthermore, let \eqn{C} denote either
#' \eqn{C_{\mathrm{cat}}}{C_cat} or \eqn{C_{\mathrm{lat}}}{C_lat}, whichever is
#' appropriate in the context where it is used (e.g., for `ref_predfun`'s
#' output, \eqn{C = C_{\mathrm{lat}}}{C = C_lat}). Similarly, let \eqn{S} denote
#' either \eqn{S_{\mathrm{ref}}}{S_ref} or \eqn{S_{\mathrm{prj}}}{S_prj},
#' whichever is appropriate in the context where it is used. Then an
#' augmented-rows matrix is a matrix with \eqn{N \cdot C}{N * C} rows in \eqn{C}
#' blocks of \eqn{N} rows, i.e., with the \eqn{N} observations nested in the
#' \eqn{C} (latent) response categories. For ordered response categories, the
#' \eqn{C} (latent) response categories (i.e., the row blocks) have to be sorted
#' increasingly. The columns of an augmented-rows matrix have to correspond to
#' the \eqn{S} parameter draws, just like for the non-augmented-data projection.
#' An augmented-rows matrix is of class `augmat` (inheriting from classes
#' `matrix` and `array`) and needs to have the value of \eqn{N} stored in an
#' attribute called `nobs_orig`. An augmented-length vector (class `augvec`) is
#' the vector resulting from subsetting an augmented-rows matrix to extract a
#' single column and thereby dropping dimensions.
#'
#' @name augdat-internals
#' @keywords internal
NULL

# Convert a 3-dimensional array to an augmented-rows matrix.
#
# @param arr If `margin_draws` is `3`, a 3-dimensional array with dimensions N x
#   C x S. If `margin_draws` is `1`, a 3-dimensional array with dimensions S x N
#   x C. See above for a definition of these dimensions.
# @param margin_draws The index of `arr`'s margin which corresponds to the
#   posterior draws (i.e., the margin of length S). Restricted to values `1` and
#   `3`.
#
# @return An augmented-rows matrix (see above for a definition).
arr2augmat <- function(arr, margin_draws = 3) {
  stopifnot(is.array(arr) && length(dim(arr)) == 3)
  stopifnot(margin_draws %in% c(1, 3))
  if (margin_draws == 1) {
    margin_obs <- 2
  } else if (margin_draws == 3) {
    margin_obs <- 1
  }
  augmat <- apply(arr, margin_draws, as.vector, simplify = FALSE)
  augmat <- do.call(cbind, augmat)
  attr(augmat, "nobs_orig") <- dim(arr)[margin_obs]
  class(augmat) <- "augmat"
  return(augmat)
}

# Convert an augmented-rows matrix (see above for a definition) to a
# 3-dimensional array.
#
# @param augmat An augmented-rows matrix.
# @param nobs_orig The number of observations (N). Usually should not have to be
#   specified manually (i.e., the default should always work).
# @param margin_draws The index of the returned array's margin which shall
#   correspond to the posterior draws (i.e., the margin which shall be of
#   length S). Restricted to values `1` and `3`.
#
# @return If `margin_draws` is `3`, a 3-dimensional array with dimensions
#   N x C x S. If `margin_draws` is `1`, a 3-dimensional array with dimensions
#   S x N x C.
augmat2arr <- function(augmat,
                       nobs_orig = attr(augmat, "nobs_orig"),
                       margin_draws = 3) {
  stopifnot(inherits(augmat, "augmat"))
  stopifnot(!is.null(dim(augmat)))
  stopifnot(!is.null(nobs_orig))
  stopifnot(margin_draws %in% c(1, 3))
  n_discr <- nrow(augmat) / nobs_orig
  stopifnot(.is.wholenumber(n_discr))
  n_discr <- as.integer(round(n_discr))
  arr <- array(augmat, dim = c(nobs_orig, n_discr, ncol(augmat)))
  if (margin_draws == 1) {
    arr <- aperm(arr, perm = c(3, 1, 2))
  }
  return(arr)
}

# A t() method for class `augmat`, dropping the class and the `nobs_orig`
# attribute. This is necessary for clustering with kmeans(), for example.
#' @noRd
#' @export
t.augmat <- function(x) {
  class(x) <- NULL
  attr(x, "nobs_orig") <- NULL
  return(t(x))
}

# A t() method for class `augvec`, dropping the class and the `nobs_orig`
# attribute. This should not be necessary, but it's probably safer to have such
# a method (to avoid that the attributes are carried around after a t() call).
#' @noRd
#' @export
t.augvec <- function(x) {
  class(x) <- NULL
  attr(x, "nobs_orig") <- NULL
  return(t(x))
}

# A method for subsetting an object of class `augmat` (mainly following
# `[.factor`). This method keeps the `nobs_orig` attribute. It also keeps the
# class, except if the result is a vector (in which case the class is changed
# from `augmat` to `augvec`).
#' @noRd
#' @export
`[.augmat` <- function(x, ..., drop = TRUE) {
  x_out <- NextMethod("[")
  attr(x_out, "nobs_orig") <- attr(x, "nobs_orig")
  cls_out <- oldClass(x)
  if (is.null(dim(x_out))) {
    cls_out <- sub("augmat", "augvec", cls_out, fixed = TRUE)
  }
  class(x_out) <- cls_out
  return(x_out)
}

# A method for subsetting an object of class `augvec` (mainly following
# `[.factor`). This method keeps the `nobs_orig` attribute and the class. It
# should not be necessary, but it's probably safer to have it.
#' @noRd
#' @export
`[.augvec` <- function(x, ..., drop = TRUE) {
  x_out <- NextMethod("[")
  attr(x_out, "nobs_orig") <- attr(x, "nobs_orig")
  class(x_out) <- oldClass(x)
  return(x_out)
}

# Convert an augmented-length vector to an augmented-rows matrix.
#
# @param augvec An augmented-length vector (see above for a definition).
#
# @return An augmented-rows matrix (see above for a definition) with a single
#   column.
augvec2augmat <- function(augvec) {
  stopifnot(inherits(augvec, "augvec"))
  return(structure(
    as.matrix(augvec),
    nobs_orig = attr(augvec, "nobs_orig"),
    class = sub("augvec", "augmat", oldClass(augvec), fixed = TRUE)
  ))
}

# Convert an augmented-rows matrix (with a single column) to an augmented-length
# vector.
#
# @param augmat An augmented-rows matrix (see above for a definition) with a
#   single column.
#
# @return An augmented-length vector (see above for a definition).
augmat2augvec <- function(augmat) {
  stopifnot(inherits(augmat, "augmat"))
  stopifnot(identical(ncol(augmat), 1L))
  return(augmat[, 1])
}

# Find the maximum-probability category for each observation (with "observation"
# meaning one of the N original observations, not one of the \eqn{N \cdot C}{N *
# C} augmented observations).
#
# @param augvec An augmented-length vector (see above for a definition)
#   containing the probabilities for the response categories.
# @param lvls The response levels (as a character vector).
#
# @return A `factor` consisting of the maximum-probability categories. The
#   levels of this `factor` are those from `lvls`.
catmaxprb <- function(augvec, lvls) {
  arr <- augmat2arr(augvec2augmat(augvec))
  idxmaxprb <- do.call(c, lapply(seq_len(dim(arr)[1]), function(i_obs) {
    which.max(arr[i_obs, , 1])
  }))
  return(factor(lvls[idxmaxprb], levels = lvls))
}

# Link and inverse-link functions with array as input and output ----------

#' Link function for augmented-data projection with binomial family
#'
#' This is the function which has to be supplied to [extend_family()]'s argument
#' `augdat_link` in case of the augmented-data projection for the [binomial()]
#' family.
#'
#' @param prb_arr An array as described in section "Augmented-data projection"
#'   of [extend_family()]'s documentation.
#' @param link The same as argument `link` of [binomial()].
#'
#' @return An array as described in section "Augmented-data projection" of
#'   [extend_family()]'s documentation.
#'
#' @export
augdat_link_binom <- function(prb_arr, link = "logit") {
  basic_link <- binomial(link = link)$linkfun
  return(basic_link(prb_arr[, , -1, drop = FALSE]))
}

#' Inverse-link function for augmented-data projection with binomial family
#'
#' This is the function which has to be supplied to [extend_family()]'s argument
#' `augdat_ilink` in case of the augmented-data projection for the [binomial()]
#' family.
#'
#' @param eta_arr An array as described in section "Augmented-data projection"
#'   of [extend_family()]'s documentation.
#' @param link The same as argument `link` of [binomial()].
#'
#' @return An array as described in section "Augmented-data projection" of
#'   [extend_family()]'s documentation.
#'
#' @export
augdat_ilink_binom <- function(eta_arr, link = "logit") {
  basic_ilink <- binomial(link = link)$linkinv
  prb_arr1 <- basic_ilink(eta_arr)
  prb_arr0 <- 1 - prb_arr1
  stopifnot(identical(dim(prb_arr0), dim(prb_arr1)))
  stopifnot(identical(dim(prb_arr1)[3], 1L))
  return(array(c(prb_arr0, prb_arr1), dim = c(dim(prb_arr1)[-3], 2L)))
}

## From brms --------------------------------------------------------------
## The functions from this (sub-)section are copied over from brms (with consent
## by Paul Buerkner) to avoid loading brms just for these special link and
## inverse-link functions. (After copying over, they have been slightly modified
## here to avoid dependencies on other brms-internal functions.)

augdat_link_cumul <- function(prb_arr, link) {
  ncat <- utils::tail(dim(prb_arr), 1)
  cumprb_arr <- apply(prb_arr[, , -ncat, drop = FALSE], c(1, 2), cumsum)
  cumprb_arr <- array(cumprb_arr,
                      dim = c(ncat - 1, utils::head(dim(prb_arr), -1)))
  cumprb_arr <- aperm(cumprb_arr, perm = c(c(1, 2) + 1, 1))
  return(linkfun_raw(cumprb_arr, link_nm = link))
}

augdat_ilink_cumul <- function(eta_arr, link) {
  cumprb_arr <- ilinkfun_raw(eta_arr, link_nm = link)
  dim_noncat <- utils::head(dim(cumprb_arr), -1)
  ones_arr <- array(1, dim = c(dim_noncat, 1))
  zeros_arr <- array(0, dim = c(dim_noncat, 1))
  return(abind::abind(cumprb_arr, ones_arr) -
           abind::abind(zeros_arr, cumprb_arr))
}
