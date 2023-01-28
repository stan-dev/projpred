#__________________________________________________________________________
# Helper functions for the latent projection
#__________________________________________________________________________

# Internal defaults for `latent_ll_oscale` and `latent_ppd_oscale` --------
# These are the functions which would have to be supplied to extend_family()'s
# arguments `latent_ll_oscale` and `latent_ppd_oscale` in certain situations of
# the latent projection. Note the "*would* have to be supplied": These functions
# are used by default (internally) in the respective situations described below.

# Situation: If `family$cats` (*after* applying extend_family(); see
# extend_family()'s argument `latent_y_unqs`) is not `NULL`.
latent_ll_oscale_cats <- function(ilpreds, y_oscale,
                                  wobs = rep(1, length(y_oscale)), cl_ref,
                                  wdraws_ref = rep(1, length(cl_ref))) {
  return(ll_cats(ilpreds, margin_draws = 1, y = y_oscale, wobs = wobs))
}
latent_ppd_oscale_cats <- function(ilpreds_resamp, wobs, cl_ref,
                                   wdraws_ref = rep(1, length(cl_ref)),
                                   idxs_prjdraws) {
  return(ppd_cats(ilpreds_resamp, margin_draws = 1, wobs = wobs))
}

# Situation: For the binomial family if `family$cats` (*after* applying
# extend_family(); see extend_family()'s argument `latent_y_unqs`) is `NULL`.
latent_ll_oscale_binom_nocats <- function(ilpreds, y_oscale,
                                          wobs = rep(1, length(y_oscale)),
                                          cl_ref,
                                          wdraws_ref = rep(1, length(cl_ref))) {
  # Ensure finite log() values:
  ilpreds[ilpreds %in% c(0)] <- .Machine$double.eps
  ilpreds[ilpreds %in% c(1)] <- 1 - .Machine$double.eps

  ilpreds <- t(ilpreds)
  ll_unw <- y_oscale * log(ilpreds) + (1 - y_oscale) * log(1 - ilpreds)
  return(t(wobs * ll_unw))
}
latent_ppd_oscale_binom_nocats <- function(ilpreds_resamp, wobs, cl_ref,
                                           wdraws_ref = rep(1, length(cl_ref)),
                                           idxs_prjdraws) {
  ilpreds_resamp <- t(ilpreds_resamp)
  ppd <- rbinom(prod(dim(ilpreds_resamp)), size = wobs, prob = ilpreds_resamp)
  ppd <- matrix(ppd, ncol = length(wobs), byrow = TRUE)
  return(ppd)
}

# Situation: For the Poisson family.
latent_ll_oscale_poiss <- function(ilpreds, y_oscale,
                                   wobs = rep(1, length(y_oscale)), cl_ref,
                                   wdraws_ref = rep(1, length(cl_ref))) {
  ll_unw <- dpois(y_oscale, lambda = t(ilpreds), log = TRUE)
  return(t(wobs * ll_unw))
}
latent_ppd_oscale_poiss <- function(ilpreds_resamp, wobs, cl_ref,
                                    wdraws_ref = rep(1, length(cl_ref)),
                                    idxs_prjdraws) {
  ppd <- rpois(prod(dim(ilpreds_resamp)), lambda = ilpreds_resamp)
  ppd <- matrix(ppd, nrow = nrow(ilpreds_resamp), ncol = ncol(ilpreds_resamp))
  return(ppd)
}

# Situation: For a family for which response-scale log predictive density (LPD)
# values cannot or should not be calculated.
latent_ll_oscale_NA <- function(ilpreds, y_oscale,
                                wobs = rep(1, length(y_oscale)), cl_ref,
                                wdraws_ref = rep(1, length(cl_ref))) {
  return(array(dim = dim(ilpreds)[1:2]))
}
latent_ppd_oscale_NA <- function(ilpreds_resamp, wobs, cl_ref,
                                 wdraws_ref = rep(1, length(cl_ref)),
                                 idxs_prjdraws) {
  return(array(dim = dim(ilpreds_resamp)[1:2]))
}

# Other -------------------------------------------------------------------

#' Weighted averaging within clusters of parameter draws
#'
#' This function aggregates \eqn{S} parameter draws that have been clustered
#' into \eqn{S_{\mathrm{cl}}}{S_cl} clusters by averaging across the draws that
#' belong to the same cluster. This averaging can be done in a weighted fashion.
#'
#' @param draws An \eqn{S \times P}{S x P} matrix of parameter draws, with
#'   \eqn{P} denoting the number of parameters.
#' @param cl A numeric vector of length \eqn{S}, giving the cluster indices for
#'   the draws. Draws that should be dropped (e.g., by thinning) need to have an
#'   `NA` in `cl`.
#' @param wdraws A numeric vector of length \eqn{S}, giving the weights of the
#'   draws. It doesn't matter whether these are normalized (i.e., sum to `1`) or
#'   not because internally, these weights are normalized to sum to `1` within
#'   each cluster. Draws that should be dropped (e.g., by thinning) can (but
#'   must not necessarily) have an `NA` in `wdraws`.
#' @param eps_wdraws A positive numeric value (typically small) which will be
#'   used to improve numerical stability: The weights of the draws within each
#'   cluster are multiplied by `1 - eps_wdraws`. The default of `0` should be
#'   fine for most cases; this argument only exists to help in those cases where
#'   numerical instabilities occur (which must be detected by the user; this
#'   function will not detect numerical instabilities itself).
#'
#' @return An \eqn{S_{\mathrm{cl}} \times P}{S_cl x P} matrix of aggregated
#'   parameter draws.
#'
#' @examples
#' set.seed(323)
#' S <- 100L
#' P <- 3L
#' draws <- matrix(rnorm(S * P), nrow = S, ncol = P)
#' # Clustering example:
#' S_cl <- 10L
#' cl_draws <- sample.int(S_cl, size = S, replace = TRUE)
#' draws_cl <- cl_agg(draws, cl = cl_draws)
#' # Clustering example with nonconstant `wdraws`:
#' w_draws <- rgamma(S, shape = 4)
#' draws_cl <- cl_agg(draws, cl = cl_draws, wdraws = w_draws)
#' # Thinning example (implying constant `wdraws`):
#' S_th <- 50L
#' idxs_thin <- round(seq(1, S, length.out = S_th))
#' th_draws <- rep(NA, S)
#' th_draws[idxs_thin] <- seq_len(S_th)
#' draws_th <- cl_agg(draws, cl = th_draws)
#'
#' @export
cl_agg <- function(draws,
                   cl = seq_len(nrow(draws)),
                   # Note: Most of the time, `wdraws` is not needed in the
                   # context of projpred because cl_agg() is meant to be applied
                   # to parameter draws from the reference model and these are
                   # usually assumed to have the same weight (see
                   # init_refmodel()), except for the `validate_search = TRUE`
                   # case of loo_varsel() where the PSIS weights are used.
                   wdraws = rep(1, nrow(draws)),
                   eps_wdraws = 0) {
  n_cl <- max(cl, na.rm = TRUE)
  return(do.call(rbind, lapply(seq_len(n_cl), function(i_cl) {
    idxs_draws_i <- which(cl == i_cl)
    wdraws_i <- wdraws[idxs_draws_i] / sum(wdraws[idxs_draws_i]) *
      (1 - eps_wdraws)
    return(wdraws_i %*% draws[idxs_draws_i, , drop = FALSE])
  })))
}
