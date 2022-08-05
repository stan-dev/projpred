#__________________________________________________________________________
# Helper functions for the latent projection
#__________________________________________________________________________

# TODO: Add tests.

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
#' # Thinning example:
#' S_th <- 50L
#' idxs_thin <- round(seq(1, S, length.out = S_th))
#' th_draws <- rep(NA, S)
#' th_draws[idxs_thin] <- seq_len(S_th)
#' draws_th <- cl_agg(draws, cl = th_draws)
#' # A thinning example with nonconstant `wdraws` doesn't make sense because in
#' # case of thinning, `wdraws` doesn't have any impact.
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
