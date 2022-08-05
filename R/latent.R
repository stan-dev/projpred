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
#'   the draws. Draws that should be dropped need to have an `NA` in `cl`.
#' @param wdraws A numeric vector of length \eqn{S}, giving the weights of the
#'   draws. Internally, these weights are normalized to sum to `1` within each
#'   cluster. Draws that should be dropped can (but must not necessarily) have
#'   an `NA` in `wdraws`.
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
#' TODO
#'
#' @export
cl_agg <- function(draws,
                   cl = seq_len(nrow(draws)),
                   # Note: Actually, the `wdraws` should not be needed in the
                   # context of projpred because cl_agg() is meant to be applied
                   # to parameter draws from the reference model and these are
                   # assumed to have all the same weight (see init_refmodel()).
                   wdraws = rep(1 / nrow(draws), nrow(draws)),
                   eps_wdraws = 0) {
  n_cl <- max(cl, na.rm = TRUE)
  return(do.call(rbind, lapply(seq_len(n_cl), function(i_cl) {
    idxs_draws_i <- which(cl == i_cl)
    wdraws_i <- wdraws[idxs_draws_i] / sum(wdraws[idxs_draws_i]) *
      (1 - eps_wdraws)
    return(wdraws_i %*% draws[idxs_draws_i, , drop = FALSE])
  })))
}
