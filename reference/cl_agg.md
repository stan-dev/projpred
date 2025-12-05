# Weighted averaging within clusters of parameter draws

This function aggregates \\S\\ parameter draws that have been clustered
into \\S\_{\mathrm{cl}}\\ clusters by averaging across the draws that
belong to the same cluster. This averaging can be done in a weighted
fashion.

## Usage

``` r
cl_agg(
  draws,
  cl = seq_len(nrow(draws)),
  wdraws = rep(1, nrow(draws)),
  eps_wdraws = 0
)
```

## Arguments

- draws:

  An \\S \times P\\ matrix of parameter draws, with \\P\\ denoting the
  number of parameters.

- cl:

  A numeric vector of length \\S\\, giving the cluster indices for the
  draws. The cluster indices need to be values from the set \\\\1, ...,
  S\_{\mathrm{cl}}\\\\, except for draws that should be dropped (e.g.,
  by thinning), in which case `NA` needs to be provided at the positions
  of `cl` corresponding to these draws.

- wdraws:

  A numeric vector of length \\S\\, giving the weights of the draws. It
  doesn't matter whether these are normalized (i.e., sum to `1`) or not
  because internally, these weights are normalized to sum to `1` within
  each cluster. Draws that should be dropped (e.g., by thinning) can
  (but must not necessarily) have an `NA` in `wdraws`.

- eps_wdraws:

  A positive numeric value (typically small) which will be used to
  improve numerical stability: The weights of the draws within each
  cluster are multiplied by `1 - eps_wdraws`. The default of `0` should
  be fine for most cases; this argument only exists to help in those
  cases where numerical instabilities occur (which must be detected by
  the user; this function will not detect numerical instabilities
  itself).

## Value

An \\S\_{\mathrm{cl}} \times P\\ matrix of aggregated parameter draws.

## Examples

``` r
set.seed(323)
S <- 100L
P <- 3L
draws <- matrix(rnorm(S * P), nrow = S, ncol = P)
# Clustering example:
S_cl <- 10L
cl_draws <- sample.int(S_cl, size = S, replace = TRUE)
draws_cl <- cl_agg(draws, cl = cl_draws)
# Clustering example with nonconstant `wdraws`:
w_draws <- rgamma(S, shape = 4)
draws_cl <- cl_agg(draws, cl = cl_draws, wdraws = w_draws)
# Thinning example (implying constant `wdraws`):
S_th <- 50L
idxs_thin <- round(seq(1, S, length.out = S_th))
th_draws <- rep(NA, S)
th_draws[idxs_thin] <- seq_len(S_th)
draws_th <- cl_agg(draws, cl = th_draws)
```
