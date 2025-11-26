# Run search and performance evaluation without cross-validation

Run the *search* part and the *evaluation* part for a projection
predictive variable selection. The search part determines the predictor
ranking (also known as solution path), i.e., the best submodel for each
submodel size (number of predictor terms). The evaluation part
determines the predictive performance of the submodels along the
predictor ranking. A special method is `varsel.vsel()` which re-uses the
search results from an earlier `varsel()` (or
[`cv_varsel()`](https://mc-stan.org/projpred/dev/reference/cv_varsel.md))
run, as illustrated in the main vignette.

## Usage

``` r
varsel(object, ...)

# Default S3 method
varsel(object, ...)

# S3 method for class 'vsel'
varsel(object, ...)

# S3 method for class 'refmodel'
varsel(
  object,
  d_test = NULL,
  method = "forward",
  ndraws = NULL,
  nclusters = 20,
  ndraws_pred = 400,
  nclusters_pred = NULL,
  refit_prj = !inherits(object, "datafit"),
  nterms_max = NULL,
  verbose = getOption("projpred.verbose", as.integer(interactive())),
  search_control = NULL,
  lambda_min_ratio = 1e-05,
  nlambda = 150,
  thresh = 1e-06,
  penalty = NULL,
  search_terms = NULL,
  search_out = NULL,
  seed = NA,
  ...
)
```

## Arguments

- object:

  An object of class `refmodel` (returned by
  [`get_refmodel()`](https://mc-stan.org/projpred/dev/reference/refmodel-init-get.md)
  or
  [`init_refmodel()`](https://mc-stan.org/projpred/dev/reference/refmodel-init-get.md))
  or an object that can be passed to argument `object` of
  [`get_refmodel()`](https://mc-stan.org/projpred/dev/reference/refmodel-init-get.md).

- ...:

  For `varsel.default()`: Arguments passed to
  [`get_refmodel()`](https://mc-stan.org/projpred/dev/reference/refmodel-init-get.md)
  as well as to `varsel.refmodel()`. For `varsel.vsel()`: Arguments
  passed to `varsel.refmodel()`. For `varsel.refmodel()`: Arguments
  passed to the divergence minimizer (see argument `div_minimizer` of
  [`init_refmodel()`](https://mc-stan.org/projpred/dev/reference/refmodel-init-get.md)
  as well as section "Draw-wise divergence minimizers" of
  [projpred-package](https://mc-stan.org/projpred/dev/reference/projpred-package.md))
  when refitting the submodels for the performance evaluation (if
  `refit_prj` is `TRUE`).

- d_test:

  A `list` of the structure outlined in section "Argument `d_test`"
  below, providing test data for evaluating the predictive performance
  of the submodels as well as of the reference model. If `NULL`, the
  training data is used.

- method:

  The method for the search part. Possible options are `"forward"` for
  forward search and `"L1"` for L1 search. See also section "Details"
  below.

- ndraws:

  Number of posterior draws used in the search part. Ignored if
  `nclusters` is not `NULL` or in case of L1 search (because L1 search
  always uses a single cluster). If both (`nclusters` and `ndraws`) are
  `NULL`, the number of posterior draws from the reference model is used
  for `ndraws`. See also section "Details" below.

- nclusters:

  Number of clusters of posterior draws used in the search part. Ignored
  in case of L1 search (because L1 search always uses a single cluster).
  For the meaning of `NULL`, see argument `ndraws`. See also section
  "Details" below.

- ndraws_pred:

  Only relevant if `refit_prj` is `TRUE`. Number of posterior draws used
  in the evaluation part. Ignored if `nclusters_pred` is not `NULL`. If
  both (`nclusters_pred` and `ndraws_pred`) are `NULL`, the number of
  posterior draws from the reference model is used for `ndraws_pred`.
  See also section "Details" below.

- nclusters_pred:

  Only relevant if `refit_prj` is `TRUE`. Number of clusters of
  posterior draws used in the evaluation part. For the meaning of
  `NULL`, see argument `ndraws_pred`. See also section "Details" below.

- refit_prj:

  For the evaluation part, should the projections onto the submodels
  along the predictor ranking be performed again using `ndraws_pred`
  draws or `nclusters_pred` clusters (`TRUE`) or should their
  projections from the search part, which used `ndraws` draws or
  `nclusters` clusters, be re-used (`FALSE`)?

- nterms_max:

  Maximum submodel size (number of predictor terms) up to which the
  search is continued. If `NULL`, then `min(19, D)` is used where `D` is
  the number of terms in the reference model (or in `search_terms`, if
  supplied). Note that `nterms_max` does not count the intercept, so use
  `nterms_max = 0` for the intercept-only model. (Correspondingly, `D`
  above does not count the intercept.)

- verbose:

  A single integer value from the set \\\\0, 1, 2, 3, 4\\\\ (for
  `varsel()`, \\3\\ and \\4\\ have the same effect), indicating how much
  information (if any) to print out during the computations. Higher
  values indicate that more information should be printed, `0`
  deactivates the verbose mode. Internally, argument `verbose` is
  coerced to integer via
  [`as.integer()`](https://rdrr.io/r/base/integer.html), so technically,
  a single logical value or a single numeric value work as well.

- search_control:

  A `list` of "control" arguments (i.e., tuning parameters) for the
  search. In case of forward search, these arguments are passed to the
  divergence minimizer (see argument `div_minimizer` of
  [`init_refmodel()`](https://mc-stan.org/projpred/dev/reference/refmodel-init-get.md)
  as well as section "Draw-wise divergence minimizers" of
  [projpred-package](https://mc-stan.org/projpred/dev/reference/projpred-package.md)).
  In case of forward search, `NULL` causes `...` to be used not only for
  the performance evaluation, but also for the search. In case of L1
  search, possible arguments are:

  - `lambda_min_ratio`: Ratio between the smallest and largest lambda in
    the L1-penalized search (default: `1e-5`). This parameter
    essentially determines how long the search is carried out, i.e., how
    large submodels are explored. No need to change this unless the
    program gives a warning about this.

  - `nlambda`: Number of values in the lambda grid for L1-penalized
    search (default: `150`). No need to change this unless the program
    gives a warning about this.

  - `thresh`: Convergence threshold when computing the L1 path (default:
    `1e-6`). Usually, there is no need to change this.

- lambda_min_ratio:

  Deprecated (please use `search_control` instead). Only relevant for L1
  search. Ratio between the smallest and largest lambda in the
  L1-penalized search. This parameter essentially determines how long
  the search is carried out, i.e., how large submodels are explored. No
  need to change this unless the program gives a warning about this.

- nlambda:

  Deprecated (please use `search_control` instead). Only relevant for L1
  search. Number of values in the lambda grid for L1-penalized search.
  No need to change this unless the program gives a warning about this.

- thresh:

  Deprecated (please use `search_control` instead). Only relevant for L1
  search. Convergence threshold when computing the L1 path. Usually,
  there is no need to change this.

- penalty:

  Only relevant for L1 search. A numeric vector determining the relative
  penalties or costs for the predictors. A value of `0` means that those
  predictors have no cost and will therefore be selected first, whereas
  `Inf` means those predictors will never be selected. If `NULL`, then
  `1` is used for each predictor.

- search_terms:

  Only relevant for forward search. A custom character vector of
  predictor term blocks to consider for the search. Section "Details"
  below describes more precisely what "predictor term block" means. The
  intercept (`"1"`) is always included internally via
  [`union()`](https://rdrr.io/r/base/sets.html), so there's no
  difference between including it explicitly or omitting it. The default
  `search_terms` considers all the terms in the reference model's
  formula.

- search_out:

  Intended for internal use.

- seed:

  Pseudorandom number generation (PRNG) seed by which the same results
  can be obtained again if needed. Passed to argument `seed` of
  [`set.seed()`](https://rdrr.io/r/base/Random.html), but can also be
  `NA` to not call [`set.seed()`](https://rdrr.io/r/base/Random.html) at
  all. If not `NA`, then the PRNG state is reset (to the state before
  calling `varsel()`) upon exiting `varsel()`. Here, `seed` is used for
  clustering the reference model's posterior draws (if
  `!is.null(nclusters)` or `!is.null(nclusters_pred)`) and for drawing
  new group-level effects when predicting from a multilevel submodel
  (however, not yet in case of a GAMM).

## Value

An object of class `vsel`. The elements of this object are not meant to
be accessed directly but instead via helper functions (see the main
vignette and
[projpred-package](https://mc-stan.org/projpred/dev/reference/projpred-package.md)).

## Details

Arguments `ndraws`, `nclusters`, `nclusters_pred`, and `ndraws_pred` are
automatically truncated at the number of posterior draws in the
reference model (which is `1` for `datafit`s). Using less draws or
clusters in `ndraws`, `nclusters`, `nclusters_pred`, or `ndraws_pred`
than posterior draws in the reference model may result in slightly
inaccurate projection performance. Increasing these arguments affects
the computation time linearly.

For argument `method`, there are some restrictions: For a reference
model with multilevel or additive formula terms or a reference model set
up for the augmented-data projection, only the forward search is
available. Furthermore, argument `search_terms` requires a forward
search to take effect.

L1 search is faster than forward search, but forward search may be more
accurate. Furthermore, forward search may find a sparser model with
comparable performance to that found by L1 search, but it may also
overfit when more predictors are added. This overfit can be detected by
running search validation (see
[`cv_varsel()`](https://mc-stan.org/projpred/dev/reference/cv_varsel.md)).

An L1 search may select an interaction term before all involved
lower-order interaction terms (including main-effect terms) have been
selected. In projpred versions \> 2.6.0, the resulting predictor ranking
is automatically modified so that the lower-order interaction terms come
before this interaction term, but if this is conceptually undesired,
choose the forward search instead.

The elements of the `search_terms` character vector don't need to be
individual predictor terms. Instead, they can be building blocks
consisting of several predictor terms connected by the `+` symbol. To
understand how these building blocks work, it is important to know how
projpred's forward search works: It starts with an empty vector `chosen`
which will later contain already selected predictor terms. Then, the
search iterates over model sizes \\j \in \\0, ..., J\\\\ (with \\J\\
denoting the maximum submodel size, not counting the intercept). The
candidate models at model size \\j\\ are constructed from those elements
from `search_terms` which yield model size \\j\\ when combined with the
`chosen` predictor terms. Note that sometimes, there may be no candidate
models for model size \\j\\. Also note that internally, `search_terms`
is expanded to include the intercept (`"1"`), so the first step of the
search (model size 0) always consists of the intercept-only model as the
only candidate.

As a `search_terms` example, consider a reference model with formula
`y ~ x1 + x2 + x3`. Then, to ensure that `x1` is always included in the
candidate models, specify
`search_terms = c("x1", "x1 + x2", "x1 + x3", "x1 + x2 + x3")` (or, in a
simpler way that leads to the same results,
`search_terms = c("x1", "x1 + x2", "x1 + x3")`, for which helper
function
[`force_search_terms()`](https://mc-stan.org/projpred/dev/reference/force_search_terms.md)
exists). This search would start with `y ~ 1` as the only candidate at
model size 0. At model size 1, `y ~ x1` would be the only candidate. At
model size 2, `y ~ x1 + x2` and `y ~ x1 + x3` would be the two
candidates. At the last model size of 3, `y ~ x1 + x2 + x3` would be the
only candidate. As another example, to exclude `x1` from the search,
specify `search_terms = c("x2", "x3", "x2 + x3")` (or, in a simpler way
that leads to the same results, `search_terms = c("x2", "x3")`).

## Argument `d_test`

If not `NULL`, then `d_test` needs to be a `list` with the following
elements:

- `data`: a `data.frame` containing the predictor variables for the test
  set. In case of the latent projection, this `data.frame` is also used
  for evaluating attribute `cens_var` of the `latent_ll_oscale`
  function, so if `cens_var` is not `NULL`, `data` also needs to contain
  the data for the variable from `cens_var`.

- `offset`: a numeric vector containing the offset values for the test
  set (if there is no offset, use a vector of zeros).

- `weights`: a numeric vector containing the observation weights for the
  test set (if there are no observation weights, use a vector of ones).

- `y`: a vector or a `factor` containing the response values for the
  test set. In case of the latent projection, this has to be a vector
  containing the *latent* response values, but it can also be a vector
  full of `NA`s if latent-scale post-processing is not needed.

- `y_oscale`: Only needs to be provided in case of the latent projection
  where this needs to be a vector or a `factor` containing the
  *original* (i.e., non-latent) response values for the test set.

## See also

[`cv_varsel()`](https://mc-stan.org/projpred/dev/reference/cv_varsel.md)

## Examples

``` r
# Data:
dat_gauss <- data.frame(y = df_gaussian$y, df_gaussian$x)

# The `stanreg` fit which will be used as the reference model (with small
# values for `chains` and `iter`, but only for technical reasons in this
# example; this is not recommended in general):
fit <- rstanarm::stan_glm(
  y ~ X1 + X2 + X3 + X4 + X5, family = gaussian(), data = dat_gauss,
  QR = TRUE, chains = 2, iter = 500, refresh = 0, seed = 9876
)

# Run varsel() (here without cross-validation, with L1 search, and with small
# values for `nterms_max` and `nclusters_pred`, but only for the sake of
# speed in this example; this is not recommended in general):
vs <- varsel(fit, method = "L1", nterms_max = 3, nclusters_pred = 10,
             seed = 5555)
# Now see, for example, `?print.vsel`, `?plot.vsel`, `?suggest_size.vsel`,
# and `?ranking` for possible post-processing functions.
```
