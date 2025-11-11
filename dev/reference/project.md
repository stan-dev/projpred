# Projection onto submodel(s)

Project the posterior of the reference model onto the parameter space of
a single submodel consisting of a specific combination of predictor
terms or (after variable selection) onto the parameter space of a single
or multiple submodels of specific sizes.

## Usage

``` r
project(
  object,
  nterms = NULL,
  solution_terms = predictor_terms,
  predictor_terms = NULL,
  refit_prj = TRUE,
  ndraws = 400,
  nclusters = NULL,
  seed = NA,
  verbose = getOption("projpred.verbose", as.integer(interactive())),
  ...
)
```

## Arguments

- object:

  An object which can be used as input to
  [`get_refmodel()`](https://mc-stan.org/projpred/dev/reference/refmodel-init-get.md)
  (in particular, objects of class `refmodel`).

- nterms:

  Only relevant if `object` is of class `vsel` (returned by
  [`varsel()`](https://mc-stan.org/projpred/dev/reference/varsel.md) or
  [`cv_varsel()`](https://mc-stan.org/projpred/dev/reference/cv_varsel.md)).
  Ignored if `!is.null(predictor_terms)`. Number of terms for the
  submodel (the corresponding combination of predictor terms is taken
  from `object`). If a numeric vector, then the projection is performed
  for each element of this vector. If `NULL` (and
  `is.null(predictor_terms)`), then the value suggested by
  [`suggest_size()`](https://mc-stan.org/projpred/dev/reference/suggest_size.md)
  is taken (with default arguments for
  [`suggest_size()`](https://mc-stan.org/projpred/dev/reference/suggest_size.md),
  implying that this suggested size is based on the ELPD). Note that
  `nterms` does not count the intercept, so use `nterms = 0` for the
  intercept-only model.

- solution_terms:

  Deprecated. Please use argument `predictor_terms` instead.

- predictor_terms:

  If not `NULL`, then this needs to be a character vector of predictor
  terms for the submodel onto which the projection will be performed.
  Argument `nterms` is ignored in that case. For an `object` which is
  not of class `vsel`, `predictor_terms` must not be `NULL`.

- refit_prj:

  A single logical value indicating whether to fit the submodels (again)
  (`TRUE`) or—if `object` is of class `vsel`—to re-use the submodel fits
  from the full-data search that was run when creating `object`
  (`FALSE`). For an `object` which is not of class `vsel`, `refit_prj`
  must be `TRUE`. See also section "Details" below.

- ndraws:

  Only relevant if `refit_prj` is `TRUE`. Number of posterior draws to
  be projected. Ignored if `nclusters` is not `NULL` or if the reference
  model is of class `datafit` (in which case one cluster is used). If
  both (`nclusters` and `ndraws`) are `NULL`, the number of posterior
  draws from the reference model is used for `ndraws`. See also section
  "Details" below.

- nclusters:

  Only relevant if `refit_prj` is `TRUE`. Number of clusters of
  posterior draws to be projected. Ignored if the reference model is of
  class `datafit` (in which case one cluster is used). For the meaning
  of `NULL`, see argument `ndraws`. See also section "Details" below.

- seed:

  Pseudorandom number generation (PRNG) seed by which the same results
  can be obtained again if needed. Passed to argument `seed` of
  [`set.seed()`](https://rdrr.io/r/base/Random.html), but can also be
  `NA` to not call [`set.seed()`](https://rdrr.io/r/base/Random.html) at
  all. If not `NA`, then the PRNG state is reset (to the state before
  calling `project()`) upon exiting `project()`. Here, `seed` is used
  for clustering the reference model's posterior draws (if
  `!is.null(nclusters)`) and for drawing new group-level effects when
  predicting from a multilevel submodel (however, not yet in case of a
  GAMM) and having global option `projpred.mlvl_pred_new` set to `TRUE`.
  (Such a prediction takes place when calculating output elements `dis`
  and `ce`.)

- verbose:

  A single integer value from the set \\\\0, 1, 2\\\\ (if
  `!is.null(predictor_terms)`, \\1\\ and \\2\\ have the same effect),
  indicating how much information (if any) to print out during the
  computations. Higher values indicate that more information should be
  printed, `0` deactivates the verbose mode. Internally, argument
  `verbose` is coerced to integer via
  [`as.integer()`](https://rdrr.io/r/base/integer.html), so technically,
  a single logical value or a single numeric value work as well.

- ...:

  Arguments passed to
  [`get_refmodel()`](https://mc-stan.org/projpred/dev/reference/refmodel-init-get.md)
  (if
  [`get_refmodel()`](https://mc-stan.org/projpred/dev/reference/refmodel-init-get.md)
  is actually used; see argument `object`) as well as to the divergence
  minimizer (if `refit_prj` is `TRUE`).

## Value

If the projection is performed onto a single submodel (i.e.,
`length(nterms) == 1 || !is.null(predictor_terms)`), an object of class
`projection` which is a `list` containing the following elements:

- `dis`:

  Projected draws for the dispersion parameter.

- `ce`:

  The cross-entropy part of the Kullback-Leibler (KL) divergence from
  the reference model to the submodel. For some families, this is not
  the actual cross-entropy, but a reduced one where terms which would
  cancel out when calculating the KL divergence have been dropped. In
  case of the Gaussian family, that reduced cross-entropy is further
  modified, yielding merely a proxy.

- `wdraws_prj`:

  Weights for the projected draws.

- `predictor_terms`:

  A character vector of the submodel's predictor terms.

- `outdmin`:

  A `list` containing the submodel fits (one fit per projected draw).
  This is the same as the return value of the `div_minimizer` function
  (see
  [`init_refmodel()`](https://mc-stan.org/projpred/dev/reference/refmodel-init-get.md)),
  except if `project()` was used with an `object` of class `vsel` based
  on an L1 search as well as with `refit_prj = FALSE`, in which case
  this is the output from an internal *L1-penalized* divergence
  minimizer.

- `cl_ref`:

  A numeric vector of length equal to the number of posterior draws in
  the reference model, containing the cluster indices of these draws.

- `wdraws_ref`:

  A numeric vector of length equal to the number of posterior draws in
  the reference model, giving the weights of these draws. These weights
  should be treated as not being normalized (i.e., they don't
  necessarily sum to `1`).

- `const_wdraws_prj`:

  A single logical value indicating whether the projected draws have
  constant weights (`TRUE`) or not (`FALSE`).

- `refmodel`:

  The reference model object.

If the projection is performed onto more than one submodel, the output
from above is returned for each submodel, giving a `list` with one
element for each submodel.

The elements of an object of class `projection` are not meant to be
accessed directly but instead via helper functions (see the main
vignette and
[projpred-package](https://mc-stan.org/projpred/dev/reference/projpred-package.md);
see also
[`as_draws_matrix.projection()`](https://mc-stan.org/projpred/dev/reference/as_draws_matrix.projection.md),
argument `return_draws_matrix` of
[`proj_linpred()`](https://mc-stan.org/projpred/dev/reference/pred-projection.md),
and argument `nresample_clusters` of
[`proj_predict()`](https://mc-stan.org/projpred/dev/reference/pred-projection.md)
for the intended use of the weights stored in element `wdraws_prj`).

## Details

Arguments `ndraws` and `nclusters` are automatically truncated at the
number of posterior draws in the reference model (which is `1` for
`datafit`s). Using less draws or clusters in `ndraws` or `nclusters`
than posterior draws in the reference model may result in slightly
inaccurate projection performance. Increasing these arguments affects
the computation time linearly.

If `refit_prj = FALSE` (which is only possible if `object` is of class
`vsel`), `project()` retrieves the submodel fits from the full-data
search that was run when creating `object`. Usually, the search relies
on a rather coarse clustering or thinning of the reference model's
posterior draws (by default,
[`varsel()`](https://mc-stan.org/projpred/dev/reference/varsel.md) and
[`cv_varsel()`](https://mc-stan.org/projpred/dev/reference/cv_varsel.md)
use `nclusters = 20`). Consequently, `project()` with
`refit_prj = FALSE` then inherits this coarse clustering or thinning.

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

# Projection onto the best submodel with 2 predictor terms (with a small
# value for `nclusters`, but only for the sake of speed in this example;
# this is not recommended in general):
prj_from_vs <- project(vs, nterms = 2, nclusters = 10, seed = 9182)

# Projection onto an arbitrary combination of predictor terms (with a small
# value for `nclusters`, but only for the sake of speed in this example;
# this is not recommended in general):
prj <- project(fit, predictor_terms = c("X1", "X3", "X5"), nclusters = 10,
               seed = 9182)
```
