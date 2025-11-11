# Predictions from a submodel (after projection)

After the projection of the reference model onto a submodel, the linear
predictors (for the original or a new dataset) based on that submodel
can be calculated by `proj_linpred()`. These linear predictors can also
be transformed to response scale and averaged across the projected
parameter draws. Furthermore, `proj_linpred()` returns the corresponding
log predictive density values if the (original or new) dataset contains
response values. The `proj_predict()` function draws from the predictive
distributions (there is one such distribution for each observation from
the original or new dataset) of the submodel that the reference model
has been projected onto. If the projection has not been performed yet,
both functions call
[`project()`](https://mc-stan.org/projpred/dev/reference/project.md)
internally to perform the projection. Both functions can also handle
multiple submodels at once (for `object`s of class `vsel` or `object`s
returned by a
[`project()`](https://mc-stan.org/projpred/dev/reference/project.md)
call to an object of class `vsel`; see
[`project()`](https://mc-stan.org/projpred/dev/reference/project.md)).

## Usage

``` r
proj_linpred(
  object,
  newdata = NULL,
  offsetnew = NULL,
  weightsnew = NULL,
  filter_nterms = NULL,
  transform = FALSE,
  integrated = FALSE,
  allow_nonconst_wdraws_prj = return_draws_matrix,
  return_draws_matrix = FALSE,
  .seed = NA,
  ...
)

proj_predict(
  object,
  newdata = NULL,
  offsetnew = NULL,
  weightsnew = NULL,
  filter_nterms = NULL,
  nresample_clusters = 1000,
  return_draws_matrix = FALSE,
  .seed = NA,
  resp_oscale = TRUE,
  ...
)
```

## Arguments

- object:

  An object returned by
  [`project()`](https://mc-stan.org/projpred/dev/reference/project.md)
  or an object that can be passed to argument `object` of
  [`project()`](https://mc-stan.org/projpred/dev/reference/project.md).

- newdata:

  Passed to argument `newdata` of the reference model's
  `extract_model_data` function (see
  [`init_refmodel()`](https://mc-stan.org/projpred/dev/reference/refmodel-init-get.md)).
  Provides the predictor (and possibly also the response) data for the
  new (or old) observations. May also be `NULL` for using the original
  dataset. If not `NULL`, any `NA`s will trigger an error.

- offsetnew:

  Passed to argument `orhs` of the reference model's
  `extract_model_data` function (see
  [`init_refmodel()`](https://mc-stan.org/projpred/dev/reference/refmodel-init-get.md)).
  Used to get the offsets for the new (or old) observations.

- weightsnew:

  Passed to argument `wrhs` of the reference model's
  `extract_model_data` function (see
  [`init_refmodel()`](https://mc-stan.org/projpred/dev/reference/refmodel-init-get.md)).
  Used to get the weights for the new (or old) observations.

- filter_nterms:

  Only applies if `object` is an object returned by
  [`project()`](https://mc-stan.org/projpred/dev/reference/project.md).
  In that case, `filter_nterms` can be used to filter `object` for only
  those elements (submodels) with a number of predictor terms in
  `filter_nterms`. Therefore, needs to be a numeric vector or `NULL`. If
  `NULL`, use all submodels.

- transform:

  For `proj_linpred()` only. A single logical value indicating whether
  the linear predictor should be transformed to response scale using the
  inverse-link function (`TRUE`) or not (`FALSE`). In case of the latent
  projection, argument `transform` is similar in spirit to argument
  `resp_oscale` from other functions and affects the scale of both
  output elements `pred` and `lpd` (see sections "Details" and "Value"
  below).

- integrated:

  For `proj_linpred()` only. A single logical value indicating whether
  the output should be averaged across the projected posterior draws
  (`TRUE`) or not (`FALSE`).

- allow_nonconst_wdraws_prj:

  Only relevant for `proj_linpred()` and only if `integrated` is
  `FALSE`. A single logical value indicating whether to allow projected
  draws with different (i.e., nonconstant) weights (`TRUE`) or not
  (`FALSE`). If `return_draws_matrix` is `TRUE`,
  `allow_nonconst_wdraws_prj` is internally set to `TRUE` as well.
  **CAUTION**: Expert use only because if set to `TRUE`, the weights of
  the projected draws are stored in attributes `wdraws_prj` and handling
  these attributes requires special care (e.g., when subsetting the
  returned matrices).

- return_draws_matrix:

  A single logical value indicating whether to return an object (in case
  of `proj_predict()`) or objects (in case of `proj_linpred()`) of class
  `draws_matrix` (see
  [`posterior::draws_matrix()`](https://mc-stan.org/posterior/reference/draws_matrix.html)).
  In case of `proj_linpred()` and projected draws with nonconstant
  weights (as well as `integrated` being `FALSE`),
  [`posterior::weight_draws()`](https://mc-stan.org/posterior/reference/weight_draws.html)
  is applied internally.

- .seed:

  Pseudorandom number generation (PRNG) seed by which the same results
  can be obtained again if needed. Passed to argument `seed` of
  [`set.seed()`](https://rdrr.io/r/base/Random.html), but can also be
  `NA` to not call [`set.seed()`](https://rdrr.io/r/base/Random.html) at
  all. If not `NA`, then the PRNG state is reset (to the state before
  calling `proj_linpred()` or `proj_predict()`) upon exiting
  `proj_linpred()` or `proj_predict()`. Here, `.seed` is used for
  drawing new group-level effects in case of a multilevel submodel
  (however, not yet in case of a GAMM) and for drawing from the
  predictive distributions of the submodel(s) in case of
  `proj_predict()`. If a clustered projection was performed, then in
  `proj_predict()`, `.seed` is also used for drawing from the set of
  projected clusters of posterior draws (see argument
  `nresample_clusters`). If
  [`project()`](https://mc-stan.org/projpred/dev/reference/project.md)
  is called internally with `seed = NA` (or with `seed` being a lazily
  evaluated expression that uses the PRNG), then `.seed` also affects
  the PRNG usage there.

- ...:

  Arguments passed to
  [`project()`](https://mc-stan.org/projpred/dev/reference/project.md)
  if `object` is not already an object returned by
  [`project()`](https://mc-stan.org/projpred/dev/reference/project.md).

- nresample_clusters:

  For `proj_predict()` with clustered projection (and nonconstant
  weights for the projected draws) only. Number of draws to return from
  the predictive distributions of the submodel(s). Not to be confused
  with argument `nclusters` of
  [`project()`](https://mc-stan.org/projpred/dev/reference/project.md):
  `nresample_clusters` gives the number of draws (*with* replacement)
  from the set of clustered posterior draws after projection (with this
  set being determined by argument `nclusters` of
  [`project()`](https://mc-stan.org/projpred/dev/reference/project.md)).

- resp_oscale:

  Only relevant for the latent projection. A single logical value
  indicating whether to draw from the posterior-projection predictive
  distributions on the original response scale (`TRUE`) or on latent
  scale (`FALSE`).

## Value

In the following, \\S\_{\mathrm{prj}}\\, \\N\\, \\C\_{\mathrm{cat}}\\,
and \\C\_{\mathrm{lat}}\\ from help topic
[refmodel-init-get](https://mc-stan.org/projpred/dev/reference/refmodel-init-get.md)
are used. (For `proj_linpred()` with `integrated = TRUE`, we have
\\S\_{\mathrm{prj}} = 1\\.) Furthermore, let \\C\\ denote either
\\C\_{\mathrm{cat}}\\ (if `transform = TRUE`) or \\C\_{\mathrm{lat}}\\
(if `transform = FALSE`). Then, if the prediction is done for one
submodel only (i.e., `length(nterms) == 1 || !is.null(predictor_terms)`
in the explicit or implicit call to
[`project()`](https://mc-stan.org/projpred/dev/reference/project.md),
see argument `object`):

- `proj_linpred()` returns a `list` with the following elements:

  - Element `pred` contains the actual predictions, i.e., the linear
    predictors, possibly transformed to response scale (depending on
    argument `transform`).

  - Element `lpd` is non-`NULL` only if `newdata` is `NULL` or if
    `newdata` contains response values in the corresponding column. In
    that case, it contains the log predictive density values
    (conditional on each of the projected parameter draws if
    `integrated = FALSE` and averaged across the projected parameter
    draws if `integrated = TRUE`).

  In case of (i) the traditional projection, (ii) the latent projection
  with `transform = FALSE`, or (iii) the latent projection with
  `transform = TRUE` and `<refmodel>$family$cats` (where `<refmodel>` is
  an object resulting from
  [`init_refmodel()`](https://mc-stan.org/projpred/dev/reference/refmodel-init-get.md);
  see also
  [`extend_family()`](https://mc-stan.org/projpred/dev/reference/extend_family.md)'s
  argument `latent_y_unqs`) being `NULL`, both elements are
  \\S\_{\mathrm{prj}} \times N\\ matrices (converted to a—possibly
  weighted—`draws_matrix` if argument `return_draws_matrix` is `TRUE`,
  see the description of this argument). In case of (i) the
  augmented-data projection or (ii) the latent projection with
  `transform = TRUE` and `<refmodel>$family$cats` being not `NULL`,
  `pred` is an \\S\_{\mathrm{prj}} \times N \times C\\ array (if
  argument `return_draws_matrix` is `TRUE`, this array is "compressed"
  to an \\S\_{\mathrm{prj}} \times (N \cdot C)\\ matrix—with the columns
  consisting of \\C\\ blocks of \\N\\ rows—and then converted to
  a—possibly weighted—`draws_matrix`) and `lpd` is an
  \\S\_{\mathrm{prj}} \times N\\ matrix (converted to a—possibly
  weighted—`draws_matrix` if argument `return_draws_matrix` is `TRUE`).
  If `return_draws_matrix` is `FALSE` and `allow_nonconst_wdraws_prj` is
  `TRUE` and `integrated` is `FALSE` and the projected draws have
  nonconstant weights, then both `list` elements have the weights of
  these draws stored in an attribute `wdraws_prj`. (If
  `return_draws_matrix`, `allow_nonconst_wdraws_prj`, and `integrated`
  are all `FALSE`, then projected draws with nonconstant weights cause
  an error.)

- `proj_predict()` returns an \\S\_{\mathrm{prj}} \times N\\ matrix of
  predictions where \\S\_{\mathrm{prj}}\\ denotes `nresample_clusters`
  in case of clustered projection (or, more generally, in case of
  projected draws with nonconstant weights). If argument
  `return_draws_matrix` is `TRUE`, the returned matrix is converted to a
  `draws_matrix` (see
  [`posterior::draws_matrix()`](https://mc-stan.org/posterior/reference/draws_matrix.html)).
  In case of (i) the augmented-data projection or (ii) the latent
  projection with `resp_oscale = TRUE` and `<refmodel>$family$cats`
  being not `NULL`, the returned matrix (or `draws_matrix`) has an
  attribute called `cats` (the character vector of response categories)
  and the values of the matrix (or `draws_matrix`) are the predicted
  indices of the response categories (these indices refer to the order
  of the response categories from attribute `cats`).

If the prediction is done for more than one submodel, the output from
above is returned for each submodel, giving a named `list` with one
element for each submodel (the names of this `list` being the numbers of
predictor terms of the submodels when counting the intercept, too).

## Details

Currently, `proj_predict()` ignores observation weights that are not
equal to `1`. A corresponding warning is thrown if this is the case.

In case of the latent projection and `transform = FALSE`:

- Output element `pred` contains the linear predictors without any
  modifications that may be due to the original response distribution
  (e.g., for a
  [`brms::cumulative()`](https://paulbuerkner.com/brms/reference/brmsfamily.html)
  model, the ordered thresholds are not taken into account).

- Output element `lpd` contains the *latent* log predictive density
  values, i.e., those corresponding to the latent Gaussian distribution.
  If `newdata` is not `NULL`, this requires the latent response values
  to be supplied in a column called `.<response_name>` of `newdata`
  where `<response_name>` needs to be replaced by the name of the
  original response variable (if `<response_name>` contained
  parentheses, these have been stripped off by
  [`init_refmodel()`](https://mc-stan.org/projpred/dev/reference/refmodel-init-get.md);
  see the left-hand side of `formula(<refmodel>)`). For technical
  reasons, the existence of column `<response_name>` in `newdata` is
  another requirement (even though `.<response_name>` is actually used).

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

# Projection onto an arbitrary combination of predictor terms (with a small
# value for `ndraws`, but only for the sake of speed in this example; this
# is not recommended in general):
prj <- project(fit, predictor_terms = c("X1", "X3", "X5"), ndraws = 21,
               seed = 9182)

# Predictions (at the training points) from the submodel onto which the
# reference model was projected:
prjl <- proj_linpred(prj)
prjp <- proj_predict(prj, .seed = 7364)
```
