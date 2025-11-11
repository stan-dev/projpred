# Reference model and more general information

Function `get_refmodel()` is a generic function whose methods usually
call `init_refmodel()` which is the underlying workhorse (and may also
be used directly without a call to `get_refmodel()`).

Both, `get_refmodel()` and `init_refmodel()`, create an object
containing information needed for the projection predictive variable
selection, namely about the reference model, the submodels, and how the
projection should be carried out. For the sake of simplicity, the
documentation may refer to the resulting object also as "reference
model" or "reference model object", even though it also contains
information about the submodels and the projection.

A "typical" reference model object is created by
`get_refmodel.stanreg()` and
[`brms::get_refmodel.brmsfit()`](https://paulbuerkner.com/brms/reference/get_refmodel.brmsfit.html),
either implicitly by a call to a top-level function such as
[`project()`](https://mc-stan.org/projpred/dev/reference/project.md),
[`varsel()`](https://mc-stan.org/projpred/dev/reference/varsel.md), and
[`cv_varsel()`](https://mc-stan.org/projpred/dev/reference/cv_varsel.md)
or explicitly by a call to `get_refmodel()`. All non-"typical" reference
model objects will be called "custom" reference model objects.

Some arguments are for \\K\\-fold cross-validation (\\K\\-fold CV) only;
see
[`cv_varsel()`](https://mc-stan.org/projpred/dev/reference/cv_varsel.md)
for the use of \\K\\-fold CV in projpred.

## Usage

``` r
get_refmodel(object, ...)

# S3 method for class 'refmodel'
get_refmodel(object, ...)

# S3 method for class 'vsel'
get_refmodel(object, ...)

# S3 method for class 'projection'
get_refmodel(object, ...)

# Default S3 method
get_refmodel(object, family = NULL, ...)

# S3 method for class 'stanreg'
get_refmodel(object, latent = FALSE, dis = NULL, ...)

init_refmodel(
  object,
  data,
  formula,
  family,
  ref_predfun = NULL,
  div_minimizer = NULL,
  proj_predfun = NULL,
  extract_model_data = NULL,
  cvfun = NULL,
  cvfits = NULL,
  dis = NULL,
  cvrefbuilder = NULL,
  called_from_cvrefbuilder = FALSE,
  ...
)
```

## Arguments

- object:

  For `init_refmodel()`, an object that the functions from arguments
  `extract_model_data` and `ref_predfun` can be applied to, with a
  `NULL` object being treated specially (see section "Value" below). For
  `get_refmodel.default()`, an object that function
  [`family()`](https://rdrr.io/r/stats/family.html) can be applied to in
  order to retrieve the family (if argument `family` is `NULL`),
  additionally to the properties required for `init_refmodel()`. For
  non-default methods of `get_refmodel()`, an object of the
  corresponding class.

- ...:

  For `get_refmodel.default()` and `get_refmodel.stanreg()`: arguments
  passed to `init_refmodel()`. For the `get_refmodel()` generic:
  arguments passed to the appropriate method. For `init_refmodel()`:
  arguments passed to
  [`extend_family()`](https://mc-stan.org/projpred/dev/reference/extend_family.md)
  (apart from `family`).

- family:

  An object of class `family` representing the observation model (i.e.,
  the distributional family for the response) of the *submodels*.
  (However, the link and the inverse-link function of this `family` are
  also used for quantities like predictions and fitted values related to
  the *reference model*.) May be `NULL` for `get_refmodel.default()` in
  which case the family is retrieved from `object`. For custom reference
  models, `family` does not have to coincide with the family of the
  reference model (if the reference model possesses a formal `family` at
  all). In typical reference models, however, these families do
  coincide. Furthermore, the latent projection is an exception where
  `family` is not the family of the submodels (in that case, the family
  of the submodels is the
  [`gaussian()`](https://rdrr.io/r/stats/family.html) family).

- latent:

  A single logical value indicating whether to use the latent projection
  (`TRUE`) or not (`FALSE`). Note that setting `latent = TRUE` causes
  all arguments starting with `augdat_` to be ignored.

- dis:

  A vector of posterior draws for the reference model's dispersion
  parameter or—more precisely—the posterior values for the reference
  model's parameter-conditional predictive variance (assuming that this
  variance is the same for all observations). May be `NULL` if the
  submodels have no dispersion parameter or if the submodels do have a
  dispersion parameter, but `object` is `NULL` (in which case `0` is
  used for `dis`). Note that for the
  [`gaussian()`](https://rdrr.io/r/stats/family.html) `family`, `dis` is
  the standard deviation, not the variance.

- data:

  A `data.frame` containing the data to use for the projection
  predictive variable selection. Any `contrasts` attributes of the
  dataset's columns are silently removed. For custom reference models,
  the columns of `data` do not necessarily have to coincide with those
  of the dataset used for fitting the reference model, but keep in mind
  that a row-subset of `data` is used for argument `newdata` of
  `ref_predfun` during \\K\\-fold CV.

- formula:

  The full formula to use for the search procedure. For custom reference
  models, this does not necessarily coincide with the reference model's
  formula. For general information about formulas in R, see
  [`formula`](https://rdrr.io/r/stats/formula.html). For information
  about possible right-hand side (i.e., predictor) terms in `formula`
  here, see the main vignette and section "Formula terms" below. For
  multilevel formulas, see also package lme4 (in particular, functions
  [`lme4::lmer()`](https://rdrr.io/pkg/lme4/man/lmer.html) and
  [`lme4::glmer()`](https://rdrr.io/pkg/lme4/man/glmer.html)). For
  additive formulas, see also packages mgcv (in particular, function
  [`mgcv::gam()`](https://rdrr.io/pkg/mgcv/man/gam.html)) and gamm4 (in
  particular, function
  [`gamm4::gamm4()`](https://rdrr.io/pkg/gamm4/man/gamm4.html)).

- ref_predfun:

  Prediction function for the linear predictor of the reference model,
  including offsets (if existing). See also section "Arguments
  `ref_predfun`, `proj_predfun`, and `div_minimizer`" below. If `object`
  is `NULL`, `ref_predfun` is ignored and an internal default is used
  instead.

- div_minimizer:

  A function for minimizing the Kullback-Leibler (KL) divergence from
  the reference model to a submodel (i.e., for performing the projection
  of the reference model onto a submodel). The output of `div_minimizer`
  is used, e.g., by `proj_predfun`'s argument `fits`. See also section
  "Arguments `ref_predfun`, `proj_predfun`, and `div_minimizer`" below.

- proj_predfun:

  Prediction function for the linear predictor of a submodel onto which
  the reference model is projected. See also section "Arguments
  `ref_predfun`, `proj_predfun`, and `div_minimizer`" below.

- extract_model_data:

  A function for fetching some variables (response, observation weights,
  offsets) from the original dataset (supplied to argument `data`) or
  from a new dataset. May be `NULL` for using an internal default that
  essentially corresponds to
  [`y_wobs_offs()`](https://mc-stan.org/projpred/dev/reference/y_wobs_offs.md).
  See also section "Argument `extract_model_data`" below.

- cvfun:

  For \\K\\-fold CV only. A function that, given a fold indices vector,
  fits the reference model separately for each fold and returns the
  \\K\\ model fits as a `list`. If `object` is `NULL`, `cvfun` may be
  `NULL` for using an internal default. Only one of `cvfits` and `cvfun`
  needs to be provided (for \\K\\-fold CV). Note that `cvfits` takes
  precedence over `cvfun`, i.e., if both are provided, `cvfits` is used.

- cvfits:

  For \\K\\-fold CV only. A `list` containing the \\K\\ reference model
  refits from which reference model objects are created. This `list`
  needs to have an attribute called `folds`, consisting of an integer
  vector giving the fold indices (one fold index per observation). Only
  one of `cvfits` and `cvfun` needs to be provided (for \\K\\-fold CV).
  Note that `cvfits` takes precedence over `cvfun`, i.e., if both are
  provided, `cvfits` is used.

- cvrefbuilder:

  For \\K\\-fold CV only. A function that, given a reference model fit
  for fold \\k \in \\1, ..., K\\\\, returns an object of the same type
  as `init_refmodel()` does. The reference model fit for fold \\k\\ is
  the \\k\\-th element of the return value of `cvfun` or the \\k\\-th
  element of the `list` supplied to `cvfits` (either here in
  `init_refmodel()` or in
  [`cv_varsel.refmodel()`](https://mc-stan.org/projpred/dev/reference/cv_varsel.md)),
  extended by elements `omitted` (containing the indices of the left-out
  observations in that fold) and `projpred_k` (containing the integer
  \\k\\) if that \\k\\-th element is a `list` itself (otherwise,
  `omitted` and `projpred_k` are appended as attributes). Argument
  `cvrefbuilder` may be `NULL` for using an internal default:
  `get_refmodel()` if `object` is not `NULL` and a function calling
  `init_refmodel()` appropriately (with the assumption `dis = 0`) if
  `object` is `NULL`.

- called_from_cvrefbuilder:

  A single logical value indicating whether `init_refmodel()` is called
  from a `cvrefbuilder` function (`TRUE`) or not (`FALSE`). Currently,
  `TRUE` only causes some warnings to be suppressed (warnings which
  don't need to be thrown for each of the \\K\\ reference model objects
  because it is sufficient to throw them for the original reference
  model object only). This argument is mainly for internal use, but may
  also be helpful for users with a custom `cvrefbuilder` function.

## Value

An object that can be passed to all the functions that take the
reference model fit as the first argument, such as
[`varsel()`](https://mc-stan.org/projpred/dev/reference/varsel.md),
[`cv_varsel()`](https://mc-stan.org/projpred/dev/reference/cv_varsel.md),
[`project()`](https://mc-stan.org/projpred/dev/reference/project.md),
[`proj_linpred()`](https://mc-stan.org/projpred/dev/reference/pred-projection.md),
and
[`proj_predict()`](https://mc-stan.org/projpred/dev/reference/pred-projection.md).
Usually, the returned object is of class `refmodel`. However, if
`object` is `NULL`, the returned object is of class `datafit` as well as
of class `refmodel` (with `datafit` being first). Objects of class
`datafit` are handled differently at several places throughout this
package.

The elements of the returned object are not meant to be accessed
directly but instead via downstream functions (see the functions
mentioned above as well as
[`predict.refmodel()`](https://mc-stan.org/projpred/dev/reference/predict.refmodel.md)).

## Formula terms

Although bad practice (in general), a reference model lacking an
intercept can be used within projpred. However, it will always be
projected onto submodels which *include* an intercept. The reason is
that even if the true intercept in the reference model is zero, this
does not need to hold for the submodels.

In multilevel (group-level) terms, function calls on the right-hand side
of the `|` character (e.g., `(1 | gr(group_variable))`, which is
possible in brms) are currently not allowed in projpred.

For additive models (still an experimental feature), only
[`mgcv::s()`](https://rdrr.io/pkg/mgcv/man/s.html) and
[`mgcv::t2()`](https://rdrr.io/pkg/mgcv/man/t2.html) are currently
supported as smooth terms. Furthermore, these need to be called without
any arguments apart from the predictor names (symbols). For example, for
smoothing the effect of a predictor `x`, only `s(x)` or `t2(x)` are
allowed. As another example, for smoothing the joint effect of two
predictors `x` and `z`, only `s(x, z)` or `t2(x, z)` are allowed (and
analogously for higher-order joint effects, e.g., of three predictors).
Note that all smooth terms need to be included in `formula` (there is no
`random` argument as in
[`rstanarm::stan_gamm4()`](https://mc-stan.org/rstanarm/reference/stan_gamm4.html),
for example).

## Arguments `ref_predfun`, `proj_predfun`, and `div_minimizer`

Arguments `ref_predfun`, `proj_predfun`, and `div_minimizer` may be
`NULL` for using an internal default (see
[projpred-package](https://mc-stan.org/projpred/dev/reference/projpred-package.md)
for the functions used by the default divergence minimizers). Otherwise,
let \\N\\ denote the number of observations (in case of CV, these may be
reduced to each fold), \\S\_{\mathrm{ref}}\\ the number of posterior
draws for the reference model's parameters, and \\S\_{\mathrm{prj}}\\
the number of draws for the parameters of a submodel that the reference
model has been projected onto (short: the number of projected draws).
For the augmented-data projection, let \\C\_{\mathrm{cat}}\\ denote the
number of response categories, \\C\_{\mathrm{lat}}\\ the number of
latent response categories (which typically equals \\C\_{\mathrm{cat}} -
1\\), and define \\N\_{\mathrm{augcat}} := N \cdot C\_{\mathrm{cat}}\\
as well as \\N\_{\mathrm{auglat}} := N \cdot C\_{\mathrm{lat}}\\. Then
the functions supplied to these arguments need to have the following
prototypes:

- `ref_predfun`: `ref_predfun(fit, newdata = NULL)` where:

  - `fit` accepts the reference model fit as given in argument `object`
    (but possibly refitted to a subset of the observations, as done in
    \\K\\-fold CV).

  - `newdata` accepts either `NULL` (for using the original dataset,
    typically stored in `fit`) or data for new observations (at least in
    the form of a `data.frame`).

- `proj_predfun`: `proj_predfun(fits, newdata)` where:

  - `fits` accepts a `list` of length \\S\_{\mathrm{prj}}\\ containing
    this number of submodel fits. This `list` is the same as that
    returned by
    [`project()`](https://mc-stan.org/projpred/dev/reference/project.md)
    in its output element `outdmin` (which in turn is the same as the
    return value of `div_minimizer`, except if
    [`project()`](https://mc-stan.org/projpred/dev/reference/project.md)
    was used with an `object` of class `vsel` based on an L1 search as
    well as with `refit_prj = FALSE`).

  - `newdata` accepts data for new observations (at least in the form of
    a `data.frame`).

- `div_minimizer` does not need to have a specific prototype, but it
  needs to be able to be called with the following arguments:

  - `formula` accepts either a standard
    [`formula`](https://rdrr.io/r/stats/formula.html) with a single
    response (if \\S\_{\mathrm{prj}} = 1\\ or in case of the
    augmented-data projection) or a
    [`formula`](https://rdrr.io/r/stats/formula.html) with
    \\S\_{\mathrm{prj}} \> 1\\ response variables
    [`cbind()`](https://rdrr.io/r/base/cbind.html)-ed on the left-hand
    side in which case the projection has to be performed for each of
    the response variables separately.

  - `data` accepts a `data.frame` to be used for the projection. In case
    of the traditional or the latent projection, this dataset has \\N\\
    rows. In case of the augmented-data projection, this dataset has
    \\N\_{\mathrm{augcat}}\\ rows.

  - `family` accepts an object of class `family`.

  - `weights` accepts either observation weights (at least in the form
    of a numeric vector) or `NULL` (for using a vector of ones as
    weights).

  - `projpred_var` accepts an \\N \times S\_{\mathrm{prj}}\\ matrix of
    predictive variances (necessary for projpred's internal GLM fitter)
    in case of the traditional or the latent projection and an
    \\N\_{\mathrm{augcat}} \times S\_{\mathrm{prj}}\\ matrix (containing
    only `NA`s) in case of the augmented-data projection.

  - `projpred_ws_aug` accepts an \\N \times S\_{\mathrm{prj}}\\ matrix
    of expected values for the response in case of the traditional or
    the latent projection and an \\N\_{\mathrm{augcat}} \times
    S\_{\mathrm{prj}}\\ matrix of probabilities for the response
    categories in case of the augmented-data projection.

  - `...` accepts further arguments specified by the user (or by
    projpred).

The return value of these functions needs to be:

- `ref_predfun`: for the traditional or the latent projection, an \\N
  \times S\_{\mathrm{ref}}\\ matrix; for the augmented-data projection,
  an \\S\_{\mathrm{ref}} \times N \times C\_{\mathrm{lat}}\\ array (the
  only exception is the augmented-data projection for the
  [`binomial()`](https://rdrr.io/r/stats/family.html) family in which
  case `ref_predfun` needs to return an \\N \times S\_{\mathrm{ref}}\\
  matrix just like for the traditional projection because the array is
  constructed by an internal wrapper function).

- `proj_predfun`: for the traditional or the latent projection, an \\N
  \times S\_{\mathrm{prj}}\\ matrix; for the augmented-data projection,
  an \\N \times C\_{\mathrm{lat}} \times S\_{\mathrm{prj}}\\ array.

- `div_minimizer`: a `list` of length \\S\_{\mathrm{prj}}\\ containing
  this number of submodel fits.

## Argument `extract_model_data`

The function supplied to argument `extract_model_data` needs to have the
prototype

    extract_model_data(object, newdata, wrhs = NULL, orhs = NULL,
                       extract_y = TRUE)

where:

- `object` accepts the reference model fit as given in argument `object`
  (but possibly refitted to a subset of the observations, as done in
  \\K\\-fold CV).

- `newdata` accepts data for new observations (at least in the form of a
  `data.frame`).

- `wrhs` accepts at least (i) a right-hand side formula consisting only
  of the variable in `newdata` containing the observation weights
  or (ii) `NULL` for using the observation weights corresponding to
  `newdata` (typically, the observation weights are stored in a column
  of `newdata`; if the model was fitted without observation weights, a
  vector of ones should be used).

- `orhs` accepts at least (i) a right-hand side formula consisting only
  of the variable in `newdata` containing the offsets or (ii) `NULL` for
  using the offsets corresponding to `newdata` (typically, the offsets
  are stored in a column of `newdata`; if the model was fitted without
  offsets, a vector of zeros should be used).

- `extract_y` accepts a single logical value indicating whether output
  element `y` (see below) shall be `NULL` (`TRUE`) or not (`FALSE`).

The return value of `extract_model_data` needs to be a `list` with
elements `y`, `weights`, and `offset`, each being a numeric vector
containing the data for the response, the observation weights, and the
offsets, respectively. An exception is that `y` may also be `NULL`
(depending on argument `extract_y`), a non-numeric vector, or a
`factor`.

The weights and offsets returned by `extract_model_data` will be assumed
to hold for the reference model as well as for the submodels.

Above, arguments `wrhs` and `orhs` were assumed to have defaults of
`NULL`. It should be possible to use defaults other than `NULL`, but we
strongly recommend to use `NULL`. If defaults other than `NULL` are
used, they need to imply the behaviors described at items "(ii)" (see
the descriptions of `wrhs` and `orhs`).

## Augmented-data projection

If a custom reference model for an augmented-data projection is needed,
see also
[`extend_family()`](https://mc-stan.org/projpred/dev/reference/extend_family.md).

For the augmented-data projection, the response vector resulting from
`extract_model_data` is internally coerced to a `factor` (using
[`as.factor()`](https://rdrr.io/r/base/factor.html)). The levels of this
`factor` have to be identical to `family$cats` (*after* applying
[`extend_family()`](https://mc-stan.org/projpred/dev/reference/extend_family.md)
internally; see
[`extend_family()`](https://mc-stan.org/projpred/dev/reference/extend_family.md)'s
argument `augdat_y_unqs`).

Note that response-specific offsets (i.e., one length-\\N\\ offset
vector per response category) are not supported by projpred yet. So far,
only offsets which are the same across all response categories are
supported. This is why in case of the
[`brms::categorical()`](https://paulbuerkner.com/brms/reference/brmsfamily.html)
family, offsets are currently not supported at all.

Currently, `object = NULL` (i.e., a `datafit`; see section "Value") is
not supported in case of the augmented-data projection.

## Latent projection

If a custom reference model for a latent projection is needed, see also
[`extend_family()`](https://mc-stan.org/projpred/dev/reference/extend_family.md).

For the latent projection, `family$cats` (*after* applying
[`extend_family()`](https://mc-stan.org/projpred/dev/reference/extend_family.md)
internally; see
[`extend_family()`](https://mc-stan.org/projpred/dev/reference/extend_family.md)'s
argument `latent_y_unqs`) currently must not be `NULL` if the original
(i.e., non-latent) response is a `factor`. Conversely, if `family$cats`
(*after* applying
[`extend_family()`](https://mc-stan.org/projpred/dev/reference/extend_family.md))
is non-`NULL`, the response vector resulting from `extract_model_data`
is internally coerced to a `factor` (using
[`as.factor()`](https://rdrr.io/r/base/factor.html)). The levels of this
`factor` have to be identical to that non-`NULL` element `family$cats`.

Currently, `object = NULL` (i.e., a `datafit`; see section "Value") is
not supported in case of the latent projection.

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

# Define the reference model object explicitly:
ref <- get_refmodel(fit)
print(class(ref)) # gives `"refmodel"`
#> [1] "refmodel"
# Now see, for example, `?varsel`, `?cv_varsel`, and `?project` for
# possible post-processing functions. Most of the post-processing functions
# call get_refmodel() internally at the beginning, so you will rarely need
# to call get_refmodel() yourself.

# A custom reference model object which may be used in a variable selection
# where the candidate predictors are not a subset of those used for the
# reference model's predictions:
ref_cust <- init_refmodel(
  fit,
  data = dat_gauss,
  formula = y ~ X6 + X7,
  family = gaussian(),
  cvfun = function(folds) {
    kfold(
      fit, K = max(folds), save_fits = TRUE, folds = folds, cores = 1
    )$fits[, "fit"]
  },
  dis = as.matrix(fit)[, "sigma"],
  cvrefbuilder = function(cvfit) {
    init_refmodel(cvfit,
                  data = dat_gauss[-cvfit$omitted, , drop = FALSE],
                  formula = y ~ X6 + X7,
                  family = gaussian(),
                  dis = as.matrix(cvfit)[, "sigma"],
                  called_from_cvrefbuilder = TRUE)
  }
)
# Now, the post-processing functions mentioned above (for example,
# varsel(), cv_varsel(), and project()) may be applied to `ref_cust`.
```
