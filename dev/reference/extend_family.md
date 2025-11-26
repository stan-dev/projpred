# Extend a family

This function adds some internally required elements to an object of
class `family` (see, e.g.,
[`family()`](https://rdrr.io/r/stats/family.html)). It is called
internally by
[`init_refmodel()`](https://mc-stan.org/projpred/dev/reference/refmodel-init-get.md),
so you will rarely need to call it yourself.

## Usage

``` r
extend_family(
  family,
  latent = FALSE,
  latent_y_unqs = NULL,
  latent_ilink = NULL,
  latent_ll_oscale = NULL,
  latent_ppd_oscale = NULL,
  augdat_y_unqs = NULL,
  augdat_link = NULL,
  augdat_ilink = NULL,
  augdat_args_link = list(),
  augdat_args_ilink = list(),
  ...
)
```

## Arguments

- family:

  An object of class `family`.

- latent:

  A single logical value indicating whether to use the latent projection
  (`TRUE`) or not (`FALSE`). Note that setting `latent = TRUE` causes
  all arguments starting with `augdat_` to be ignored.

- latent_y_unqs:

  Only relevant for a latent projection where the original response
  space has finite support (i.e., the original response values may be
  regarded as categories), in which case this needs to be the character
  vector of unique response values (which will be assigned to
  `family$cats` internally) or may be left at `NULL` (so that projpred
  will try to infer it from `family$cats`). See also section "Latent
  projection" below.

- latent_ilink:

  Only relevant for the latent projection, in which case this needs to
  be the inverse-link function. If the original response family was the
  [`binomial()`](https://rdrr.io/r/stats/family.html) or the
  [`poisson()`](https://rdrr.io/r/stats/family.html) family, then
  `latent_ilink` can be `NULL`, in which case an internal default will
  be used. Can also be `NULL` in all other cases, but then an internal
  default based on `family$linkinv` will be used which might not work
  for all families. See also section "Latent projection" below.

- latent_ll_oscale:

  Only relevant for the latent projection, in which case this needs to
  be the function computing response-scale (not latent-scale)
  log-likelihood values. If `!is.null(family$cats)` (after taking
  `latent_y_unqs` into account) or if the original response family was
  the [`binomial()`](https://rdrr.io/r/stats/family.html) or the
  [`poisson()`](https://rdrr.io/r/stats/family.html) family, then
  `latent_ll_oscale` can be `NULL`, in which case an internal default
  will be used. Can also be `NULL` in all other cases, but then
  downstream functions will have limited functionality (a message thrown
  by `extend_family()` will state what exactly won't be available). See
  also section "Latent projection" below.

- latent_ppd_oscale:

  Only relevant for the latent projection, in which case this needs to
  be the function sampling response values given latent predictors that
  have been transformed to response scale using `latent_ilink`. If
  `!is.null(family$cats)` (after taking `latent_y_unqs` into account) or
  if the original response family was the
  [`binomial()`](https://rdrr.io/r/stats/family.html) or the
  [`poisson()`](https://rdrr.io/r/stats/family.html) family, then
  `latent_ppd_oscale` can be `NULL`, in which case an internal default
  will be used. Can also be `NULL` in all other cases, but then
  downstream functions will have limited functionality (a message thrown
  by `extend_family()` will state what exactly won't be available). See
  also section "Latent projection" below. Note that although this
  function has the abbreviation "PPD" in its name (which stands for
  "posterior predictive distribution"), projpred currently only uses it
  in
  [`proj_predict()`](https://mc-stan.org/projpred/dev/reference/pred-projection.md),
  i.e., for sampling from what would better be termed
  posterior-projection predictive distribution (PPPD).

- augdat_y_unqs:

  Only relevant for augmented-data projection, in which case this needs
  to be the character vector of unique response values (which will be
  assigned to `family$cats` internally) or may be left at `NULL` if
  `family$cats` is already non-`NULL`. See also section "Augmented-data
  projection" below.

- augdat_link:

  Only relevant for augmented-data projection, in which case this needs
  to be the link function. Use `NULL` for the traditional projection.
  See also section "Augmented-data projection" below.

- augdat_ilink:

  Only relevant for augmented-data projection, in which case this needs
  to be the inverse-link function. Use `NULL` for the traditional
  projection. See also section "Augmented-data projection" below.

- augdat_args_link:

  Only relevant for augmented-data projection, in which case this may be
  a named `list` of arguments to pass to the function supplied to
  `augdat_link`.

- augdat_args_ilink:

  Only relevant for augmented-data projection, in which case this may be
  a named `list` of arguments to pass to the function supplied to
  `augdat_ilink`.

- ...:

  Ignored (exists only to swallow up further arguments which might be
  passed to this function).

## Value

The `family` object extended in the way needed by projpred.

## Details

In the following, \\N\\, \\C\_{\mathrm{cat}}\\, \\C\_{\mathrm{lat}}\\,
\\S\_{\mathrm{ref}}\\, and \\S\_{\mathrm{prj}}\\ from help topic
[refmodel-init-get](https://mc-stan.org/projpred/dev/reference/refmodel-init-get.md)
are used. Note that \\N\\ does not necessarily denote the number of
original observations; it can also refer to new observations.
Furthermore, let \\S\\ denote either \\S\_{\mathrm{ref}}\\ or
\\S\_{\mathrm{prj}}\\, whichever is appropriate in the context where it
is used.

## Augmented-data projection

As their first input, the functions supplied to arguments `augdat_link`
and `augdat_ilink` have to accept:

- For `augdat_link`: an \\S \times N \times C\_{\mathrm{cat}}\\ array
  containing the probabilities for the response categories. The order of
  the response categories is the same as in `family$cats` (see argument
  `augdat_y_unqs`).

- For `augdat_ilink`: an \\S \times N \times C\_{\mathrm{lat}}\\ array
  containing the linear predictors.

The return value of these functions needs to be:

- For `augdat_link`: an \\S \times N \times C\_{\mathrm{lat}}\\ array
  containing the linear predictors.

- For `augdat_ilink`: an \\S \times N \times C\_{\mathrm{cat}}\\ array
  containing the probabilities for the response categories. The order of
  the response categories has to be the same as in `family$cats` (see
  argument `augdat_y_unqs`).

For the augmented-data projection, the response vector resulting from
`extract_model_data` (see
[`init_refmodel()`](https://mc-stan.org/projpred/dev/reference/refmodel-init-get.md))
is coerced to a `factor` (using
[`as.factor()`](https://rdrr.io/r/base/factor.html)) at multiple places
throughout this package. Inside of
[`init_refmodel()`](https://mc-stan.org/projpred/dev/reference/refmodel-init-get.md),
the levels of this `factor` have to be identical to `family$cats`
(*after* applying `extend_family()` inside of
[`init_refmodel()`](https://mc-stan.org/projpred/dev/reference/refmodel-init-get.md)).
Everywhere else, these levels have to be a subset of
`<refmodel>$family$cats` (where `<refmodel>` is an object resulting from
[`init_refmodel()`](https://mc-stan.org/projpred/dev/reference/refmodel-init-get.md)).
See argument `augdat_y_unqs` for how to control `family$cats`.

For ordinal brms families, be aware that the submodels (onto which the
reference model is projected) currently have the following restrictions:

- The discrimination parameter `disc` is not supported (i.e., it is a
  constant with value 1).

- The thresholds are `"flexible"` (see
  [`brms::brmsfamily()`](https://paulbuerkner.com/brms/reference/brmsfamily.html)).

- The thresholds do not vary across the levels of a `factor`-like
  variable (see argument `gr` of
  [`brms::resp_thres()`](https://paulbuerkner.com/brms/reference/addition-terms.html)).

- The `"probit_approx"` link is replaced by `"probit"`.

For the
[`brms::categorical()`](https://paulbuerkner.com/brms/reference/brmsfamily.html)
family, be aware that:

- For multilevel submodels, the group-level effects are allowed to be
  correlated between different response categories.

- For multilevel submodels, mclogit versions \< 0.9.4 may throw the
  error `'a' (<number> x 1) must be square`. Updating mclogit to a
  version \>= 0.9.4 should fix this.

## Latent projection

The function supplied to argument `latent_ilink` needs to have the
prototype

    latent_ilink(lpreds, cl_ref, wdraws_ref = rep(1, length(cl_ref)))

where:

- `lpreds` accepts an \\S \times N\\ matrix containing the linear
  predictors.

- `cl_ref` accepts a numeric vector of length \\S\_{\mathrm{ref}}\\,
  containing projpred's internal cluster indices for these draws.

- `wdraws_ref` accepts a numeric vector of length \\S\_{\mathrm{ref}}\\,
  containing weights for these draws. These weights should be treated as
  not being normalized (i.e., they don't necessarily sum to `1`).

The return value of `latent_ilink` needs to contain the linear
predictors transformed to the original response space, with the
following structure:

- If `is.null(family$cats)` (after taking `latent_y_unqs` into account):
  an \\S \times N\\ matrix.

- If `!is.null(family$cats)` (after taking `latent_y_unqs` into
  account): an \\S \times N \times C\_{\mathrm{cat}}\\ array. In that
  case, `latent_ilink` needs to return *probabilities* (for the response
  categories given in `family$cats`, after taking `latent_y_unqs` into
  account).

The function supplied to argument `latent_ll_oscale` needs to have the
prototype

    latent_ll_oscale(ilpreds, dis, y_oscale, wobs = rep(1, ncol(ilpreds)),
                     cens, cl_ref, wdraws_ref = rep(1, length(cl_ref)))

where:

- `ilpreds` accepts the return value from `latent_ilink`.

- `dis` accepts a vector of length \\S\\ containing dispersion parameter
  draws.

- `y_oscale` accepts a vector of length \\N\\ containing response values
  on the original response scale.

- `wobs` accepts a numeric vector of length \\N\\ containing observation
  weights.

- `cens` accepts a vector containing censoring indicators for the
  observations for which to calculate the response-scale log-likelihood
  values (i.e., for the observations from the second dimension of
  `ilpreds`). When calling `latent_ll_oscale`, projpred always specifies
  argument `cens` (with value `NULL` if attribute `cens_var` of
  `latent_ll_oscale` does not exist or is `NULL`), so a default value of
  `cens` can be defined, but will not be used.

- `cl_ref` accepts the same input as argument `cl_ref` of
  `latent_ilink`.

- `wdraws_ref` accepts the same input as argument `wdraws_ref` of
  `latent_ilink`.

In case of censoring (in the response values, i.e., survival or
time-to-event analysis), the latent projection (with response-scale
analyses) can be used by setting an attribute `cens_var` of the
`latent_ll_oscale` function to a right-hand side formula with the name
of the variable containing the censoring indicators (e.g., `0` =
uncensored, `1` = censored) on its right-hand side. This variable named
in the `cens_var` attribute is then retrieved (internally, whenever
calling the `latent_ll_oscale` function) from the original dataset
(possibly subsetted to the observations corresponding to the second
dimension of `ilpreds`), `newdata`, or element `data` from
[`varsel()`](https://mc-stan.org/projpred/dev/reference/varsel.md)'s
argument `d_test`, whichever is applicable. The content of the retrieved
variable is passed to argument `cens` of the `latent_ll_oscale`
function. Note that only the performance statistics `"elpd"`, `"mlpd"`,
and `"gmpd"` take censoring into account (on response scale).

The return value of `latent_ll_oscale` needs to be an \\S \times N\\
matrix containing the response-scale (not latent-scale) log-likelihood
values for the \\N\\ observations from its inputs.

The function supplied to argument `latent_ppd_oscale` needs to have the
prototype

    latent_ppd_oscale(ilpreds_resamp, dis_resamp,
                      wobs = rep(1, ncol(ilpreds_resamp)), cl_ref,
                      wdraws_ref = rep(1, length(cl_ref)), idxs_prjdraws)

where:

- `ilpreds_resamp` accepts the return value from `latent_ilink`, but
  possibly with resampled (clustered) draws (see argument
  `nresample_clusters` of
  [`proj_predict()`](https://mc-stan.org/projpred/dev/reference/pred-projection.md)).

- `dis_resamp` accepts a vector of length `dim(ilpreds_resamp)[1]`
  containing dispersion parameter draws, possibly resampled (in the same
  way as the draws in `ilpreds_resamp`, see also argument
  `idxs_prjdraws`).

- `wobs` accepts a numeric vector of length \\N\\ containing observation
  weights.

- `cl_ref` accepts the same input as argument `cl_ref` of
  `latent_ilink`.

- `wdraws_ref` accepts the same input as argument `wdraws_ref` of
  `latent_ilink`.

- `idxs_prjdraws` accepts a numeric vector of length
  `dim(ilpreds_resamp)[1]` containing the resampled indices of the
  projected draws (i.e., these indices are values from the set \\\\1,
  ..., \texttt{dim(ilpreds)\[1\]}\\\\ where `ilpreds` denotes the return
  value of `latent_ilink`).

The return value of `latent_ppd_oscale` needs to be a
\\\texttt{dim(ilpreds\\resamp)\[1\]} \times N\\ matrix containing the
response-scale (not latent-scale) draws from the posterior(-projection)
predictive distributions for the \\N\\ observations from its inputs.

If the bodies of these three functions involve parameter draws from the
reference model which have not been projected (e.g., for `latent_ilink`,
the thresholds in an ordinal model),
[`cl_agg()`](https://mc-stan.org/projpred/dev/reference/cl_agg.md) is
provided as a helper function for aggregating these reference model
draws in the same way as the draws have been aggregated for the first
argument of these functions (e.g., `lpreds` in case of `latent_ilink`).

In fact, the weights passed to argument `wdraws_ref` are nonconstant
only in case of
[`cv_varsel()`](https://mc-stan.org/projpred/dev/reference/cv_varsel.md)
with `cv_method = "LOO"` and `validate_search = TRUE`. In that case, the
weights passed to this argument are the PSIS-LOO CV weights for one
observation. Note that although argument `wdraws_ref` has the suffix
`_ref`, `wdraws_ref` does not necessarily obtain weights for the
*initial* reference model's posterior draws: In case of
[`cv_varsel()`](https://mc-stan.org/projpred/dev/reference/cv_varsel.md)
with `cv_method = "kfold"`, these weights may refer to one of the \\K\\
reference model refits (but in that case, they are constant anyway).

If `family$cats` is not `NULL` (after taking `latent_y_unqs` into
account), then the response vector resulting from `extract_model_data`
(see
[`init_refmodel()`](https://mc-stan.org/projpred/dev/reference/refmodel-init-get.md))
is coerced to a `factor` (using
[`as.factor()`](https://rdrr.io/r/base/factor.html)) at multiple places
throughout this package. Inside of
[`init_refmodel()`](https://mc-stan.org/projpred/dev/reference/refmodel-init-get.md),
the levels of this `factor` have to be identical to `family$cats`
(*after* applying `extend_family()` inside of
[`init_refmodel()`](https://mc-stan.org/projpred/dev/reference/refmodel-init-get.md)).
Everywhere else, these levels have to be a subset of
`<refmodel>$family$cats` (where `<refmodel>` is an object resulting from
[`init_refmodel()`](https://mc-stan.org/projpred/dev/reference/refmodel-init-get.md)).
