# Predictions or log posterior predictive densities from a reference model

This is the [`predict()`](https://rdrr.io/r/stats/predict.html) method
for `refmodel` objects (returned by
[`get_refmodel()`](https://mc-stan.org/projpred/dev/reference/refmodel-init-get.md)
or
[`init_refmodel()`](https://mc-stan.org/projpred/dev/reference/refmodel-init-get.md)).
It offers three types of output which are all based on the reference
model and new (or old) observations: Either the linear predictor on link
scale, the linear predictor transformed to response scale, or the log
posterior predictive density.

## Usage

``` r
# S3 method for class 'refmodel'
predict(
  object,
  newdata = NULL,
  ynew = NULL,
  offsetnew = NULL,
  weightsnew = NULL,
  type = "response",
  ...
)
```

## Arguments

- object:

  An object of class `refmodel` (returned by
  [`get_refmodel()`](https://mc-stan.org/projpred/dev/reference/refmodel-init-get.md)
  or
  [`init_refmodel()`](https://mc-stan.org/projpred/dev/reference/refmodel-init-get.md)).

- newdata:

  Passed to argument `newdata` of the reference model's
  `extract_model_data` function (see
  [`init_refmodel()`](https://mc-stan.org/projpred/dev/reference/refmodel-init-get.md)).
  Provides the predictor (and possibly also the response) data for the
  new (or old) observations. May also be `NULL` for using the original
  dataset. If not `NULL`, any `NA`s will trigger an error.

- ynew:

  If not `NULL`, then this needs to be a vector of new (or old) response
  values. See also section "Value" below. In case of (i) the
  augmented-data projection or (ii) the latent projection with
  `type = "response"` and `object$family$cats` being not `NULL`, `ynew`
  is internally coerced to a `factor` (using
  [`as.factor()`](https://rdrr.io/r/base/factor.html)). The levels of
  this `factor` have to be a subset of `object$family$cats` (see
  [`extend_family()`](https://mc-stan.org/projpred/dev/reference/extend_family.md)'s
  arguments `augdat_y_unqs` and `latent_y_unqs`, respectively).

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

- type:

  Usually only relevant if `is.null(ynew)`, but for the latent
  projection, this also affects the `!is.null(ynew)` case (see below).
  The scale on which the predictions are returned, either `"link"` or
  `"response"` (see
  [`predict.glm()`](https://rdrr.io/r/stats/predict.glm.html) but note
  that `predict.refmodel()` does not adhere to the typical R convention
  of a default prediction on link scale). For both scales, the
  predictions are averaged across the posterior draws. In case of the
  latent projection, argument `type` is similar in spirit to argument
  `resp_oscale` from other functions: If (i) `is.null(ynew)`, then
  argument `type` affects the predictions as described above. In that
  case, note that `type = "link"` yields the linear predictors without
  any modifications that may be due to the original response
  distribution (e.g., for a
  [`brms::cumulative()`](https://paulbuerkner.com/brms/reference/brmsfamily.html)
  model, the ordered thresholds are not taken into account). If (ii)
  `!is.null(ynew)`, then argument `type` also affects the scale of the
  log posterior predictive densities (`type = "response"` for the
  original response scale, `type = "link"` for the latent Gaussian
  scale).

- ...:

  Currently ignored.

## Value

In the following, \\N\\, \\C\_{\mathrm{cat}}\\, and
\\C\_{\mathrm{lat}}\\ from help topic
[refmodel-init-get](https://mc-stan.org/projpred/dev/reference/refmodel-init-get.md)
are used. Furthermore, let \\C\\ denote either \\C\_{\mathrm{cat}}\\ (if
`type = "response"`) or \\C\_{\mathrm{lat}}\\ (if `type = "link"`).
Then, if `is.null(ynew)`, the returned object contains the reference
model's predictions (with the scale depending on argument `type`) as:

- a length-\\N\\ vector in case of (i) the traditional projection, (ii)
  the latent projection with `type = "link"`, or (iii) the latent
  projection with `type = "response"` and `object$family$cats` being
  `NULL`;

- an \\N \times C\\ matrix in case of (i) the augmented-data projection
  or (ii) the latent projection with `type = "response"` and
  `object$family$cats` being not `NULL`.

If `!is.null(ynew)`, the returned object is a length-\\N\\ vector of log
posterior predictive densities evaluated at `ynew`.

## Details

Argument `weightsnew` is only relevant if `!is.null(ynew)`.

In case of a multilevel reference model, group-level effects for new
group levels are drawn randomly from a (multivariate) Gaussian
distribution. When setting `projpred.mlvl_pred_new` to `TRUE`, all group
levels from `newdata` (even those that already exist in the original
dataset) are treated as new group levels (if `is.null(newdata)`, all
group levels from the original dataset are considered as new group
levels in that case).
