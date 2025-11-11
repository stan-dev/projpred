# Summary of a [`varsel()`](https://mc-stan.org/projpred/dev/reference/varsel.md) or [`cv_varsel()`](https://mc-stan.org/projpred/dev/reference/cv_varsel.md) run

This is the [`summary()`](https://rdrr.io/r/base/summary.html) method
for `vsel` objects (returned by
[`varsel()`](https://mc-stan.org/projpred/dev/reference/varsel.md) or
[`cv_varsel()`](https://mc-stan.org/projpred/dev/reference/cv_varsel.md)).
Apart from some general information about the
[`varsel()`](https://mc-stan.org/projpred/dev/reference/varsel.md) or
[`cv_varsel()`](https://mc-stan.org/projpred/dev/reference/cv_varsel.md)
run, it shows the full-data predictor ranking, basic information about
the (CV) variability in the ranking of the predictors (if available;
inferred from
[`cv_proportions()`](https://mc-stan.org/projpred/dev/reference/cv_proportions.md)),
and estimates for user-specified predictive performance statistics. For
a graphical representation, see
[`plot.vsel()`](https://mc-stan.org/projpred/dev/reference/plot.vsel.md).
For extracting the predictive performance results printed at the bottom
of the output created by this
[`summary()`](https://rdrr.io/r/base/summary.html) method, see
[`performances()`](https://mc-stan.org/projpred/dev/reference/performances.md).

## Usage

``` r
# S3 method for class 'vsel'
summary(
  object,
  nterms_max = NULL,
  stats = "elpd",
  type = c("mean", "se", "diff", "diff.se"),
  deltas = FALSE,
  alpha = 2 * pnorm(-1),
  baseline = if (!inherits(object$refmodel, "datafit")) "ref" else "best",
  resp_oscale = TRUE,
  cumulate = FALSE,
  ...
)
```

## Arguments

- object:

  An object of class `vsel` (returned by
  [`varsel()`](https://mc-stan.org/projpred/dev/reference/varsel.md) or
  [`cv_varsel()`](https://mc-stan.org/projpred/dev/reference/cv_varsel.md)).

- nterms_max:

  Maximum submodel size (number of predictor terms) for which the
  performance statistics are calculated. Using `NULL` is effectively the
  same as `length(ranking(object)$fulldata)`. Note that `nterms_max`
  does not count the intercept, so use `nterms_max = 0` for the
  intercept-only model. For
  [`plot.vsel()`](https://mc-stan.org/projpred/dev/reference/plot.vsel.md),
  `nterms_max` must be at least `1`.

- stats:

  One or more character strings determining which performance statistics
  (i.e., utilities or losses) to estimate based on the observations in
  the evaluation (or "test") set (in case of cross-validation, these are
  all observations because they are partitioned into multiple test sets;
  in case of
  [`varsel()`](https://mc-stan.org/projpred/dev/reference/varsel.md)
  with `d_test = NULL`, these are again all observations because the
  test set is the same as the training set). Available statistics are:

  - `"elpd"`: expected log (pointwise) predictive density (for a new
    dataset) (ELPD). Estimated by the sum of the observation-specific
    log predictive density values (with each of these predictive density
    values being a—possibly weighted—average across the parameter
    draws). For the corresponding uncertainty interval, a normal
    approximation is used.

  - `"mlpd"`: mean log predictive density (MLPD), that is, the ELPD
    divided by the number of observations. For the corresponding
    uncertainty interval, a normal approximation is used.

  - `"gmpd"`: geometric mean predictive density (GMPD), that is,
    [`exp()`](https://rdrr.io/r/base/Log.html) of the MLPD. The GMPD is
    especially helpful for discrete response families (because there,
    the GMPD is bounded by zero and one). For the corresponding standard
    error, the delta method is used. The corresponding uncertainty
    interval type is "exponentiated normal approximation" because the
    uncertainty interval bounds are the exponentiated uncertainty
    interval bounds of the MLPD.

  - `"mse"`: mean squared error (only available in the situations
    mentioned in section "Details" below). For the corresponding
    uncertainty interval, a log-normal approximation is used if `deltas`
    is `FALSE` and a normal approximation is used if `deltas` is `TRUE`
    (or `"mixed"`, in case of
    [`plot.vsel()`](https://mc-stan.org/projpred/dev/reference/plot.vsel.md)).

  - `"rmse"`: root mean squared error (only available in the situations
    mentioned in section "Details" below). For the corresponding
    standard error, the delta method is used. For the corresponding
    uncertainty interval, a log-normal approximation is used if `deltas`
    is `FALSE` and a normal approximation is used if `deltas` is `TRUE`
    (or `"mixed"`, in case of
    [`plot.vsel()`](https://mc-stan.org/projpred/dev/reference/plot.vsel.md)).

  - `"R2"`: R-squared, i.e., coefficient of determination (only
    available in the situations mentioned in section "Details" below).
    For the corresponding standard error, the delta method is used. For
    the corresponding uncertainty interval, a normal approximation is
    used.

  - `"acc"` (or its alias, `"pctcorr"`): classification accuracy (only
    available in the situations mentioned in section "Details" below).
    By "classification accuracy", we mean the proportion of correctly
    classified observations. For this, the response category ("class")
    with highest probability (the probabilities are model-based) is
    taken as the prediction ("classification") for an observation. For
    the corresponding uncertainty interval, a normal approximation is
    used.

  - `"auc"`: area under the ROC curve (only available in the situations
    mentioned in section "Details" below). For the corresponding
    standard error and lower and upper uncertainty interval bounds,
    bootstrapping is used. Not supported in case of subsampled LOO-CV
    (see argument `nloo` of
    [`cv_varsel()`](https://mc-stan.org/projpred/dev/reference/cv_varsel.md)).

- type:

  One or more items from `"mean"`, `"se"`, `"lower"`, `"upper"`,
  `"diff"`, `"diff.lower"`, `"diff.upper"`, and `"diff.se"` indicating
  which of these to compute for each item from `stats` (mean, standard
  error, lower and upper uncertainty interval bounds, mean difference to
  the corresponding statistic of the reference model, lower and upper
  uncertainty interval bound for this difference, and standard error of
  this difference, respectively; note that for the GMPD, `"diff"`,
  `"diff.lower"`, `"diff.upper"`, and `"diff.se"` actually refer to the
  ratio vs. the reference model, not the difference). The uncertainty
  interval bounds belong to uncertainty intervals with (nominal)
  coverage `1 - alpha`. Items `"diff"`, `"diff.lower"`, `"diff.upper"`,
  and `"diff.se"` are only supported if `deltas` is `FALSE`.

- deltas:

  May be set to `FALSE` or `TRUE`. If `FALSE`, the submodel performance
  statistics are estimated on their actual scale. If `TRUE`, the
  submodel statistics are estimated relatively to the baseline model
  (see argument `baseline`). For the GMPD, the term "relatively" refers
  to the *ratio* vs. the baseline model (i.e., the submodel statistic
  divided by the baseline model statistic). For all other `stats`,
  "relatively" refers to the *difference* from the baseline model (i.e.,
  the submodel statistic minus the baseline model statistic).

- alpha:

  A number determining the (nominal) coverage `1 - alpha` of the
  uncertainty intervals. For example, in case of a normal-approximation
  uncertainty interval, `alpha = 2 * pnorm(-1)` corresponds to a
  uncertainty interval stretching by one standard error on either side
  of the point estimate.

- baseline:

  For `summary.vsel()`: Only relevant if `deltas` is `TRUE`. For
  [`plot.vsel()`](https://mc-stan.org/projpred/dev/reference/plot.vsel.md):
  Always relevant. Either `"ref"` or `"best"`, indicating whether the
  baseline is the reference model or the best submodel found (in terms
  of `stats[1]`), respectively. In case of subsampled LOO-CV,
  `baseline = "best"` is not supported.

- resp_oscale:

  Only relevant for the latent projection. A single logical value
  indicating whether to calculate the performance statistics on the
  original response scale (`TRUE`) or on latent scale (`FALSE`).

- cumulate:

  Passed to argument `cumulate` of
  [`cv_proportions()`](https://mc-stan.org/projpred/dev/reference/cv_proportions.md).
  Affects column `cv_proportions_diag` of the summary table.

- ...:

  Arguments passed to the internal function which is used for
  bootstrapping (if applicable; see argument `stats`). Currently,
  relevant arguments are `B` (the number of bootstrap samples,
  defaulting to `2000`) and `seed` (see
  [`set.seed()`](https://rdrr.io/r/base/Random.html), but defaulting to
  `NA` so that [`set.seed()`](https://rdrr.io/r/base/Random.html) is not
  called within that function at all).

## Value

An object of class `vselsummary`. The elements of this object are not
meant to be accessed directly but instead via helper functions
([`print.vselsummary()`](https://mc-stan.org/projpred/dev/reference/print.vselsummary.md)
and
[`performances.vselsummary()`](https://mc-stan.org/projpred/dev/reference/performances.md)).

## Details

The `stats` options `"mse"`, `"rmse"`, and `"R2"` are only available
for:

- the traditional projection,

- the latent projection with `resp_oscale = FALSE`,

- the latent projection with `resp_oscale = TRUE` in combination with
  `<refmodel>$family$cats` being `NULL`.

The `stats` option `"acc"` (= `"pctcorr"`) is only available for:

- the [`binomial()`](https://rdrr.io/r/stats/family.html) family in case
  of the traditional projection,

- all families in case of the augmented-data projection,

- the [`binomial()`](https://rdrr.io/r/stats/family.html) family (on the
  original response scale) in case of the latent projection with
  `resp_oscale = TRUE` in combination with `<refmodel>$family$cats`
  being `NULL`,

- all families (on the original response scale) in case of the latent
  projection with `resp_oscale = TRUE` in combination with
  `<refmodel>$family$cats` being not `NULL`.

The `stats` option `"auc"` is only available for:

- the [`binomial()`](https://rdrr.io/r/stats/family.html) family in case
  of the traditional projection,

- the [`binomial()`](https://rdrr.io/r/stats/family.html) family (on the
  original response scale) in case of the latent projection with
  `resp_oscale = TRUE` in combination with `<refmodel>$family$cats`
  being `NULL`.

Note that the `stats` option `"auc"` is not supported in case of
subsampled LOO-CV (see argument `nloo` of
[`cv_varsel()`](https://mc-stan.org/projpred/dev/reference/cv_varsel.md)).

## See also

[`print.vselsummary()`](https://mc-stan.org/projpred/dev/reference/print.vselsummary.md),
[`performances.vselsummary()`](https://mc-stan.org/projpred/dev/reference/performances.md)

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
print(summary(vs), digits = 1)
#> 
#> Family: gaussian 
#> Link function: identity 
#> 
#> Formula: y ~ X1 + X2 + X3 + X4 + X5
#> Observations: 100
#> Projection method: traditional
#> Search method: L1
#> Maximum submodel size for the search: 3
#> Number of projected draws in the search: 1 (from clustered projection)
#> Number of projected draws in the performance evaluation: 10 (from clustered projection)
#> Argument `refit_prj`: TRUE
#> 
#> Submodel performance evaluation summary with `deltas = FALSE` and `cumulate = FALSE`:
#>  size ranking_fulldata cv_proportions_diag elpd elpd.se elpd.diff elpd.diff.se
#>     0      (Intercept)                  NA -249       5     -39.1          5.8
#>     1               X1                  NA -231       6     -20.5          4.5
#>     2               X5                  NA -220       6      -9.7          3.4
#>     3               X3                  NA -211       7      -0.4          0.9
#> 
#> Reference model performance evaluation summary with `deltas = FALSE`:
#>    elpd elpd.se 
#>    -210       7 
```
