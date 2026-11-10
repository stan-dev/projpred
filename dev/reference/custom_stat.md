# Custom predictive performance statistic

Defines a predictive performance statistic (a utility or a loss) from a
user function. Use the result in argument `stats` of
[`summary.vsel()`](https://mc-stan.org/projpred/dev/reference/summary.vsel.md),
[`plot.vsel()`](https://mc-stan.org/projpred/dev/reference/plot.vsel.md)
and
[`performances()`](https://mc-stan.org/projpred/dev/reference/performances.md),
or in argument `stat` of
[`suggest_size()`](https://mc-stan.org/projpred/dev/reference/suggest_size.md).
The statistic affects the evaluation of the submodels and the suggested
size, but not the search.

## Usage

``` r
custom_stat(
  fun,
  name,
  type = c("pointwise", "aggregate"),
  higher_better = TRUE,
  label = name
)
```

## Arguments

- fun:

  A function with arguments `y`, `mu` and `wobs` (see section
  "Details").

- name:

  A single character string. The output of
  [`summary.vsel()`](https://mc-stan.org/projpred/dev/reference/summary.vsel.md)
  and
  [`performances()`](https://mc-stan.org/projpred/dev/reference/performances.md)
  uses it in the column names and adds suffixes such as `.se`, so it
  must start with a letter and contain only letters, digits and
  underscores. It must not be the name of a built-in statistic or of
  another column of these tables (`"size"`, `"ranking_fulldata"`,
  `"cv_proportions_diag"`).

- type:

  Either `"pointwise"` or `"aggregate"` (see section "Details").

- higher_better:

  `TRUE` for a utility (higher values are better), `FALSE` for a loss.
  `baseline = "best"` and
  [`suggest_size()`](https://mc-stan.org/projpred/dev/reference/suggest_size.md)
  use it.

- label:

  A single character string for the facet label in
  [`plot.vsel()`](https://mc-stan.org/projpred/dev/reference/plot.vsel.md).

## Value

An object of class `custom_stat`.

## Details

Function `fun` gets these arguments for each model (the reference model
and each submodel size):

- `y`: the observed response values of the evaluation (test) set. For
  the [`binomial()`](https://rdrr.io/r/stats/family.html) family, `y`
  holds the proportions of successes (the number of successes divided by
  the number of trials). For the latent projection, `y` holds the
  response values on the original response scale if `resp_oscale = TRUE`
  and the latent response values if `resp_oscale = FALSE`.

- `mu`: the predictive means for the observations in `y`, averaged over
  the parameter draws of the reference model or over the projected draws
  of the submodel. For the
  [`binomial()`](https://rdrr.io/r/stats/family.html) family, these are
  the predicted probabilities of success. With cross-validation, `mu`
  holds the cross-validated predictions. For the latent projection, `mu`
  is on the same scale as `y`. `mu` can contain `NA`s. With subsampled
  LOO CV (see argument `nloo` of
  [`cv_varsel()`](https://mc-stan.org/projpred/dev/reference/cv_varsel.md)),
  the submodel predictions are `NA` for the observations outside the
  subsample. For a reference model of class `datafit`, all reference
  model predictions are `NA`.

- `wobs`: the observation weights (for the
  [`binomial()`](https://rdrr.io/r/stats/family.html) family: the number
  of trials).

For `type = "pointwise"`, `fun` returns a numeric (or logical) vector of
the same length as `y` with the utility (or loss) of each observation.
Element `i` must depend only on `y[i]`, `mu[i]` and `wobs[i]`. Do not
apply `wobs` in `fun`: the statistic is the mean of these values,
weighted by `wobs`. As for `"mse"`, the standard error treats each row
of the evaluation set as one observation. For the
[`binomial()`](https://rdrr.io/r/stats/family.html) family, a row with
several trials thus counts once, whereas `"acc"` splits it into single
trials. The difference from the baseline model pairs the two models
observation by observation. The uncertainty intervals use a normal
approximation. Subsampled LOO CV works with this type.

For `type = "aggregate"`, `fun` returns a single number for the whole
evaluation set. As for `"auc"`, the standard error and the uncertainty
intervals come from a bootstrap over the rows of the evaluation set (see
argument `...` of
[`summary.vsel()`](https://mc-stan.org/projpred/dev/reference/summary.vsel.md)
for the number of bootstrap samples and the seed). The difference from
the baseline model uses the same bootstrap samples for both models. For
the [`binomial()`](https://rdrr.io/r/stats/family.html) family, a row
with several trials counts once, whereas `"auc"` splits it into single
trials. This type does not support subsampled LOO CV.

Custom statistics are not available for the augmented-data projection,
or for the latent projection with `resp_oscale = TRUE` and a non-`NULL`
element `family$cats`.

## Examples

``` r
# Net benefit at the decision threshold `thres` (a utility). Each observation
# contributes one pointwise term:
net_benefit <- function(thres) {
  custom_stat(
    function(y, mu, wobs) {
      treat <- mu >= thres
      treat * (y - (1 - y) * thres / (1 - thres))
    },
    name = paste0("nb", 100 * thres),
    label = paste0("Net benefit \u2265 ", thres)
  )
}

# Data:
dat_binom <- data.frame(y = df_binom$y, df_binom$x)

# The `stanreg` fit which will be used as the reference model (with small
# values for `chains` and `iter`, but only for technical reasons in this
# example; this is not recommended in general):
fit <- rstanarm::stan_glm(
  y ~ X1 + X2 + X3 + X4 + X5, family = binomial(), data = dat_binom,
  chains = 2, iter = 500, refresh = 0, seed = 9876
)

# Run varsel() (here without cross-validation and with small values for
# `nterms_max` and `nclusters_pred`, but only for the sake of speed in this
# example; this is not recommended in general):
vs <- varsel(fit, nterms_max = 3, nclusters_pred = 10, seed = 5555)
print(summary(vs, stats = list("elpd", net_benefit(0.2))), digits = 2)
#> 
#> Family: binomial 
#> Link function: logit 
#> 
#> Formula: y ~ X1 + X2 + X3 + X4 + X5
#> Observations: 100
#> Projection method: traditional
#> Search method: forward
#> Maximum submodel size for the search: 3
#> Number of projected draws in the search: 20 (from clustered projection)
#> Number of projected draws in the performance evaluation: 10 (from clustered projection)
#> Argument `refit_prj`: TRUE
#> 
#> Submodel performance evaluation summary with `deltas = FALSE` and `cumulate = FALSE`:
#>  size ranking_fulldata cv_proportions_diag elpd elpd.se elpd.diff elpd.diff.se
#>     0      (Intercept)                  NA  -69     1.2     -20.3          5.1
#>     1               X4                  NA  -55     4.6      -6.8          3.6
#>     2               X5                  NA  -52     4.9      -3.3          2.7
#>     3               X2                  NA  -50     5.0      -1.7          1.9
#>  nb20 nb20.se nb20.diff nb20.diff.se
#>  0.45   0.062   -0.0050        0.020
#>  0.47   0.060    0.0100        0.015
#>  0.45   0.060   -0.0025        0.012
#>  0.46   0.059    0.0050        0.012
#> 
#> Reference model performance evaluation summary with `deltas = FALSE`:
#>    elpd elpd.se    nb20 nb20.se 
#> -48.327   5.162   0.455   0.059 
```
