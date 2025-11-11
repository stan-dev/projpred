# Predictive performance results

Retrieves the predictive performance summaries after running
[`varsel()`](https://mc-stan.org/projpred/dev/reference/varsel.md) or
[`cv_varsel()`](https://mc-stan.org/projpred/dev/reference/cv_varsel.md).
These summaries are computed by
[`summary.vsel()`](https://mc-stan.org/projpred/dev/reference/summary.vsel.md),
so the main method of `performances()` is `performances.vselsummary()`
(objects of class `vselsummary` are returned by
[`summary.vsel()`](https://mc-stan.org/projpred/dev/reference/summary.vsel.md)).
As a shortcut method, `performances.vsel()` is provided as well (objects
of class `vsel` are returned by
[`varsel()`](https://mc-stan.org/projpred/dev/reference/varsel.md) and
[`cv_varsel()`](https://mc-stan.org/projpred/dev/reference/cv_varsel.md)).
For a graphical representation, see
[`plot.vsel()`](https://mc-stan.org/projpred/dev/reference/plot.vsel.md).

## Usage

``` r
performances(object, ...)

# S3 method for class 'vselsummary'
performances(object, ...)

# S3 method for class 'vsel'
performances(object, ...)
```

## Arguments

- object:

  The object from which to retrieve the predictive performance results.
  Possible classes may be inferred from the names of the corresponding
  methods (see also the description).

- ...:

  For `performances.vsel()`: arguments passed to
  [`summary.vsel()`](https://mc-stan.org/projpred/dev/reference/summary.vsel.md).
  For `performances.vselsummary()`: currently ignored.

## Value

An object of class `performances` which is a `list` with the following
elements:

- `submodels`: The predictive performance results for the submodels, as
  a `data.frame`.

- `reference_model`: The predictive performance results for the
  reference model, as a named vector.

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
print(performances(vs))
#> $submodels
#>   size      elpd  elpd.se   elpd.diff elpd.diff.se
#> 1    0 -249.1981 5.256908 -39.1341918     5.759373
#> 2    1 -230.5763 5.621175 -20.5123784     4.479675
#> 3    2 -219.8008 6.029368  -9.7369536     3.363230
#> 4    3 -210.5013 6.559977  -0.4374522     0.896227
#> 
#> $reference_model
#>        elpd     elpd.se 
#> -210.063880    6.562731 
#> 
#> attr(,"class")
#> [1] "performances"
```
