# Print results (summary) of a [`varsel()`](https://mc-stan.org/projpred/dev/reference/varsel.md) or [`cv_varsel()`](https://mc-stan.org/projpred/dev/reference/cv_varsel.md) run

This is the [`print()`](https://rdrr.io/r/base/print.html) method for
`vsel` objects (returned by
[`varsel()`](https://mc-stan.org/projpred/dev/reference/varsel.md) or
[`cv_varsel()`](https://mc-stan.org/projpred/dev/reference/cv_varsel.md)).
It displays a summary of a
[`varsel()`](https://mc-stan.org/projpred/dev/reference/varsel.md) or
[`cv_varsel()`](https://mc-stan.org/projpred/dev/reference/cv_varsel.md)
run by first calling
[`summary.vsel()`](https://mc-stan.org/projpred/dev/reference/summary.vsel.md)
and then
[`print.vselsummary()`](https://mc-stan.org/projpred/dev/reference/print.vselsummary.md).

## Usage

``` r
# S3 method for class 'vsel'
print(x, digits = getOption("projpred.digits", 2), ...)
```

## Arguments

- x:

  An object of class `vsel` (returned by
  [`varsel()`](https://mc-stan.org/projpred/dev/reference/varsel.md) or
  [`cv_varsel()`](https://mc-stan.org/projpred/dev/reference/cv_varsel.md)).

- digits:

  Passed to argument `digits` of
  [`print.vselsummary()`](https://mc-stan.org/projpred/dev/reference/print.vselsummary.md).

- ...:

  Arguments passed to
  [`summary.vsel()`](https://mc-stan.org/projpred/dev/reference/summary.vsel.md).

## Value

The output of
[`summary.vsel()`](https://mc-stan.org/projpred/dev/reference/summary.vsel.md)
(invisible).
