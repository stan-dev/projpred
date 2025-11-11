# Print summary of a [`varsel()`](https://mc-stan.org/projpred/dev/reference/varsel.md) or [`cv_varsel()`](https://mc-stan.org/projpred/dev/reference/cv_varsel.md) run

This is the [`print()`](https://rdrr.io/r/base/print.html) method for
summary objects created by
[`summary.vsel()`](https://mc-stan.org/projpred/dev/reference/summary.vsel.md).
It displays a summary of the results from a
[`varsel()`](https://mc-stan.org/projpred/dev/reference/varsel.md) or
[`cv_varsel()`](https://mc-stan.org/projpred/dev/reference/cv_varsel.md)
run.

## Usage

``` r
# S3 method for class 'vselsummary'
print(x, digits = getOption("projpred.digits", 2), ...)
```

## Arguments

- x:

  An object of class `vselsummary`.

- digits:

  Passed to
  [`print.data.frame()`](https://rdrr.io/r/base/print.dataframe.html)
  (for the table containing the submodel performance evaluation results)
  and [`print.default()`](https://rdrr.io/r/base/print.default.html)
  (for the vector containing the reference model performance evaluation
  results).

- ...:

  Arguments passed to
  [`print.data.frame()`](https://rdrr.io/r/base/print.dataframe.html)
  (for the table containing the submodel performance evaluation results)
  and [`print.default()`](https://rdrr.io/r/base/print.default.html)
  (for the vector containing the reference model performance evaluation
  results).

## Value

The output of
[`summary.vsel()`](https://mc-stan.org/projpred/dev/reference/summary.vsel.md)
(invisible).

## Details

In the submodel predictive performance table printed at (or towards) the
bottom, column `ranking_fulldata` contains the full-data predictor
ranking and column `cv_proportions_diag` contains the main diagonal of
the matrix returned by
[`cv_proportions()`](https://mc-stan.org/projpred/dev/reference/cv_proportions.md)
(with `cumulate` as set in the
[`summary.vsel()`](https://mc-stan.org/projpred/dev/reference/summary.vsel.md)
call that created `x`). To retrieve the fold-wise predictor rankings,
use the
[`ranking()`](https://mc-stan.org/projpred/dev/reference/ranking.md)
function, possibly followed by
[`cv_proportions()`](https://mc-stan.org/projpred/dev/reference/cv_proportions.md)
for computing the ranking proportions (which can be visualized by
[`plot.cv_proportions()`](https://mc-stan.org/projpred/dev/reference/plot.cv_proportions.md)).
