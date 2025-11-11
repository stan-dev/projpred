# Ranking proportions from fold-wise predictor rankings

Calculates the *ranking proportions* from the fold-wise predictor
rankings in a cross-validation (CV) with fold-wise searches. For a given
predictor \\x\\ and a given submodel size \\j\\, the ranking proportion
is the proportion of CV folds which have predictor \\x\\ at position
\\j\\ of their predictor ranking. While these ranking proportions are
helpful for investigating variability in the predictor ranking, they can
also be *cumulated* across submodel sizes. The cumulated ranking
proportions are more helpful when it comes to model selection.

## Usage

``` r
cv_proportions(object, ...)

# S3 method for class 'ranking'
cv_proportions(object, cumulate = FALSE, ...)

# S3 method for class 'vsel'
cv_proportions(object, ...)
```

## Arguments

- object:

  For `cv_proportions.ranking()`: an object of class `ranking` (returned
  by
  [`ranking()`](https://mc-stan.org/projpred/dev/reference/ranking.md)).
  For `cv_proportions.vsel()`: an object of class `vsel` (returned by
  [`varsel()`](https://mc-stan.org/projpred/dev/reference/varsel.md) or
  [`cv_varsel()`](https://mc-stan.org/projpred/dev/reference/cv_varsel.md))
  that
  [`ranking()`](https://mc-stan.org/projpred/dev/reference/ranking.md)
  will be applied to internally before then calling
  `cv_proportions.ranking()`.

- ...:

  For `cv_proportions.vsel()`: arguments passed to
  [`ranking.vsel()`](https://mc-stan.org/projpred/dev/reference/ranking.md)
  and `cv_proportions.ranking()`. For `cv_proportions.ranking()`:
  currently ignored.

- cumulate:

  A single logical value indicating whether the ranking proportions
  should be cumulated across increasing submodel sizes (`TRUE`) or not
  (`FALSE`).

## Value

A numeric matrix containing the ranking proportions. This matrix has
`nterms_max` rows and `nterms_max` columns, with `nterms_max` as
specified in the (possibly implicit)
[`ranking()`](https://mc-stan.org/projpred/dev/reference/ranking.md)
call. The rows correspond to the submodel sizes and the columns to the
predictor terms (sorted according to the full-data predictor ranking).
If `cumulate` is `FALSE`, then the returned matrix is of class
`cv_proportions`. If `cumulate` is `TRUE`, then the returned matrix is
of classes `cv_proportions_cumul` and `cv_proportions` (in this order).

Note that if `cumulate` is `FALSE`, then the values in the returned
matrix only need to sum to 1 (column-wise and row-wise) if `nterms_max`
(see above) is equal to the full model size. Likewise, if `cumulate` is
`TRUE`, then the value `1` only needs to occur in each column of the
returned matrix if `nterms_max` is equal to the full model size.

The `cv_proportions()` function is only applicable if the `ranking`
object includes fold-wise predictor rankings (i.e., if it is based on a
`vsel` object created by
[`cv_varsel()`](https://mc-stan.org/projpred/dev/reference/cv_varsel.md)
with `validate_search = TRUE`). If the `ranking` object contains only a
full-data predictor ranking (i.e., if it is based on a `vsel` object
created by
[`varsel()`](https://mc-stan.org/projpred/dev/reference/varsel.md) or by
[`cv_varsel()`](https://mc-stan.org/projpred/dev/reference/cv_varsel.md),
but the latter with `validate_search = FALSE`), then an error is thrown
because in that case, there are no fold-wise predictor rankings from
which to calculate ranking proportions.

## See also

[`plot.cv_proportions()`](https://mc-stan.org/projpred/dev/reference/plot.cv_proportions.md)

## Examples

``` r
# For an example, see `?plot.cv_proportions`.
```
