# Predictor ranking(s)

Extracts the *predictor ranking(s)* from an object of class `vsel`
(returned by
[`varsel()`](https://mc-stan.org/projpred/dev/reference/varsel.md) or
[`cv_varsel()`](https://mc-stan.org/projpred/dev/reference/cv_varsel.md)).
A predictor ranking is simply a character vector of predictor terms
ranked by predictive relevance (with the most relevant term first). In
any case, objects of class `vsel` contain the predictor ranking based on
the *full-data* search. If an object of class `vsel` is based on a
cross-validation (CV) with fold-wise searches (i.e., if it was created
by
[`cv_varsel()`](https://mc-stan.org/projpred/dev/reference/cv_varsel.md)
with `validate_search = TRUE`), then it also contains *fold-wise*
predictor rankings.

## Usage

``` r
ranking(object, ...)

# S3 method for class 'vsel'
ranking(object, nterms_max = NULL, ...)
```

## Arguments

- object:

  The object from which to retrieve the predictor ranking(s). Possible
  classes may be inferred from the names of the corresponding methods
  (see also the description).

- ...:

  Currently ignored.

- nterms_max:

  Maximum submodel size (number of predictor terms) for the predictor
  ranking(s), i.e., the submodel size at which to cut off the predictor
  ranking(s). Using `NULL` is effectively the same as setting
  `nterms_max` to the full model size, i.e., this means to not cut off
  the predictor ranking(s) at all. Note that `nterms_max` does not count
  the intercept, so `nterms_max = 1` corresponds to the submodel
  consisting of the first (non-intercept) predictor term.

## Value

An object of class `ranking` which is a `list` with the following
elements:

- `fulldata`: The predictor ranking from the full-data search.

- `foldwise`: The predictor rankings from the fold-wise searches in the
  form of a character matrix (only available if `object` is based on a
  CV with fold-wise searches, otherwise element `foldwise` is `NULL`).
  The rows of this matrix correspond to the CV folds and the columns to
  the submodel sizes. Each row contains the predictor ranking from the
  search of that CV fold.

## See also

[`cv_proportions()`](https://mc-stan.org/projpred/dev/reference/cv_proportions.md)

## Examples

``` r
# For an example, see `?plot.cv_proportions`.
```
