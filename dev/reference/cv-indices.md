# Create cross-validation folds

These are helper functions to create cross-validation (CV) folds, i.e.,
to split up the indices from 1 to `n` into `K` subsets ("folds") for
\\K\\-fold CV. These functions are potentially useful when creating the
input for arguments `cvfits` and `cvfun` of
[`init_refmodel()`](https://mc-stan.org/projpred/dev/reference/refmodel-init-get.md)
(or argument `cvfits` of
[`cv_varsel.refmodel()`](https://mc-stan.org/projpred/dev/reference/cv_varsel.md)).
Function `cvfolds()` is deprecated; please use `cv_folds()` instead
(apart from the name, they are the same). The return value of
`cv_folds()` and `cv_ids()` is different, see below for details.

## Usage

``` r
cv_folds(n, K, seed = NA)

cvfolds(n, K, seed = NA)

cv_ids(n, K, out = c("foldwise", "indices"), seed = NA)
```

## Arguments

- n:

  Number of observations.

- K:

  Number of folds. Must be at least 2 and not exceed `n`.

- seed:

  Pseudorandom number generation (PRNG) seed by which the same results
  can be obtained again if needed. Passed to argument `seed` of
  [`set.seed()`](https://rdrr.io/r/base/Random.html), but can also be
  `NA` to not call [`set.seed()`](https://rdrr.io/r/base/Random.html) at
  all. If not `NA`, then the PRNG state is reset (to the state before
  calling `cv_folds()` or `cv_ids()`) upon exiting `cv_folds()` or
  `cv_ids()`.

- out:

  Format of the output, either `"foldwise"` or `"indices"`. See below
  for details.

## Value

`cv_folds()` returns a vector of length `n` such that each element is an
integer between 1 and `K` denoting which fold the corresponding data
point belongs to. The return value of `cv_ids()` depends on the `out`
argument. If `out = "foldwise"`, the return value is a `list` with `K`
elements, each being a `list` with elements `tr` and `ts` giving the
training and test indices, respectively, for the corresponding fold. If
`out = "indices"`, the return value is a `list` with elements `tr` and
`ts` each being a `list` with `K` elements giving the training and test
indices, respectively, for each fold.

## Examples

``` r
n <- 100
set.seed(1234)
y <- rnorm(n)
cv <- cv_ids(n, K = 5)
# Mean within the test set of each fold:
cvmeans <- sapply(cv, function(fold) mean(y[fold$ts]))
```
