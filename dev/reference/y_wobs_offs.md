# Extract response values, observation weights, and offsets

A helper function for extracting response values, observation weights,
and offsets from a dataset. It is designed for use in the
`extract_model_data` function of custom reference model objects (see
[`init_refmodel()`](https://mc-stan.org/projpred/dev/reference/refmodel-init-get.md)).

## Usage

``` r
y_wobs_offs(newdata, wrhs = NULL, orhs = NULL, resp_form)
```

## Arguments

- newdata:

  The `data.frame` from which at least the response values should be
  extracted.

- wrhs:

  Either a right-hand side formula consisting only of the variable in
  `newdata` containing the weights, `NULL` (for using a vector of ones),
  or directly the numeric vector of observation weights.

- orhs:

  Either a right-hand side formula consisting only of the variable in
  `newdata` containing the offsets, `NULL` (for using a vector of
  zeros), or directly the numeric vector of offsets.

- resp_form:

  If this is a formula, then the second element of this formula (if the
  formula is a standard formula with both left-hand and right-hand side,
  then its second element is the left-hand side; if the formula is a
  right-hand side formula, then its second element is the right-hand
  side) will be extracted from `newdata` (so `resp_form` may be either a
  standard formula or a right-hand side formula, but in the latter case,
  the right-hand side should consist only of the response variable). In
  all other cases, `NULL` will be returned for element `y` of the output
  `list`.

## Value

A `list` with elements `y`, `weights`, and `offset`, each being a
numeric vector containing the data for the response, the observation
weights, and the offsets, respectively. An exception is that `y` may
also be `NULL` (depending on argument `resp_form`), a non-numeric
vector, or a `factor`.

## See also

[`init_refmodel()`](https://mc-stan.org/projpred/dev/reference/refmodel-init-get.md)
