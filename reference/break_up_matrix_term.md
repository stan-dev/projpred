# Break up matrix terms

Sometimes there can be terms in a formula that refer to a matrix instead
of a single predictor. This function breaks up the matrix term into
individual predictors to handle separately, as that is probably the
intention of the user.

## Usage

``` r
break_up_matrix_term(formula, data)
```

## Arguments

- formula:

  A [`formula`](https://rdrr.io/r/stats/formula.html) for a valid model.

- data:

  The original `data.frame` with a matrix as predictor.

## Value

A `list` containing the expanded
[`formula`](https://rdrr.io/r/stats/formula.html) and the expanded
`data.frame`.
