# Retrieve the full-data solution path from a [`varsel()`](https://mc-stan.org/projpred/dev/reference/varsel.md) or [`cv_varsel()`](https://mc-stan.org/projpred/dev/reference/cv_varsel.md) run or the predictor combination from a [`project()`](https://mc-stan.org/projpred/dev/reference/project.md) run

The `solution_terms.vsel()` method retrieves the solution path from a
full-data search (`vsel` objects are returned by
[`varsel()`](https://mc-stan.org/projpred/dev/reference/varsel.md) or
[`cv_varsel()`](https://mc-stan.org/projpred/dev/reference/cv_varsel.md)).
The `solution_terms.projection()` method retrieves the predictor
combination onto which a projection was performed (`projection` objects
are returned by
[`project()`](https://mc-stan.org/projpred/dev/reference/project.md),
possibly as elements of a `list`). Both methods (and hence also the
`solution_terms()` generic) are deprecated and will be removed in a
future release. Please use
[`ranking()`](https://mc-stan.org/projpred/dev/reference/ranking.md)
instead of `solution_terms.vsel()`
([`ranking()`](https://mc-stan.org/projpred/dev/reference/ranking.md)'s
output element `fulldata` contains the full-data predictor ranking that
is extracted by `solution_terms.vsel()`;
[`ranking()`](https://mc-stan.org/projpred/dev/reference/ranking.md)'s
output element `foldwise` contains the fold-wise predictor rankings—if
available—which were previously not accessible via a built-in function)
and
[`predictor_terms()`](https://mc-stan.org/projpred/dev/reference/predictor_terms.md)
instead of `solution_terms.projection()`.

## Usage

``` r
solution_terms(object, ...)

# S3 method for class 'vsel'
solution_terms(object, ...)

# S3 method for class 'projection'
solution_terms(object, ...)
```

## Arguments

- object:

  The object from which to retrieve the predictor terms. Possible
  classes may be inferred from the names of the corresponding methods
  (see also the description).

- ...:

  Currently ignored.

## Value

A character vector of predictor terms.
