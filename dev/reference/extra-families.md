# Extra family objects

Family objects not in the set of default
[`family`](https://rdrr.io/r/stats/family.html) objects.

## Usage

``` r
Student_t(link = "identity", nu = 3)
```

## Arguments

- link:

  Name of the link function. In contrast to the default
  [`family`](https://rdrr.io/r/stats/family.html) objects, this has to
  be a character string here.

- nu:

  Degrees of freedom for the Student-\\t\\ distribution.

## Value

A family object analogous to those described in
[`family`](https://rdrr.io/r/stats/family.html).

## Note

Support for the `Student_t()` family is still experimental.
