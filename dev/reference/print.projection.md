# Print information about [`project()`](https://mc-stan.org/projpred/dev/reference/project.md) output

This is the [`print()`](https://rdrr.io/r/base/print.html) method for
objects of class `projection`. This method mainly exists to avoid
cluttering the console when printing such objects accidentally.

## Usage

``` r
# S3 method for class 'projection'
print(x, ...)
```

## Arguments

- x:

  An object of class `projection` (returned by
  [`project()`](https://mc-stan.org/projpred/dev/reference/project.md),
  possibly as elements of a `list`).

- ...:

  Currently ignored.

## Value

The input object `x` (invisible).
