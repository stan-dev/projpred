# Print information about a reference model object

This is the [`print()`](https://rdrr.io/r/base/print.html) method for
reference model objects (objects of class `refmodel`). This method
mainly exists to avoid cluttering the console when printing such objects
accidentally.

## Usage

``` r
# S3 method for class 'refmodel'
print(x, ...)
```

## Arguments

- x:

  An object of class `refmodel` (returned by
  [`get_refmodel()`](https://mc-stan.org/projpred/dev/reference/refmodel-init-get.md)
  or
  [`init_refmodel()`](https://mc-stan.org/projpred/dev/reference/refmodel-init-get.md)).

- ...:

  Currently ignored.

## Value

The input object `x` (invisible).
