# Execute a function call

Execute a function call similar to
[`do.call()`](https://rdrr.io/r/base/do.call.html), but without
deparsing function arguments.

## Usage

``` r
do_call(what, args, pkg = NULL)
```

## Arguments

- what:

  Either a function or a non-empty character string naming the function
  to be called.

- args:

  A `list` of arguments to the function call. The
  [`names`](https://rdrr.io/r/base/names.html) attribute of `args` gives
  the argument names.

- pkg:

  Optional name of the package in which to search for the function if
  `what` is a character string.

## Value

The result of the (evaluated) function call.
