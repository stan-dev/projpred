# Predictor terms used in a [`project()`](https://mc-stan.org/projpred/dev/reference/project.md) run

For a `projection` object (returned by
[`project()`](https://mc-stan.org/projpred/dev/reference/project.md),
possibly as elements of a `list`), this function extracts the
combination of predictor terms onto which the projection was performed.

## Usage

``` r
predictor_terms(object, ...)

# S3 method for class 'projection'
predictor_terms(object, ...)
```

## Arguments

- object:

  An object of class `projection` (returned by
  [`project()`](https://mc-stan.org/projpred/dev/reference/project.md),
  possibly as elements of a `list`) from which to retrieve the predictor
  terms.

- ...:

  Currently ignored.

## Value

A character vector of predictor terms.

## Examples

``` r
# Data:
dat_gauss <- data.frame(y = df_gaussian$y, df_gaussian$x)

# The `stanreg` fit which will be used as the reference model (with small
# values for `chains` and `iter`, but only for technical reasons in this
# example; this is not recommended in general):
fit <- rstanarm::stan_glm(
  y ~ X1 + X2 + X3 + X4 + X5, family = gaussian(), data = dat_gauss,
  QR = TRUE, chains = 2, iter = 500, refresh = 0, seed = 9876
)

# Projection onto an arbitrary combination of predictor terms (with a small
# value for `nclusters`, but only for the sake of speed in this example;
# this is not recommended in general):
prj <- project(fit, predictor_terms = c("X1", "X3", "X5"), nclusters = 10,
               seed = 9182)
print(predictor_terms(prj)) # gives `c("X1", "X3", "X5")`
#> [1] "X1" "X3" "X5"
```
