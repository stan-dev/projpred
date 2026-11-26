# Formula of the submodel of a [`project()`](https://mc-stan.org/projpred/dev/reference/project.md) run

For a `projection` object (returned by
[`project()`](https://mc-stan.org/projpred/dev/reference/project.md)),
this method gives the formula of the submodel of the projection. The
formula has the form that projpred uses for the submodel fits. Use it,
for example, to set up the submodel in another package. For a `list` of
`projection` objects, use
[`lapply()`](https://rdrr.io/r/base/lapply.html).

## Usage

``` r
# S3 method for class 'projection'
formula(x, ...)
```

## Arguments

- x:

  An object of class `projection` (returned by
  [`project()`](https://mc-stan.org/projpred/dev/reference/project.md)).

- ...:

  Currently ignored.

## Value

A [`formula`](https://rdrr.io/r/stats/formula.html).

## Details

The right-hand side has an intercept and the predictor terms of `x` (see
[`predictor_terms()`](https://mc-stan.org/projpred/dev/reference/predictor_terms.md)),
with two changes that projpred applies to its submodel fits:

- The group-level terms of one grouping factor are merged into one term.
  For example, `(1 | g) + (x1 | g)` becomes `(x1 | g)`.

- A linear main-effect term is dropped if the submodel also has a smooth
  term of the same variable. For example, `x1 + s(x1)` becomes `s(x1)`.

The term labels of the formula can differ from
[`predictor_terms()`](https://mc-stan.org/projpred/dev/reference/predictor_terms.md)
in the order of the variables of an interaction (for example, `x2:x1`
instead of `x1:x2`).

The left-hand side is the response of the formula of the reference model
object (element `formula` of the output of
[`get_refmodel()`](https://mc-stan.org/projpred/dev/reference/refmodel-init-get.md)).
It has only the response variable, so it can differ from the left-hand
side of the reference model fit:

- For a two-column binomial response (`cbind(successes, failures)`), the
  left-hand side is only the first column. The projection still uses the
  numbers of trials (`successes + failures`): projpred keeps them in the
  reference model object (element `wobs`), not in the formula.

- For a brms reference model, the left-hand side is only the response
  variable, without addition terms such as `trials()` or
  [`weights()`](https://rdrr.io/r/stats/weights.html). The projection
  still uses them: projpred keeps the numbers of trials (or the weights)
  in the reference model object (element `wobs`), not in the formula.

- For the latent projection, the left-hand side is the latent response.

Observation weights and offsets are not part of the formula. projpred
fits the submodels without the offsets and adds the offsets of the
reference model to the predictions. If you fit the submodel in another
package, add the offsets there.

The environment of the formula is that of the reference model formula.

With `refit_prj = FALSE` in
[`project()`](https://mc-stan.org/projpred/dev/reference/project.md),
the submodel fits come from the search. After an L1 search, these fits
are the L1-penalized fits of the search path. A fit of this formula in
another package gives other coefficients.

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
print(formula(prj)) # gives `y ~ X1 + X3 + X5`
#> y ~ X1 + X3 + X5
#> <environment: 0x5c35b6d30b90>
```
