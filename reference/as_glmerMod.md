# Convert a submodel fit of class `projpred_mlvl` to class `glmerMod`

For binomial and Poisson multilevel submodels with random intercepts,
projpred uses its own fit, which has class `projpred_mlvl` (see section
"Draw-wise divergence minimizers" of
[projpred-package](https://mc-stan.org/projpred/reference/projpred-package.md)).
This function builds the corresponding lme4 object, for use with the
methods of lme4. It does not fit the submodel again: it evaluates the
objective of [`lme4::glmer()`](https://rdrr.io/pkg/lme4/man/glmer.html)
(with `nAGQ = 1`) once, at the parameter estimates of `object`.

## Usage

``` r
as_glmerMod(object, data, weights = NULL)
```

## Arguments

- object:

  An object of class `projpred_mlvl`, for example an element of element
  `outdmin` of an object returned by
  [`project()`](https://mc-stan.org/projpred/reference/project.md).

- data:

  A `data.frame` with the data that `object` was fitted to: the
  predictors, the grouping variables and the response, with the names
  used in `formula(object)`. The object does not store these data. In
  the traditional projection, the response of a submodel fit is the
  vector of the reference model's predicted means (success probabilities
  for the binomial family) for the projected draw or cluster. In the
  augmented-data projection, the data are the augmented rows (one row
  per observation and response category), and the predicted
  probabilities are in the weights.

- weights:

  A numeric vector of the observation weights that `object` was fitted
  with, or `NULL` for weights of 1. Rows of `data` with weight 0 are
  removed, as projpred does before it fits a multilevel submodel.

## Value

An object of class `glmerMod` (see
[`lme4::glmer()`](https://rdrr.io/pkg/lme4/man/glmer.html)) with the
fixed effects, the group-level standard deviations and the group-level
intercepts of `object`.

## Details

The group-level intercepts of the returned object are computed by lme4
from `data` and `weights`. If they or the value of the objective
function at the parameter estimates differ from those of `object`, then
`data` or `weights` are not those of the fit, and this function throws a
warning.

For non-integer responses, lme4 versions up to 2.1-0 round binomial
responses in the log-likelihood and give a log-likelihood of `-Inf` for
the Poisson family. With these versions, the log-likelihood of the
returned object (and what depends on it, for example
[`AIC()`](https://rdrr.io/r/stats/AIC.html)) is not that of the fit; the
parameter estimates and the predictions are not affected.

The call stored in the returned object is the
[`lme4::glmer()`](https://rdrr.io/pkg/lme4/man/glmer.html) call for
`data` and `weights`, with the default settings of
[`lme4::glmer()`](https://rdrr.io/pkg/lme4/man/glmer.html).

## Examples

``` r
# Data:
set.seed(1)
dat <- data.frame(x1 = rnorm(120), x2 = rnorm(120),
                  g = factor(rep(1:8, each = 15)))
dat$y <- rbinom(120, 1, plogis(0.8 * dat$x1 + rnorm(8)[dat$g]))

# The `stanreg` fit which will be used as the reference model (with small
# values for `chains` and `iter`, but only for technical reasons in this
# example; this is not recommended in general):
fit <- rstanarm::stan_glmer(
  y ~ x1 + x2 + (1 | g), family = binomial(), data = dat,
  chains = 2, iter = 500, refresh = 0, seed = 9876
)
#> Warning: Bulk Effective Samples Size (ESS) is too low, indicating posterior means and medians may be unreliable.
#> Running the chains for more iterations may help. See
#> https://mc-stan.org/misc/warnings.html#bulk-ess
#> Warning: Tail Effective Samples Size (ESS) is too low, indicating posterior variances and tail quantiles may be unreliable.
#> Running the chains for more iterations may help. See
#> https://mc-stan.org/misc/warnings.html#tail-ess

# Projection onto a multilevel submodel, with a single cluster (so the
# response of the submodel fit is the mean of the reference model's
# predicted success probabilities over the posterior draws):
refm <- suppressWarnings(get_refmodel(fit))
prj <- project(refm, predictor_terms = c("x1", "(1 | g)"), nclusters = 1)
fit_sub <- prj$outdmin[[1]]

# Conversion (the name of the response is taken from the formula of the
# submodel fit):
dat_sub <- dat
dat_sub[[all.vars(formula(fit_sub))[1]]] <- rowMeans(refm$mu)
fit_lme4 <- as_glmerMod(fit_sub, data = dat_sub, weights = refm$wobs)
lme4::ranef(fit_lme4)
#> $g
#>   (Intercept)
#> 1  0.20921382
#> 2  0.67702021
#> 3  0.08850761
#> 4 -0.02651275
#> 5  0.45002772
#> 6 -0.49389935
#> 7 -0.26729176
#> 8 -0.59042784
#> 
#> with conditional variances for “g” 
```
