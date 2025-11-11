# Create `cvfits` from `cvfun`

A helper function that can be used to create input for
[`cv_varsel.refmodel()`](https://mc-stan.org/projpred/dev/reference/cv_varsel.md)'s
argument `cvfits` by running first
[`cv_folds()`](https://mc-stan.org/projpred/dev/reference/cv-indices.md)
and then the reference model object's `cvfun` (see
[`init_refmodel()`](https://mc-stan.org/projpred/dev/reference/refmodel-init-get.md)).
This is helpful if \\K\\-fold CV is run multiple times based on the same
\\K\\ reference model refits.

## Usage

``` r
run_cvfun(object, ...)

# Default S3 method
run_cvfun(object, ...)

# S3 method for class 'refmodel'
run_cvfun(
  object,
  K = if (!inherits(object, "datafit")) 5 else 10,
  folds = NULL,
  seed = NA,
  ...
)
```

## Arguments

- object:

  An object of class `refmodel` (returned by
  [`get_refmodel()`](https://mc-stan.org/projpred/dev/reference/refmodel-init-get.md)
  or
  [`init_refmodel()`](https://mc-stan.org/projpred/dev/reference/refmodel-init-get.md))
  or an object that can be passed to argument `object` of
  [`get_refmodel()`](https://mc-stan.org/projpred/dev/reference/refmodel-init-get.md).

- ...:

  For `run_cvfun.default()`: Arguments passed to
  [`get_refmodel()`](https://mc-stan.org/projpred/dev/reference/refmodel-init-get.md).
  For `run_cvfun.refmodel()`: Currently ignored.

- K:

  Number of folds. Must be at least 2 and not exceed the number of
  observations. Ignored if `folds` is not `NULL`.

- folds:

  Either `NULL` for determining the CV folds automatically via
  [`cv_folds()`](https://mc-stan.org/projpred/dev/reference/cv-indices.md)
  (using argument `K`) or a numeric (in fact, integer) vector giving the
  fold index for each observation. In the latter case, argument `K` is
  ignored.

- seed:

  Pseudorandom number generation (PRNG) seed by which the same results
  can be obtained again if needed. Passed to argument `seed` of
  [`set.seed()`](https://rdrr.io/r/base/Random.html), but can also be
  `NA` to not call [`set.seed()`](https://rdrr.io/r/base/Random.html) at
  all. If not `NA`, then the PRNG state is reset (to the state before
  calling `run_cvfun()`) upon exiting `run_cvfun()`.

## Value

An object that can be used as input for
[`cv_varsel.refmodel()`](https://mc-stan.org/projpred/dev/reference/cv_varsel.md)'s
argument `cvfits`.

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

# Define the reference model object explicitly (not really necessary here
# because the get_refmodel() call is quite fast in this example, but in
# general, this approach is faster than defining the reference model object
# multiple times implicitly):
ref <- get_refmodel(fit)

# Run the reference model object's `cvfun` (with a small value for `K`, but
# only for the sake of speed in this example; this is not recommended in
# general):
cv_fits <- run_cvfun(ref, K = 2, seed = 184)
#> Fitting model 1 out of 2
#> Warning: Bulk Effective Samples Size (ESS) is too low, indicating posterior means and medians may be unreliable.
#> Running the chains for more iterations may help. See
#> https://mc-stan.org/misc/warnings.html#bulk-ess
#> Fitting model 2 out of 2
#> Warning: Bulk Effective Samples Size (ESS) is too low, indicating posterior means and medians may be unreliable.
#> Running the chains for more iterations may help. See
#> https://mc-stan.org/misc/warnings.html#bulk-ess

# Run cv_varsel() (with L1 search and small values for `nterms_max` and
# `nclusters_pred`, but only for the sake of speed in this example; this is
# not recommended in general) and use `cv_fits` there:
cvvs_L1 <- cv_varsel(ref, method = "L1", cv_method = "kfold",
                     cvfits = cv_fits, nterms_max = 3, nclusters_pred = 10,
                     seed = 5555)
# Now see, for example, `?print.vsel`, `?plot.vsel`, `?suggest_size.vsel`,
# and `?ranking` for possible post-processing functions.

# The purpose of run_cvfun() is to create an object that can be used in
# multiple cv_varsel() calls, e.g., to check the sensitivity to the search
# method (L1 or forward):
cvvs_fw <- cv_varsel(ref, method = "forward", cv_method = "kfold",
                     cvfits = cv_fits, nterms_max = 3, nclusters = 5,
                     nclusters_pred = 10, seed = 5555)

# Stratified K-fold CV is straightforward:
n_strat <- 3L
set.seed(692)
# Some example strata:
strat_fac <- sample(paste0("lvl", seq_len(n_strat)), size = nrow(dat_gauss),
                    replace = TRUE,
                    prob = diff(c(0, pnorm(seq_len(n_strat - 1L) - 0.5), 1)))
table(strat_fac)
#> strat_fac
#> lvl1 lvl2 lvl3 
#>   70   24    6 
# Use loo::kfold_split_stratified() to create the folds vector:
folds_strat <- loo::kfold_split_stratified(K = 2, x = strat_fac)
table(folds_strat, strat_fac)
#>            strat_fac
#> folds_strat lvl1 lvl2 lvl3
#>           1   35   12    3
#>           2   35   12    3
# Call run_cvfun(), but this time with argument `folds` instead of `K` (here,
# specifying argument `seed` would not be necessary because of the set.seed()
# call above, but we specify it nonetheless for the sake of generality):
cv_fits_strat <- run_cvfun(ref, folds = folds_strat, seed = 391)
#> Fitting model 1 out of 2
#> Fitting model 2 out of 2
# Now use `cv_fits_strat` analogously to `cv_fits` from above.
```
