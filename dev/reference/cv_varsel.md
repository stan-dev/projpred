# Run search and performance evaluation with cross-validation

Run the *search* part and the *evaluation* part for a projection
predictive variable selection. The search part determines the predictor
ranking (also known as solution path), i.e., the best submodel for each
submodel size (number of predictor terms). The evaluation part
determines the predictive performance of the submodels along the
predictor ranking. In contrast to
[`varsel()`](https://mc-stan.org/projpred/dev/reference/varsel.md),
`cv_varsel()` performs a cross-validation (CV) by running the search
part with the training data of each CV fold separately (an exception is
explained in section "Note" below) and by running the evaluation part on
the corresponding test set of each CV fold. A special method is
`cv_varsel.vsel()` because it re-uses the search results from an earlier
`cv_varsel()` (or
[`varsel()`](https://mc-stan.org/projpred/dev/reference/varsel.md)) run,
as illustrated in the main vignette.

## Usage

``` r
cv_varsel(object, ...)

# Default S3 method
cv_varsel(object, ...)

# S3 method for class 'vsel'
cv_varsel(
  object,
  cv_method = object$cv_method %||% "LOO",
  nloo = object$nloo,
  K = object$K %||% if (!inherits(object, "datafit")) 5 else 10,
  cvfits = object$cvfits,
  validate_search = object$validate_search %||% TRUE,
  ...
)

# S3 method for class 'refmodel'
cv_varsel(
  object,
  method = "forward",
  cv_method = if (!inherits(object, "datafit")) "LOO" else "kfold",
  ndraws = NULL,
  nclusters = 20,
  ndraws_pred = 400,
  nclusters_pred = NULL,
  refit_prj = !inherits(object, "datafit"),
  nterms_max = NULL,
  penalty = NULL,
  verbose = getOption("projpred.verbose", as.integer(interactive())),
  nloo = if (cv_method == "LOO") object$nobs else NULL,
  K = if (!inherits(object, "datafit")) 5 else 10,
  cvfits = object$cvfits,
  search_control = NULL,
  lambda_min_ratio = 1e-05,
  nlambda = 150,
  thresh = 1e-06,
  validate_search = TRUE,
  seed = NA,
  search_terms = NULL,
  search_out = NULL,
  parallel = getOption("projpred.parallel_cv", FALSE),
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

  For `cv_varsel.default()`: Arguments passed to
  [`get_refmodel()`](https://mc-stan.org/projpred/dev/reference/refmodel-init-get.md)
  as well as to `cv_varsel.refmodel()`. For `cv_varsel.vsel()`:
  Arguments passed to `cv_varsel.refmodel()`. For
  `cv_varsel.refmodel()`: Arguments passed to the divergence minimizer
  (see argument `div_minimizer` of
  [`init_refmodel()`](https://mc-stan.org/projpred/dev/reference/refmodel-init-get.md)
  as well as section "Draw-wise divergence minimizers" of
  [projpred-package](https://mc-stan.org/projpred/dev/reference/projpred-package.md))
  when refitting the submodels for the performance evaluation (if
  `refit_prj` is `TRUE`).

- cv_method:

  The CV method, either `"LOO"` or `"kfold"`. In the `"LOO"` case, a
  Pareto-smoothed importance sampling leave-one-out CV (PSIS-LOO CV) is
  performed, which avoids refitting the reference model `nloo` times (in
  contrast to a standard LOO-CV). In the `"kfold"` case, a \\K\\-fold CV
  is performed. See also section "Note" below.

- nloo:

  Only relevant if `cv_method = "LOO"` and `validate_search = TRUE`. If
  `nloo > 0` is smaller than the number of all observations, full LOO-CV
  (i.e., PSIS-LOO CV with `validate_search = TRUE` and with `nloo = n`
  where `n` denotes the number of all observations) is approximated by
  subsampled LOO-CV, i.e., by combining the fast (i.e.,
  `validate_search = FALSE`) LOO result for the selected models and
  `nloo` leave-one-out searches using the difference estimator with
  simple random sampling (SRS) without replacement (WOR) (Magnusson et
  al., 2020). Smaller `nloo` values lead to faster computation, but
  higher uncertainty in the evaluation part. If `NULL`, all observations
  are used (as by default). Note that performance statistic `"auc"` (see
  argument `stats` of
  [`summary.vsel()`](https://mc-stan.org/projpred/dev/reference/summary.vsel.md)
  and
  [`plot.vsel()`](https://mc-stan.org/projpred/dev/reference/plot.vsel.md))
  is not supported in case of subsampled LOO-CV. Furthermore, option
  `"best"` for argument `baseline` of
  [`summary.vsel()`](https://mc-stan.org/projpred/dev/reference/summary.vsel.md)
  and
  [`plot.vsel()`](https://mc-stan.org/projpred/dev/reference/plot.vsel.md)
  is not supported in case of subsampled LOO-CV.

- K:

  Only relevant if `cv_method = "kfold"` and if `cvfits` is `NULL`
  (which is the case for reference model objects created by
  [`get_refmodel.stanreg()`](https://mc-stan.org/projpred/dev/reference/refmodel-init-get.md)
  or
  [`brms::get_refmodel.brmsfit()`](https://paulbuerkner.com/brms/reference/get_refmodel.brmsfit.html)).
  Number of folds in \\K\\-fold CV.

- cvfits:

  Only relevant if `cv_method = "kfold"`. The same as argument `cvfits`
  of
  [`init_refmodel()`](https://mc-stan.org/projpred/dev/reference/refmodel-init-get.md),
  but repeated here so that output from
  [`run_cvfun()`](https://mc-stan.org/projpred/dev/reference/run_cvfun.md)
  can be inserted here straightforwardly.

- validate_search:

  A single logical value indicating whether to cross-validate also the
  search part, i.e., whether to run the search separately for each CV
  fold (`TRUE`) or not (`FALSE`). With `FALSE`, the computation is
  faster, but the predictive performance estimates of the selected
  submodels are optimistically biased. However, these fast biased
  estimates can be useful to obtain initial information on the
  usefulness of projection predictive variable selection.

- method:

  The method for the search part. Possible options are `"forward"` for
  forward search and `"L1"` for L1 search. See also section "Details"
  below.

- ndraws:

  Number of posterior draws used in the search part. Ignored if
  `nclusters` is not `NULL` or in case of L1 search (because L1 search
  always uses a single cluster). If both (`nclusters` and `ndraws`) are
  `NULL`, the number of posterior draws from the reference model is used
  for `ndraws`. See also section "Details" below.

- nclusters:

  Number of clusters of posterior draws used in the search part. Ignored
  in case of L1 search (because L1 search always uses a single cluster).
  For the meaning of `NULL`, see argument `ndraws`. See also section
  "Details" below.

- ndraws_pred:

  Only relevant if `refit_prj` is `TRUE`. Number of posterior draws used
  in the evaluation part. Ignored if `nclusters_pred` is not `NULL`. If
  both (`nclusters_pred` and `ndraws_pred`) are `NULL`, the number of
  posterior draws from the reference model is used for `ndraws_pred`.
  See also section "Details" below.

- nclusters_pred:

  Only relevant if `refit_prj` is `TRUE`. Number of clusters of
  posterior draws used in the evaluation part. For the meaning of
  `NULL`, see argument `ndraws_pred`. See also section "Details" below.

- refit_prj:

  For the evaluation part, should the projections onto the submodels
  along the predictor ranking be performed again using `ndraws_pred`
  draws or `nclusters_pred` clusters (`TRUE`) or should their
  projections from the search part, which used `ndraws` draws or
  `nclusters` clusters, be re-used (`FALSE`)?

- nterms_max:

  Maximum submodel size (number of predictor terms) up to which the
  search is continued. If `NULL`, then `min(19, D)` is used where `D` is
  the number of terms in the reference model (or in `search_terms`, if
  supplied). Note that `nterms_max` does not count the intercept, so use
  `nterms_max = 0` for the intercept-only model. (Correspondingly, `D`
  above does not count the intercept.)

- penalty:

  Only relevant for L1 search. A numeric vector determining the relative
  penalties or costs for the predictors. A value of `0` means that those
  predictors have no cost and will therefore be selected first, whereas
  `Inf` means those predictors will never be selected. If `NULL`, then
  `1` is used for each predictor.

- verbose:

  A single integer value from the set \\\\0, 1, 2, 3, 4\\\\ (for
  [`varsel()`](https://mc-stan.org/projpred/dev/reference/varsel.md),
  \\3\\ and \\4\\ have the same effect), indicating how much information
  (if any) to print out during the computations. Higher values indicate
  that more information should be printed, `0` deactivates the verbose
  mode. Internally, argument `verbose` is coerced to integer via
  [`as.integer()`](https://rdrr.io/r/base/integer.html), so technically,
  a single logical value or a single numeric value work as well.

- search_control:

  A `list` of "control" arguments (i.e., tuning parameters) for the
  search. In case of forward search, these arguments are passed to the
  divergence minimizer (see argument `div_minimizer` of
  [`init_refmodel()`](https://mc-stan.org/projpred/dev/reference/refmodel-init-get.md)
  as well as section "Draw-wise divergence minimizers" of
  [projpred-package](https://mc-stan.org/projpred/dev/reference/projpred-package.md)).
  In case of forward search, `NULL` causes `...` to be used not only for
  the performance evaluation, but also for the search. In case of L1
  search, possible arguments are:

  - `lambda_min_ratio`: Ratio between the smallest and largest lambda in
    the L1-penalized search (default: `1e-5`). This parameter
    essentially determines how long the search is carried out, i.e., how
    large submodels are explored. No need to change this unless the
    program gives a warning about this.

  - `nlambda`: Number of values in the lambda grid for L1-penalized
    search (default: `150`). No need to change this unless the program
    gives a warning about this.

  - `thresh`: Convergence threshold when computing the L1 path (default:
    `1e-6`). Usually, there is no need to change this.

- lambda_min_ratio:

  Deprecated (please use `search_control` instead). Only relevant for L1
  search. Ratio between the smallest and largest lambda in the
  L1-penalized search. This parameter essentially determines how long
  the search is carried out, i.e., how large submodels are explored. No
  need to change this unless the program gives a warning about this.

- nlambda:

  Deprecated (please use `search_control` instead). Only relevant for L1
  search. Number of values in the lambda grid for L1-penalized search.
  No need to change this unless the program gives a warning about this.

- thresh:

  Deprecated (please use `search_control` instead). Only relevant for L1
  search. Convergence threshold when computing the L1 path. Usually,
  there is no need to change this.

- seed:

  Pseudorandom number generation (PRNG) seed by which the same results
  can be obtained again if needed. Passed to argument `seed` of
  [`set.seed()`](https://rdrr.io/r/base/Random.html), but can also be
  `NA` to not call [`set.seed()`](https://rdrr.io/r/base/Random.html) at
  all. If not `NA`, then the PRNG state is reset (to the state before
  calling `cv_varsel()`) upon exiting `cv_varsel()`. Here, `seed` is
  used for clustering the reference model's posterior draws (if
  `!is.null(nclusters)` or `!is.null(nclusters_pred)`), for subsampling
  PSIS-LOO CV folds (if `nloo` is smaller than the number of
  observations), for sampling the folds in \\K\\-fold CV, and for
  drawing new group-level effects when predicting from a multilevel
  submodel (however, not yet in case of a GAMM).

- search_terms:

  Only relevant for forward search. A custom character vector of
  predictor term blocks to consider for the search. Section "Details"
  below describes more precisely what "predictor term block" means. The
  intercept (`"1"`) is always included internally via
  [`union()`](https://rdrr.io/r/base/sets.html), so there's no
  difference between including it explicitly or omitting it. The default
  `search_terms` considers all the terms in the reference model's
  formula.

- search_out:

  Intended for internal use.

- parallel:

  A single logical value indicating whether to run costly parts of the
  CV in parallel (`TRUE`) or not (`FALSE`). See also section "Note"
  below as well as section "Parallelization" in
  [projpred-package](https://mc-stan.org/projpred/dev/reference/projpred-package.md).

## Value

An object of class `vsel`. The elements of this object are not meant to
be accessed directly but instead via helper functions (see the main
vignette and
[projpred-package](https://mc-stan.org/projpred/dev/reference/projpred-package.md)).

## Details

Arguments `ndraws`, `nclusters`, `nclusters_pred`, and `ndraws_pred` are
automatically truncated at the number of posterior draws in the
reference model (which is `1` for `datafit`s). Using less draws or
clusters in `ndraws`, `nclusters`, `nclusters_pred`, or `ndraws_pred`
than posterior draws in the reference model may result in slightly
inaccurate projection performance. Increasing these arguments affects
the computation time linearly.

For argument `method`, there are some restrictions: For a reference
model with multilevel or additive formula terms or a reference model set
up for the augmented-data projection, only the forward search is
available. Furthermore, argument `search_terms` requires a forward
search to take effect.

L1 search is faster than forward search, but forward search may be more
accurate. Furthermore, forward search may find a sparser model with
comparable performance to that found by L1 search, but it may also
overfit when more predictors are added. This overfit can be detected by
running search validation (see `cv_varsel()`).

An L1 search may select an interaction term before all involved
lower-order interaction terms (including main-effect terms) have been
selected. In projpred versions \> 2.6.0, the resulting predictor ranking
is automatically modified so that the lower-order interaction terms come
before this interaction term, but if this is conceptually undesired,
choose the forward search instead.

The elements of the `search_terms` character vector don't need to be
individual predictor terms. Instead, they can be building blocks
consisting of several predictor terms connected by the `+` symbol. To
understand how these building blocks work, it is important to know how
projpred's forward search works: It starts with an empty vector `chosen`
which will later contain already selected predictor terms. Then, the
search iterates over model sizes \\j \in \\0, ..., J\\\\ (with \\J\\
denoting the maximum submodel size, not counting the intercept). The
candidate models at model size \\j\\ are constructed from those elements
from `search_terms` which yield model size \\j\\ when combined with the
`chosen` predictor terms. Note that sometimes, there may be no candidate
models for model size \\j\\. Also note that internally, `search_terms`
is expanded to include the intercept (`"1"`), so the first step of the
search (model size 0) always consists of the intercept-only model as the
only candidate.

As a `search_terms` example, consider a reference model with formula
`y ~ x1 + x2 + x3`. Then, to ensure that `x1` is always included in the
candidate models, specify
`search_terms = c("x1", "x1 + x2", "x1 + x3", "x1 + x2 + x3")` (or, in a
simpler way that leads to the same results,
`search_terms = c("x1", "x1 + x2", "x1 + x3")`, for which helper
function
[`force_search_terms()`](https://mc-stan.org/projpred/dev/reference/force_search_terms.md)
exists). This search would start with `y ~ 1` as the only candidate at
model size 0. At model size 1, `y ~ x1` would be the only candidate. At
model size 2, `y ~ x1 + x2` and `y ~ x1 + x3` would be the two
candidates. At the last model size of 3, `y ~ x1 + x2 + x3` would be the
only candidate. As another example, to exclude `x1` from the search,
specify `search_terms = c("x2", "x3", "x2 + x3")` (or, in a simpler way
that leads to the same results, `search_terms = c("x2", "x3")`).

## Note

If `validate_search` is `FALSE`, the search is not included in the CV so
that only a single full-data search is run. If the number of
observations is large, the fast PSIS-LOO CV along the full-data search
path is likely to be accurate. If the number of observations is small or
moderate, the fast PSIS-LOO CV along the full-data search path is likely
to have optimistic bias in the middle of the search path. This result
can be used to guide further actions and the optimistic bias can be
greatly reduced by using `validate_search = TRUE`.

PSIS uses the Pareto-\\\hat{k}\\ diagnostic to assess the reliability of
PSIS-LOO CV. Global option `projpred.warn_psis` (default `TRUE`)
controls whether the Pareto-\\\hat{k}\\ diagnostics may result in
warnings. See
[loo::loo-glossary](https://mc-stan.org/loo/reference/loo-glossary.html)
for how to interpret the Pareto-\\\hat{k}\\ values and the warning
thresholds. projpred does not support the usually recommended
moment-matching (see
[`loo::loo_moment_match()`](https://mc-stan.org/loo/reference/loo_moment_match.html)
and
[`brms::loo_moment_match()`](https://paulbuerkner.com/brms/reference/loo_moment_match.brmsfit.html)),
mixture importance sampling
([`vignette("loo2-mixis", package="loo")`](https://mc-stan.org/loo/articles/loo2-mixis.html)),
or `reloo`-ing
([`brms::reloo()`](https://paulbuerkner.com/brms/reference/reloo.brmsfit.html)).
If the reference model PSIS-LOO CV Pareto-\\\hat{k}\\ values are good,
but there are high Pareto-\\\hat{k}\\ values for the projected models,
you can try increasing the number of draws used for the PSIS-LOO CV
(`ndraws` in case of `refit_prj = FALSE`; `ndraws_pred` in case of
`refit_prj = TRUE`). If increasing the number of draws does not help and
if the reference model PSIS-LOO CV Pareto-\\\hat{k}\\ values are high,
and the reference model PSIS-LOO CV results change substantially when
using moment-matching, mixture importance sampling, or `reloo`-ing, we
recommend to use \\K\\-fold CV within `projpred`.

For PSIS-LOO CV, projpred calls
[`loo::psis()`](https://mc-stan.org/loo/reference/psis.html) (or,
exceptionally,
[`loo::sis()`](https://mc-stan.org/loo/reference/sis.html), see below)
with `r_eff = NA`. This is only a problem if there was extreme
autocorrelation between the MCMC iterations when the reference model was
built. In those cases however, the reference model should not have been
used anyway, so we don't expect projpred's `r_eff = NA` to be a problem.

PSIS cannot be used if the number of draws or clusters is too small. In
such cases, projpred resorts to standard importance sampling (SIS) and
shows a message about this. Throughout the documentation, the term
"PSIS" is used even though in fact, projpred resorts to SIS in these
special cases. If SIS is used, check that the reference model PSIS-LOO
CV Pareto-\\\hat{k}\\ values are good.

With `parallel = TRUE`, costly parts of projpred's CV can be run in
parallel. Costly parts are the fold-wise searches and performance
evaluations in case of `validate_search = TRUE`. (Note that in case of
\\K\\-fold CV, the \\K\\ reference model refits are not affected by
argument `parallel`; only projpred's CV is affected.) The
parallelization is powered by the foreach package. Thus, any parallel
(or sequential) backend compatible with foreach can be used, e.g., the
backends from packages doParallel, doMPI, or doFuture. For GLMs, this CV
parallelization should work reliably, but for other models (such as
GLMMs), it may lead to excessive memory usage which in turn may crash
the R session (on Unix systems, setting an appropriate memory limit via
[`unix::rlimit_as()`](https://jeroen.r-universe.dev/unix/reference/rlimit.html)
may avoid crashing the whole machine). However, the problem of excessive
memory usage is less pronounced for the CV parallelization than for the
projection parallelization described in
[projpred-package](https://mc-stan.org/projpred/dev/reference/projpred-package.md).
In that regard, the CV parallelization is recommended over the
projection parallelization.

## References

Magnusson, Måns, Michael Riis Andersen, Johan Jonasson, Aki Vehtari.
2020. "Leave-One-Out Cross-Validation for Bayesian Model Comparison in
Large Data." In *Proceedings of the 23rd International Conference on
Artificial Intelligence and Statistics*, edited by Silvia Chiappa and
Roberto Calandra, 108:341–351. Proceedings of Machine Learning Research.
PMLR. <https://proceedings.mlr.press/v108/magnusson20a.html>.

McLatchie, Yann, Sölvi Rögnvaldsson, Frank Weber, and Aki Vehtari. 2025.
"Advances in Projection Predictive Inference." *Statistical Science*, 40
(1):128–147. [doi:10.1214/24-STS949](https://doi.org/10.1214/24-STS949)
.

Piironen, Juho, Markus Paasiniemi, and Aki Vehtari. 2020. "Projective
Inference in High-Dimensional Problems: Prediction and Feature
Selection." *Electronic Journal of Statistics*, 14 (1):2155–2197.
[doi:10.1214/20-EJS1711](https://doi.org/10.1214/20-EJS1711) .

Vehtari, Aki, Andrew Gelman, and Jonah Gabry. 2017. "Practical Bayesian
Model Evaluation Using Leave-One-Out Cross-Validation and WAIC."
*Statistics and Computing*, 27 (5):1413–32.
[doi:10.1007/s11222-016-9696-4](https://doi.org/10.1007/s11222-016-9696-4)
.

Vehtari, Aki, Daniel Simpson, Andrew Gelman, Yuling Yao, and Jonah
Gabry. 2024. "Pareto Smoothed Importance Sampling." *Journal of Machine
Learning Research*, 25 (72):1–58.
<https://jmlr.org/papers/v25/19-556.html>.

## See also

[`varsel()`](https://mc-stan.org/projpred/dev/reference/varsel.md)

## Examples

``` r
# Data:
dat_gauss <- data.frame(y = df_gaussian$y, df_gaussian$x)

# The `stanreg` fit which will be used as the reference model (with small
# values for `chains` and `iter`, but only for technical reasons in this
# example; this is not recommended in general):
fit <- rstanarm::stan_glm(
  y ~ X1 + X2 + X3 + X4 + X5, family = gaussian(), data = dat_gauss,
  QR = TRUE, chains = 2, iter = 1000, refresh = 0, seed = 9876
)

# Run cv_varsel() (with L1 search and small values for `K`, `nterms_max`, and
# `nclusters_pred`, but only for the sake of speed in this example; this is
# not recommended in general):
cvvs <- cv_varsel(fit, method = "L1", cv_method = "kfold", K = 2,
                  nterms_max = 3, nclusters_pred = 10, seed = 5555)
#> Fitting model 1 out of 2
#> Fitting model 2 out of 2
# Now see, for example, `?print.vsel`, `?plot.vsel`, `?suggest_size.vsel`,
# and `?ranking` for possible post-processing functions.
```
