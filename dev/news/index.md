# Changelog

## projpred 2.10.0

CRAN release: 2025-12-06

### Major changes

- Added support for censored observations when using the latent
  projection (with response-scale analyses). This makes it possible,
  e.g., to use the latent projection for time-to-event models (also
  known as models for survival analysis). Due to this new feature, the
  function passed to argument `latent_ll_oscale` of
  [`extend_family()`](https://mc-stan.org/projpred/dev/reference/extend_family.md)
  now needs to have an argument `cens`. See
  [`?extend_family`](https://mc-stan.org/projpred/dev/reference/extend_family.md)
  (section “Latent projection”) as well as the latent-projection
  vignette (section [“Censored observations (survival
  analysis)”](https://mc-stan.org/projpred/articles/latent.html#cens))
  for more information and examples. Note that only the performance
  statistics `"elpd"`, `"mlpd"`, and `"gmpd"` take censoring into
  account (on response scale). (GitHub:
  [\#528](https://github.com/stan-dev/projpred/issues/528))

### Minor changes

- [`suggest_size()`](https://mc-stan.org/projpred/dev/reference/suggest_size.md)
  returns `nterms_max` instead of NA when no model satisfies the
  criterion
- Use `reformulas` package directly for formula processing functions.
  This package is already a dependency of `lme4` \>= v1.1-38.

### Bug fixes

## projpred 2.9.1

CRAN release: 2025-10-28

### Major changes

- Fixed a major bug in forward search with 3-way or higher-order
  interactions. See “Bug fixes” below for details.

### Bug fixes

- Fixed a bug that caused forward search not to place lower-order terms
  of a 3-way or higher-order interaction term before that interaction
  term. (GitHub:
  [\#531](https://github.com/stan-dev/projpred/issues/531),
  [\#532](https://github.com/stan-dev/projpred/issues/532))
- Fixed a bug that caused
  [`project()`](https://mc-stan.org/projpred/dev/reference/project.md)
  to construct an incorrect internal table of predictor terms (from
  which argument `predictor_terms` can select one or several terms from)
  in case of a group-level term that does not contain a group-level
  intercept. (GitHub:
  [\#533](https://github.com/stan-dev/projpred/issues/533),
  [\#534](https://github.com/stan-dev/projpred/issues/534))
- Relaxed test to prevent spurious failures and unblock reverse
  dependencies. No functional changes.

## projpred 2.9.0

CRAN release: 2025-07-08

### Major changes

- Subsampled PSIS-LOO CV (usable via argument `nloo` of
  [`cv_varsel()`](https://mc-stan.org/projpred/dev/reference/cv_varsel.md))
  has been fixed and is not experimental anymore. There are a few
  restrictions: Performance statistic `"auc"` (see argument `stats` of
  [`summary.vsel()`](https://mc-stan.org/projpred/dev/reference/summary.vsel.md)
  and
  [`plot.vsel()`](https://mc-stan.org/projpred/dev/reference/plot.vsel.md);
  argument `stat` of
  [`suggest_size()`](https://mc-stan.org/projpred/dev/reference/suggest_size.md)
  is concerned as well) is not supported in case of subsampled
  PSIS-LOO CV. Furthermore, `baseline = "best"` (in
  [`summary.vsel()`](https://mc-stan.org/projpred/dev/reference/summary.vsel.md)
  and
  [`plot.vsel()`](https://mc-stan.org/projpred/dev/reference/plot.vsel.md))
  is not supported in case of subsampled PSIS-LOO CV either. (GitHub:
  [\#94](https://github.com/stan-dev/projpred/issues/94),
  [\#496](https://github.com/stan-dev/projpred/issues/496))

- The uncertainty interval for performance statistic `"mse"` is now
  based on a log-normal approximation (instead of a normal
  approximation) if argument `deltas` of
  [`summary.vsel()`](https://mc-stan.org/projpred/dev/reference/summary.vsel.md)
  or
  [`plot.vsel()`](https://mc-stan.org/projpred/dev/reference/plot.vsel.md)
  is `FALSE`. (GitHub:
  [\#496](https://github.com/stan-dev/projpred/issues/496))

- The standard error for performance statistic `"rmse"` is now computed
  via the delta method (instead of bootstrapping). The uncertainty
  interval for `"rmse"` is now based on a log-normal approximation
  (instead of bootstrapping) if argument `deltas` of
  [`summary.vsel()`](https://mc-stan.org/projpred/dev/reference/summary.vsel.md)
  or
  [`plot.vsel()`](https://mc-stan.org/projpred/dev/reference/plot.vsel.md)
  is `FALSE` and based on a normal approximation (instead of
  bootstrapping) if `deltas` is `TRUE`. (GitHub:
  [\#496](https://github.com/stan-dev/projpred/issues/496))

- Performance statistic `"R2"` (R-squared) has been added, see argument
  `stats` of
  [`summary.vsel()`](https://mc-stan.org/projpred/dev/reference/summary.vsel.md)
  and
  [`plot.vsel()`](https://mc-stan.org/projpred/dev/reference/plot.vsel.md);
  argument `stat` of
  [`suggest_size()`](https://mc-stan.org/projpred/dev/reference/suggest_size.md)
  supports it as well. (GitHub:
  [\#483](https://github.com/stan-dev/projpred/issues/483),
  [\#496](https://github.com/stan-dev/projpred/issues/496))

- The performance evaluation part of
  [`cv_varsel()`](https://mc-stan.org/projpred/dev/reference/cv_varsel.md)
  with `cv_method = "LOO"` and `validate_search = FALSE` now always
  applies Pareto smoothing when computing the importance sampling
  weights (as long as the number of importance ratios in the tail is
  large enough; otherwise, no Pareto smoothing is applied). Previously,
  in case of projected draws with nonconstant weights (i.e., in case of
  clustering), no Pareto smoothing had been applied. (GitHub:
  [\#496](https://github.com/stan-dev/projpred/issues/496),
  [\#507](https://github.com/stan-dev/projpred/issues/507))

- The threshold for high Pareto-\\\hat{k}\\ values was updated to the
  one presented by Vehtari et al. (2024, “Pareto smoothed importance
  sampling”, *Journal of Machine Learning Research*, 25(72):1-58,
  <https://www.jmlr.org/papers/v25/19-556.html>). This threshold depends
  on the Monte Carlo sample size and is often close to the former fixed
  threshold of 0.7 (a short introduction may also be found in the [LOO
  glossary](https://mc-stan.org/loo/reference/loo-glossary.html)).
  Correspondingly, the former “secondary” threshold of 0.5 is not used
  anymore either. (GitHub:
  [\#490](https://github.com/stan-dev/projpred/issues/490),
  [\#498](https://github.com/stan-dev/projpred/issues/498))

- Argument `type` of
  [`summary.vsel()`](https://mc-stan.org/projpred/dev/reference/summary.vsel.md)
  has gained options `"diff.lower"` and `"diff.upper"` (see the
  documentation for details). (GitHub:
  [\#511](https://github.com/stan-dev/projpred/issues/511))

- Argument `deltas` of
  [`plot.vsel()`](https://mc-stan.org/projpred/dev/reference/plot.vsel.md)
  has gained option `"mixed"` which combines the point estimates from
  `deltas = FALSE` with the uncertainty bars from `deltas = TRUE`.
  (GitHub: [\#511](https://github.com/stan-dev/projpred/issues/511))

- For the latent projection, the function passed to argument
  `latent_ll_oscale` of
  [`extend_family()`](https://mc-stan.org/projpred/dev/reference/extend_family.md)
  now needs to have an argument `dis` (at the second position).
  Similarly, the function passed to argument `latent_ppd_oscale` of
  [`extend_family()`](https://mc-stan.org/projpred/dev/reference/extend_family.md)
  now needs to have an argument `dis_resamp` (again at the second
  position). This makes it possible, e.g., to use the latent projection
  for a log-normal response family. (GitHub:
  [\#513](https://github.com/stan-dev/projpred/issues/513))

- Argument `verbose` of
  [`project()`](https://mc-stan.org/projpred/dev/reference/project.md),
  [`varsel()`](https://mc-stan.org/projpred/dev/reference/varsel.md),
  and
  [`cv_varsel()`](https://mc-stan.org/projpred/dev/reference/cv_varsel.md)
  has been changed from logical to integer. However, logical values
  continue to work (since
  [`as.integer()`](https://rdrr.io/r/base/integer.html) is applied
  internally). Global options `projpred.extra_verbose` and
  `projpred.verbose_project` are now deprecated because additional
  verbosity can be achieved via higher integer values for argument
  `verbose`. The new global option `projpred.verbose` may be used to set
  argument `verbose` of
  [`project()`](https://mc-stan.org/projpred/dev/reference/project.md),
  [`varsel()`](https://mc-stan.org/projpred/dev/reference/varsel.md),
  and
  [`cv_varsel()`](https://mc-stan.org/projpred/dev/reference/cv_varsel.md)
  globally. (GitHub:
  [\#519](https://github.com/stan-dev/projpred/issues/519))

- Some global options have been renamed, so please use their new names
  from now on (although the old names will continue to work for a while)
  (GitHub: [\#500](https://github.com/stan-dev/projpred/issues/500),
  [\#521](https://github.com/stan-dev/projpred/issues/521)):

  - Global option `projpred.prll_cv` has been renamed to
    `projpred.parallel_cv`.
  - Global option `projpred.warn_prj_drawwise` has been renamed to
    `projpred.warn_proj_drawwise`.
  - Global option `projpred.check_conv` has been renamed to
    `projpred.check_convergence`.
  - Global option `projpred.prll_prj_trigger` has been renamed to
    `projpred.parallel_proj_trigger`.

- In
  [`plot.vsel()`](https://mc-stan.org/projpred/dev/reference/plot.vsel.md),
  several defaults have been changed (GitHub:
  [\#517](https://github.com/stan-dev/projpred/issues/517),
  [\#522](https://github.com/stan-dev/projpred/issues/522)):

  - Argument `text_angle` now defaults to `45` (previously, the default
    was `NULL`).
  - Argument `size_position` now defaults to `"primary_x_top"`
    (previously, the default was `"primary_x_bottom"`).
  - Argument `show_cv_proportions` now defaults to `FALSE` (previously,
    the default was `TRUE`).

  These arguments can now also be controlled via global options, see
  section “Usage” of
  [`?plot.vsel`](https://mc-stan.org/projpred/dev/reference/plot.vsel.md)
  for their names and the main vignette for an illustration.

- The changelog for version 2.6.0 did not contain a notification that
  [`cvfolds()`](https://mc-stan.org/projpred/dev/reference/cv-indices.md)
  had been deprecated and that the new name
  [`cv_folds()`](https://mc-stan.org/projpred/dev/reference/cv-indices.md)
  should be used instead. This changelog entry has been added now (see
  below), but is also mentioned here to make users aware of it (although
  a deprecation warning was already added in version 2.6.0 and will be
  kept until
  [`cvfolds()`](https://mc-stan.org/projpred/dev/reference/cv-indices.md)
  is eventually removed in a future release). (GitHub:
  [\#411](https://github.com/stan-dev/projpred/issues/411))

### Minor changes

- When using the **doFuture** backend for parallelization, progression
  updates can now be received via the **progressr** package, see
  [`` ?`projpred-package` ``](https://mc-stan.org/projpred/dev/reference/projpred-package.md)
  (section “Parallelization”). (GitHub:
  [\#504](https://github.com/stan-dev/projpred/issues/504))
- Several enhancements concerning verbosity, e.g., the number of
  projected draws (resulting from clustering or thinning) is now printed
  out during the different steps of the computations and verbose-mode
  output is redirected to
  [`stderr()`](https://rdrr.io/r/base/showConnections.html) instead of
  [`stdout()`](https://rdrr.io/r/base/showConnections.html). (GitHub:
  [\#506](https://github.com/stan-dev/projpred/issues/506),
  [\#518](https://github.com/stan-dev/projpred/issues/518))
- For the CV parallelization (see argument `parallel` of
  [`cv_varsel()`](https://mc-stan.org/projpred/dev/reference/cv_varsel.md)),
  a new global option `projpred.export_to_workers` may be set to a
  character vector of names of objects to export from the global
  environment to the parallel workers. (GitHub:
  [\#497](https://github.com/stan-dev/projpred/issues/497),
  [\#510](https://github.com/stan-dev/projpred/issues/510))
- Added global options `projpred.foreach_errorhandling` and
  `projpred.foreach_verbose` whose values are passed to
  [`foreach::foreach()`](https://rdrr.io/pkg/foreach/man/foreach.html)’s
  arguments `.errorhandling` and `.verbose`, respectively. The defaults
  for these new global options are the same as those for the respective
  [`foreach::foreach()`](https://rdrr.io/pkg/foreach/man/foreach.html)
  arguments: `"stop"` for global option `projpred.foreach_errorhandling`
  and `FALSE` for global option `projpred.foreach_verbose`. (GitHub:
  commit 3231d13)
- Added global options to control several arguments of
  [`plot.vsel()`](https://mc-stan.org/projpred/dev/reference/plot.vsel.md)
  and
  [`plot.cv_proportions()`](https://mc-stan.org/projpred/dev/reference/plot.cv_proportions.md)
  (see section “Usage” of the help pages of these two functions).
  (GitHub: commit 3333043)
- Changed the maintainer to Osvaldo Martin.

### Bug fixes

- Fixed a bug that caused an error when using the augmented-data or
  latent projection in combination with a single projected draw for
  performance evaluation in
  [`cv_varsel()`](https://mc-stan.org/projpred/dev/reference/cv_varsel.md)
  with `cv_method = "LOO"` and `validate_search = FALSE`. (GitHub:
  [\#512](https://github.com/stan-dev/projpred/issues/512))
- Previously, in case of PSIS-LOO CV with `validate_search = TRUE` and
  thinned posterior draws for projection (i.e., argument(s) `ndraws` or
  `ndraws_pred` being used, not `nclusters` or `nclusters_pred`),
  [`print.vselsummary()`](https://mc-stan.org/projpred/dev/reference/print.vselsummary.md)
  incorrectly reported that the posterior draws had been clustered. This
  has now been fixed, so thinning is reported in such cases. (GitHub:
  [\#516](https://github.com/stan-dev/projpred/issues/516))
- Fixed the internal default `extract_model_data` function when using
  the latent projection for a custom reference model object. (GitHub:
  [\#523](https://github.com/stan-dev/projpred/issues/523))

## projpred 2.8.0

CRAN release: 2023-12-14

### Major changes

- Search results generated in an earlier
  [`varsel()`](https://mc-stan.org/projpred/dev/reference/varsel.md) or
  [`cv_varsel()`](https://mc-stan.org/projpred/dev/reference/cv_varsel.md)
  call can now be re-used by the help of the new
  [`varsel.vsel()`](https://mc-stan.org/projpred/dev/reference/varsel.md)
  and
  [`cv_varsel.vsel()`](https://mc-stan.org/projpred/dev/reference/cv_varsel.md)
  methods (i.e., by applying
  [`varsel()`](https://mc-stan.org/projpred/dev/reference/varsel.md) or
  [`cv_varsel()`](https://mc-stan.org/projpred/dev/reference/cv_varsel.md)
  to the output of the earlier
  [`varsel()`](https://mc-stan.org/projpred/dev/reference/varsel.md) or
  [`cv_varsel()`](https://mc-stan.org/projpred/dev/reference/cv_varsel.md)
  call). This can save a lot of time when re-running the predictive
  performance evaluation part multiple times based on the same search
  results. An illustration may be found in the updated main vignette
  (section [“Preliminary `cv_varsel()`
  run”](https://mc-stan.org/projpred/articles/projpred.html#preliminary-cv_varsel-run);
  a more general description may also be found in section
  [“Speed”](https://mc-stan.org/projpred/articles/projpred.html#speed)).
  (GitHub: [\#461](https://github.com/stan-dev/projpred/issues/461),
  [\#463](https://github.com/stan-dev/projpred/issues/463),
  [\#465](https://github.com/stan-dev/projpred/issues/465),
  [\#466](https://github.com/stan-dev/projpred/issues/466))
- K-fold CV can now be combined with `validate_search = FALSE`. Related
  to this is an internal change which may cause subsampled PSIS-LOO CV
  (an experimental feature controlled by argument `nloo` of
  [`cv_varsel()`](https://mc-stan.org/projpred/dev/reference/cv_varsel.md))
  with clustered projection during the search (i.e.,
  `1 < nclusters && nclusters < S`, where `S` denotes the number of
  posterior draws in the reference model) to yield slightly different
  results due to different internal pseudorandom number generator (PRNG)
  states. Furthermore, if `is.na(seed)`, then the PRNG state for code
  downstream of such a
  [`cv_varsel()`](https://mc-stan.org/projpred/dev/reference/cv_varsel.md)
  call will be different due to this internal change. (GitHub:
  [\#464](https://github.com/stan-dev/projpred/issues/464))
- [`print.vselsummary()`](https://mc-stan.org/projpred/dev/reference/print.vselsummary.md)
  (and hence also
  [`print.vsel()`](https://mc-stan.org/projpred/dev/reference/print.vsel.md))
  now prints the reference model’s performance evaluation results as
  well (not just those of the submodels). Correspondingly, a new helper
  function
  [`performances()`](https://mc-stan.org/projpred/dev/reference/performances.md)
  has been added which allows to access the reference model’s (as well
  as the submodels’) performance evaluation results. (GitHub:
  [\#471](https://github.com/stan-dev/projpred/issues/471))
- Argument `solution_terms` of
  [`project()`](https://mc-stan.org/projpred/dev/reference/project.md)
  has been deprecated. Please use the new argument `predictor_terms`
  instead. (GitHub:
  [\#472](https://github.com/stan-dev/projpred/issues/472))
- For expert users of the augmented-data projection only: Objects of
  class `augmat` or `augvec` do not need to have an attribute called
  `nobs_orig` anymore, but a new attribute called `ndiscrete`, giving
  the number of (possibly latent) response categories instead of the
  number of observations (see
  [`` ?`augdat-internals` ``](https://mc-stan.org/projpred/dev/reference/augdat-internals.md)).
  This simplifies the subsetting of such objects. (GitHub:
  [\#473](https://github.com/stan-dev/projpred/issues/473))
- By default, **projpred** now catches messages and warnings from the
  draw-wise divergence minimizers and throws their unique collection
  after performing all draw-wise divergence minimizations (i.e.,
  draw-wise projections). This can be deactivated by setting global
  option `projpred.warn_prj_drawwise` to `FALSE`. Previously,
  **projpred** suppressed such messages and warnings. (GitHub:
  [\#478](https://github.com/stan-dev/projpred/issues/478))
- By default, **projpred** now checks the convergence of the draw-wise
  divergence minimizers and throws a warning in case of potential
  convergence problems. This can be deactivated by setting global option
  `projpred.check_conv` to `FALSE`. (GitHub:
  [\#478](https://github.com/stan-dev/projpred/issues/478))

### Minor changes

- In
  [`as.matrix.projection()`](https://mc-stan.org/projpred/dev/reference/as.matrix.projection.md),
  `nm_scheme = "auto"` is deprecated. Please use `nm_scheme = NULL`
  instead.
- The plot produced by
  [`plot.vsel()`](https://mc-stan.org/projpred/dev/reference/plot.vsel.md)
  now includes a title and a subtitle, with the subtitle mentioning the
  nominal coverage as well as the type of the confidence intervals (CIs)
  explicitly. However, in case of a facetted plot (i.e., in case of
  multiple `stats`) and some `stats` implying a different CI type than
  other `stats`, the CI types are omitted (because mentioning them would
  make the subtitle too complicated). Note that title and subtitle can
  always be omitted with
  `<plot.vsel() output object> + ggplot2::labs(title = NULL, subtitle = NULL)`.
  (GitHub: [\#468](https://github.com/stan-dev/projpred/issues/468))
- [`plot.vsel()`](https://mc-stan.org/projpred/dev/reference/plot.vsel.md)
  has gained a new argument `show_cv_proportions`, allowing to omit the
  CV ranking proportions. (GitHub:
  [\#470](https://github.com/stan-dev/projpred/issues/470))
- Renamed
  [`summary.vsel()`](https://mc-stan.org/projpred/dev/reference/summary.vsel.md)’s
  output element `selection` to `perf_sub` and made the names of this
  `data.frame`’s columns more consistent so that it is easier to handle
  that `data.frame` programmatically. This should not be a breaking
  change because elements of `vselsummary` objects (i.e., elements of
  objects returned by
  [`summary.vsel()`](https://mc-stan.org/projpred/dev/reference/summary.vsel.md))
  are not meant to be accessed directly (for elements `perf_sub` and
  `perf_ref`, the new helper function
  [`performances()`](https://mc-stan.org/projpred/dev/reference/performances.md)
  has been added, see “Major changes” above). (GitHub:
  [\#471](https://github.com/stan-dev/projpred/issues/471))
- In
  [`summary.vsel()`](https://mc-stan.org/projpred/dev/reference/summary.vsel.md)
  and
  [`plot.vsel()`](https://mc-stan.org/projpred/dev/reference/plot.vsel.md),
  the `NA_character_` “string” (which was previously used as a
  placeholder for the predictor term of the intercept-only model at size
  `0`) was replaced by the string `"(Intercept)"`. (GitHub:
  [\#471](https://github.com/stan-dev/projpred/issues/471))
- Renamed
  [`project()`](https://mc-stan.org/projpred/dev/reference/project.md)’s
  output element `solution_terms` to `predictor_terms`. This should not
  be a breaking change because that element is meant to be accessed via
  [`predictor_terms()`](https://mc-stan.org/projpred/dev/reference/predictor_terms.md).
  (GitHub: [\#472](https://github.com/stan-dev/projpred/issues/472))
- Renamed elements `solution_terms` and `solution_terms_cv` of `vsel`
  objects (returned by
  [`varsel()`](https://mc-stan.org/projpred/dev/reference/varsel.md) and
  [`cv_varsel()`](https://mc-stan.org/projpred/dev/reference/cv_varsel.md))
  to `predictor_ranking` and `predictor_ranking_cv`, respectively. This
  should not be a breaking change because those elements are meant to be
  accessed via
  [`ranking()`](https://mc-stan.org/projpred/dev/reference/ranking.md).
  (GitHub: [\#472](https://github.com/stan-dev/projpred/issues/472))
- Global option `projpred.verbose_project` now affects the verbosity of
  *all* projections performed by the built-in divergence minimizers
  (except for the built-in L1-projection divergence minimizer). In
  particular, the divergence minimizer (no matter whether built-in or
  user-specified) is also employed when calling
  [`varsel()`](https://mc-stan.org/projpred/dev/reference/varsel.md) or
  [`cv_varsel()`](https://mc-stan.org/projpred/dev/reference/cv_varsel.md),
  so setting option `projpred.verbose_project` to `TRUE` now shows the
  progress of the projections during a
  [`varsel()`](https://mc-stan.org/projpred/dev/reference/varsel.md) or
  [`cv_varsel()`](https://mc-stan.org/projpred/dev/reference/cv_varsel.md)
  call. Previously, that option only affected the projections performed
  through
  [`project()`](https://mc-stan.org/projpred/dev/reference/project.md)
  (see the default for
  [`project()`](https://mc-stan.org/projpred/dev/reference/project.md)’s
  argument `verbose`). Usually, setting `projpred.verbose_project` to
  `TRUE` only makes sense when setting global option
  `projpred.extra_verbose` and argument `verbose` (of
  [`varsel()`](https://mc-stan.org/projpred/dev/reference/varsel.md) or
  [`cv_varsel()`](https://mc-stan.org/projpred/dev/reference/cv_varsel.md))
  to `TRUE` as well.
- Added [`print()`](https://rdrr.io/r/base/print.html) methods for
  objects of class `refmodel` and `projection`, mainly to avoid
  cluttering the console when printing such objects accidentally.
- Argument `extract_model_data` of
  [`init_refmodel()`](https://mc-stan.org/projpred/dev/reference/refmodel-init-get.md)
  is now allowed to be `NULL` for using an internal default.
- [`print.vselsummary()`](https://mc-stan.org/projpred/dev/reference/print.vselsummary.md)
  and
  [`print.vsel()`](https://mc-stan.org/projpred/dev/reference/print.vsel.md)
  now use a minimum number of significant digits of `2` by default. The
  previous behavior can be restored by setting
  `options(projpred.digits = getOption("digits"))`.
- Added a new performance statistic, the geometric mean predictive
  density (GMPD). This is particularly useful for discrete outcomes
  because there, the GMPD is a geometric mean of probabilities and hence
  bounded by zero and one. For details, see argument `stats` of the
  [`?summary.vsel`](https://mc-stan.org/projpred/dev/reference/summary.vsel.md)
  help. (GitHub:
  [\#476](https://github.com/stan-dev/projpred/issues/476))
- [`project()`](https://mc-stan.org/projpred/dev/reference/project.md)’s
  argument `verbose` now gets passed to argument `verbose_divmin` (not
  `projpred_verbose`) of the divergence minimizer function (see argument
  `div_minimizer` of
  [`init_refmodel()`](https://mc-stan.org/projpred/dev/reference/refmodel-init-get.md)).
- Arguments `lambda_min_ratio`, `nlambda`, and `thresh` of
  [`varsel()`](https://mc-stan.org/projpred/dev/reference/varsel.md) and
  [`cv_varsel()`](https://mc-stan.org/projpred/dev/reference/cv_varsel.md)
  have been deprecated. Instead,
  [`varsel()`](https://mc-stan.org/projpred/dev/reference/varsel.md) and
  [`cv_varsel()`](https://mc-stan.org/projpred/dev/reference/cv_varsel.md)
  have gained a new argument called `search_control` which accepts
  control arguments for the search as a `list`. Thus, former arguments
  `lambda_min_ratio`, `nlambda`, and `thresh` should now be specified
  via `search_control` (but note that `search_control` is more general
  because it also accepts control arguments for a *forward* search).
  (GitHub: [\#477](https://github.com/stan-dev/projpred/issues/477))
- [`run_cvfun()`](https://mc-stan.org/projpred/dev/reference/run_cvfun.md)
  has gained a new argument `folds`, accepting a vector of fold indices
  (the default is `NULL`, meaning that the folds are constructed
  internally, as before). This new argument is helpful, for example, to
  perform a stratified K-fold CV in a convenient manner (an example of
  this has been added to the
  [`?run_cvfun`](https://mc-stan.org/projpred/dev/reference/run_cvfun.md)
  help). (GitHub:
  [\#480](https://github.com/stan-dev/projpred/issues/480))
- [`plot.vsel()`](https://mc-stan.org/projpred/dev/reference/plot.vsel.md)
  has gained a new argument `size_position`. Setting it to
  `"primary_x_top"` moves the text for the submodel sizes *above* the
  x-axis. Setting it to `"secondary_x"` moves that text into a secondary
  x-axis located at the top of the plot. (GitHub:
  [\#484](https://github.com/stan-dev/projpred/issues/484))
- The default alignments of the x-axis text in
  [`plot.vsel()`](https://mc-stan.org/projpred/dev/reference/plot.vsel.md)
  have been changed: x-axis text is now right-aligned (left-aligned) for
  `text_angle > 0` (`< 0`) and also top-aligned for
  `-90 < text_angle && text_angle < 90 && text_angle != 0`. We emphasize
  that alignments can always be customized with
  `<plot.vsel() output object> + ggplot2::theme(axis.text.x.bottom = ggplot2::element_text(hjust = <hjust_value>, vjust = <vjust_value>))`.
  (GitHub: [\#484](https://github.com/stan-dev/projpred/issues/484))

### Bug fixes

- Fixed a bug sometimes causing
  [`plot.vsel()`](https://mc-stan.org/projpred/dev/reference/plot.vsel.md)
  to produce extra (“empty”) ticks on the x-axis. (GitHub:
  [\#462](https://github.com/stan-dev/projpred/issues/462))
- Fixed a bug in
  [`summary.vsel()`](https://mc-stan.org/projpred/dev/reference/summary.vsel.md)
  and
  [`plot.vsel()`](https://mc-stan.org/projpred/dev/reference/plot.vsel.md)
  causing bootstrap results (i.e., standard error and confidence
  interval for RMSE and AUC) to be incorrect if `deltas = TRUE`.
  (GitHub: [\#474](https://github.com/stan-dev/projpred/issues/474))
- Fixed several bugs in
  [`summary.vsel()`](https://mc-stan.org/projpred/dev/reference/summary.vsel.md)
  and
  [`plot.vsel()`](https://mc-stan.org/projpred/dev/reference/plot.vsel.md)
  sometimes causing incorrect predictive performance results in case of
  subsampled PSIS-LOO CV (an experimental feature controlled by argument
  `nloo` of
  [`cv_varsel()`](https://mc-stan.org/projpred/dev/reference/cv_varsel.md)).
  (GitHub: [\#475](https://github.com/stan-dev/projpred/issues/475))
- Fixed backward compatibility for the legacy structure of `cvfits` (the
  new structure was introduced by version 2.7.0, see GitHub pull request
  [\#456](https://github.com/stan-dev/projpred/issues/456)).

## projpred 2.7.0

CRAN release: 2023-09-30

### Major changes

- The default search `method` is now `"forward"` search for all kinds of
  models (previously, `"L1"` search was used by default where
  available). The reason for this change is that in general, forward
  search is more favorable compared to L1 search (see section “Details”
  in [`?varsel`](https://mc-stan.org/projpred/dev/reference/varsel.md)
  or
  [`?cv_varsel`](https://mc-stan.org/projpred/dev/reference/cv_varsel.md)).
  (GitHub: [\#453](https://github.com/stan-dev/projpred/issues/453),
  [\#459](https://github.com/stan-dev/projpred/issues/459))

- Several enhancements with respect to projected draws with different
  (i.e., nonconstant) weights, which typically occurs in case of
  clustered projection (GitHub:
  [\#206](https://github.com/stan-dev/projpred/issues/206),
  [\#439](https://github.com/stan-dev/projpred/issues/439)):

  - [`as.matrix.projection()`](https://mc-stan.org/projpred/dev/reference/as.matrix.projection.md)
    now throws an error if the projected draws have nonconstant weights.
    (This error is the default behavior; it can be avoided by setting
    the new argument `allow_nonconst_wdraws_prj` to `TRUE`, but this is
    for expert use only because in that case, the weights of the
    projected draws are stored in an attribute `wdraws_prj` and handling
    this attribute requires special care, e.g., when subsetting the
    returned matrix.) Instead, a
    [`posterior::as_draws_matrix()`](https://mc-stan.org/posterior/reference/draws_matrix.html)
    method
    ([`as_draws_matrix.projection()`](https://mc-stan.org/projpred/dev/reference/as_draws_matrix.projection.md))
    has been added which allows for a safer handling of these weights
    (e.g., with the help of
    [`posterior::resample_draws()`](https://mc-stan.org/posterior/reference/resample_draws.html),
    see section “Examples” of the
    [`?as_draws_matrix.projection`](https://mc-stan.org/projpred/dev/reference/as_draws_matrix.projection.md)
    help). Just like
    [`as.matrix.projection()`](https://mc-stan.org/projpred/dev/reference/as.matrix.projection.md),
    [`as_draws_matrix.projection()`](https://mc-stan.org/projpred/dev/reference/as_draws_matrix.projection.md)
    also works for the more common case of projected draws with constant
    weights. A
    [`posterior::as_draws()`](https://mc-stan.org/posterior/reference/draws.html)
    method
    ([`as_draws.projection()`](https://mc-stan.org/projpred/dev/reference/as_draws_matrix.projection.md))
    has also been added, but this is merely a wrapper for
    [`as_draws_matrix.projection()`](https://mc-stan.org/projpred/dev/reference/as_draws_matrix.projection.md).
  - [`proj_linpred()`](https://mc-stan.org/projpred/dev/reference/pred-projection.md)
    now also throws an error (by default) if the projected draws have
    nonconstant weights and has gained the new arguments
    `allow_nonconst_wdraws_prj` and `return_draws_matrix`. As in
    [`as.matrix.projection()`](https://mc-stan.org/projpred/dev/reference/as.matrix.projection.md),
    argument `allow_nonconst_wdraws_prj` is for expert use only.
    Instead, `return_draws_matrix` is the intended argument in case of
    projected draws with nonconstant weights. Similarly to
    [`as_draws_matrix.projection()`](https://mc-stan.org/projpred/dev/reference/as_draws_matrix.projection.md),
    it requires the **posterior** package and returns a `draws_matrix`
    (with weighted draws if the projected draws have nonconstant weights
    and `integrated` is `FALSE`).
  - [`proj_predict()`](https://mc-stan.org/projpred/dev/reference/pred-projection.md)
    has gained an argument `return_draws_matrix` for converting the
    returned matrix to a `draws_matrix` (which again requires the
    **posterior** package). For
    [`proj_predict()`](https://mc-stan.org/projpred/dev/reference/pred-projection.md),
    no further modifications were necessary because its argument
    `nresample_clusters` already takes weights of projected draws
    appropriately into account.

- Added helper function
  [`run_cvfun()`](https://mc-stan.org/projpred/dev/reference/run_cvfun.md)
  which can be used to create input for
  [`cv_varsel.refmodel()`](https://mc-stan.org/projpred/dev/reference/cv_varsel.md)’s
  new argument `cvfits` (which is the same as
  [`init_refmodel()`](https://mc-stan.org/projpred/dev/reference/refmodel-init-get.md)’s
  argument `cvfits`, but avoids having to call
  [`init_refmodel()`](https://mc-stan.org/projpred/dev/reference/refmodel-init-get.md)
  or
  [`get_refmodel()`](https://mc-stan.org/projpred/dev/reference/refmodel-init-get.md)
  twice). See the documentation of
  [`run_cvfun()`](https://mc-stan.org/projpred/dev/reference/run_cvfun.md)
  for details. (GitHub:
  [\#458](https://github.com/stan-dev/projpred/issues/458))

- Users applying
  [`varsel()`](https://mc-stan.org/projpred/dev/reference/varsel.md) or
  [`cv_varsel()`](https://mc-stan.org/projpred/dev/reference/cv_varsel.md)
  to an object of class `vsel` now need to use
  `varsel(get_refmodel(<vsel_object>), <...>)` or
  `cv_varsel(get_refmodel(<vsel_object>), <...>)` instead of
  `varsel(<vsel_object>, <...>)` and `cv_varsel(<vsel_object>, <...>)`,
  respectively. The reason is that new methods
  [`varsel.vsel()`](https://mc-stan.org/projpred/dev/reference/varsel.md)
  and
  [`cv_varsel.vsel()`](https://mc-stan.org/projpred/dev/reference/cv_varsel.md)
  have been added. Currently, these are only placeholders, but in a
  future release, they will offer new functionality.

### Minor changes

- If an L1 search selects an interaction term before all involved
  lower-order interaction terms (including main-effect terms) have been
  selected, the predictor ranking is now automatically modified so that
  the lower-order interaction terms come before this interaction term. A
  corresponding warning is thrown, which may be deactivated by setting
  the global option `projpred.warn_L1_interactions` to `FALSE`.
  Previously, beginning with version 2.5.0, only a warning was thrown
  and this only if an L1 search selected an interaction term before all
  involved *main-effect* terms had been selected. (GitHub:
  [\#420](https://github.com/stan-dev/projpred/issues/420))

- Added a progress bar for
  [`project()`](https://mc-stan.org/projpred/dev/reference/project.md)
  (when using the built-in divergence minimizers). For this,
  [`project()`](https://mc-stan.org/projpred/dev/reference/project.md)
  has gained a new argument `verbose` which can also be controlled via
  the global option `projpred.verbose_project`. By default, the new
  progress bar is activated. (GitHub:
  [\#421](https://github.com/stan-dev/projpred/issues/421))

- Added a new argument `parallel` to
  [`cv_varsel()`](https://mc-stan.org/projpred/dev/reference/cv_varsel.md).
  With `parallel = TRUE`, costly parts of **projpred**’s
  cross-validation (CV) can be run in parallel. See the documentation of
  that new argument (and section “Note” of
  [`cv_varsel()`](https://mc-stan.org/projpred/dev/reference/cv_varsel.md)’s
  documentation) for details. (GitHub:
  [\#422](https://github.com/stan-dev/projpred/issues/422))

- Added a warning for issue
  [\#323](https://github.com/stan-dev/projpred/issues/323) (for
  multilevel Gaussian models, the projection onto the full model can be
  instable). (GitHub:
  [\#426](https://github.com/stan-dev/projpred/issues/426))

- [`plot.vsel()`](https://mc-stan.org/projpred/dev/reference/plot.vsel.md)
  has gained the new arguments `point_size` and `bar_thickness` which
  control the size of the points and the thickness of the uncertainty
  bars, respectively. By default, the points are slightly larger now and
  the uncertainty bars slightly thicker than before. The previous
  appearance can be achieved by setting `point_size = 1.5` and
  `bar_thickness = 0.5`. (GitHub:
  [\#429](https://github.com/stan-dev/projpred/issues/429),
  [\#443](https://github.com/stan-dev/projpred/issues/443))

- [`plot.vsel()`](https://mc-stan.org/projpred/dev/reference/plot.vsel.md):
  Added argument `ranking_colored` for coloring the points and the
  uncertainty bars according to the magnitude of the (possibly
  cumulated) CV ranking proportions. (GitHub:
  [\#430](https://github.com/stan-dev/projpred/issues/430); thanks to
  [@yannmclatchie](https://github.com/yannmclatchie) for the suggestion)

- Added warnings for most of the problems described in section
  [“Troubleshooting”](https://mc-stan.org/projpred/articles/projpred.html#troubleshooting)
  of the main vignette. (GitHub:
  [\#431](https://github.com/stan-dev/projpred/issues/431))

- Output element `p_type` of
  [`project()`](https://mc-stan.org/projpred/dev/reference/project.md)
  has been removed. Instead, output element `const_wdraws_prj` has been
  added, but its definition is essentially the inverse of former element
  `p_type` (see the updated documentation of
  [`project()`](https://mc-stan.org/projpred/dev/reference/project.md)’s
  output). This should not be a breaking change for users (as `p_type`
  was mainly intended for internal use and the new element
  `const_wdraws_prj` is so, too) but this slightly enhances the cases
  where
  [`as.matrix.projection()`](https://mc-stan.org/projpred/dev/reference/as.matrix.projection.md)
  used to throw a warning (and now throws an error; see “Major changes”
  above) concerning the weights of the projected draws and the cases
  where
  [`proj_predict()`](https://mc-stan.org/projpred/dev/reference/pred-projection.md)
  resamples from the projected draws using argument
  `nresample_clusters`. (GitHub:
  [\#432](https://github.com/stan-dev/projpred/issues/432))

- Improved handling of PSIS-LOO CV warnings. (GitHub:
  [\#438](https://github.com/stan-dev/projpred/issues/438),
  [\#451](https://github.com/stan-dev/projpred/issues/451))

- Reduced peak memory usage during forward search. A global option
  `projpred.run_gc` has also been added, see the general package
  documentation (available
  [online](https://mc-stan.org/projpred/reference/projpred-package.html)
  or by typing
  [`` ?`projpred-package` ``](https://mc-stan.org/projpred/dev/reference/projpred-package.md)).
  (GitHub: [\#442](https://github.com/stan-dev/projpred/issues/442))

- Slightly improved efficiency in K-fold and PSIS-LOO CV, especially in
  case of a large number of observations. Under very special conditions
  (`refit_prj = FALSE`, `1 < nclusters && nclusters < S`, and
  `1 < nclusters_pred && nclusters_pred < S`; note that `1 < nclusters`
  requires forward search, `S` denotes the number of posterior draws in
  the reference model, and `nclusters_pred` is essentially unused if
  `refit_prj = FALSE`), this change might affect K-fold CV results, due
  to a different pseudorandom number generator (PRNG) state in folds
  other than the first one. Under similarly special conditions
  (`refit_prj = FALSE` and `1 < nclusters_pred && nclusters_pred < S`),
  the PRNG state for LOO subsampling (see argument `nloo`) is affected.
  Furthermore, if `is.na(seed)`, then the PRNG state for code downstream
  of such
  [`cv_varsel()`](https://mc-stan.org/projpred/dev/reference/cv_varsel.md)
  calls will be different due to this change. (GitHub:
  [\#446](https://github.com/stan-dev/projpred/issues/446))

- Slightly improved efficiency at the end of
  [`cv_varsel()`](https://mc-stan.org/projpred/dev/reference/cv_varsel.md),
  especially in case of a large number of observations. If
  `is.na(seed)`, then the PRNG state for code downstream of a
  [`cv_varsel()`](https://mc-stan.org/projpred/dev/reference/cv_varsel.md)
  call with `refit_prj = TRUE` and
  `1 < nclusters_pred && nclusters_pred < S` (where `S` denotes the
  number of posterior draws in the reference model) will be different
  due to this change. (GitHub:
  [\#447](https://github.com/stan-dev/projpred/issues/447))

- Slightly improved memory usage in
  [`varsel()`](https://mc-stan.org/projpred/dev/reference/varsel.md),
  [`cv_varsel()`](https://mc-stan.org/projpred/dev/reference/cv_varsel.md),
  and
  [`project()`](https://mc-stan.org/projpred/dev/reference/project.md).
  In case of LOO subsampling (see argument `nloo`) with clustered
  projection (i.e., `1 < nclusters && nclusters < S` or
  `1 < nclusters_pred && nclusters_pred < S`, where `S` denotes the
  number of posterior draws in the reference model), this change may
  lead to slightly different results due to different internal PRNG
  states. Furthermore, if `is.na(seed)`, then the PRNG state for code
  downstream of such a
  [`cv_varsel()`](https://mc-stan.org/projpred/dev/reference/cv_varsel.md)
  call will be different due to this change. (GitHub:
  [\#448](https://github.com/stan-dev/projpred/issues/448))

- The internal function `.extract_model_data` has been removed. As an
  alternative (with some differences compared to `.extract_model_data`),
  the new function
  [`y_wobs_offs()`](https://mc-stan.org/projpred/dev/reference/y_wobs_offs.md)
  is exported.

- Fixes/enhancements with respect to observation weights and offsets
  (GitHub: [\#449](https://github.com/stan-dev/projpred/issues/449)):

  - In case of an **rstanarm** reference model, the defaults for
    arguments `weightsnew` and `offsetnew` (see
    [`proj_linpred()`](https://mc-stan.org/projpred/dev/reference/pred-projection.md),
    [`proj_predict()`](https://mc-stan.org/projpred/dev/reference/pred-projection.md),
    and
    [`predict.refmodel()`](https://mc-stan.org/projpred/dev/reference/predict.refmodel.md))
    now cause the original observation weights and offsets to be used if
    possible (instead of ones and zeros, respectively, which could even
    be considered to have been a bug—hence why this is mentioned under
    “Bug fixes” as well). For **brms** reference models, this behavior
    had already been implemented before.
  - An error is now thrown if a length-zero element `weights` or
    `offset` is returned by the function supplied to argument
    `extract_model_data` of
    [`init_refmodel()`](https://mc-stan.org/projpred/dev/reference/refmodel-init-get.md)
    (before, a vector of ones or zeros was used silently for the
    observation weights and offsets, respectively).

- Added the helper function
  [`force_search_terms()`](https://mc-stan.org/projpred/dev/reference/force_search_terms.md)
  which allows to construct `search_terms` where certain predictor terms
  are forced to be included (i.e., they are forced to be selected first)
  whereas other predictor terms are optional (i.e., they are subject to
  the variable selection, but only after the inclusion of the “forced”
  terms). (GitHub:
  [\#346](https://github.com/stan-dev/projpred/issues/346))

- Reduced peak memory usage during performance evaluation (more
  precisely, during the re-projections done for the performance
  evaluation). This reduction is considerable especially for multilevel
  submodels, but possibly also for additive submodels. (GitHub:
  [\#440](https://github.com/stan-dev/projpred/issues/440),
  [\#450](https://github.com/stan-dev/projpred/issues/450))

- A message is now thrown when cutting off the search at `nterms_max`’s
  internal default of (currently) `19`. (GitHub:
  [\#452](https://github.com/stan-dev/projpred/issues/452))

- Added sub-section
  [“Speed”](https://mc-stan.org/projpred/articles/projpred.html#speed)
  to the main vignette’s
  [“Troubleshooting”](https://mc-stan.org/projpred/articles/projpred.html#troubleshooting)
  section. (GitHub:
  [\#455](https://github.com/stan-dev/projpred/issues/455))

- In case of K-fold CV, the `list` passed to argument `cvfits` of
  [`init_refmodel()`](https://mc-stan.org/projpred/dev/reference/refmodel-init-get.md)
  should not have a sub-`list` called `fits` anymore. Instead, the
  content of this former sub-`list` called `fits` should be moved one
  level up, i.e., should be placed directly in the `list` passed to
  `cvfits` (the empty element `fits` should then be removed). For some
  time, the old structure will continue to work, but this possibility is
  deprecated and will be removed in the future. (GitHub:
  [\#456](https://github.com/stan-dev/projpred/issues/456))

- In case of K-fold CV, the `K` reference model fits (i.e., the elements
  of the return value of the function passed to argument `cvfun` of
  [`init_refmodel()`](https://mc-stan.org/projpred/dev/reference/refmodel-init-get.md)
  or the elements of the `list` supplied to argument `cvfits` of
  [`init_refmodel()`](https://mc-stan.org/projpred/dev/reference/refmodel-init-get.md))
  do not need to be `list`s anymore (see the documentation for argument
  `cvrefbuilder` of
  [`init_refmodel()`](https://mc-stan.org/projpred/dev/reference/refmodel-init-get.md)).
  (GitHub: [\#457](https://github.com/stan-dev/projpred/issues/457))

### Bug fixes

- Fixed a bug in the printed number of projected draws for the
  performance evaluation when calling
  [`print.vselsummary()`](https://mc-stan.org/projpred/dev/reference/print.vselsummary.md)
  based on output from
  [`varsel()`](https://mc-stan.org/projpred/dev/reference/varsel.md)
  with `refit_prj = FALSE`.
- Fixed a bug sometimes causing an error when predicting from a submodel
  that is a GLM and has interactions. (GitHub:
  [\#420](https://github.com/stan-dev/projpred/issues/420))
- Fixed a bug introduced in version 2.6.0, causing an incompatibility of
  K-fold CV with R versions \< 4.2.0. (GitHub:
  [\#423](https://github.com/stan-dev/projpred/issues/423),
  [\#427](https://github.com/stan-dev/projpred/issues/427))
- Fixed a bug for the augmented-data projection in combination with
  subsampled PSIS-LOO CV. (GitHub:
  [\#433](https://github.com/stan-dev/projpred/issues/433))
- [`cv_varsel()`](https://mc-stan.org/projpred/dev/reference/cv_varsel.md)
  with `validate_search = FALSE` used to call
  [`loo::psis()`](https://mc-stan.org/loo/reference/psis.html) (for the
  submodel performance evaluation PSIS-LOO CV) even in case of draws
  with different (i.e., nonconstant) weights. In such cases,
  [`loo::sis()`](https://mc-stan.org/loo/reference/sis.html) is called
  now (with a warning). (GitHub:
  [\#438](https://github.com/stan-dev/projpred/issues/438))
- Fixed a bug for **rstanarm** (and custom) multilevel reference models
  with interactions (`:` syntax) between grouping variables, caused by
  missing columns in the reference model’s `data.frame` (for **brms**
  reference models, this was already done correctly). (GitHub:
  [\#445](https://github.com/stan-dev/projpred/issues/445))
- In case of an **rstanarm** reference model, the defaults for arguments
  `weightsnew` and `offsetnew` (see
  [`proj_linpred()`](https://mc-stan.org/projpred/dev/reference/pred-projection.md),
  [`proj_predict()`](https://mc-stan.org/projpred/dev/reference/pred-projection.md),
  and
  [`predict.refmodel()`](https://mc-stan.org/projpred/dev/reference/predict.refmodel.md))
  now cause the original observation weights and offsets to be used if
  possible (instead of ones and zeros, respectively, which could be
  considered to have been a bug). For **brms** reference models, this
  behavior had already been implemented before. (GitHub:
  [\#449](https://github.com/stan-dev/projpred/issues/449))
- Fixed a bug causing PSIS-LOO CV with `validate_search = FALSE` to fail
  in case of a single projected draw. (GitHub:
  [\#454](https://github.com/stan-dev/projpred/issues/454))

## projpred 2.6.0

CRAN release: 2023-06-01

### Major changes

- In anticipation of a larger overhaul of the **projpred** user
  interface, this release comes with several new functions for accessing
  and investigating solution paths (which are now termed *predictor
  rankings* by these new functions, a term that is hopefully easier to
  grasp for new users):

  - Added a new function called
    [`ranking()`](https://mc-stan.org/projpred/dev/reference/ranking.md)
    which returns the predictor ranking from the full-data search and
    possibly also the predictor rankings from fold-wise searches in case
    of cross-validation (CV). (More precisely,
    [`ranking()`](https://mc-stan.org/projpred/dev/reference/ranking.md)
    is a generic. The only method is
    [`ranking.vsel()`](https://mc-stan.org/projpred/dev/reference/ranking.md),
    applicable to objects returned by
    [`varsel()`](https://mc-stan.org/projpred/dev/reference/varsel.md)
    or
    [`cv_varsel()`](https://mc-stan.org/projpred/dev/reference/cv_varsel.md).
    The output is of class `ranking`.)
  - Added a new function called
    [`cv_proportions()`](https://mc-stan.org/projpred/dev/reference/cv_proportions.md)
    which computes ranking proportions (across CV folds, see
    [`?cv_proportions`](https://mc-stan.org/projpred/dev/reference/cv_proportions.md)
    for details) from fold-wise predictor rankings. (More precisely,
    [`cv_proportions()`](https://mc-stan.org/projpred/dev/reference/cv_proportions.md)
    is a generic. The main method is
    [`cv_proportions.ranking()`](https://mc-stan.org/projpred/dev/reference/cv_proportions.md),
    but as a shortcut,
    [`cv_proportions.vsel()`](https://mc-stan.org/projpred/dev/reference/cv_proportions.md)
    has also been added. The output is of class `cv_proportions`.)
  - Added a new [`plot()`](https://rdrr.io/r/graphics/plot.default.html)
    method called
    [`plot.cv_proportions()`](https://mc-stan.org/projpred/dev/reference/plot.cv_proportions.md)
    for plotting ranking proportions from fold-wise predictor rankings.
    (As a shortcut,
    [`plot.ranking()`](https://mc-stan.org/projpred/dev/reference/plot.cv_proportions.md)
    has also been added.)

  Because of these new functions, a message has been added to
  [`print.vselsummary()`](https://mc-stan.org/projpred/dev/reference/print.vselsummary.md),
  mentioning how to access and investigate the fold-wise predictor
  rankings (if they exist). Furthermore, due to these changes, element
  `pct_solution_terms_cv` of `vsel` objects has been replaced with
  element `solution_terms_cv` which contains the fold-wise predictor
  rankings instead of the corresponding ranking proportions. However,
  elements of `vsel` objects are not meant to be accessed directly, so
  this replacement should not be a breaking change for most users.
  Finally, method
  [`solution_terms.vsel()`](https://mc-stan.org/projpred/dev/reference/solution_terms.md)
  (which—until now—was the only possibility to extract the full-data
  predictor ranking) has now been deprecated and will be removed in a
  future release. Please use the new function
  [`ranking()`](https://mc-stan.org/projpred/dev/reference/ranking.md)
  instead (more precisely,
  [`ranking()`](https://mc-stan.org/projpred/dev/reference/ranking.md)’s
  output element `fulldata` contains the full-data predictor ranking
  that is also extracted by
  [`solution_terms.vsel()`](https://mc-stan.org/projpred/dev/reference/solution_terms.md);
  [`ranking()`](https://mc-stan.org/projpred/dev/reference/ranking.md)’s
  output element `foldwise` contains the fold-wise predictor rankings—if
  available—which were previously not accessible via a built-in
  function). (GitHub:
  [\#289](https://github.com/stan-dev/projpred/issues/289),
  [\#406](https://github.com/stan-dev/projpred/issues/406),
  [\#411](https://github.com/stan-dev/projpred/issues/411))

- Added function
  [`predictor_terms()`](https://mc-stan.org/projpred/dev/reference/predictor_terms.md)
  which retrieves the predictor terms used in a
  [`project()`](https://mc-stan.org/projpred/dev/reference/project.md)
  run. Correspondingly, method
  [`solution_terms.projection()`](https://mc-stan.org/projpred/dev/reference/solution_terms.md)
  has now been deprecated and will be removed in a future release.
  Please use
  [`predictor_terms()`](https://mc-stan.org/projpred/dev/reference/predictor_terms.md)
  instead. (GitHub:
  [\#411](https://github.com/stan-dev/projpred/issues/411))

- Renamed function
  [`cvfolds()`](https://mc-stan.org/projpred/dev/reference/cv-indices.md)
  to
  [`cv_folds()`](https://mc-stan.org/projpred/dev/reference/cv-indices.md)
  (more precisely, the former variant still exists, but is deprecated
  and will be removed in a future release). (GitHub:
  [\#411](https://github.com/stan-dev/projpred/issues/411))

- `seed` (and `.seed`) arguments now have a default of `NA` instead of
  `sample.int(.Machine$integer.max, 1)` and the pseudorandom number
  generator (PRNG) state is reset only if the user-supplied seed is not
  `NA`. This allows setting a seed once at the beginning of any
  **projpred**-related code and then leaving all `seed` (and `.seed`)
  arguments at their default. Previously, such practice could lead to
  results which were “less random” than they should have been because
  the former default of `sample.int(.Machine$integer.max, 1)` caused
  **projpred** functions with a `seed` (or `.seed`) argument to reset
  the PRNG state upon exit, meaning that two repeated calls to
  [`cv_varsel()`](https://mc-stan.org/projpred/dev/reference/cv_varsel.md)
  (for example) with no PRNG-using code between them would use the same
  seed internally. (GitHub:
  [\#412](https://github.com/stan-dev/projpred/issues/412))

- Added the main diagonal of the matrix returned by
  [`cv_proportions()`](https://mc-stan.org/projpred/dev/reference/cv_proportions.md)
  to a new column called `cv_proportions_diag` of the summary table
  computed by
  [`summary.vsel()`](https://mc-stan.org/projpred/dev/reference/summary.vsel.md).
  The purpose of this new column is to give a basic sense for the (CV)
  variability in the ranking of the predictors. Argument `cumulate` of
  [`cv_proportions()`](https://mc-stan.org/projpred/dev/reference/cv_proportions.md)
  has been added to
  [`summary.vsel()`](https://mc-stan.org/projpred/dev/reference/summary.vsel.md)
  as well (to allow the ranking proportions in the newly added column to
  be *cumulated* ranking proportions, if desired). (GitHub:
  [\#289](https://github.com/stan-dev/projpred/issues/289),
  [\#413](https://github.com/stan-dev/projpred/issues/413))

- Added the full-data predictor ranking and the main diagonal of the
  matrix returned by
  [`cv_proportions()`](https://mc-stan.org/projpred/dev/reference/cv_proportions.md)
  to the plot created by
  [`plot.vsel()`](https://mc-stan.org/projpred/dev/reference/plot.vsel.md).
  These new elements can be omitted by setting
  [`plot.vsel()`](https://mc-stan.org/projpred/dev/reference/plot.vsel.md)’s
  new argument `ranking_nterms_max` to `NA` (setting it to some specific
  submodel size causes the full-data predictor ranking and the
  corresponding ranking proportions to be omitted after that size).
  Argument `cumulate` of
  [`cv_proportions()`](https://mc-stan.org/projpred/dev/reference/cv_proportions.md)
  has been added to
  [`plot.vsel()`](https://mc-stan.org/projpred/dev/reference/plot.vsel.md)
  as well (to allow the ranking proportions to be *cumulated* ranking
  proportions, if desired). Other new arguments are `ranking_abbreviate`
  (together with `ranking_abbreviate_args`), `ranking_repel` (together
  with `ranking_repel_args`), and `text_angle` (see the
  [`plot.vsel()`](https://mc-stan.org/projpred/dev/reference/plot.vsel.md)
  documentation for details). (GitHub:
  [\#289](https://github.com/stan-dev/projpred/issues/289),
  [\#414](https://github.com/stan-dev/projpred/issues/414),
  [\#416](https://github.com/stan-dev/projpred/issues/416),
  [\#417](https://github.com/stan-dev/projpred/issues/417))

### Minor changes

- Enhancements in the vignettes. In particular, the new functions
  [`ranking()`](https://mc-stan.org/projpred/dev/reference/ranking.md),
  [`cv_proportions()`](https://mc-stan.org/projpred/dev/reference/cv_proportions.md),
  and
  [`plot.cv_proportions()`](https://mc-stan.org/projpred/dev/reference/plot.cv_proportions.md)
  (see “Major changes” above) are now illustrated in the main vignette.
  (GitHub: [\#407](https://github.com/stan-dev/projpred/issues/407),
  [\#411](https://github.com/stan-dev/projpred/issues/411))
- Reduced the peak memory usage of
  [`cv_varsel()`](https://mc-stan.org/projpred/dev/reference/cv_varsel.md)
  with `cv_method = "kfold"`. This may slightly change results from such
  a
  [`cv_varsel()`](https://mc-stan.org/projpred/dev/reference/cv_varsel.md)
  run compared to older **projpred** versions due to different
  pseudorandom number generator (PRNG) states when clustering posterior
  draws. (GitHub:
  [\#419](https://github.com/stan-dev/projpred/issues/419))
- The `cvfits` list (see
  [`init_refmodel()`](https://mc-stan.org/projpred/dev/reference/refmodel-init-get.md))
  does not need to have an attribute called `K` anymore.

### Bug fixes

- Fixed a bug causing L1 search to throw an error in case of some
  [`I()`](https://rdrr.io/r/base/AsIs.html) terms. (GitHub:
  [\#404](https://github.com/stan-dev/projpred/issues/404),
  [\#408](https://github.com/stan-dev/projpred/issues/408))
- Fixed a bug causing L1 search to throw an error in case of
  [`poly()`](https://rdrr.io/r/stats/poly.html) or
  [`polym()`](https://rdrr.io/r/stats/poly.html) terms. Note that just
  like [`step()`](https://rdrr.io/r/stats/step.html) and
  [`MASS::stepAIC()`](https://rdrr.io/pkg/MASS/man/stepAIC.html),
  **projpred**’s search algorithms do not split up a
  [`poly()`](https://rdrr.io/r/stats/poly.html) or
  [`polym()`](https://rdrr.io/r/stats/poly.html) term into its
  lower-degree polynomial terms (which would be helpful, for example, if
  the linear part of a [`poly()`](https://rdrr.io/r/stats/poly.html)
  term with `degrees = 2` was relevant but the quadratic part not). Such
  a split-up of a [`poly()`](https://rdrr.io/r/stats/poly.html) or
  [`polym()`](https://rdrr.io/r/stats/poly.html) term needs to be
  performed manually (if desired). (GitHub:
  [\#183](https://github.com/stan-dev/projpred/issues/183),
  [\#409](https://github.com/stan-dev/projpred/issues/409))
- Fixed a bug causing some non-smooth predictor terms to be treated as
  smooth terms. (GitHub:
  [\#182](https://github.com/stan-dev/projpred/issues/182),
  [\#410](https://github.com/stan-dev/projpred/issues/410))
- See “Major changes” above: Fixed a bug causing **projpred** functions
  with a `seed` (or `.seed`) argument to use the same seed internally
  when users set a seed once at the beginning (via
  [`set.seed()`](https://rdrr.io/r/base/Random.html)) and then had two
  or more calls to such **projpred** functions with their `seed` (or
  `.seed`) argument being at its default and no PRNG-using code between
  those calls. (GitHub:
  [\#412](https://github.com/stan-dev/projpred/issues/412))

## projpred 2.5.0

CRAN release: 2023-04-05

### Minor changes

- Setting the new global option `projpred.extra_verbose` to `TRUE` will
  print out which submodel **projpred** is currently projecting onto.
  Furthermore, if `method = "forward"` and `verbose = TRUE` in
  [`varsel()`](https://mc-stan.org/projpred/dev/reference/varsel.md) or
  [`cv_varsel()`](https://mc-stan.org/projpred/dev/reference/cv_varsel.md),
  this new option will also make **projpred** print out which submodel
  has been selected at those steps of the forward search for which a
  percentage is printed (the percentage refers to the maximum submodel
  size that the search is run up to). In general, however, we cannot
  recommend setting this new global option to `TRUE` for
  [`cv_varsel()`](https://mc-stan.org/projpred/dev/reference/cv_varsel.md)
  with `validate_search = TRUE` (simply due to the amount of information
  that will be printed, but also due to the progress bar which will not
  work anymore as intended). (GitHub:
  [\#363](https://github.com/stan-dev/projpred/issues/363); thanks to
  [@jtimonen](https://github.com/jtimonen))

- Enhanced `verbose` output. In particular,
  [`varsel()`](https://mc-stan.org/projpred/dev/reference/varsel.md) is
  now more verbose, similarly to how
  [`cv_varsel()`](https://mc-stan.org/projpred/dev/reference/cv_varsel.md)
  has already been for a long time. The `verbose` output for
  [`cv_varsel()`](https://mc-stan.org/projpred/dev/reference/cv_varsel.md)
  has also been updated, with the aim to give users a better
  understanding of the methodology behind **projpred**. (GitHub:
  [\#382](https://github.com/stan-dev/projpred/issues/382))

- Slightly improved the calculation of predictive variances to make them
  less prone to numerical inaccuracies. (GitHub:
  [\#199](https://github.com/stan-dev/projpred/issues/199))

- Improved computational efficiency by avoiding an unnecessary final
  full-data performance evaluation (including costly re-projections if
  `refit_prj = TRUE`, which is the default for non-`datafit` reference
  models) in
  [`cv_varsel()`](https://mc-stan.org/projpred/dev/reference/cv_varsel.md)
  with `validate_search = TRUE`. Due to this change, results from
  [`cv_varsel()`](https://mc-stan.org/projpred/dev/reference/cv_varsel.md)
  (with `validate_search = TRUE`) may slightly change due to a different
  pseudorandom number generator (PRNG) state when clustering posterior
  draws. The different PRNG state was necessary to make the PRNG state
  for the full-data search in the `validate_search = TRUE` case
  consistent to the PRNG state for the full-data search in the
  `validate_search = FALSE` case. (GitHub:
  [\#385](https://github.com/stan-dev/projpred/issues/385))

- Reduced dependencies. (GitHub:
  [\#388](https://github.com/stan-dev/projpred/issues/388))

- Argument `digits` of
  [`print.vselsummary()`](https://mc-stan.org/projpred/dev/reference/print.vselsummary.md)
  which used to be passed to an internal
  [`round()`](https://rdrr.io/r/base/Round.html) call was removed.
  Instead, `digits` can now be passed to
  [`print.data.frame()`](https://rdrr.io/r/base/print.dataframe.html)
  via `...`, thereby determining the minimum number of *significant
  digits* to be printed. (GitHub:
  [\#389](https://github.com/stan-dev/projpred/issues/389))

- Although bad practice (in general), a reference model lacking an
  intercept can now be used within **projpred**. However, it will always
  be projected onto submodels which *include* an intercept. The reason
  is that even if the true intercept in the reference model is zero,
  this does not need to hold for the submodels. An informational message
  mentioning the projection onto intercept-including submodels is thrown
  when **projpred** encounters a reference model lacking an intercept.
  (GitHub: [\#96](https://github.com/stan-dev/projpred/issues/96),
  [\#391](https://github.com/stan-dev/projpred/issues/391))

- In case of non-predictor arguments of `s()` or `t2()`, **projpred**
  now throws an error. (This had already been documented before, but a
  suitable error message was missing.) (GitHub:
  [\#393](https://github.com/stan-dev/projpred/issues/393), based on
  [\#156](https://github.com/stan-dev/projpred/issues/156) and
  [\#269](https://github.com/stan-dev/projpred/issues/269))

- In case of the
  [`brms::categorical()`](https://paulbuerkner.com/brms/reference/brmsfamily.html)
  family (supported since version 2.4.0), **projpred** now strips
  underscores from response category names in
  [`as.matrix.projection()`](https://mc-stan.org/projpred/dev/reference/as.matrix.projection.md)
  output, as done by **brms**. (GitHub:
  [\#394](https://github.com/stan-dev/projpred/issues/394))

- L1 search now throws a warning if an interaction term is selected
  before all involved main-effect terms have been selected. (GitHub:
  [\#395](https://github.com/stan-dev/projpred/issues/395))

- Documented that in multilevel (group-level) terms, function calls on
  the right-hand side of the `|` character (e.g.,
  `(1 | gr(group_variable))`, which is possible in **brms**) are
  currently not allowed in **projpred**. A corresponding error message
  has also been added. (GitHub:
  [\#319](https://github.com/stan-dev/projpred/issues/319))

- Due to internal refactoring:

  - [`project()`](https://mc-stan.org/projpred/dev/reference/project.md)’s
    output elements `submodl` and `weights` have been renamed to
    `outdmin` and `wdraws_prj`, respectively.
  - [`varsel()`](https://mc-stan.org/projpred/dev/reference/varsel.md)’s
    and
    [`cv_varsel()`](https://mc-stan.org/projpred/dev/reference/cv_varsel.md)’s
    output element `d_test` has been replaced with new output elements
    `type_test` and `y_wobs_test`.

  Apart from
  [`project()`](https://mc-stan.org/projpred/dev/reference/project.md)’s
  output element `wdraws_prj`, these elements are not meant to be
  accessed manually, so changes are mentioned here only for the sake of
  completeness. Output element `wdraws_prj` of
  [`project()`](https://mc-stan.org/projpred/dev/reference/project.md)
  is only needed if
  [`project()`](https://mc-stan.org/projpred/dev/reference/project.md)
  was used for a clustered projection, which is not the default (and
  discouraged in most applied cases, at least with a small number of
  clusters). Thus, these renamings are breaking changes only in very
  rare cases.

- [`print.vselsummary()`](https://mc-stan.org/projpred/dev/reference/print.vselsummary.md)
  now also prints `K` in case of K-fold CV.

- The
  [`print.vselsummary()`](https://mc-stan.org/projpred/dev/reference/print.vselsummary.md)
  output has been slightly improved, e.g., adding a remark what “search
  included” or “search not included” means.

- [`print.vselsummary()`](https://mc-stan.org/projpred/dev/reference/print.vselsummary.md)
  now also prints whether `deltas = TRUE` or `deltas = FALSE` was used.

- Output element `pct_solution_terms_cv` has now also been added to
  `vsel` objects returned by
  [`varsel()`](https://mc-stan.org/projpred/dev/reference/varsel.md),
  but in that case, it is simply `NULL`. This (`pct_solution_terms_cv`
  being `NULL`) is now also the case if `validate_search = FALSE` was
  used in
  [`cv_varsel()`](https://mc-stan.org/projpred/dev/reference/cv_varsel.md).

- Minor enhancements in the documentation.

- Enhancements in the vignettes. In particular, section
  [“Troubleshooting”](https://mc-stan.org/projpred/articles/projpred.html#troubleshooting)
  of the main vignette has been revised.

- If
  [`proj_predict()`](https://mc-stan.org/projpred/dev/reference/pred-projection.md)
  is used with observation weights that are not all equal to `1`, a
  warning is now thrown. (GitHub: starts to address
  [\#402](https://github.com/stan-dev/projpred/issues/402))

### Bug fixes

- Fixed a long-standing bug (existing at least from version 2.0.2 on)
  causing
  [`predict.refmodel()`](https://mc-stan.org/projpred/dev/reference/predict.refmodel.md)
  to require `newdata` to contain the response variable in case of a
  **brms** reference model. This is similar to
  [paul-buerkner/brms#1457](https://github.com/paul-buerkner/brms/issues/1457),
  but concerns
  [`predict.refmodel()`](https://mc-stan.org/projpred/dev/reference/predict.refmodel.md)
  ([paul-buerkner/brms#1457](https://github.com/paul-buerkner/brms/issues/1457)
  referred to predictions from the *submodels*). In order to make this
  [`predict.refmodel()`](https://mc-stan.org/projpred/dev/reference/predict.refmodel.md)
  fix work, **brms** version 2.19.0 or later is needed. (GitHub:
  [\#381](https://github.com/stan-dev/projpred/issues/381))
- Fixed a long-standing bug (existing from version 2.1.0 on) causing
  output element `p_type` of
  [`project()`](https://mc-stan.org/projpred/dev/reference/project.md)
  to be incorrect in case of `refit_prj = FALSE`, `!is.null(nclusters)`,
  and an `object` of class `vsel` that was created with a non-clustered
  (thinned) projection during the search phase. The fix comes with a
  slightly different behavior of
  [`proj_predict()`](https://mc-stan.org/projpred/dev/reference/pred-projection.md)
  for `datafit`s: It will not draw `nresample_clusters` times from the
  posterior-projection predictive distribution (which is based on the
  same single projected draw), but only once. (GitHub:
  [\#211](https://github.com/stan-dev/projpred/issues/211),
  [\#401](https://github.com/stan-dev/projpred/issues/401))
- When performing predictions from submodels which are GLMs (or from
  submodels which are L1-penalized GLMs, which is only possible in case
  of `refit_prj = FALSE` after an L1 search), a new dataset containing a
  `character` predictor variable with only a single unique value (or a
  new dataset containing a `factor` predictor variable with a single
  level) used to cause an error. The case of a `character` (not
  `factor`) predictor variable with only a single unique value occurred,
  e.g., during the performance evaluation in a LOO CV if a `character`
  predictor got selected into a fold’s solution path. The `character`
  issue existed from version 2.1.0 on (in earlier versions, however,
  there were other issues which caused `character` predictors to throw
  an error). Now, all issues with respect to `character` predictor
  variables should be resolved. The issue with single-level `factor`
  predictor variables is resolved now as well. (GitHub:
  [\#403](https://github.com/stan-dev/projpred/issues/403))
- When performing predictions from submodels which are GLMs (or from
  submodels which are L1-penalized GLMs, which is only possible in case
  of `refit_prj = FALSE` after an L1 search), a new dataset containing a
  `factor` predictor with re-ordered levels (compared to this same
  `factor` in the original dataset) used to lead to incorrect
  predictions. This bug existed at least from version 2.0.2 on (possibly
  even in earlier versions), but has been resolved now. (GitHub:
  [\#403](https://github.com/stan-dev/projpred/issues/403))
- Fixed an error thrown by **projpred**’s internal GLM submodel fitter
  in case of unused levels of a `factor`. This issue existed at least
  from version 2.0.2 on (possibly even in earlier versions), but should
  have only affected **rstanarm** reference model fits (**brms**
  reference model fits were only affected in case of a
  [`brms::brm()`](https://paulbuerkner.com/brms/reference/brm.html) call
  with `drop_unused_levels = FALSE`, which is not the default). (GitHub:
  [\#403](https://github.com/stan-dev/projpred/issues/403))
- Fixed a bug that caused an L1 search combined with `refit_prj = FALSE`
  (which is the default only for `datafit`s, not for the reference model
  objects of class `refmodel` that are usually employed in practice) to
  lead to incorrect predictions from the L1-searched submodels (which
  are L1-penalized GLMs) if the solution path had a main effect ranked
  after an interaction term. This bug existed at least from version
  2.0.2 on (possibly even in earlier versions). The mentioned submodel
  predictions did not only affect the performance evaluation, but also
  the projected dispersion parameter and the returned Kullback-Leibler
  divergence (and the corresponding cross-entropy). (GitHub:
  [\#403](https://github.com/stan-dev/projpred/issues/403))

## projpred 2.4.0

CRAN release: 2023-02-12

### Major changes

- Introduction of the augmented-data projection [(Weber et al.,
  2023)](https://doi.org/10.48550/arXiv.2301.01660) (see section
  [“Supported types of
  models”](https://mc-stan.org/projpred/articles/projpred.html#modtypes)
  of the main vignette for details). (GitHub:
  [\#70](https://github.com/stan-dev/projpred/issues/70),
  [\#322](https://github.com/stan-dev/projpred/issues/322))
- Introduction of the latent projection [(Catalina et al.,
  2021)](https://doi.org/10.48550/arXiv.2109.04702) (see section
  [“Supported types of
  models”](https://mc-stan.org/projpred/articles/projpred.html#modtypes)
  of the main vignette and the new [latent-projection
  vignette](https://mc-stan.org/projpred/articles/latent.html) for
  details). A consequence of the latent projection (more precisely, of
  the `resp_oscale = TRUE` default in
  [`summary.vsel()`](https://mc-stan.org/projpred/dev/reference/summary.vsel.md))
  is that
  [`varsel()`](https://mc-stan.org/projpred/dev/reference/varsel.md) and
  [`cv_varsel()`](https://mc-stan.org/projpred/dev/reference/cv_varsel.md)
  no longer call
  [`suggest_size()`](https://mc-stan.org/projpred/dev/reference/suggest_size.md)
  internally at the end. Thus,
  [`print()`](https://rdrr.io/r/base/print.html)-ing an object of class
  `vsel` no longer includes the suggested projection size in the output
  (the `stat` for this suggested size was fixed to `"elpd"` anyway, a
  fact that many users were probably not aware of). (GitHub:
  [\#372](https://github.com/stan-dev/projpred/issues/372))
- In case of multilevel models, **projpred** now has two global options
  for “integrating out” group-level effects: `projpred.mlvl_pred_new`
  and `projpred.mlvl_proj_ref_new`. These are explained in detail in the
  general package documentation (available
  [online](https://mc-stan.org/projpred/reference/projpred-package.html)
  or by typing
  [`` ?`projpred-package` ``](https://mc-stan.org/projpred/dev/reference/projpred-package.md)).
  (GitHub: [\#379](https://github.com/stan-dev/projpred/issues/379))

### Minor changes

- Improvements in the numerical stability of internal link and
  inverse-link functions. (GitHub:
  [\#376](https://github.com/stan-dev/projpred/issues/376))

### Bug fixes

- Fix a bug for offsets in cases where `family` (see
  [`init_refmodel()`](https://mc-stan.org/projpred/dev/reference/refmodel-init-get.md))
  has a non-identity link function: After clustering the reference
  model’s posterior draws, we need to aggregate (within a given cluster)
  the reference model’s fitted values which already take the offsets
  into account instead of taking the offsets into account after
  aggregating the fitted values which do *not* take the offsets into
  account. This fix should affect results only in a very slight manner.
  Due to **projpred**’s internal adjustment for numerical stability when
  averaging a quantity across the draws within a given cluster, this
  also changes the projected residual standard deviations in Gaussian
  models in the order of `1e-10`. (GitHub:
  [\#374](https://github.com/stan-dev/projpred/issues/374))

## projpred 2.3.0

CRAN release: 2023-01-10

### Major changes

- In
  [`plot.vsel()`](https://mc-stan.org/projpred/dev/reference/plot.vsel.md)
  and
  [`summary.vsel()`](https://mc-stan.org/projpred/dev/reference/summary.vsel.md),
  the default of `alpha = 0.32` is replaced by `alpha = 2 * pnorm(-1)`
  (= `1 - diff(pnorm(c(-1, 1)))`, which is only *approximately* 0.32) so
  that now, a normal-approximation confidence interval with default
  `alpha` stretches by exactly one standard error on either side of the
  point estimate. Typically, this changes results only slightly. In some
  cases, however, the new default may lead to a different suggested
  size, explaining why this is regarded as a major change. (GitHub:
  [\#371](https://github.com/stan-dev/projpred/issues/371))

### Minor changes

- The deprecated function
  [`ggplot2::aes_string()`](https://ggplot2.tidyverse.org/reference/aes_.html)
  is not used anymore, thereby avoiding an occasional soft-deprecation
  warning thrown by **ggplot2** 3.4.0. (GitHub:
  [\#367](https://github.com/stan-dev/projpred/issues/367))
- The KL divergence from the reference model to a submodel is simplified
  to the corresponding cross-entropy (i.e., the reference model’s
  entropy is dropped), with some caveats described in the documentation
  for output element `ce` of
  [`project()`](https://mc-stan.org/projpred/dev/reference/project.md).
  The reason for this change is that the former KL divergence assumed
  the reference model’s family to be the same as the submodel’s family,
  which does not need to be the case for custom reference models. This
  should not be a user-facing change as users are discouraged to make
  use of specific output elements (like the former element `kl` of
  objects of class `projection` or `vsel`) directly. (GitHub:
  [\#369](https://github.com/stan-dev/projpred/issues/369))
- Improvements in the documentation (especially for argument `family` of
  [`init_refmodel()`](https://mc-stan.org/projpred/dev/reference/refmodel-init-get.md)
  and
  [`get_refmodel.default()`](https://mc-stan.org/projpred/dev/reference/refmodel-init-get.md)).

## projpred 2.2.2

CRAN release: 2022-11-09

### Major changes

- Several important bug fixes (see below).

### Minor changes

- Improvements in documentation and vignette, especially to emphasize
  the generality of the reference model object resulting from
  [`get_refmodel()`](https://mc-stan.org/projpred/dev/reference/refmodel-init-get.md)
  and
  [`init_refmodel()`](https://mc-stan.org/projpred/dev/reference/refmodel-init-get.md)
  (thereby also distinguishing more clearly between “typical” and
  “custom” reference model objects) in (i) the description and several
  arguments of
  [`get_refmodel()`](https://mc-stan.org/projpred/dev/reference/refmodel-init-get.md)
  and
  [`init_refmodel()`](https://mc-stan.org/projpred/dev/reference/refmodel-init-get.md), (ii)
  sections [“Reference
  model”](https://mc-stan.org/projpred/articles/projpred.html#refmod)
  and [“Supported types of
  models”](https://mc-stan.org/projpred/articles/projpred.html#modtypes)
  of the vignette. (GitHub:
  [\#357](https://github.com/stan-dev/projpred/issues/357),
  [\#359](https://github.com/stan-dev/projpred/issues/359),
  [\#364](https://github.com/stan-dev/projpred/issues/364),
  [\#365](https://github.com/stan-dev/projpred/issues/365),
  [\#366](https://github.com/stan-dev/projpred/issues/366))
- Minor improvement in terms of efficiency in the
  `validate_search = FALSE` case of
  [`cv_varsel()`](https://mc-stan.org/projpred/dev/reference/cv_varsel.md).
- Improvement in terms of efficiency in case of a forward search with
  custom `search_terms` (at least in some instances), also affecting the
  output of `solution_terms(<vsel_object>)` in those cases. (GitHub:
  [\#360](https://github.com/stan-dev/projpred/issues/360); thanks to
  [@sor16](https://github.com/sor16))
- Update [Catalina et
  al. (2020)](https://doi.org/10.48550/arXiv.2010.06994) to [Catalina et
  al. (2022)](https://proceedings.mlr.press/v151/catalina22a.html).
  (GitHub: [\#364](https://github.com/stan-dev/projpred/issues/364))

### Bug fixes

- Fix a bug causing offsets not to be taken into account appropriately
  when calculating the PSIS weights (those used for the submodels) in
  the `validate_search = FALSE` case of
  [`cv_varsel()`](https://mc-stan.org/projpred/dev/reference/cv_varsel.md).
  This bug was introduced in v2.2.0 (and existed up
  to—including—v2.2.1).
- Fix a (long-standing) bug causing offsets not to be taken into account
  appropriately when calculating the predictive variances for a
  reference model that has a dispersion parameter and a non-identity
  link function. (GitHub:
  [\#186](https://github.com/stan-dev/projpred/issues/186) (partly),
  [\#355](https://github.com/stan-dev/projpred/issues/355))
- Fix a (long-standing) bug causing offsets not to be taken into account
  appropriately when calculating the reference model’s summary
  statistics in case of
  [`cv_varsel()`](https://mc-stan.org/projpred/dev/reference/cv_varsel.md)
  with `cv_method = "LOO"` (more precisely, only the LOO posterior
  predictive expected values `<vsel_object>$summaries$ref$mu` were
  affected, not the (pointwise) LOO log posterior predictive density
  values `<vsel_object>$summaries$ref$lppd`). (GitHub:
  [\#186](https://github.com/stan-dev/projpred/issues/186) (partly),
  [\#356](https://github.com/stan-dev/projpred/issues/356))
- Fix a (long-standing) bug leading to an error when trying to use
  [`cv_varsel()`](https://mc-stan.org/projpred/dev/reference/cv_varsel.md)
  with custom `search_terms` (in some instances). (GitHub:
  [\#345](https://github.com/stan-dev/projpred/issues/345),
  [\#360](https://github.com/stan-dev/projpred/issues/360); thanks to
  [@sor16](https://github.com/sor16))

## projpred 2.2.1

CRAN release: 2022-09-20

### Minor changes

- Several improvements in the documentation.
- For the RMSE as well as the AUC (see argument `stats` of
  [`summary.vsel()`](https://mc-stan.org/projpred/dev/reference/summary.vsel.md)),
  the bootstrapping results are now also used for inferring the lower
  and upper confidence interval bounds. (GitHub:
  [\#318](https://github.com/stan-dev/projpred/issues/318),
  [\#347](https://github.com/stan-dev/projpred/issues/347); thanks to
  [@awd97](https://github.com/awd97) and
  [@VisionResearchBlog](https://github.com/VisionResearchBlog))
- For `datafit`s, offsets are not supported anymore. (GitHub:
  [\#186](https://github.com/stan-dev/projpred/issues/186) (partly),
  [\#351](https://github.com/stan-dev/projpred/issues/351))

### Bug fixes

- Fix GitHub issue
  [\#348](https://github.com/stan-dev/projpred/issues/348) (L1 search in
  the presence of interaction terms). This bug was introduced in v2.1.0
  (and existed up to—including—v2.2.0).
- Fix incorrectly thrown messages in case of `datafit`s (and
  other—unlikely—cases where `nclusters == S` and `S <= 20`, with `S`
  denoting the number of draws in the reference model).
- Fix GitHub issue
  [\#349](https://github.com/stan-dev/projpred/issues/349) (only
  concerned `datafit`s). (GitHub:
  [\#350](https://github.com/stan-dev/projpred/issues/350))

## projpred 2.2.0

CRAN release: 2022-08-19

### Major changes

- In the `validate_search = FALSE` case of
  [`cv_varsel()`](https://mc-stan.org/projpred/dev/reference/cv_varsel.md)
  (with `cv_method = "LOO"`), the PSIS weights are now calculated based
  on the reference model (they used to be calculated based on the
  submodels which is incorrect). (GitHub:
  [\#325](https://github.com/stan-dev/projpred/issues/325))
- Some long-standing severe bugs (GitHub issues
  [\#329](https://github.com/stan-dev/projpred/issues/329),
  [\#330](https://github.com/stan-dev/projpred/issues/330), and
  [\#342](https://github.com/stan-dev/projpred/issues/342)) have been
  fixed, concerning the performance evaluation of models with nontrivial
  observation weights (i.e., models where at least one observation had a
  weight differing from 1). Concerned performance statistics were
  `"mse"`, `"rmse"`, `"acc"` (= `"pctcorr"`), and `"auc"` (i.e., all
  performance statistics except for `"elpd"` and `"mlpd"`).
- [`plot.vsel()`](https://mc-stan.org/projpred/dev/reference/plot.vsel.md)
  and
  [`suggest_size()`](https://mc-stan.org/projpred/dev/reference/suggest_size.md)
  gain a new argument `thres_elpd`. By default, this argument doesn’t
  have any impact, but a non-`NA` value can be used for a customized
  model size selection rule (see
  [`?suggest_size`](https://mc-stan.org/projpred/dev/reference/suggest_size.md)
  for details). (GitHub:
  [\#335](https://github.com/stan-dev/projpred/issues/335))

### Minor changes

- Several improvements in the documentation (especially in the
  explanation of the
  [`suggest_size()`](https://mc-stan.org/projpred/dev/reference/suggest_size.md)
  heuristic).
- Improvement of the numerical stability for some link functions,
  achieved by avoiding unnecessary back-and-forth transformations
  between latent space and response space. (GitHub:
  [\#337](https://github.com/stan-dev/projpred/issues/337),
  [\#338](https://github.com/stan-dev/projpred/issues/338))
- All arguments `seed` and `.seed` are now allowed to be `NA` for not
  calling [`set.seed()`](https://rdrr.io/r/base/Random.html) internally
  at all.
- Argument `d_test` of
  [`varsel()`](https://mc-stan.org/projpred/dev/reference/varsel.md) is
  not considered as an internal feature anymore. This was possible after
  fixing a bug for `d_test` (see below). (GitHub:
  [\#341](https://github.com/stan-dev/projpred/issues/341))
- The order of the observations in the sub-elements of
  `<vsel_object>$summaries` and `<vsel_object>$d_test` now corresponds
  to the order of the observations in the original dataset if
  `<vsel_object>` was created by a call to
  `cv_varsel(<...>, cv_method = "kfold")` (formerly, in that case, the
  observations in those sub-elements were ordered by fold). Thereby, the
  order of the observations in those sub-elements now always corresponds
  to the order of the observations in the original dataset, except if
  `<vsel_object>` was created by a call to
  `varsel(<...>, d_test = <non-NULL_d_test_object>)`, in which case the
  order of the observations in those sub-elements corresponds to the
  order of the observations in `<non-NULL_d_test_object>`. (GitHub:
  [\#341](https://github.com/stan-dev/projpred/issues/341))

### Bug fixes

- Fix GitHub issue
  [\#324](https://github.com/stan-dev/projpred/issues/324) (large
  `search_terms` caused the R session to crash).
- Fix GitHub issue
  [\#204](https://github.com/stan-dev/projpred/issues/204). (GitHub:
  [\#325](https://github.com/stan-dev/projpred/issues/325))
- Fix the `validate_search = FALSE` bug described above in “Major
  changes”: The PSIS weights are now calculated based on the reference
  model (they used to be calculated based on the submodels which is
  incorrect). (GitHub:
  [\#325](https://github.com/stan-dev/projpred/issues/325))
- Fix `\mbox{}` commands displayed incorrectly in the HTML help from R
  version 4.2.0 on. (GitHub:
  [\#326](https://github.com/stan-dev/projpred/issues/326))
- Fix GitHub issue
  [\#329](https://github.com/stan-dev/projpred/issues/329) (see also
  “Major changes” above).
- Fix GitHub issue
  [\#331](https://github.com/stan-dev/projpred/issues/331).
- [`plot.vsel()`](https://mc-stan.org/projpred/dev/reference/plot.vsel.md)
  now draws the dashed red horizontal line for the reference model
  (and—if present—the dotted black horizontal line for the baseline
  model) first (i.e., before the submodel-specific graphical elements),
  to avoid overplotting.
- Fix GitHub issue
  [\#339](https://github.com/stan-dev/projpred/issues/339). (GitHub:
  [\#340](https://github.com/stan-dev/projpred/issues/340))
- Fix argument `d_test` of
  [`varsel()`](https://mc-stan.org/projpred/dev/reference/varsel.md):
  Not only the predictive performance of the *reference model* needs to
  be evaluated on the test data, but also the predictive performance of
  the *submodels*. (GitHub:
  [\#341](https://github.com/stan-dev/projpred/issues/341))
- Fix GitHub issue
  [\#342](https://github.com/stan-dev/projpred/issues/342) (see also
  “Major changes” above).
- Fix GitHub issue
  [\#330](https://github.com/stan-dev/projpred/issues/330) (see also
  “Major changes” above). (GitHub:
  [\#344](https://github.com/stan-dev/projpred/issues/344), commit
  23e7101)

## projpred 2.1.2

CRAN release: 2022-05-13

### Minor changes

- Account for changes concerning the handling of offsets in **rstanarm**
  version 2.21.3. In particular, issue
  [stan-dev/rstanarm#542](https://github.com/stan-dev/rstanarm/issues/542)
  was fixed in **rstanarm** 2.21.3.
- Show the output of the vignette on CRAN.
- In the vignette, use
  [`cv_varsel()`](https://mc-stan.org/projpred/dev/reference/cv_varsel.md)
  with LOO CV and `validate_search = FALSE` instead of K-fold CV.
  (GitHub: [\#305](https://github.com/stan-dev/projpred/issues/305))
- Improve the documentation for argument `search_terms` of
  [`varsel()`](https://mc-stan.org/projpred/dev/reference/varsel.md) and
  [`cv_varsel()`](https://mc-stan.org/projpred/dev/reference/cv_varsel.md).
  (GitHub: [\#155](https://github.com/stan-dev/projpred/issues/155),
  [\#308](https://github.com/stan-dev/projpred/issues/308))
- In case of user-specified (non-`NULL`) `search_terms`, `method = NULL`
  is internally changed to `method = "forward"` and `method = "L1"`
  throws a warning. This is done because `search_terms` only takes
  effect in case of a forward search. (GitHub:
  [\#155](https://github.com/stan-dev/projpred/issues/155),
  [\#308](https://github.com/stan-dev/projpred/issues/308))
- Internally, the intercept is now always included in `search_terms`.
  This is necessary to prevent a bug described below. (GitHub:
  [\#308](https://github.com/stan-dev/projpred/issues/308))
- When fitting multilevel submodels via **lme4**, **projpred** now tries
  to handle `PIRLS loop resulted in NaN value` errors automatically.
  (GitHub: [\#314](https://github.com/stan-dev/projpred/issues/314))
- The fix for GitHub issue
  [\#320](https://github.com/stan-dev/projpred/issues/320) (see below)
  required to rename argument `b` of `projpred:::bootstrap()` to `B`.

### Bug fixes

- Throw a more informative error message in case of special group-level
  terms which are currently not supported (in particular, nested ones).
- Previously, using a `search_terms` vector which excluded the intercept
  in conjunction with `refit_prj = FALSE` (the latter in
  [`project()`](https://mc-stan.org/projpred/dev/reference/project.md),
  [`varsel()`](https://mc-stan.org/projpred/dev/reference/varsel.md), or
  [`cv_varsel()`](https://mc-stan.org/projpred/dev/reference/cv_varsel.md))
  led to incorrect submodels being fetched from the search or to an
  error while doing so. This has been fixed now by internally forcing
  the inclusion of the intercept in `search_terms`. (GitHub:
  [\#308](https://github.com/stan-dev/projpred/issues/308))
- Fix GitHub issues
  [\#147](https://github.com/stan-dev/projpred/issues/147) and
  [\#202](https://github.com/stan-dev/projpred/issues/202). (GitHub:
  [\#312](https://github.com/stan-dev/projpred/issues/312))
- Fix GitHub issue
  [\#320](https://github.com/stan-dev/projpred/issues/320). (GitHub:
  [\#321](https://github.com/stan-dev/projpred/issues/321))

## projpred 2.1.1

CRAN release: 2022-04-03

### Bug fixes

- Fix the order of the package authors.
- Fix failing CRAN checks.
- Add an input check for argument `solution_terms` of
  [`project()`](https://mc-stan.org/projpred/dev/reference/project.md)
  to fix a test failure in R versions \>= 4.2.

## projpred 2.1.0

CRAN release: 2022-04-01

### Major changes

- Added support for weighted LOO proportional-to-size subsampling based
  on [Magnusson et
  al. (2019)](https://proceedings.mlr.press/v97/magnusson19a.html).
  However, subsampled PSIS-LOO CV is currently regarded as experimental.
  Therefore, a corresponding warning is thrown when calling
  [`cv_varsel()`](https://mc-stan.org/projpred/dev/reference/cv_varsel.md)
  with `nloo < n` where `n` denotes the number of observations. (GitHub:
  [\#94](https://github.com/stan-dev/projpred/issues/94),
  [\#252](https://github.com/stan-dev/projpred/issues/252), commit
  feea39e)
- Automatically explore both linear and smooths components in GAM
  models. This allows the user to gauge the impact of the smooth term
  against its linear counterpart.
- Fast approximate LOO computation for `validate_search = FALSE` in
  [`cv_varsel()`](https://mc-stan.org/projpred/dev/reference/cv_varsel.md).
- Formerly, the defaults for arguments `nclusters` (= `1`) and
  `nclusters_pred` (= `5`) of
  [`varsel()`](https://mc-stan.org/projpred/dev/reference/varsel.md) and
  [`cv_varsel()`](https://mc-stan.org/projpred/dev/reference/cv_varsel.md)
  were set internally (the user-visible defaults were `NULL`). Now,
  `nclusters` and `ndraws_pred` (note the `ndraws_pred`, not
  `nclusters_pred`) have non-`NULL` user-visible defaults of `20` and
  `400`, respectively. In general, this increases the runtime of these
  functions a lot. With respect to
  [`cv_varsel()`](https://mc-stan.org/projpred/dev/reference/cv_varsel.md),
  the new vignette (see
  [vignettes](https://mc-stan.org/projpred/articles/)) mentions two ways
  to quickly obtain some rough preliminary results which in general
  should not be used as final results, though: (i)
  [`varsel()`](https://mc-stan.org/projpred/dev/reference/varsel.md)
  and (ii)
  [`cv_varsel()`](https://mc-stan.org/projpred/dev/reference/cv_varsel.md)
  with `validate_search = FALSE` (which only takes effect for
  `cv_method = "LOO"`). (GitHub:
  [\#291](https://github.com/stan-dev/projpred/issues/291) and several
  commits beforehand, in particular bbd0f0a, babe031, 4ef95d3, and
  ce7d1e0)
- For
  [`proj_linpred()`](https://mc-stan.org/projpred/dev/reference/pred-projection.md)
  and
  [`proj_predict()`](https://mc-stan.org/projpred/dev/reference/pred-projection.md),
  arguments `nterms`, `ndraws`, and `seed` have been removed to allow
  the user to pass them to
  [`project()`](https://mc-stan.org/projpred/dev/reference/project.md).
  New arguments `filter_nterms`, `nresample_clusters`, and `.seed` have
  been introduced (see the documentation for details). (GitHub:
  [\#92](https://github.com/stan-dev/projpred/issues/92),
  [\#135](https://github.com/stan-dev/projpred/issues/135))
- Reference models lacking an intercept are not supported anymore
  (actually, the previous implementation for such models was
  incomplete). Support might be re-introduced in the future (when
  fixed), but for now it is withdrawn as it requires some larger
  changes. (GitHub:
  [\#124](https://github.com/stan-dev/projpred/issues/124), but see also
  [\#96](https://github.com/stan-dev/projpred/issues/96) and
  [\#100](https://github.com/stan-dev/projpred/issues/100))
- In the output of
  [`proj_linpred()`](https://mc-stan.org/projpred/dev/reference/pred-projection.md),
  dimensions are not dropped anymore (i.e., output elements `pred` and
  `lpd` are always S x N matrices now). (GitHub:
  [\#143](https://github.com/stan-dev/projpred/issues/143))
- In case of `integrated = TRUE`,
  [`proj_linpred()`](https://mc-stan.org/projpred/dev/reference/pred-projection.md)
  now averages the LPD (across the projected posterior draws) instead of
  taking the LPD at the averaged linear predictors. (GitHub:
  [\#143](https://github.com/stan-dev/projpred/issues/143))
- If `newdata` does not contain the response variable,
  [`proj_linpred()`](https://mc-stan.org/projpred/dev/reference/pred-projection.md)
  now returns `NULL` for output element `lpd`. (GitHub:
  [\#143](https://github.com/stan-dev/projpred/issues/143))
- The fix for the offset issues (listed below under “Bug fixes”)
  requires reference model fits of class `stanreg` (from package
  **rstanarm**) with offsets to have these offsets specified via an
  [`offset()`](https://rdrr.io/r/stats/offset.html) term in the model
  formula (and not via argument `offset`).
- Improved handling of errors when fitting multilevel submodels.
  (GitHub: [\#201](https://github.com/stan-dev/projpred/issues/201))
- Some defaults have been changed from `NULL` to a user-visible value
  (and `NULL` is not allowed anymore).
- Argument `data` of
  [`get_refmodel.stanreg()`](https://mc-stan.org/projpred/dev/reference/refmodel-init-get.md)
  has been removed. (GitHub:
  [\#219](https://github.com/stan-dev/projpred/issues/219))
- The function passed to argument `div_minimizer` of
  [`init_refmodel()`](https://mc-stan.org/projpred/dev/reference/refmodel-init-get.md)
  now always needs to return a `list` of submodels (see the
  documentation for details). Correspondingly, the function passed to
  argument `proj_predfun` of
  [`init_refmodel()`](https://mc-stan.org/projpred/dev/reference/refmodel-init-get.md)
  can now always expect a `list` as input for argument `fits` (see the
  documentation for details). (GitHub:
  [\#230](https://github.com/stan-dev/projpred/issues/230))
- The function passed to argument `proj_predfun` of
  [`init_refmodel()`](https://mc-stan.org/projpred/dev/reference/refmodel-init-get.md)
  now always needs to return a matrix (see the documentation for
  details). (GitHub:
  [\#230](https://github.com/stan-dev/projpred/issues/230))
- The projection can be run in parallel now. However, we cannot
  recommend this for all kinds of platforms and all kinds of models. For
  more information, see the general package documentation available at
  [`` ?`projpred-package` ``](https://mc-stan.org/projpred/dev/reference/projpred-package.md).
  (GitHub: [\#235](https://github.com/stan-dev/projpred/issues/235))
- Support for the
  [`Student_t()`](https://mc-stan.org/projpred/dev/reference/extra-families.md)
  family is regarded as experimental. Therefore, a corresponding warning
  is thrown when creating the reference model. (GitHub:
  [\#233](https://github.com/stan-dev/projpred/issues/233),
  [\#252](https://github.com/stan-dev/projpred/issues/252))
- Support for additive models (i.e., GAMs and GAMMs) is regarded as
  experimental. Therefore, a corresponding warning is thrown when
  creating the reference model. (GitHub:
  [\#237](https://github.com/stan-dev/projpred/issues/237),
  [\#252](https://github.com/stan-dev/projpred/issues/252))
- Support for the [`Gamma()`](https://rdrr.io/r/stats/family.html)
  family is regarded as experimental. Therefore, a corresponding warning
  is thrown when creating the reference model. (GitHub:
  [paul-buerkner/brms#1255](https://github.com/paul-buerkner/brms/issues/1255),
  [\#240](https://github.com/stan-dev/projpred/issues/240),
  [\#252](https://github.com/stan-dev/projpred/issues/252))
- The previous behavior of
  [`init_refmodel()`](https://mc-stan.org/projpred/dev/reference/refmodel-init-get.md)
  in case of argument `dis` being `NULL` (the default) was dangerous for
  custom reference models with a `family` having a dispersion parameter
  (in that case, `dis` values of all-zeros were used silently). The new
  behavior now requires a non-`NULL` argument `dis` in that case.
  (GitHub: [\#254](https://github.com/stan-dev/projpred/issues/254))
- Argument `cv_search` has been renamed to `refit_prj`. (GitHub:
  [\#154](https://github.com/stan-dev/projpred/issues/154),
  [\#265](https://github.com/stan-dev/projpred/issues/265))
- [`as.matrix.projection()`](https://mc-stan.org/projpred/dev/reference/as.matrix.projection.md)
  has gained a new argument `nm_scheme` which allows to choose the
  naming scheme for the column names of the returned matrix. The default
  (`"auto"`) follows the naming scheme of the reference model fit (and
  uses the `"rstanarm"` naming scheme if the reference model fit is of
  an unknown class). (GitHub:
  [\#82](https://github.com/stan-dev/projpred/issues/82),
  [\#279](https://github.com/stan-dev/projpred/issues/279))
- `seed` (and `.seed`) arguments now have a default of
  `sample.int(.Machine$integer.max, 1)` instead of `NULL`. Furthermore,
  the value supplied to these arguments is now used to generate new
  seeds internally on-the-fly. In many cases, this will change results
  compared to older **projpred** versions. Also note that now, the
  internal seeds are never fixed to a specific value if `seed` (and
  `.seed`) arguments are set to `NULL`. (GitHub:
  [\#84](https://github.com/stan-dev/projpred/issues/84),
  [\#286](https://github.com/stan-dev/projpred/issues/286))

### Minor changes

- Improved summary output with important details.
- For group-level effects, the
  [`as.matrix.projection()`](https://mc-stan.org/projpred/dev/reference/as.matrix.projection.md)
  method now also returns the estimated group-level effects themselves.
  (GitHub: [\#75](https://github.com/stan-dev/projpred/issues/75))
- For group-level effects, the
  [`as.matrix.projection()`](https://mc-stan.org/projpred/dev/reference/as.matrix.projection.md)
  method now returns the variance components (population SD(s) and
  population correlation(s)) instead of the empirical SD(s) of the
  group-level effects. (GitHub:
  [\#74](https://github.com/stan-dev/projpred/issues/74))
- Improved documentation. (GitHub: especially
  [\#233](https://github.com/stan-dev/projpred/issues/233))
- Replaced the two vignettes by a single one which also has new content.
  (GitHub: [\#237](https://github.com/stan-dev/projpred/issues/237))
- Updated the `README` file. (GitHub:
  [\#245](https://github.com/stan-dev/projpred/issues/245))
- Some error and warning messages have been improved and added. (GitHub:
  especially [\#219](https://github.com/stan-dev/projpred/issues/219),
  [\#221](https://github.com/stan-dev/projpred/issues/221),
  [\#223](https://github.com/stan-dev/projpred/issues/223),
  [\#252](https://github.com/stan-dev/projpred/issues/252),
  [\#263](https://github.com/stan-dev/projpred/issues/263))
- For K-fold cross-validation, an internally hard-coded value of 5 for
  `nclusters_pred` was removed. (GitHub: commit 5062f2f)
- Throw a proper error message for unsupported families. (GitHub:
  [\#140](https://github.com/stan-dev/projpred/issues/140))
- Show the README also on the CRAN website. (GitHub:
  [\#140](https://github.com/stan-dev/projpred/issues/140))
- [`project()`](https://mc-stan.org/projpred/dev/reference/project.md):
  Warn if elements of `solution_terms` are not found in the reference
  model (and therefore ignored). (GitHub:
  [\#140](https://github.com/stan-dev/projpred/issues/140))
- [`get_refmodel.default()`](https://mc-stan.org/projpred/dev/reference/refmodel-init-get.md)
  now passes arguments via the ellipsis (`...`) to
  [`init_refmodel()`](https://mc-stan.org/projpred/dev/reference/refmodel-init-get.md).
  (GitHub: [\#153](https://github.com/stan-dev/projpred/issues/153),
  commit dd3716e)
- Remove dependency on package **rngtools** (version 2.0.0 of
  **projpred** re-introduced this dependency after it was already
  removed in version 1.1.2). (GitHub:
  [\#189](https://github.com/stan-dev/projpred/issues/189))
- [`init_refmodel()`](https://mc-stan.org/projpred/dev/reference/refmodel-init-get.md):
  The default (`NULL`) for argument `extract_model_data` has been
  removed as it wasn’t meaningful anyway. (GitHub:
  [\#219](https://github.com/stan-dev/projpred/issues/219))
- Argument `folds` of
  [`init_refmodel()`](https://mc-stan.org/projpred/dev/reference/refmodel-init-get.md)
  has been removed as it was effectively unused. (GitHub:
  [\#220](https://github.com/stan-dev/projpred/issues/220))
- Use the S3 system for
  [`solution_terms()`](https://mc-stan.org/projpred/dev/reference/solution_terms.md).
  This allowed the introduction of a
  [`solution_terms.projection()`](https://mc-stan.org/projpred/dev/reference/solution_terms.md)
  method. (GitHub:
  [\#223](https://github.com/stan-dev/projpred/issues/223))
- [`predict.refmodel()`](https://mc-stan.org/projpred/dev/reference/predict.refmodel.md)
  now uses a default of `newdata = NULL`. (GitHub:
  [\#223](https://github.com/stan-dev/projpred/issues/223))
- Argument `weights` of
  [`init_refmodel()`](https://mc-stan.org/projpred/dev/reference/refmodel-init-get.md)’s
  argument `proj_predfun` has been removed. (GitHub:
  [\#163](https://github.com/stan-dev/projpred/issues/163),
  [\#224](https://github.com/stan-dev/projpred/issues/224))
- **projpred**’s internal `div_minimizer` functions have been unified
  into a single `div_minimizer` which chooses an appropriate submodel
  fitter based on the formula of the submodel, not based on that of the
  reference model. Furthermore, the automatic handling of errors in the
  submodel fitters has been improved. (GitHub:
  [\#230](https://github.com/stan-dev/projpred/issues/230))
- Improve the axis labels in
  [`plot.vsel()`](https://mc-stan.org/projpred/dev/reference/plot.vsel.md).
  (GitHub: [\#234](https://github.com/stan-dev/projpred/issues/234),
  [\#270](https://github.com/stan-dev/projpred/issues/270))
- Handle **rstanarm**’s GitHub issue
  [\#551](https://github.com/stan-dev/projpred/issues/551). This implies
  that **projpred**’s default `cvfun` for `stanreg` fits will now always
  use *inner* parallelization in
  [`rstanarm::kfold.stanreg()`](https://mc-stan.org/rstanarm/reference/kfold.stanreg.html)
  (i.e., across chains, not across CV folds), with
  `getOption("mc.cores", 1)` cores. We do so on all systems (not only
  Windows). (GitHub:
  [\#249](https://github.com/stan-dev/projpred/issues/249))
- Argument `fit` of
  [`init_refmodel()`](https://mc-stan.org/projpred/dev/reference/refmodel-init-get.md)’s
  argument `proj_predfun` was renamed to `fits`. This is a non-breaking
  change since all calls to `proj_predfun` in **projpred** have that
  argument unnamed. However, this cannot be guaranteed in the future, so
  we strongly encourage users with a custom `proj_predfun` to rename
  argument `fit` to `fits`. (GitHub:
  [\#263](https://github.com/stan-dev/projpred/issues/263))
- [`init_refmodel()`](https://mc-stan.org/projpred/dev/reference/refmodel-init-get.md)
  has gained argument `cvrefbuilder` which may be a custom function for
  constructing the K reference models in a K-fold CV. (GitHub:
  [\#271](https://github.com/stan-dev/projpred/issues/271))
- Allow arguments to be passed from
  [`project()`](https://mc-stan.org/projpred/dev/reference/project.md),
  [`varsel()`](https://mc-stan.org/projpred/dev/reference/varsel.md),
  and
  [`cv_varsel()`](https://mc-stan.org/projpred/dev/reference/cv_varsel.md)
  to the divergence minimizer. (GitHub:
  [\#278](https://github.com/stan-dev/projpred/issues/278))
- In
  [`init_refmodel()`](https://mc-stan.org/projpred/dev/reference/refmodel-init-get.md),
  any `contrasts` attributes of the dataset’s columns are silently
  removed. (GitHub:
  [\#284](https://github.com/stan-dev/projpred/issues/284))
- `NA`s in data supplied to `newdata` arguments now trigger an error.
  (GitHub: [\#285](https://github.com/stan-dev/projpred/issues/285))

### Bug fixes

- Fixed a bug in
  [`as.matrix.projection()`](https://mc-stan.org/projpred/dev/reference/as.matrix.projection.md)
  (causing incorrect column names for the returned matrix). (GitHub:
  [\#72](https://github.com/stan-dev/projpred/issues/72),
  [\#73](https://github.com/stan-dev/projpred/issues/73))
- Fixed a bug raising an error when not projecting from a `vsel` object.
  (GitHub: [\#79](https://github.com/stan-dev/projpred/issues/79),
  [\#80](https://github.com/stan-dev/projpred/issues/80))
- Fixed a bug in the calculation of the Gaussian deviance. (GitHub:
  [\#81](https://github.com/stan-dev/projpred/issues/81))
- Fixed a bug in the calculation of the predictive statistics of the
  reference model on test data in
  [`varsel()`](https://mc-stan.org/projpred/dev/reference/varsel.md).
  (GitHub [\#90](https://github.com/stan-dev/projpred/issues/90))
- Fixed a bug in an input check for argument `nloo` of
  [`cv_varsel()`](https://mc-stan.org/projpred/dev/reference/cv_varsel.md).
  (GitHub: [\#93](https://github.com/stan-dev/projpred/issues/93))
- Fixed a bug in
  [`cv_varsel()`](https://mc-stan.org/projpred/dev/reference/cv_varsel.md),
  causing an error in case of `!validate_search && cv_method != "LOO"`.
  (GitHub: [\#95](https://github.com/stan-dev/projpred/issues/95))
- Fixed bugs related to the setting of the seed. (GitHub: commit
  02cd50d)
- Fixed a bug causing
  [`proj_linpred()`](https://mc-stan.org/projpred/dev/reference/pred-projection.md)
  to raise an error if argument `newdata` was `NULL`. (GitHub:
  [\#97](https://github.com/stan-dev/projpred/issues/97))
- Fixed an incorrect usage of the dispersion parameter values when
  calculating output element `lpd` in
  [`proj_linpred()`](https://mc-stan.org/projpred/dev/reference/pred-projection.md)
  (for `integrated = TRUE` as well as for `integrated = FALSE`).
  (GitHub: [\#105](https://github.com/stan-dev/projpred/issues/105))
- Fixed bugs in
  [`proj_linpred()`](https://mc-stan.org/projpred/dev/reference/pred-projection.md)’s
  calculation of output element `lpd` (for `integrated = TRUE`).
  (GitHub: [\#106](https://github.com/stan-dev/projpred/issues/106),
  [\#112](https://github.com/stan-dev/projpred/issues/112))
- Fixed an inconsistency in the dimensions of
  [`proj_linpred()`](https://mc-stan.org/projpred/dev/reference/pred-projection.md)’s
  output elements `pred` and `lpd` (for `integrated = FALSE`): Now, they
  are both S x N matrices, with S denoting the number of (possibly
  clustered) posterior draws and N denoting the number of observations.
  (GitHub: [\#107](https://github.com/stan-dev/projpred/issues/107),
  [\#112](https://github.com/stan-dev/projpred/issues/112))
- Fixed a bug causing
  [`proj_predict()`](https://mc-stan.org/projpred/dev/reference/pred-projection.md)’s
  output matrix to be transposed in case of `nrow(newdata) == 1`.
  (GitHub: [\#112](https://github.com/stan-dev/projpred/issues/112))
- Fixed a bug when using weights or offsets e.g. in
  [`proj_linpred()`](https://mc-stan.org/projpred/dev/reference/pred-projection.md).
  (GitHub: [\#114](https://github.com/stan-dev/projpred/issues/114))
- Fixed a bug causing
  [`varsel()`](https://mc-stan.org/projpred/dev/reference/varsel.md)/`make_formula`
  to fail with multidimensional interaction terms. (GitHub:
  [\#102](https://github.com/stan-dev/projpred/issues/102),
  [\#103](https://github.com/stan-dev/projpred/issues/103))
- Fixed an indexing bug in
  [`cv_varsel()`](https://mc-stan.org/projpred/dev/reference/cv_varsel.md)
  for models with a single predictor. (GitHub:
  [\#115](https://github.com/stan-dev/projpred/issues/115))
- Fixed bugs for argument `nterms` of
  [`proj_linpred()`](https://mc-stan.org/projpred/dev/reference/pred-projection.md)
  and
  [`proj_predict()`](https://mc-stan.org/projpred/dev/reference/pred-projection.md).
  (GitHub: [\#110](https://github.com/stan-dev/projpred/issues/110))
- Fixed an inconsistency for some intercept-only submodels. (GitHub:
  [\#119](https://github.com/stan-dev/projpred/issues/119))
- Fix a bug for
  [`as.matrix.projection()`](https://mc-stan.org/projpred/dev/reference/as.matrix.projection.md)
  in case of 1 (clustered) draw after projection. (GitHub:
  [\#130](https://github.com/stan-dev/projpred/issues/130))
- For submodels of class `subfit`, make the column names of
  [`as.matrix.projection()`](https://mc-stan.org/projpred/dev/reference/as.matrix.projection.md)’s
  output matrix consistent with other classes of submodels. (GitHub:
  [\#132](https://github.com/stan-dev/projpred/issues/132))
- Fix a bug for argument `nterms_max` of
  [`plot.vsel()`](https://mc-stan.org/projpred/dev/reference/plot.vsel.md)
  if there is just the intercept-only submodel. (GitHub:
  [\#138](https://github.com/stan-dev/projpred/issues/138))
- Throw an appropriate error message when trying to apply an L1 search
  to an empty (i.e. intercept-only) reference model. (GitHub:
  [\#139](https://github.com/stan-dev/projpred/issues/139))
- Fix the list names of element `search_path` in, e.g.,
  [`varsel()`](https://mc-stan.org/projpred/dev/reference/varsel.md)’s
  output. (GitHub:
  [\#140](https://github.com/stan-dev/projpred/issues/140))
- Fix a bug (error `unused argument`) when initializing the K reference
  models in a K-fold CV with CV fits not of class `brmsfit` or
  `stanreg`. (GitHub:
  [\#140](https://github.com/stan-dev/projpred/issues/140))
- In
  [`get_refmodel.default()`](https://mc-stan.org/projpred/dev/reference/refmodel-init-get.md),
  remove old defunct arguments `fetch_data`, `wobs`, and `offset`.
  (GitHub: [\#140](https://github.com/stan-dev/projpred/issues/140))
- Fix a bug in
  [`get_refmodel.stanreg()`](https://mc-stan.org/projpred/dev/reference/refmodel-init-get.md).
  (GitHub: [\#142](https://github.com/stan-dev/projpred/issues/142),
  [\#184](https://github.com/stan-dev/projpred/issues/184))
- Fix a possible bug related to `extract_model_data()`’s argument
  `extract_y` in
  [`get_refmodel.default()`](https://mc-stan.org/projpred/dev/reference/refmodel-init-get.md).
  (GitHub: [\#153](https://github.com/stan-dev/projpred/issues/153),
  commit 39fece8)
- Fix a possible bug related to `extract_model_data()` in K-fold CV.
  (GitHub: [\#153](https://github.com/stan-dev/projpred/issues/153),
  commit 4f32195)
- Fix GitHub issue
  [\#161](https://github.com/stan-dev/projpred/issues/161).
- Fix GitHub issue
  [\#162](https://github.com/stan-dev/projpred/issues/162).
- Fix GitHub issue
  [\#164](https://github.com/stan-dev/projpred/issues/164).
- Fix GitHub issue
  [\#160](https://github.com/stan-dev/projpred/issues/160).
- Fix GitHub issue
  [\#159](https://github.com/stan-dev/projpred/issues/159).
- Fix GitHub issue
  [\#158](https://github.com/stan-dev/projpred/issues/158).
- Fix GitHub issue
  [\#157](https://github.com/stan-dev/projpred/issues/157).
- Fix GitHub issue
  [\#144](https://github.com/stan-dev/projpred/issues/144).
- Fix GitHub issue
  [\#146](https://github.com/stan-dev/projpred/issues/146).
- Fix GitHub issue
  [\#169](https://github.com/stan-dev/projpred/issues/169).
- Fix GitHub issue
  [\#167](https://github.com/stan-dev/projpred/issues/167).
- Fix a bug in the default `proj_predfun()` for GLMMs. (GitHub:
  [\#174](https://github.com/stan-dev/projpred/issues/174))
- Fix GitHub issue
  [\#171](https://github.com/stan-dev/projpred/issues/171).
- Fix GitHub issue
  [\#172](https://github.com/stan-dev/projpred/issues/172).
- Fix a bug in the default `proj_predfun()` for `datafit`s. (GitHub:
  [\#177](https://github.com/stan-dev/projpred/issues/177))
- Fix the names of `summary.vsel()$selection` for objects of class
  `vsel` created by
  [`varsel()`](https://mc-stan.org/projpred/dev/reference/varsel.md).
  (GitHub: [\#179](https://github.com/stan-dev/projpred/issues/179))
- Fix forward search when `search_terms` are not consecutive in size.
  (GitHub: commit 34e24de)
- Fix a bug in `cv_varsel()$pct_solution_terms_cv`. (GitHub:
  [\#188](https://github.com/stan-dev/projpred/issues/188), commit
  e529ec1)
- Fix GitHub issue
  [\#185](https://github.com/stan-dev/projpred/issues/185). (GitHub:
  [\#193](https://github.com/stan-dev/projpred/issues/193),
  [\#194](https://github.com/stan-dev/projpred/issues/194))
- Fix a bug in forward searches with interaction terms. (GitHub:
  [\#191](https://github.com/stan-dev/projpred/issues/191))
- Fix offset issues. (GitHub:
  [\#196](https://github.com/stan-dev/projpred/issues/196),
  [\#203](https://github.com/stan-dev/projpred/issues/203),
  [\#228](https://github.com/stan-dev/projpred/issues/228))
- Fix a bug in `glm_elnet()` (the workhorse for L1 search), causing the
  grid for lambda to be constructed without taking observation weights
  into account. (GitHub:
  [\#198](https://github.com/stan-dev/projpred/issues/198); note that
  the second part of
  [\#198](https://github.com/stan-dev/projpred/issues/198) did not have
  any consequences for users)
- Fix GitHub issue
  [\#136](https://github.com/stan-dev/projpred/issues/136). (GitHub:
  [\#221](https://github.com/stan-dev/projpred/issues/221))
- Fix a bug in
  [`print.vsel()`](https://mc-stan.org/projpred/dev/reference/print.vsel.md)
  causing argument `digits` to be ignored. (GitHub:
  [\#222](https://github.com/stan-dev/projpred/issues/222))
- Fix a bug causing the default of argument `cv_search` in
  [`varsel()`](https://mc-stan.org/projpred/dev/reference/varsel.md) and
  [`cv_varsel()`](https://mc-stan.org/projpred/dev/reference/cv_varsel.md)
  to be `TRUE` for `datafit`s, although it should be `FALSE` in that
  case. (GitHub:
  [\#223](https://github.com/stan-dev/projpred/issues/223))
- Fix a bug
  (`Error: Levels '<...>' of grouping factor '<...>' cannot be found in the fitted model. Consider setting argument 'allow_new_levels' to TRUE.`)
  when predicting from submodels which are GLMMs for `newdata`
  containing new levels for grouping factors. (GitHub:
  [\#223](https://github.com/stan-dev/projpred/issues/223))
- [`predict.refmodel()`](https://mc-stan.org/projpred/dev/reference/predict.refmodel.md):
  Fix a bug for integer `ynew`. (GitHub:
  [\#223](https://github.com/stan-dev/projpred/issues/223))
- [`predict.refmodel()`](https://mc-stan.org/projpred/dev/reference/predict.refmodel.md):
  Fix input checks for `offsetnew` and `weightsnew`. (GitHub:
  [\#223](https://github.com/stan-dev/projpred/issues/223))
- After all calls to `extract_model_data()`, the weights and offsets are
  now checked if they are of length 0 (and if yes, then they are set to
  vectors of ones and zeros, respectively). This is important for
  `extract_model_data()` functions which return weights and offsets of
  length 0 (see, e.g., `brms` version \<= 2.16.1). (GitHub:
  [\#223](https://github.com/stan-dev/projpred/issues/223))
- Handle **rstanarm**’s GitHub issue
  [\#546](https://github.com/stan-dev/projpred/issues/546). (GitHub:
  [\#227](https://github.com/stan-dev/projpred/issues/227))
- Fix a bug causing the internal submodel fitter for GLMMs to not pass
  arguments `var` (the predictive variances) and `regul` (amount of
  ridge regularization) to the internal submodel fitter for GLMs.
  (GitHub: [\#230](https://github.com/stan-dev/projpred/issues/230))
- Fix GitHub issue
  [\#210](https://github.com/stan-dev/projpred/issues/210). (GitHub:
  [\#234](https://github.com/stan-dev/projpred/issues/234))
- Fix GitHub issue
  [\#242](https://github.com/stan-dev/projpred/issues/242). (GitHub:
  [\#253](https://github.com/stan-dev/projpred/issues/253))
- Fix GitHub issue
  [\#244](https://github.com/stan-dev/projpred/issues/244). (GitHub:
  [\#255](https://github.com/stan-dev/projpred/issues/255))
- Fix GitHub issue
  [\#243](https://github.com/stan-dev/projpred/issues/243). (GitHub:
  [\#262](https://github.com/stan-dev/projpred/issues/262))
- Fix GitHub issue
  [\#213](https://github.com/stan-dev/projpred/issues/213). (GitHub:
  [\#264](https://github.com/stan-dev/projpred/issues/264))
- Fix GitHub issue
  [\#215](https://github.com/stan-dev/projpred/issues/215). (GitHub:
  [\#266](https://github.com/stan-dev/projpred/issues/266))
- Fix GitHub issue
  [\#212](https://github.com/stan-dev/projpred/issues/212). (GitHub:
  [\#267](https://github.com/stan-dev/projpred/issues/267))
- Fix GitHub issue
  [\#156](https://github.com/stan-dev/projpred/issues/156). (GitHub:
  [\#269](https://github.com/stan-dev/projpred/issues/269))
- If the data used for the reference model contains `NA`s, an
  appropriate error is now thrown. Previously, the reference model was
  created successfully, but this caused opaque errors in downstream code
  such as
  [`project()`](https://mc-stan.org/projpred/dev/reference/project.md).
  (GitHub: [\#274](https://github.com/stan-dev/projpred/issues/274))
- Fix GitHub issue
  [\#268](https://github.com/stan-dev/projpred/issues/268). (GitHub:
  [\#287](https://github.com/stan-dev/projpred/issues/287))
- Fix GitHub issue
  [\#149](https://github.com/stan-dev/projpred/issues/149). (GitHub:
  [\#288](https://github.com/stan-dev/projpred/issues/288))

## projpred 2.0.2

CRAN release: 2020-10-28

We have fully rewritten the internals in several ways. Most importantly,
we now leverage maximum likelihood estimation to third parties depending
on the reference model’s family. This allows a lot of flexibility and
extensibility for various models. Functionality wise, the major updates
since the last release are:

- Added support for GLMMs and GAMMs via **lme4** and **gamm4**.
- Formula syntax support internally that allows for easier building upon
  projections.
- Thanks to the above point, we save some computation by only
  considering sensible projections during forward search instead of
  fitting every possible submodel.
- We have added a new argument `search_terms` that allows the user to
  specify custom unit building blocks of the projections. New vignette
  coming up.
- We have fully changed the way to define custom reference models. The
  user now provides projection fitting and prediction functions (more
  information in a new upcoming vignette).

## projpred 1.1.4

CRAN release: 2019-10-02

Better validation of function arguments.

## projpred 1.1.3

Added print methods for vsel and cvsel objects. Added AUC statistics for
binomial family. A few additional minor patches.

## projpred 1.1.2

CRAN release: 2019-05-31

Removed the dependency on the **rngtools** package.

## projpred 1.1.1

CRAN release: 2019-03-12

This version contains only a few patches, no new features to the user.

## projpred 1.1.0

CRAN release: 2018-10-23

### New features

- Added support for [**brms**](https://paulbuerkner.com/brms/) models.

### Bug fixes

- The program crashed with [**rstanarm**](https://mc-stan.org/rstanarm/)
  models fitted with syntax like `stan_glm(log(y) ~ log(x), ...)`, that
  is, it did not allow transformation for `y`.

## projpred 1.0.0

CRAN release: 2018-09-18

### New features and improvements

- Changed the internals so that now all fit objects (such as rstanarm
  fits) are converted to `refmodel`-objects using the generic
  `get_refmodel`-function, and all the functions use only this object.
  This makes it much easier to use projpred with other reference models
  by writing them a new `get_refmodel`-function. The syntax is now
  changed so that `varsel` and `cv_varsel` both return an object that
  has similar structure always, and the reference model is stored into
  this object.
- Added more examples to the vignette.
- Added possibility to change the baseline in `plot/summary`. Now it is
  possible to compare also to the best submodel found, not only to the
  reference model.
- Bug fix: RMSE was previously computed wrong, this is now fixed.
- Small changes: `nloo = n` by default in `cv_varsel`. `regul=1e-4` now
  by default in all functions.

## projpred 0.9.0

### New features and improvements

- Added the `cv_search` argument for the main functions
  (`varsel`,`cv_varsel`,`project` and the prediction functions). Now it
  is possible to make predictions also with those parameter estimates
  that were computed during the L1-penalized search. This change also
  allows the user to compute the Lasso-solution by providing the
  observed data as the ‘reference fit’ for init_refmodel. An example
  will be added to the vignette.

### Bug fixes

- The projection with a nonzero regularization parameter value did not
  produce exactly correct result, although the difference to the correct
  result was often so small that user would not see the difference.
  Fixed this.

## projpred 0.8.0 and earlier

CRAN release: 2018-04-16

Until this version, we did not keep record of the changes between
different versions. Started to do this from version 0.9.0 onwards.
