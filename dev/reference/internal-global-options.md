# Internal global options

The following global options are for internal use:

- `projpred.mssg_ndraws`, `projpred.mssg_cut_search`,
  `projpred.mssg_time`, `projpred.warn_wobs_ppd`,
  `projpred.warn_additive_experimental`, `projpred.warn_allrandom_dis`,
  `projpred.warn_instable_projections`,
  `projpred.warn_cvrefbuilder_NULL`, `projpred.warn_kfold_refits`: A
  single logical value indicating whether to throw certain messages or
  warnings (depending on the midfix `mssg` or `warn`, respectively). For
  the exact meaning of these global options, see their occurrences in
  the codebase. With the exception of `projpred.warn_allrandom_dis`,
  these global options are currently used in the unit tests to
  deactivate these messages and warnings. Global option
  `projpred.warn_instable_projections` is also used (invisibly) in the
  latent vignette to suppress the corresponding warnings while
  illustrating the underlying issue (instable projections).

- `projpred.additional_checks`: A single logical value indicating
  whether to run some additional checks that are not necessary to be run
  when users call the corresponding projpred functions. Currently, these
  checks are activated during the unit tests.

- `projpred.glm_fitter`: A character string naming the function to be
  used as the submodel fitter for non-multilevel, non-additive
  projections. Currently, this is an experimental feature and allowed
  values are `"fit_glm_ridge_callback"` (the default) and
  `"fit_glm_callback"`.

- `projpred.gaussian_not_as_generalized`: A single logical value
  indicating whether to treat the
  [`gaussian()`](https://rdrr.io/r/stats/family.html) family not as a
  family for a *generalized linear* model (i.e., for which
  [`glm()`](https://rdrr.io/r/stats/glm.html) would typically be used as
  a model fitting function outside of projpred), but as the family for
  an explicit *linear* model (i.e., for which
  [`lm()`](https://rdrr.io/r/stats/lm.html) would typically be used as a
  model fitting function outside of projpred). This also holds for
  models with multilevel terms (because lme4 offers both
  [`lme4::glmer()`](https://rdrr.io/pkg/lme4/man/glmer.html) and
  [`lme4::lmer()`](https://rdrr.io/pkg/lme4/man/lmer.html)). Currently,
  this is an experimental feature.

- `projpred.PQL`: A single logical value indicating whether to use
  [`MASS::glmmPQL()`](https://rdrr.io/pkg/MASS/man/glmmPQL.html) as the
  submodel fitter for multilevel (non-additive) projections (see GitHub
  issue [\#207](https://github.com/stan-dev/projpred/issues/207) and
  GitHub pull request
  [\#353](https://github.com/stan-dev/projpred/pull/353)). Currently,
  this is an experimental feature.
