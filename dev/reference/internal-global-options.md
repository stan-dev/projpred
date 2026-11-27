# Internal global options

The following global options are for internal use:

- `projpred.mssg_ndraws`, `projpred.mssg_cut_search`,
  `projpred.mssg_time`, `projpred.mssg_smooth_linear`,
  `projpred.warn_wobs_ppd`, `projpred.warn_additive_experimental`,
  `projpred.warn_allrandom_dis`, `projpred.warn_instable_projections`,
  `projpred.warn_cvrefbuilder_NULL`, `projpred.warn_kfold_refits`: A
  single logical value indicating whether to throw certain messages or
  warnings (depending on the midfix `mssg` or `warn`, respectively). For
  the exact meaning of these global options, see their occurrences in
  the codebase. With the exception of `projpred.warn_allrandom_dis`,
  these global options are currently used in the unit tests to
  deactivate these messages and warnings.

- `projpred.additional_checks`: A single logical value indicating
  whether to run some additional checks that are not necessary to be run
  when users call the corresponding projpred functions. Currently, these
  checks are activated during the unit tests.

- `projpred.glm_fitter`: A character string naming the function to be
  used as the submodel fitter for non-multilevel, non-additive
  projections. Allowed values are `"fit_glm_ridge_callback"` (the
  default) and `"fit_glm_callback"`. For a reference model with offsets,
  projpred accepts only `"fit_glm_ridge_callback"` (except for the
  Gaussian family with the identity link, the latent projection and a
  custom `div_minimizer` in
  [`init_refmodel()`](https://mc-stan.org/projpred/dev/reference/refmodel-init-get.md)).

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
  [`lme4::lmer()`](https://rdrr.io/pkg/lme4/man/lmer.html)).

- `projpred.PQL`: A single logical value indicating whether to use
  [`MASS::glmmPQL()`](https://rdrr.io/pkg/MASS/man/glmmPQL.html) as the
  submodel fitter for multilevel (non-additive) projections. For a
  reference model with offsets, projpred accepts only `FALSE` (with the
  same exceptions as for `projpred.glm_fitter`).

- `projpred.search_gauss_crossprod`: A single logical value indicating
  whether to compute the forward search and the evaluation in the LOO CV
  with validated search for Gaussian submodels with the identity link
  (including the latent projection) from cross-products of the model
  matrix (default `TRUE`). With `FALSE`, projpred projects onto each
  candidate submodel draw by draw. The cross-products are used only with
  a positive `regul`.

- `projpred.categ_mlvl_mblogit`: A single logical value indicating
  whether to fit the multilevel submodels of the augmented-data
  projection for the
  [`brms::categorical()`](https://paulbuerkner.com/brms/reference/brmsfamily.html)
  family with
  [`mclogit::mblogit()`](https://melff.github.io/mclogit/reference/mblogit.html)
  (`TRUE`) or with projpred's own implementation of the PQL algorithm of
  [`mclogit::mblogit()`](https://melff.github.io/mclogit/reference/mblogit.html)
  (`FALSE`, the default).

- `projpred.nb_gsum`: A character string, `"quadrature"` (the default)
  or `"exact"`. It selects how projpred computes the aggregated
  predictive probability mass function of the reference model for the
  projection of the shape of the negative binomial family.
  `"quadrature"` sums the probabilities exactly for y \< 100 and uses a
  quadrature for larger y. `"exact"` sums the probabilities of all y up
  to an upper tail probability of 1e-17.
