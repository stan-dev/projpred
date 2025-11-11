# Package index

## All functions

- [`as.matrix(`*`<projection>`*`)`](https://mc-stan.org/projpred/dev/reference/as.matrix.projection.md)
  : Extract projected parameter draws and coerce to matrix

- [`as_draws_matrix(`*`<projection>`*`)`](https://mc-stan.org/projpred/dev/reference/as_draws_matrix.projection.md)
  [`as_draws(`*`<projection>`*`)`](https://mc-stan.org/projpred/dev/reference/as_draws_matrix.projection.md)
  :

  Extract projected parameter draws and coerce to `draws_matrix` (see
  package posterior)

- [`augdat_ilink_binom()`](https://mc-stan.org/projpred/dev/reference/augdat_ilink_binom.md)
  : Inverse-link function for augmented-data projection with binomial
  family

- [`augdat_link_binom()`](https://mc-stan.org/projpred/dev/reference/augdat_link_binom.md)
  : Link function for augmented-data projection with binomial family

- [`break_up_matrix_term()`](https://mc-stan.org/projpred/dev/reference/break_up_matrix_term.md)
  : Break up matrix terms

- [`cl_agg()`](https://mc-stan.org/projpred/dev/reference/cl_agg.md) :
  Weighted averaging within clusters of parameter draws

- [`cv_folds()`](https://mc-stan.org/projpred/dev/reference/cv-indices.md)
  [`cvfolds()`](https://mc-stan.org/projpred/dev/reference/cv-indices.md)
  [`cv_ids()`](https://mc-stan.org/projpred/dev/reference/cv-indices.md)
  : Create cross-validation folds

- [`cv_proportions()`](https://mc-stan.org/projpred/dev/reference/cv_proportions.md)
  : Ranking proportions from fold-wise predictor rankings

- [`cv_varsel()`](https://mc-stan.org/projpred/dev/reference/cv_varsel.md)
  : Run search and performance evaluation with cross-validation

- [`df_binom`](https://mc-stan.org/projpred/dev/reference/df_binom.md) :
  Binomial toy example

- [`df_gaussian`](https://mc-stan.org/projpred/dev/reference/df_gaussian.md)
  : Gaussian toy example

- [`extend_family()`](https://mc-stan.org/projpred/dev/reference/extend_family.md)
  : Extend a family

- [`Student_t()`](https://mc-stan.org/projpred/dev/reference/extra-families.md)
  : Extra family objects

- [`force_search_terms()`](https://mc-stan.org/projpred/dev/reference/force_search_terms.md)
  : Force search terms

- [`mesquite`](https://mc-stan.org/projpred/dev/reference/mesquite.md) :
  Mesquite data set

- [`performances()`](https://mc-stan.org/projpred/dev/reference/performances.md)
  : Predictive performance results

- [`plot(`*`<cv_proportions>`*`)`](https://mc-stan.org/projpred/dev/reference/plot.cv_proportions.md)
  [`plot(`*`<ranking>`*`)`](https://mc-stan.org/projpred/dev/reference/plot.cv_proportions.md)
  : Plot ranking proportions from fold-wise predictor rankings

- [`plot(`*`<vsel>`*`)`](https://mc-stan.org/projpred/dev/reference/plot.vsel.md)
  : Plot predictive performance

- [`proj_linpred()`](https://mc-stan.org/projpred/dev/reference/pred-projection.md)
  [`proj_predict()`](https://mc-stan.org/projpred/dev/reference/pred-projection.md)
  : Predictions from a submodel (after projection)

- [`predict(`*`<refmodel>`*`)`](https://mc-stan.org/projpred/dev/reference/predict.refmodel.md)
  : Predictions or log posterior predictive densities from a reference
  model

- [`predictor_terms()`](https://mc-stan.org/projpred/dev/reference/predictor_terms.md)
  :

  Predictor terms used in a
  [`project()`](https://mc-stan.org/projpred/dev/reference/project.md)
  run

- [`print(`*`<projection>`*`)`](https://mc-stan.org/projpred/dev/reference/print.projection.md)
  :

  Print information about
  [`project()`](https://mc-stan.org/projpred/dev/reference/project.md)
  output

- [`print(`*`<refmodel>`*`)`](https://mc-stan.org/projpred/dev/reference/print.refmodel.md)
  : Print information about a reference model object

- [`print(`*`<vsel>`*`)`](https://mc-stan.org/projpred/dev/reference/print.vsel.md)
  :

  Print results (summary) of a
  [`varsel()`](https://mc-stan.org/projpred/dev/reference/varsel.md) or
  [`cv_varsel()`](https://mc-stan.org/projpred/dev/reference/cv_varsel.md)
  run

- [`print(`*`<vselsummary>`*`)`](https://mc-stan.org/projpred/dev/reference/print.vselsummary.md)
  :

  Print summary of a
  [`varsel()`](https://mc-stan.org/projpred/dev/reference/varsel.md) or
  [`cv_varsel()`](https://mc-stan.org/projpred/dev/reference/cv_varsel.md)
  run

- [`project()`](https://mc-stan.org/projpred/dev/reference/project.md) :
  Projection onto submodel(s)

- [`projpred`](https://mc-stan.org/projpred/dev/reference/projpred-package.md)
  [`projpred-package`](https://mc-stan.org/projpred/dev/reference/projpred-package.md)
  : Projection predictive feature selection

- [`ranking()`](https://mc-stan.org/projpred/dev/reference/ranking.md) :
  Predictor ranking(s)

- [`get_refmodel()`](https://mc-stan.org/projpred/dev/reference/refmodel-init-get.md)
  [`init_refmodel()`](https://mc-stan.org/projpred/dev/reference/refmodel-init-get.md)
  : Reference model and more general information

- [`run_cvfun()`](https://mc-stan.org/projpred/dev/reference/run_cvfun.md)
  :

  Create `cvfits` from `cvfun`

- [`solution_terms()`](https://mc-stan.org/projpred/dev/reference/solution_terms.md)
  :

  Retrieve the full-data solution path from a
  [`varsel()`](https://mc-stan.org/projpred/dev/reference/varsel.md) or
  [`cv_varsel()`](https://mc-stan.org/projpred/dev/reference/cv_varsel.md)
  run or the predictor combination from a
  [`project()`](https://mc-stan.org/projpred/dev/reference/project.md)
  run

- [`suggest_size()`](https://mc-stan.org/projpred/dev/reference/suggest_size.md)
  : Suggest submodel size

- [`summary(`*`<vsel>`*`)`](https://mc-stan.org/projpred/dev/reference/summary.vsel.md)
  :

  Summary of a
  [`varsel()`](https://mc-stan.org/projpred/dev/reference/varsel.md) or
  [`cv_varsel()`](https://mc-stan.org/projpred/dev/reference/cv_varsel.md)
  run

- [`varsel()`](https://mc-stan.org/projpred/dev/reference/varsel.md) :
  Run search and performance evaluation without cross-validation

- [`y_wobs_offs()`](https://mc-stan.org/projpred/dev/reference/y_wobs_offs.md)
  : Extract response values, observation weights, and offsets
