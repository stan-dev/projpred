# Note

If you read this from a place other than <https://mc-stan.org/projpred/news/index.html>, please consider switching to that website since it features better formatting and cross-linking.

# projpred 2.6.0.9000

## Minor changes

* If an L1 search selects an interaction term before all involved lower-order interaction terms (including main-effect terms) have been selected, the predictor ranking is now automatically modified so that the lower-order interaction terms come before this interaction term. A corresponding warning is thrown, which may be deactivated by setting the global option `projpred.warn_L1_interactions` to `FALSE`. Previously, beginning with version 2.5.0, only a warning was thrown and this only if an L1 search selected an interaction term before all involved *main-effect* terms had been selected. (GitHub: #420)
* Added a progress bar for `project()` (when using the built-in divergence minimizers). For this, `project()` has gained a new argument `verbose` which can also be controlled via the global option `projpred.verbose_project`. By default, the new progress bar is activated. (GitHub: #421)
* Added a new argument `parallel` to `cv_varsel()`. With `parallel = TRUE`, costly parts of **projpred**'s cross-validation (CV) can be run in parallel. See the documentation of that new argument (and section "Note" of `cv_varsel()`'s documentation) for details. (GitHub: #422)
* Added a warning for issue #323 (for multilevel Gaussian models, the projection onto the full model can be instable). (GitHub: #426)
* `plot.vsel()`: The uncertainty bars are slightly thicker now and the points slightly larger. (GitHub: #429)
* `plot.vsel()`: Added argument `ranking_colored` for coloring the points and the uncertainty bars according to the magnitude of the (possibly cumulated) CV ranking proportions. (GitHub: #430; thanks to @yannmclatchie for the suggestion)
* Added warnings for most of the problems described in section ["Troubleshooting"](https://mc-stan.org/projpred/articles/projpred.html#troubleshooting) of the main vignette. (GitHub: #431)
* Output element `p_type` of `project()` has been removed. Instead, output element `const_wdraws_prj` has been added, but its definition is essentially the inverse of former element `p_type` (see the updated documentation of `project()`'s output). This should not be a breaking change for users (as `p_type` was mainly intended for internal use and the new element `const_wdraws_prj` is so, too) but this slightly enhances the cases where `as.matrix.projection()` throws a warning concerning the weights of the projected draws and the cases where `proj_predict()` resamples from the projected draws using argument `nresample_clusters`. (GitHub: #432)

## Bug fixes

* Fixed a bug in the printed number of projected draws for the performance evaluation when calling `print.vselsummary()` based on output from `varsel()` with `refit_prj = FALSE`.
* Fixed a bug sometimes causing an error when predicting from a submodel that is a GLM and has interactions. (GitHub: #420)
* Fixed a bug introduced in version 2.6.0, causing an incompatibility of K-fold CV with R versions < 4.2.0. (GitHub: #423, #427)

# projpred 2.6.0

## Major changes

* In anticipation of a larger overhaul of the **projpred** user interface, this release comes with several new functions for accessing and investigating solution paths (which are now termed *predictor rankings* by these new functions, a term that is hopefully easier to grasp for new users):

    + Added a new function called `ranking()` which returns the predictor ranking from the full-data search and possibly also the predictor rankings from fold-wise searches in case of cross-validation (CV). (More precisely, `ranking()` is a generic. The only method is `ranking.vsel()`, applicable to objects returned by `varsel()` or `cv_varsel()`. The output is of class `ranking`.)
    + Added a new function called `cv_proportions()` which computes ranking proportions (across CV folds, see `?cv_proportions` for details) from fold-wise predictor rankings. (More precisely, `cv_proportions()` is a generic. The main method is `cv_proportions.ranking()`, but as a shortcut, `cv_proportions.vsel()` has also been added. The output is of class `cv_proportions`.)
    + Added a new `plot()` method called `plot.cv_proportions()` for plotting ranking proportions from fold-wise predictor rankings. (As a shortcut, `plot.ranking()` has also been added.)
    
    Because of these new functions, a message has been added to `print.vselsummary()`, mentioning how to access and investigate the fold-wise predictor rankings (if they exist). Furthermore, due to these changes, element `pct_solution_terms_cv` of `vsel` objects has been replaced with element `solution_terms_cv` which contains the fold-wise predictor rankings instead of the corresponding ranking proportions. However, elements of `vsel` objects are not meant to be accessed directly, so this replacement should not be a breaking change for most users. Finally, method `solution_terms.vsel()` (which---until now---was the only possibility to extract the full-data predictor ranking) has now been deprecated and will be removed in a future release. Please use the new function `ranking()` instead (more precisely, `ranking()`'s output element `fulldata` contains the full-data predictor ranking that is also extracted by `solution_terms.vsel()`; `ranking()`'s output element `foldwise` contains the fold-wise predictor rankings---if available---which were previously not accessible via a built-in function). (GitHub: #289, #406, #411)
* Added function `predictor_terms()` which retrieves the predictor terms used in a `project()` run. Correspondingly, method `solution_terms.projection()` has now been deprecated and will be removed in a future release. Please use `predictor_terms()` instead. (GitHub: #411)
* `seed` (and `.seed`) arguments now have a default of `NA` instead of `sample.int(.Machine$integer.max, 1)` and the pseudorandom number generator (PRNG) state is reset only if the user-supplied seed is not `NA`. This allows setting a seed once at the beginning of any **projpred**-related code and then leaving all `seed` (and `.seed`) arguments at their default. Previously, such practice could lead to results which were "less random" than they should have been because the former default of `sample.int(.Machine$integer.max, 1)` caused **projpred** functions with a `seed` (or `.seed`) argument to reset the PRNG state upon exit, meaning that two repeated calls to `cv_varsel()` (for example) with no PRNG-using code between them would use the same seed internally. (GitHub: #412)
* Added the main diagonal of the matrix returned by `cv_proportions()` to a new column called `cv_proportions_diag` of the summary table computed by `summary.vsel()`. The purpose of this new column is to give a basic sense for the (CV) variability in the ranking of the predictors. Argument `cumulate` of `cv_proportions()` has been added to `summary.vsel()` as well (to allow the ranking proportions in the newly added column to be *cumulated* ranking proportions, if desired). (GitHub: #289, #413)
* Added the full-data predictor ranking and the main diagonal of the matrix returned by `cv_proportions()` to the plot created by `plot.vsel()`. These new elements can be omitted by setting `plot.vsel()`'s new argument `ranking_nterms_max` to `NA` (setting it to some specific submodel size causes the full-data predictor ranking and the corresponding ranking proportions to be omitted after that size). Argument `cumulate` of `cv_proportions()` has been added to `plot.vsel()` as well (to allow the ranking proportions to be *cumulated* ranking proportions, if desired). Other new arguments are `ranking_abbreviate` (together with `ranking_abbreviate_args`), `ranking_repel` (together with `ranking_repel_args`), and `text_angle` (see the `plot.vsel()` documentation for details). (GitHub: #289, #414, #416, #417)

## Minor changes

* Enhancements in the vignettes. In particular, the new functions `ranking()`, `cv_proportions()`, and `plot.cv_proportions()` (see "Major changes" above) are now illustrated in the main vignette. (GitHub: #407, #411)
* Reduced the peak memory usage of `cv_varsel()` with `cv_method = "kfold"`. This may slightly change results from such a `cv_varsel()` run compared to older **projpred** versions due to different pseudorandom number generator (PRNG) states when clustering posterior draws. (GitHub: #419)

## Bug fixes

* Fixed a bug causing L1 search to throw an error in case of some `I()` terms. (GitHub: #404, #408)
* Fixed a bug causing L1 search to throw an error in case of `poly()` or `polym()` terms. Note that just like `step()` and `MASS::stepAIC()`, **projpred**'s search algorithms do not split up a `poly()` or `polym()` term into its lower-degree polynomial terms (which would be helpful, for example, if the linear part of a `poly()` term with `degrees = 2` was relevant but the quadratic part not). Such a split-up of a `poly()` or `polym()` term needs to be performed manually (if desired). (GitHub: #183, #409)
* Fixed a bug causing some non-smooth predictor terms to be treated as smooth terms. (GitHub: #182, #410)
* See "Major changes" above: Fixed a bug causing **projpred** functions with a `seed` (or `.seed`) argument to use the same seed internally when users set a seed once at the beginning (via `set.seed()`) and then had two or more calls to such **projpred** functions with their `seed` (or `.seed`) argument being at its default and no PRNG-using code between those calls. (GitHub: #412)

# projpred 2.5.0

## Minor changes

* Setting the new global option `projpred.extra_verbose` to `TRUE` will print out which submodel **projpred** is currently projecting onto. Furthermore, if `method = "forward"` and `verbose = TRUE` in `varsel()` or `cv_varsel()`, this new option will also make **projpred** print out which submodel has been selected at those steps of the forward search for which a percentage is printed (the percentage refers to the maximum submodel size that the search is run up to). In general, however, we cannot recommend setting this new global option to `TRUE` for `cv_varsel()` with `validate_search = TRUE` (simply due to the amount of information that will be printed, but also due to the progress bar which will not work anymore as intended). (GitHub: #363; thanks to @jtimonen)
* Enhanced `verbose` output. In particular, `varsel()` is now more verbose, similarly to how `cv_varsel()` has already been for a long time. The  `verbose` output for `cv_varsel()` has also been updated, with the aim to give users a better understanding of the methodology behind **projpred**. (GitHub: #382)
* Slightly improved the calculation of predictive variances to make them less prone to numerical inaccuracies. (GitHub: #199)
* Improved computational efficiency by avoiding an unnecessary final full-data performance evaluation (including costly re-projections if `refit_prj = TRUE`, which is the default for non-`datafit` reference models) in `cv_varsel()` with `validate_search = TRUE`. Due to this change, results from `cv_varsel()` (with `validate_search = TRUE`) may slightly change due to a different pseudorandom number generator (PRNG) state when clustering posterior draws. The different PRNG state was necessary to make the PRNG state for the full-data search in the `validate_search = TRUE` case consistent to the PRNG state for the full-data search in the `validate_search = FALSE` case. (GitHub: #385)
* Reduced dependencies. (GitHub: #388)
* Argument `digits` of `print.vselsummary()` which used to be passed to an internal `round()` call was removed. Instead, `digits` can now be passed to `print.data.frame()` via `...`, thereby determining the minimum number of *significant digits* to be printed. (GitHub: #389)
* Although bad practice (in general), a reference model lacking an intercept can now be used within **projpred**. However, it will always be projected onto submodels which *include* an intercept. The reason is that even if the true intercept in the reference model is zero, this does not need to hold for the submodels. An informational message mentioning the projection onto intercept-including submodels is thrown when **projpred** encounters a reference model lacking an intercept. (GitHub: #96, #391)
* In case of non-predictor arguments of `s()` or `t2()`, **projpred** now throws an error. (This had already been documented before, but a suitable error message was missing.) (GitHub: #393, based on #156 and #269)
* In case of the `brms::categorical()` family (supported since version 2.4.0), **projpred** now strips underscores from response category names in `as.matrix.projection()` output, as done by **brms**. (GitHub: #394)
* L1 search now throws a warning if an interaction term is selected before all involved main-effect terms have been selected. (GitHub: #395)
* Documented that in multilevel (group-level) terms, function calls on the right-hand side of the `|` character (e.g., `(1 | gr(group_variable))`, which is possible in **brms**) are currently not allowed in **projpred**. A corresponding error message has also been added. (GitHub: #319)
* Due to internal refactoring:
    
    + `project()`'s output elements `submodl` and `weights` have been renamed to `outdmin` and `wdraws_prj`, respectively.
    + `varsel()`'s and `cv_varsel()`'s output element `d_test` has been replaced with new output elements `type_test` and `y_wobs_test`.
    
    Apart from `project()`'s output element `wdraws_prj`, these elements are not meant to be accessed manually, so changes are mentioned here only for the sake of completeness. Output element `wdraws_prj` of `project()` is only needed if `project()` was used for a clustered projection, which is not the default (and discouraged in most applied cases, at least with a small number of clusters). Thus, these renamings are breaking changes only in very rare cases.
* `print.vselsummary()` now also prints `K` in case of K-fold CV.
* The `print.vselsummary()` output has been slightly improved, e.g., adding a remark what "search included" or "search not included" means.
* `print.vselsummary()` now also prints whether `deltas = TRUE` or `deltas = FALSE` was used.
* Output element `pct_solution_terms_cv` has now also been added to `vsel` objects returned by `varsel()`, but in that case, it is simply `NULL`. This (`pct_solution_terms_cv` being `NULL`) is now also the case if `validate_search = FALSE` was used in `cv_varsel()`.
* Minor enhancements in the documentation.
* Enhancements in the vignettes. In particular, section ["Troubleshooting"](https://mc-stan.org/projpred/articles/projpred.html#troubleshooting) of the main vignette has been revised.
* If `proj_predict()` is used with observation weights that are not all equal to `1`, a warning is now thrown. (GitHub: starts to address #402)

## Bug fixes

* Fixed a long-standing bug (existing at least from version 2.0.2 on) causing `predict.refmodel()` to require `newdata` to contain the response variable in case of a **brms** reference model. This is similar to [paul-buerkner/brms#1457](https://github.com/paul-buerkner/brms/issues/1457), but concerns `predict.refmodel()` ([paul-buerkner/brms#1457](https://github.com/paul-buerkner/brms/issues/1457) referred to predictions from the *submodels*). In order to make this `predict.refmodel()` fix work, **brms** version 2.19.0 or later is needed. (GitHub: #381)
* Fixed a long-standing bug (existing from version 2.1.0 on) causing output element `p_type` of `project()` to be incorrect in case of `refit_prj = FALSE`, `!is.null(nclusters)`, and an `object` of class `vsel` that was created with a non-clustered (thinned) projection during the search phase. The fix comes with a slightly different behavior of `proj_predict()` for `datafit`s: It will not draw `nresample_clusters` times from the posterior-projection predictive distribution (which is based on the same single projected draw), but only once. (GitHub: #211, #401)
* When performing predictions from submodels which are GLMs (or from submodels which are L1-penalized GLMs, which is only possible in case of `refit_prj = FALSE` after an L1 search), a new dataset containing a `character` predictor variable with only a single unique value (or a new dataset containing a `factor` predictor variable with a single level) used to cause an error. The case of a `character` (not `factor`) predictor variable with only a single unique value occurred, e.g., during the performance evaluation in a LOO CV if a `character` predictor got selected into a fold's solution path. The `character` issue existed from version 2.1.0 on (in earlier versions, however, there were other issues which caused `character` predictors to throw an error). Now, all issues with respect to `character` predictor variables should be resolved. The issue with single-level `factor` predictor variables is resolved now as well. (GitHub: #403)
* When performing predictions from submodels which are GLMs (or from submodels which are L1-penalized GLMs, which is only possible in case of `refit_prj = FALSE` after an L1 search), a new dataset containing a `factor` predictor with re-ordered levels (compared to this same `factor` in the original dataset) used to lead to incorrect predictions. This bug existed at least from version 2.0.2 on (possibly even in earlier versions), but has been resolved now. (GitHub: #403)
* Fixed an error thrown by **projpred**'s internal GLM submodel fitter in case of unused levels of a `factor`. This issue existed at least from version 2.0.2 on (possibly even in earlier versions), but should have only affected **rstanarm** reference model fits (**brms** reference model fits were only affected in case of a `brms::brm()` call with `drop_unused_levels = FALSE`, which is not the default). (GitHub: #403)
* Fixed a bug that caused an L1 search combined with `refit_prj = FALSE` (which is the default only for `datafit`s, not for the reference model objects of class `refmodel` that are usually employed in practice) to lead to incorrect predictions from the L1-searched submodels (which are L1-penalized GLMs) if the solution path had a main effect ranked after an interaction term. This bug existed at least from version 2.0.2 on (possibly even in earlier versions). The mentioned submodel predictions did not only affect the performance evaluation, but also the projected dispersion parameter and the returned Kullback-Leibler divergence (and the corresponding cross-entropy). (GitHub: #403)

# projpred 2.4.0

## Major changes

* Introduction of the augmented-data projection [(Weber and Vehtari, 2023)](https://doi.org/10.48550/arXiv.2301.01660) (see section ["Supported types of models"](https://mc-stan.org/projpred/articles/projpred.html#modtypes) of the main vignette for details). (GitHub: #70, #322)
* Introduction of the latent projection [(Catalina et al., 2021)](https://doi.org/10.48550/arXiv.2109.04702) (see section ["Supported types of models"](https://mc-stan.org/projpred/articles/projpred.html#modtypes) of the main vignette and the new [latent-projection vignette](https://mc-stan.org/projpred/articles/latent.html) for details). A consequence of the latent projection (more precisely, of the `resp_oscale = TRUE` default in `summary.vsel()`) is that `varsel()` and `cv_varsel()` no longer call `suggest_size()` internally at the end. Thus, `print()`-ing an object of class `vsel` no longer includes the suggested projection size in the output (the `stat` for this suggested size was fixed to `"elpd"` anyway, a fact that many users were probably not aware of). (GitHub: #372)
* In case of multilevel models, **projpred** now has two global options for "integrating out" group-level effects: `projpred.mlvl_pred_new` and `projpred.mlvl_proj_ref_new`. These are explained in detail in the general package documentation (available [online](https://mc-stan.org/projpred/reference/projpred-package.html) or by typing `` ?`projpred-package` ``). (GitHub: #379)

## Minor changes

* Improvements in the numerical stability of internal link and inverse-link functions. (GitHub: #376)

## Bug fixes

* Fix a bug for offsets in cases where `family` (see `init_refmodel()`) has a non-identity link function: After clustering the reference model's posterior draws, we need to aggregate (within a given cluster) the reference model's fitted values which already take the offsets into account instead of taking the offsets into account after aggregating the fitted values which do *not* take the offsets into account. This fix should affect results only in a very slight manner. Due to **projpred**'s internal adjustment for numerical stability when averaging a quantity across the draws within a given cluster, this also changes the projected residual standard deviations in Gaussian models in the order of `1e-10`. (GitHub: #374)

# projpred 2.3.0

## Major changes

* In `plot.vsel()` and `summary.vsel()`, the default of `alpha = 0.32` is replaced by `alpha = 2 * pnorm(-1)` (= `1 - diff(pnorm(c(-1, 1)))`, which is only *approximately* 0.32) so that now, a normal-approximation confidence interval with default `alpha` stretches by exactly one standard error on either side of the point estimate. Typically, this changes results only slightly. In some cases, however, the new default may lead to a different suggested size, explaining why this is regarded as a major change. (GitHub: #371)

## Minor changes

* The deprecated function `ggplot2::aes_string()` is not used anymore, thereby avoiding an occasional soft-deprecation warning thrown by **ggplot2** 3.4.0. (GitHub: #367)
* The KL divergence from the reference model to a submodel is simplified to the corresponding cross-entropy (i.e., the reference model's entropy is dropped), with some caveats described in the documentation for output element `ce` of `project()`. The reason for this change is that the former KL divergence assumed the reference model's family to be the same as the submodel's family, which does not need to be the case for custom reference models. This should not be a user-facing change as users are discouraged to make use of specific output elements (like the former element `kl` of objects of class `projection` or `vsel`) directly. (GitHub: #369)
* Improvements in the documentation (especially for argument `family` of `init_refmodel()` and `get_refmodel.default()`).

# projpred 2.2.2

## Major changes

* Several important bug fixes (see below).

## Minor changes

* Improvements in documentation and vignette, especially to emphasize the generality of the reference model object resulting from `get_refmodel()` and `init_refmodel()` (thereby also distinguishing more clearly between "typical" and "custom" reference model objects) in (i) the description and several arguments of `get_refmodel()` and `init_refmodel()`, (ii) sections ["Reference model"](https://mc-stan.org/projpred/articles/projpred.html#refmod) and ["Supported types of models"](https://mc-stan.org/projpred/articles/projpred.html#modtypes) of the vignette. (GitHub: #357, #359, #364, #365, #366)
* Minor improvement in terms of efficiency in the `validate_search = FALSE` case of `cv_varsel()`.
* Improvement in terms of efficiency in case of a forward search with custom `search_terms` (at least in some instances), also affecting the output of `solution_terms(<vsel_object>)` in those cases. (GitHub: #360; thanks to @sor16)
* Update [Catalina et al. (2020)](https://doi.org/10.48550/arXiv.2010.06994) to [Catalina et al. (2022)](https://proceedings.mlr.press/v151/catalina22a.html). (GitHub: #364)

## Bug fixes

* Fix a bug causing offsets not to be taken into account appropriately when calculating the PSIS weights (those used for the submodels) in the `validate_search = FALSE` case of `cv_varsel()`. This bug was introduced in v2.2.0 (and existed up to---including---v2.2.1).
* Fix a (long-standing) bug causing offsets not to be taken into account appropriately when calculating the predictive variances for a reference model that has a dispersion parameter and a non-identity link function. (GitHub: #186 (partly), #355)
* Fix a (long-standing) bug causing offsets not to be taken into account appropriately when calculating the reference model's summary statistics in case of `cv_varsel()` with `cv_method = "LOO"` (more precisely, only the LOO posterior predictive expected values `<vsel_object>$summaries$ref$mu` were affected, not the (pointwise) LOO log posterior predictive density values `<vsel_object>$summaries$ref$lppd`). (GitHub: #186 (partly), #356)
* Fix a (long-standing) bug leading to an error when trying to use `cv_varsel()` with custom `search_terms` (in some instances). (GitHub: #345, #360; thanks to @sor16)

# projpred 2.2.1

## Minor changes

* Several improvements in the documentation.
* For the RMSE as well as the AUC (see argument `stats` of `summary.vsel()`), the bootstrapping results are now also used for inferring the lower and upper confidence interval bounds. (GitHub: #318, #347; thanks to @awd97 and @VisionResearchBlog)
* For `datafit`s, offsets are not supported anymore. (GitHub: #186 (partly), #351)

## Bug fixes

* Fix GitHub issue #348 (L1 search in the presence of interaction terms). This bug was introduced in v2.1.0 (and existed up to---including---v2.2.0).
* Fix incorrectly thrown messages in case of `datafit`s (and other---unlikely---cases where `nclusters == S` and `S <= 20`, with `S` denoting the number of draws in the reference model).
* Fix GitHub issue #349 (only concerned `datafit`s). (GitHub: #350)

# projpred 2.2.0

## Major changes

* In the `validate_search = FALSE` case of `cv_varsel()` (with `cv_method = "LOO"`), the PSIS weights are now calculated based on the reference model (they used to be calculated based on the submodels which is incorrect). (GitHub: #325)
* Some long-standing severe bugs (GitHub issues #329, #330, and #342) have been fixed, concerning the performance evaluation of models with nontrivial observation weights (i.e., models where at least one observation had a weight differing from 1). Concerned performance statistics were `"mse"`, `"rmse"`, `"acc"` (= `"pctcorr"`), and `"auc"` (i.e., all performance statistics except for `"elpd"` and `"mlpd"`).
* `plot.vsel()` and `suggest_size()` gain a new argument `thres_elpd`. By default, this argument doesn't have any impact, but a non-`NA` value can be used for a customized model size selection rule (see `?suggest_size` for details). (GitHub: #335)

## Minor changes

* Several improvements in the documentation (especially in the explanation of the `suggest_size()` heuristic).
* Improvement of the numerical stability for some link functions, achieved by avoiding unnecessary back-and-forth transformations between latent space and response space. (GitHub: #337, #338)
* All arguments `seed` and `.seed` are now allowed to be `NA` for not calling `set.seed()` internally at all.
* Argument `d_test` of `varsel()` is not considered as an internal feature anymore. This was possible after fixing a bug for `d_test` (see below). (GitHub: #341)
* The order of the observations in the sub-elements of `<vsel_object>$summaries` and `<vsel_object>$d_test` now corresponds to the order of the observations in the original dataset if `<vsel_object>` was created by a call to `cv_varsel([...], cv_method = "kfold")` (formerly, in that case, the observations in those sub-elements were ordered by fold). Thereby, the order of the observations in those sub-elements now always corresponds to the order of the observations in the original dataset, except if `<vsel_object>` was created by a call to `varsel([...], d_test = <non-NULL_d_test_object>)`, in which case the order of the observations in those sub-elements corresponds to the order of the observations in `<non-NULL_d_test_object>`. (GitHub: #341)

## Bug fixes

* Fix GitHub issue #324 (large `search_terms` caused the R session to crash).
* Fix GitHub issue #204. (GitHub: #325)
* Fix the `validate_search = FALSE` bug described above in "Major changes": The PSIS weights are now calculated based on the reference model (they used to be calculated based on the submodels which is incorrect). (GitHub: #325)
* Fix `\mbox{}` commands displayed incorrectly in the HTML help from R version 4.2.0 on. (GitHub: #326)
* Fix GitHub issue #329 (see also "Major changes" above).
* Fix GitHub issue #331.
* `plot.vsel()` now draws the dashed red horizontal line for the reference model (and---if present---the dotted black horizontal line for the baseline model) first (i.e., before the submodel-specific graphical elements), to avoid overplotting.
* Fix GitHub issue #339. (GitHub: #340)
* Fix argument `d_test` of `varsel()`: Not only the predictive performance of the *reference model* needs to be evaluated on the test data, but also the predictive performance of the *submodels*. (GitHub: #341)
* Fix GitHub issue #342 (see also "Major changes" above).
* Fix GitHub issue #330 (see also "Major changes" above). (GitHub: #344, commit 23e7101)

# projpred 2.1.2

## Minor changes

* Account for changes concerning the handling of offsets in **rstanarm** version 2.21.3. In particular, issue [stan-dev/rstanarm#542](https://github.com/stan-dev/rstanarm/issues/542) was fixed in **rstanarm** 2.21.3.
* Show the output of the vignette on CRAN.
* In the vignette, use `cv_varsel()` with LOO CV and `validate_search = FALSE` instead of K-fold CV. (GitHub: #305)
* Improve the documentation for argument `search_terms` of `varsel()` and `cv_varsel()`. (GitHub: #155, #308)
* In case of user-specified (non-`NULL`) `search_terms`, `method = NULL` is internally changed to `method = "forward"` and `method = "L1"` throws a warning. This is done because `search_terms` only takes effect in case of a forward search. (GitHub: #155, #308)
* Internally, the intercept is now always included in `search_terms`. This is necessary to prevent a bug described below. (GitHub: #308)
* When fitting multilevel submodels via **lme4**, **projpred** now tries to handle `PIRLS loop resulted in NaN value` errors automatically. (GitHub: #314)
* The fix for GitHub issue #320 (see below) required to rename argument `b` of `projpred:::bootstrap()` to `B`.

## Bug fixes

* Throw a more informative error message in case of special group-level terms which are currently not supported (in particular, nested ones).
* Previously, using a `search_terms` vector which excluded the intercept in conjunction with `refit_prj = FALSE` (the latter in `project()`, `varsel()`, or `cv_varsel()`) led to incorrect submodels being fetched from the search or to an error while doing so. This has been fixed now by internally forcing the inclusion of the intercept in `search_terms`. (GitHub: #308)
* Fix GitHub issues #147 and #202. (GitHub: #312)
* Fix GitHub issue #320. (GitHub: #321)

# projpred 2.1.1

## Bug fixes

* Fix the order of the package authors.
* Fix failing CRAN checks.
* Add an input check for argument `solution_terms` of `project()` to fix a test failure in R versions >= 4.2.

# projpred 2.1.0

## Major changes

* Added support for weighted LOO proportional-to-size subsampling based on [Magnusson et al. (2019)](https://proceedings.mlr.press/v97/magnusson19a.html). However, subsampled LOO CV is now regarded as experimental. Therefore, a corresponding warning is thrown when calling `cv_varsel()` with `nloo < n` where `n` denotes the number of observations. (GitHub: #94, #252, commit feea39e)
* Automatically explore both linear and smooths components in GAM models. This allows the user to gauge the impact of the smooth term against its linear counterpart.
* Fast approximate LOO computation for `validate_search = FALSE` in `cv_varsel()`.
* Formerly, the defaults for arguments `nclusters` (= `1`) and `nclusters_pred` (= `5`) of `varsel()` and `cv_varsel()` were set internally (the user-visible defaults were `NULL`). Now, `nclusters` and `ndraws_pred` (note the `ndraws_pred`, not `nclusters_pred`) have non-`NULL` user-visible defaults of `20` and `400`, respectively. In general, this increases the runtime of these functions a lot. With respect to `cv_varsel()`, the new vignette (see [vignettes](https://mc-stan.org/projpred/articles/)) mentions two ways to quickly obtain some rough preliminary results which in general should not be used as final results, though: (i) `varsel()` and (ii) `cv_varsel()` with `validate_search = FALSE` (which only takes effect for `cv_method = "LOO"`). (GitHub: #291 and several commits beforehand, in particular bbd0f0a, babe031, 4ef95d3, and ce7d1e0)
* For `proj_linpred()` and `proj_predict()`, arguments `nterms`, `ndraws`, and `seed` have been removed to allow the user to pass them to `project()`. New arguments `filter_nterms`, `nresample_clusters`, and `.seed` have been introduced (see the documentation for details). (GitHub: #92, #135)
* Reference models lacking an intercept are not supported anymore (actually, the previous implementation for such models was incomplete). Support might be re-introduced in the future (when fixed), but for now it is withdrawn as it requires some larger changes. (GitHub: #124, but see also #96 and #100)
* In the output of `proj_linpred()`, dimensions are not dropped anymore (i.e., output elements `pred` and `lpd` are always S x N matrices now). (GitHub: #143)
* In case of `integrated = TRUE`, `proj_linpred()` now averages the LPD (across the projected posterior draws) instead of taking the LPD at the averaged linear predictors. (GitHub: #143)
* If `newdata` does not contain the response variable, `proj_linpred()` now returns `NULL` for output element `lpd`. (GitHub: #143)
* The fix for the offset issues (listed below under "Bug fixes") requires reference model fits of class `stanreg` (from package **rstanarm**) with offsets to have these offsets specified via an `offset()` term in the model formula (and not via argument `offset`).
* Improved handling of errors when fitting multilevel submodels. (GitHub: #201)
* Some defaults have been changed from `NULL` to a user-visible value (and `NULL` is not allowed anymore).
* Argument `data` of `get_refmodel.stanreg()` has been removed. (GitHub: #219)
* The function passed to argument `div_minimizer` of `init_refmodel()` now always needs to return a `list` of submodels (see the documentation for details). Correspondingly, the function passed to argument `proj_predfun` of `init_refmodel()` can now always expect a `list` as input for argument `fits` (see the documentation for details). (GitHub: #230)
* The function passed to argument `proj_predfun` of `init_refmodel()` now always needs to return a matrix (see the documentation for details). (GitHub: #230)
* The projection can be run in parallel now. However, we cannot recommend this for all kinds of platforms and all kinds of models. For more information, see the general package documentation available at ``?`projpred-package` ``. (GitHub: #235)
* Support for the `Student_t()` family is regarded as experimental. Therefore, a corresponding warning is thrown when creating the reference model. (GitHub: #233, #252)
* Support for additive models (i.e., GAMs and GAMMs) is regarded as experimental. Therefore, a corresponding warning is thrown when creating the reference model. (GitHub: #237, #252)
* Support for the `Gamma()` family is regarded as experimental. Therefore, a corresponding warning is thrown when creating the reference model. (GitHub: [paul-buerkner/brms#1255](https://github.com/paul-buerkner/brms/issues/1255), #240, #252)
* The previous behavior of `init_refmodel()` in case of argument `dis` being `NULL` (the default) was dangerous for custom reference models with a `family` having a dispersion parameter (in that case, `dis` values of all-zeros were used silently). The new behavior now requires a non-`NULL` argument `dis` in that case. (GitHub: #254)
* Argument `cv_search` has been renamed to `refit_prj`. (GitHub: #154, #265)
* `as.matrix.projection()` has gained a new argument `nm_scheme` which allows to choose the naming scheme for the column names of the returned matrix. The default (`"auto"`) follows the naming scheme of the reference model fit (and uses the `"rstanarm"` naming scheme if the reference model fit is of an unknown class). (GitHub: #82, #279)
* `seed` (and `.seed`) arguments now have a default of `sample.int(.Machine$integer.max, 1)` instead of `NULL`. Furthermore, the value supplied to these arguments is now used to generate new seeds internally on-the-fly. In many cases, this will change results compared to older **projpred** versions. Also note that now, the internal seeds are never fixed to a specific value if `seed` (and `.seed`) arguments are set to `NULL`. (GitHub: #84, #286)

## Minor changes

* Improved summary output with important details.
* For group-level effects, the `as.matrix.projection()` method now also returns the estimated group-level effects themselves. (GitHub: #75)
* For group-level effects, the `as.matrix.projection()` method now returns the variance components (population SD(s) and population correlation(s)) instead of the empirical SD(s) of the group-level effects. (GitHub: #74)
* Improved documentation. (GitHub: especially #233)
* Replaced the two vignettes by a single one which also has new content. (GitHub: #237)
* Updated the `README` file. (GitHub: #245)
* Some error and warning messages have been improved and added. (GitHub: especially #219, #221, #223, #252, #263)
* For K-fold cross-validation, an internally hard-coded value of 5 for `nclusters_pred` was removed. (GitHub: commit 5062f2f)
* Throw a proper error message for unsupported families. (GitHub: #140)
* Show the README also on the CRAN website. (GitHub: #140)
* `project()`: Warn if elements of `solution_terms` are not found in the reference model (and therefore ignored). (GitHub: #140)
* `get_refmodel.default()` now passes arguments via the ellipsis (`...`) to `init_refmodel()`. (GitHub: #153, commit dd3716e)
* Remove dependency on package **rngtools** (version 2.0.0 of **projpred** re-introduced this dependency after it was already removed in version 1.1.2). (GitHub: #189)
* `init_refmodel()`: The default (`NULL`) for argument `extract_model_data` has been removed as it wasn't meaningful anyway. (GitHub: #219)
* Argument `folds` of `init_refmodel()` has been removed as it was effectively unused. (GitHub: #220)
* Use the S3 system for `solution_terms()`. This allowed the introduction of a `solution_terms.projection()` method. (GitHub: #223)
* `predict.refmodel()` now uses a default of `newdata = NULL`. (GitHub: #223)
* Argument `weights` of `init_refmodel()`'s argument `proj_predfun` has been removed. (GitHub: #163, #224)
* **projpred**'s internal `div_minimizer` functions have been unified into a single `div_minimizer` which chooses an appropriate submodel fitter based on the formula of the submodel, not based on that of the reference model. Furthermore, the automatic handling of errors in the submodel fitters has been improved. (GitHub: #230)
* Improve the axis labels in `plot.vsel()`. (GitHub: #234, #270)
* Handle **rstanarm**'s GitHub issue #551. This implies that **projpred**'s default `cvfun` for `stanreg` fits will now always use *inner* parallelization in `rstanarm::kfold.stanreg()` (i.e., across chains, not across CV folds), with `getOption("mc.cores", 1)` cores. We do so on all systems (not only Windows). (GitHub: #249)
* Argument `fit` of `init_refmodel()`'s argument `proj_predfun` was renamed to `fits`. This is a non-breaking change since all calls to `proj_predfun` in **projpred** have that argument unnamed. However, this cannot be guaranteed in the future, so we strongly encourage users with a custom `proj_predfun` to rename argument `fit` to `fits`. (GitHub: #263)
* `init_refmodel()` has gained argument `cvrefbuilder` which may be a custom function for constructing the K reference models in a K-fold CV. (GitHub: #271)
* Allow arguments to be passed from `project()`, `varsel()`, and `cv_varsel()` to the divergence minimizer. (GitHub: #278)
* In `init_refmodel()`, any `contrasts` attributes of the dataset's columns are silently removed. (GitHub: #284)
* `NA`s in data supplied to `newdata` arguments now trigger an error. (GitHub: #285)

## Bug fixes

* Fixed a bug in `as.matrix.projection()` (causing incorrect column names for the returned matrix). (GitHub: #72, #73)
* Fixed a bug raising an error when not projecting from a `vsel` object. (GitHub: #79, #80)
* Fixed a bug in the calculation of the Gaussian deviance. (GitHub: #81)
* Fixed a bug in the calculation of the predictive statistics of the reference model on test data in `varsel()`. (GitHub #90)
* Fixed a bug in an input check for argument `nloo` of `cv_varsel()`. (GitHub: #93)
* Fixed a bug in `cv_varsel()`, causing an error in case of `!validate_search && cv_method != "LOO"`. (GitHub: #95)
* Fixed bugs related to the setting of the seed. (GitHub: commit 02cd50d)
* Fixed a bug causing `proj_linpred()` to raise an error if argument `newdata` was `NULL`. (GitHub: #97)
* Fixed an incorrect usage of the dispersion parameter values when calculating output element `lpd` in `proj_linpred()` (for `integrated = TRUE` as well as for `integrated = FALSE`). (GitHub: #105)
* Fixed bugs in `proj_linpred()`'s calculation of output element `lpd` (for `integrated = TRUE`). (GitHub: #106, #112)
* Fixed an inconsistency in the dimensions of `proj_linpred()`'s output elements `pred` and `lpd` (for `integrated = FALSE`): Now, they are both S x N matrices, with S denoting the number of (possibly clustered) posterior draws and N denoting the number of observations. (GitHub: #107, #112)
* Fixed a bug causing `proj_predict()`'s output matrix to be transposed in case of `nrow(newdata) == 1`. (GitHub: #112)
* Fixed a bug when using weights or offsets e.g. in `proj_linpred()`. (GitHub: #114)
* Fixed a bug causing `varsel()`/`make_formula` to fail with multidimensional interaction terms. (GitHub: #102, #103)
* Fixed an indexing bug in `cv_varsel()` for models with a single predictor. (GitHub: #115)
* Fixed bugs for argument `nterms` of `proj_linpred()` and `proj_predict()`. (GitHub: #110)
* Fixed an inconsistency for some intercept-only submodels. (GitHub: #119)
* Fix a bug for `as.matrix.projection()` in case of 1 (clustered) draw after projection. (GitHub: #130)
* For submodels of class `subfit`, make the column names of `as.matrix.projection()`'s output matrix consistent with other classes of submodels. (GitHub: #132)
* Fix a bug for argument `nterms_max` of `plot.vsel()` if there is just the intercept-only submodel. (GitHub: #138)
* Throw an appropriate error message when trying to apply an L1 search to an empty (i.e. intercept-only) reference model. (GitHub: #139)
* Fix the list names of element `search_path` in, e.g., `varsel()`'s output. (GitHub: #140)
* Fix a bug (error `unused argument`) when initializing the K reference models in a K-fold CV with CV fits not of class `brmsfit` or `stanreg`. (GitHub: #140)
* In `get_refmodel.default()`, remove old defunct arguments `fetch_data`, `wobs`, and `offset`. (GitHub: #140)
* Fix a bug in `get_refmodel.stanreg()`. (GitHub: #142, #184)
* Fix a possible bug related to `extract_model_data()`'s argument `extract_y` in `get_refmodel.default()`. (GitHub: #153, commit 39fece8)
* Fix a possible bug related to `extract_model_data()` in K-fold CV. (GitHub: #153, commit 4f32195)
* Fix GitHub issue #161.
* Fix GitHub issue #162.
* Fix GitHub issue #164.
* Fix GitHub issue #160.
* Fix GitHub issue #159.
* Fix GitHub issue #158.
* Fix GitHub issue #157.
* Fix GitHub issue #144.
* Fix GitHub issue #146.
* Fix GitHub issue #169.
* Fix GitHub issue #167.
* Fix a bug in the default `proj_predfun()` for GLMMs. (GitHub: #174)
* Fix GitHub issue #171.
* Fix GitHub issue #172.
* Fix a bug in the default `proj_predfun()` for `datafit`s. (GitHub: #177)
* Fix the names of `summary.vsel()$selection` for objects of class `vsel` created by `varsel()`. (GitHub: #179)
* Fix forward search when `search_terms` are not consecutive in size. (GitHub: commit 34e24de)
* Fix a bug in `cv_varsel()$pct_solution_terms_cv`. (GitHub: #188, commit e529ec1)
* Fix GitHub issue #185. (GitHub: #193, #194)
* Fix a bug in forward searches with interaction terms. (GitHub: #191)
* Fix offset issues. (GitHub: #196, #203, #228)
* Fix a bug in `glm_elnet()` (the workhorse for L1 search), causing the grid for lambda to be constructed without taking observation weights into account. (GitHub: #198; note that the second part of #198 did not have any consequences for users)
* Fix GitHub issue #136. (GitHub: #221)
* Fix a bug in `print.vsel()` causing argument `digits` to be ignored. (GitHub: #222)
* Fix a bug causing the default of argument `cv_search` in `varsel()` and `cv_varsel()` to be `TRUE` for `datafit`s, although it should be `FALSE` in that case. (GitHub: #223)
* Fix a bug (`Error: Levels '<...>' of grouping factor '<...>' cannot be found in the fitted model. Consider setting argument 'allow_new_levels' to TRUE.`) when predicting from submodels which are GLMMs for `newdata` containing new levels for grouping factors. (GitHub: #223)
* `predict.refmodel()`: Fix a bug for integer `ynew`. (GitHub: #223)
* `predict.refmodel()`: Fix input checks for `offsetnew` and `weightsnew`. (GitHub: #223)
* After all calls to `extract_model_data()`, the weights and offsets are now checked if they are of length 0 (and if yes, then they are set to vectors of ones and zeros, respectively). This is important for `extract_model_data()` functions which return weights and offsets of length 0 (see, e.g., `brms` version <= 2.16.1). (GitHub: #223)
* Handle **rstanarm**'s GitHub issue #546. (GitHub: #227)
* Fix a bug causing the internal submodel fitter for GLMMs to not pass arguments `var` (the predictive variances) and `regul` (amount of ridge regularization) to the internal submodel fitter for GLMs. (GitHub: #230)
* Fix GitHub issue #210. (GitHub: #234)
* Fix GitHub issue #242. (GitHub: #253)
* Fix GitHub issue #244. (GitHub: #255)
* Fix GitHub issue #243. (GitHub: #262)
* Fix GitHub issue #213. (GitHub: #264)
* Fix GitHub issue #215. (GitHub: #266)
* Fix GitHub issue #212. (GitHub: #267)
* Fix GitHub issue #156. (GitHub: #269)
* If the data used for the reference model contains `NA`s, an appropriate error is now thrown. Previously, the reference model was created successfully, but this caused opaque errors in downstream code such as `project()`. (GitHub: #274)
* Fix GitHub issue #268. (GitHub: #287)
* Fix GitHub issue #149. (GitHub: #288)

# projpred 2.0.2

We have fully rewritten the internals in several ways. Most importantly, we now leverage maximum likelihood estimation to third parties depending on the reference model's family. This allows a lot of flexibility and extensibility for various models. Functionality wise, the major updates since the last release are:

* Added support for GLMMs and GAMMs via **lme4** and **gamm4**.
* Formula syntax support internally that allows for easier building upon projections.
* Thanks to the above point, we save some computation by only considering sensible projections during forward search instead of fitting every possible submodel.
* We have added a new argument `search_terms` that allows the user to specify custom unit building blocks of the projections. New vignette coming up.
* We have fully changed the way to define custom reference models. The user now provides projection fitting and prediction functions (more information in a new upcoming vignette).

# projpred 1.1.4

Better validation of function arguments.

# projpred 1.1.3

Added print methods for vsel and cvsel objects. Added AUC statistics for binomial family. A few additional minor patches.

# projpred 1.1.2

Removed the dependency on the **rngtools** package.

# projpred 1.1.1

This version contains only a few patches, no new features to the user.

# projpred 1.1.0

## New features 

* Added support for [**brms**](https://paul-buerkner.github.io/brms/) models. 

## Bug fixes

* The program crashed with [**rstanarm**](https://mc-stan.org/rstanarm/) models fitted with syntax like `stan_glm(log(y) ~ log(x), ...)`, that is, it did not allow transformation for `y`.

# projpred 1.0.0

## New features and improvements

* Changed the internals so that now all fit objects (such as rstanarm fits) are converted to `refmodel`-objects using the generic `get_refmodel`-function, and all the functions use only this object. This makes it much easier to use projpred with other reference models by writing them a new `get_refmodel`-function. The syntax is now changed so that  `varsel` and `cv_varsel` both return an object that has similar structure always, and the reference model is stored into this object.
* Added more examples to the vignette.
* Added possibility to change the baseline in `plot/summary`. Now it is possible to compare also to the best submodel found, not only to the reference model.
* Bug fix: RMSE was previously computed wrong, this is now fixed.
* Small changes: `nloo = n` by default in `cv_varsel`. `regul=1e-4` now by default in all functions.

# projpred 0.9.0

## New features and improvements

* Added the `cv_search` argument for the main functions (`varsel`,`cv_varsel`,`project` and the prediction functions). Now it is possible to make predictions also with those parameter estimates that were computed during the L1-penalized search. This change also allows the user to compute the Lasso-solution by providing the observed data as the 'reference fit' for init_refmodel. An example will be added to the vignette.

## Bug fixes

* The projection with a nonzero regularization parameter value did not produce exactly correct result, although the difference to the correct result was often so small that user would not see the difference. Fixed this.

# projpred 0.8.0 and earlier

Until this version, we did not keep record of the changes between different versions. Started to do this from version 0.9.0 onwards.
