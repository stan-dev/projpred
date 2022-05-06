# projpred 2.1.1.9000

## Minor changes

* Account for changes concerning the handling of offsets in **rstanarm** version 2.21.3. In particular, issue stan-dev/rstanarm#542 was fixed in **rstanarm** 2.21.3.
* Show the output of the vignette on CRAN.
* In the vignette, use `cv_varsel()` with LOO CV and `validate_search = FALSE` instead of K-fold CV. (GitHub: #305)
* Improve the documentation for argument `search_terms` of `varsel()` and `cv_varsel()`. (GitHub: #155, #308)
* In case of user-specified (non-`NULL`) `search_terms`, `method = NULL` is internally changed to `method = "forward"` and `method = "L1"` throws a warning. This is done because `search_terms` only takes effect in case of a forward search. (GitHub: #155, #308)
* Internally, the intercept is now always included in `search_terms`. This is necessary to prevent a bug described below. (GitHub: #308)
* When fitting multilevel submodels via **lme4**, **projpred** now tries to handle `PIRLS loop resulted in NaN value` errors automatically.

## Bug fixes

* Throw a more informative error message in case of special group-level terms which are currently not supported (in particular, nested ones).
* Previously, using a `search_terms` vector which excluded the intercept in conjunction with `refit_prj = FALSE` (the latter in `project()`, `varsel()`, or `cv_varsel()`) led to incorrect submodels being fetched from the search or to an error while doing so. This has been fixed now by internally forcing the inclusion of the intercept in `search_terms`. (GitHub: #308)
* Fix GitHub issues #147 and #202. (GitHub: #312)

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
* Support for the `Gamma()` family is regarded as experimental. Therefore, a corresponding warning is thrown when creating the reference model. (GitHub: paul-buerkner/brms#1255, #240, #252)
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
* Handle **rstanarm**'s GitHub issue #551. This implies that **projpred**'s default `cvfun` for `stanreg` fits will now always use *inner* parallelization in `rstanarm::kfold()` (i.e., across chains, not across CV folds), with `getOption("mc.cores", 1)` cores. We do so on all systems (not only Windows). (GitHub: #249)
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
