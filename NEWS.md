

# News

## projpred 2.0.5.9000

### Major changes

* The behavior of arguments `ndraws`, `nclusters`, `ndraws_pred`, and `nclusters_pred` in `varsel()`, `cv_varsel()`, and `project()` has been changed: Now, `ndraws` and `ndraws_pred` have non-`NULL` defaults and for `ndraws <= 20` or `ndraws_pred <= 20`, the value of `ndraws` or `ndraws_pred` is passed to `nclusters` or `nclusters_pred`, respectively (so that clustering is used). (GitHub: commits babe031db7732e0d81dd2591938551d02dcf374d, 4ef95d3b4ab85eaaa5177c4d40f33b2943bff37c, and ce7d1e001fd76830c4379cbbe0dfe730cba8d9e5)
* For `proj_linpred()` and `proj_predict()`, arguments `nterms`, `ndraws`, and `seed` have been removed to allow the user to pass them to `project()`. New arguments `filter_nterms`, `nresample_clusters`, and `ppd_seed` have been introduced (see the documentation for details). (GitHub: #92, #135)
* Reference models lacking an intercept are not supported anymore (actually, the previous implementation for such models was incomplete). Support might be re-introduced in the future (when fixed), but for now it is withdrawn as it requires some larger changes. (GitHub: #124, but see also #96 and #100)

## Minor changes

* Minor documentation improvements.
* Minor improvements of error messages.
* For K-fold cross-validation, an internally hard-coded value of 5 for `nclusters_pred` was removed. (GitHub: commit 5062f2ff6f981ab0e4be06b9aaf694dcaa27afa8)
* Throw a proper error message for nonsupported families. (GitHub: #140)
* Show the README also on the CRAN website. (GitHub: #140)
* `project()`: Warn in case of `solution_terms` not being found in the reference model (and therefore getting ignored). (GitHub: #140)
* Set defaults for `get_refmodel.default()`'s arguments `ref_predfun`, `proj_predfun`, and `div_minimizer`. (GitHub: #140)
* Explicitly throw an error if argument `extract_model_data` of `get_refmodel.default()` is `NULL`. (GitHub: #140)

### Bug fixes

* Fixed a bug when using weights or offsets e.g. in `proj_linpred()`. (GitHub: #114)
* Fixed a bug causing `varsel()`/`make_formula` to fail with multidimensional interaction terms. (GitHub: #102, #103)
* Fixed an indexing bug in `cv_varsel()` for models with a single predictor. (GitHub: #115)
* Fixed bugs for argument `nterms` of `proj_linpred()` and `proj_predict()`. (GitHub: #110)
* Fixed an inconsistency for some intercept-only submodels. (GitHub: #119)
* Minor documentation fixes.
* Fix a bug for `as.matrix.projection()` in case of 1 (clustered) draw after projection. (GitHub: #130)
* For submodels of class `"subfit"`, make the column names of `as.matrix.projection()`'s output matrix consistent with other classes of submodels. (GitHub: #132)
* Fix a bug for argument `nterms_max` of `plot.vsel()` if there is just the intercept-only submodel. (GitHub: #138)
* Throw an appropriate error message when trying to apply an L1 search to an empty (i.e. intercept-only) reference model. (GitHub: #139)
* Fix the list names of element `search_path` in, e.g., `varsel()`'s output. (GitHub: #140)
* Fix a bug (error `unused argument`) when initializing the K reference models in a K-fold CV with CV fits not of class `"brmsfit"` or `"stanreg"`. (GitHub: #140)
* In `get_refmodel.default()`, remove old defunct arguments `fetch_data`, `wobs`, and `offset`. (GitHub: #140)

## projpred 2.0.5

### Minor changes

* For GLMMs, the column names of the matrix returned by the `as.matrix.projection()` method follow [**brms**](https://paul-buerkner.github.io/brms/)'s naming convention, also for the new columns introduced by **projpred** version 2.0.4 (see below). (GitHub: #82)
* Internally, the seed is not fixed to a specific value when `NULL`. (GitHub: #84)

### Bug fixes

* Fixed a bug raising an error when not projecting from a `vsel` object. (GitHub: #79, #80)
* Fixed a bug in the calculation of the Gaussian deviance. (GitHub: #81)
* Fixed a bug in the calculation of the predictive statistics of the reference model on test data in `varsel()`. (GitHub #90)
* In `init_refmodel()`: Raise an error if neither `cvfun` nor `cvfits` is provided (in cases where at least one of them is necessary). (GitHub: #91)
* Fixed a bug in an input check for argument `nloo` of `cv_varsel()`. (GitHub: #93)
* Fixed a bug in `cv_varsel()`, causing an error in case of `!validate_search && cv_method != "LOO"`. (GitHub: #95)
* Fixed bugs related to the setting of the seed. (GitHub: commit 02cd50db76b0f2d835ce8f8f39cbe94353540d64)
* Fixed a bug causing `proj_linpred()` to raise an error if argument `newdata` was `NULL`. (GitHub: #97)
* Fixed an incorrect usage of the dispersion parameter values when calculating output element `lpd` in `proj_linpred()` (for `integrated = TRUE` as well as for `integrated = FALSE`). (GitHub: #105)
* Fixed bugs in `proj_linpred()`'s calculation of output element `lpd` (for `integrated = TRUE`). (GitHub: #106, #112)
* Fixed an inconsistency in the dimensions of `proj_linpred()`'s output elements `pred` and `lpd` (for `integrated = FALSE`): Now, they are both S x N matrices, with S denoting the number of (possibly clustered) posterior draws and N denoting the number of observations. (GitHub: #107, #112)
* Fixed a bug causing `proj_predict()`'s output matrix to be transposed in case of `nrow(newdata) == 1`. (GitHub: #112)

## projpred 2.0.4

* Added support for weighted LOO proportional-to-size subsampling based on Magnusson, M., Riis Andersen, M., Jonasson, J. and Vehtari, A. (2019). Leave-One-Out Cross-Validation for Large Data. In International Conference on Machine Learning.
* Automatically explore both linear and smooths components in GAM models. This allows the user to gauge the impact of the smooth term against its linear counterpart. 
* Fast approximate LOO computation for `validate_search = FALSE` calls in `cv_varsel(...)`.
* Improved summary output with important details.
* The (internally set) default for argument `nclusters` of `varsel()` and `cv_varsel()` was increased from 10 to 20.
* For group-level effects, the `as.matrix.projection()` method now also returns the estimated group-level effects themselves. (GitHub: #75)
* For group-level effects, the `as.matrix.projection()` method now returns the variance components (population SD(s) and population correlation(s)) instead of the empirical SD(s) of the group-level effects. (GitHub: #74)

### Bug fixes

* Fixed a bug in the handling of arguments `ndraws` and `nclusters` in `varsel()` and `cv_varsel()`. (GitHub: commit bbd0f0ad8041cc938e57b2064773fac7d0d2ce9b)
* Fixed a bug in `as.matrix.projection()` (causing incorrect column names for the returned matrix). (GitHub: #72, #73)

## projpred 2.0.3

* Minor fixes for stability, no new features.
* The (internally set) default for argument `nclusters` of `varsel()` and `cv_varsel()` was increased from 1 to 10.
* The (internally set) default for argument `nclusters_pred` of `varsel()` and `cv_varsel()` was increased from 5 to 400.
* In `varsel()` and `cv_varsel()`, always perform clustering of the posterior draws (argument `ndraws` is effectively ignored). (GitHub: starting with commit 80b10dca483ad90edd7d3a1aa8b881351b9673ac)

## projpred 2.0.2

We have fully rewritten the internals in several ways. Most importantly, we now leverage maximum likelihood estimation to third parties depending on the reference model's family. This allows a lot of flexibility and extensibility for various models. Functionality wise, the major updates since the last release are:

* Added support for GLMMs and GAMMs via ```lme4``` and ```gamm4```.
* Formula syntax support internally that allows for easier building upon projections.
* Thanks to the above point, we save some computation by only considering sensible projections during forward search instead of fitting every possible submodel.
* We have added a new argument ```search_terms``` that allows the user to specify custom unit building blocks of the projections. New vignette coming up.
* We have fully changed the way to define custom reference models. The user now provides projection fitting and prediction functions (more information in a new upcoming vignette).

## projpred 1.1.4

Better validation of function arguments.

## projpred 1.1.3

Added print methods for vsel and cvsel objects. Added AUC statistics for binomial family. A few additional minor patches.

## projpred 1.1.2

Removed the dependency on the ```rngtools``` package.

## projpred 1.1.1

This version contains only a few patches, no new features to the user.

## projpred 1.1.0

### New features 

* Added support for ```brms``` models. 

### Bug fixes
* The program crashed with ```rstanarm``` models fitted with syntax like ```stan_glm(log(y) ~ log(x), ...)```, that is, it did not allow transformation for ```y```.


## projpred 1.0.0

### New features and improvements ###

* Changed the internals so that now all fit objects (such as rstanarm fits) are converted to ```refmodel```-objects using the generict ```get_refmodel```-function, and all the functions use only this object. This makes it much easier to use projpred with other reference models by writing them a new ```get_refmodel```-function. The syntax is now changed so that  ```varsel``` and ```cv_varsel``` both return an object that has similar structure always, and the reference model is stored into this object.

* Added more examples to the vignette.

* Added possibility to change the baseline in ```plot/summary```. Now it is possible to compare also to the best submodel found, not only to the reference model.

* Bug fix: RMSE was previously computed wrong, this is now fixed.

* Small changes: ```nloo = n``` by default in ```cv_varsel```. ```regul=1e-4``` now by default in all functions.


## projpred 0.9.0

### New features and improvements

* Added the ```cv_search``` argument for the main functions (```varsel```,```cv_varsel```,```project``` and the prediction functions). Now it is possible to make predictions also with those parameter estimates that were computed during the L1-penalized search. This change also allows the user to compute the Lasso-solution by providing the observed data as the 'reference fit' for init_refmodel. An example will be added to the vignette.

### Bug fixes

* The projection with a nonzero regularization parameter value did not produce exactly correct result, although the difference to the correct result was often so small that user would not see the difference. Fixed this.


## projpred 0.8.0 and earlier

Until this version, we did not keep record of the changes between different versions. Started to do this from version 0.9.0 onwards.
