

# News

## projpred 2.0.4.9000

### Major changes

### Minor changes

* For GLMMs, the column names of the matrix returned by the `as.matrix.projection()` method follow [**brms**](https://paul-buerkner.github.io/brms/)'s naming convention, also for the new columns introduced by **projpred** version 2.0.4 (see below). (GitHub: #82)
* Internally, the seed is not fixed to a specific value when `NULL`. (GitHub: #84)

### Bug fixes

* Fixed a bug raising an error when not projecting from a `vsel` object. (GitHub: #79, #80)
* Fixed a bug in the calculation of the Gaussian deviance. (GitHub: #81)
* Fixed a bug in the calculation of the predictive statistics of the reference model on test data in `varsel()`. (GitHub #90)
* In `init_refmodel()`: Raise an error if neither `cvfun` nor `cvfits` is provided (in cases where at least one of them necessary). (GitHub: #91)
* Fixed a bug in an input check for argument `nloo` of `cv_varsel()`. (GitHub: #93)
* Fixed a bug in `cv_varsel()`, causing an error in case of `!validate_search && cv_method != "LOO"`. (GitHub: #95)
* Fixed bugs related to the setting of the seed. (GitHub: commit 02cd50db76b0f2d835ce8f8f39cbe94353540d64)

## projpred 2.0.4

* Added support for weighted LOO proportional-to-size subsampling based on Magnusson, M., Riis Andersen, M., Jonasson, J. and Vehtari, A. (2019). Leave-One-Out  Cross-Validation for Large Data. In International Conference on Machine Learning. 
* Automatically explore both linear and smooths components in GAM models. This allows the user to gauge the impact of the smooth term against its linear counterpart. 
* Fast approximate LOO computation for `validate_search = FALSE` calls in `cv_varsel(...)`.
* Improved summary output with important details.
* For group-level effects, the `as.matrix.projection()` method now also returns the estimated group-level effects themselves. (GitHub: #75)
* For group-level effects, the `as.matrix.projection()` method now returns the variance components (population SD(s) and population correlation(s)) instead of the empirical SD(s) of the group-level effects. (GitHub: #74)

### Bug fixes

* Fixed a bug in `as.matrix.projection()` (causing incorrect column names for the returned matrix). (GitHub: #72, #73)

## projpred 2.0.3

Minor fixes for stability, no new features.

## projpred 2.0.2

We have fully rewritten the internals in several ways. Most importantly, we now leverage maximum likelihood estimation to third parties depending on the reference model's family. This allows a lot of flexibility and extensibility for various models. Functionality wise, the major updates since the last release are:

* Added support for GLMMs and GAMMs via ```lme4``` and ```gamm4```.
* Formula syntax support internally that allows for easier building upon projections.
* Thanks to the above point, we save some computation by only considering sensible projections during forward search instead of fitting every possible submodel.
* We have added a new argument ```search_terms``` that allows the user to specify custom unit building blocks of the projections. This can be used to include _fixed terms_ across all projections, for instance. New vignette coming up.
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
