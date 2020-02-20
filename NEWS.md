

# News

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

* Added possibility to change the baseline in ```varsel_plot/varsel_stats```. Now it is possible to compare also to the best submodel found, not only to the reference model.

* Bug fix: RMSE was previously computed wrong, this is now fixed.

* Small changes: ```nloo = n``` by default in ```cv_varsel```. ```regul=1e-4``` now by default in all functions.


## projpred 0.9.0

### New features and improvements

* Added the ```cv_search``` argument for the main functions (```varsel```,```cv_varsel```,```project``` and the prediction functions). Now it is possible to make predictions also with those parameter estimates that were computed during the L1-penalized search. This change also allows the user to compute the Lasso-solution by providing the observed data as the 'reference fit' for init_refmodel. An example will be added to the vignette.

### Bug fixes

* The projection with a nonzero regularization parameter value did not produce exactly correct result, although the difference to the correct result was often so small that user would not see the difference. Fixed this.


## projpred 0.8.0 and earlier

Until this version, we did not keep record of the changes between different versions. Started to do this from version 0.9.0 onwards.
