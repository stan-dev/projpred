#' Projection predictive feature selection
#'
#' @useDynLib projpred
#' @importFrom Rcpp sourceCpp
#'
#' @import stats
#' @import ggplot2
#' @importFrom rstantools posterior_linpred
#' @importFrom loo kfold
#'
#' @description
#'
#' The \R package \pkg{projpred} performs the projection predictive variable (or
#' "feature") selection for various regression models. We recommend to read the
#' `README` file (available with enhanced formatting
#' [online](https://mc-stan.org/projpred/)) and the main vignette (`topic =
#' "projpred"`, but also available
#' [online](https://mc-stan.org/projpred/articles/projpred.html)) before
#' continuing here.
#'
#' @details
#'
#' # Terminology
#'
#' Throughout the whole package documentation, we use the term "submodel" for
#' all kinds of candidate models onto which the reference model is projected.
#' For custom reference models, the candidate models don't need to be actual
#' *sub*models of the reference model, but in any case (even for custom
#' reference models), the candidate models are always actual *sub*models of the
#' full [`formula`] used by the search procedure. In this regard, it is correct
#' to speak of *sub*models, even in case of a custom reference model.
#'
#' The following model type abbreviations will be used at multiple places
#' throughout the documentation: GLM (generalized linear model), GLMM
#' (generalized linear multilevel---or "mixed"---model), GAM (generalized
#' additive model), and GAMM (generalized additive multilevel---or
#' "mixed"---model). Note that the term "generalized" includes the Gaussian
#' family as well.
#'
#' # Draw-wise divergence minimizers
#'
#' For the projection of the reference model onto a submodel, \pkg{projpred}
#' currently relies on the following functions as draw-wise divergence
#' minimizers (in other words, these are the workhorse functions employed by
#' \pkg{projpred}'s internal default `div_minimizer` functions, see
#' [init_refmodel()]):
#' * Submodel without multilevel or additive terms:
#'     + For the traditional (or latent) projection (or the augmented-data
#'     projection in case of the [binomial()] or [brms::bernoulli()] family): An
#'     internal C++ function which basically serves the same purpose as [lm()]
#'     for the [gaussian()] family and [glm()] for all other families. The
#'     returned object inherits from class `subfit`. Possible tuning parameters
#'     for this internal C++ function are: `regul` (amount of ridge
#'     regularization; default: `1e-4`), `thresh_conv` (convergence threshold;
#'     default: `1e-7`, but not used in case of the Gaussian family with
#'     identity link), `qa_updates_max` (maximum number of quadratic
#'     approximation updates; default: `100`, but fixed to `1` in case of the
#'     Gaussian family with identity link), `ls_iter_max` (maximum number of
#'     line search iterations; default: `30`, but fixed to `1` in case of the
#'     Gaussian family with identity link), `normalize` (single logical value
#'     indicating whether to scale the predictors internally with the returned
#'     regression coefficient estimates being back-adjusted appropriately;
#'     default: `TRUE`), `beta0_init` (single numeric value giving the starting
#'     value for the intercept at centered predictors; default: `0`), and
#'     `beta_init` (numeric vector giving the starting values for the regression
#'     coefficients; default: vector of `0`s).
#'     + For the augmented-data projection: [MASS::polr()] (the returned object
#'     inherits from class `polr`) for the [brms::cumulative()] family or
#'     [rstanarm::stan_polr()] fits, [nnet::multinom()] (the returned object
#'     inherits from class `multinom`) for the [brms::categorical()] family.
#' * Submodel with multilevel but no additive terms:
#'     + For the traditional (or latent) projection (or the augmented-data
#'     projection in case of the [binomial()] or [brms::bernoulli()] family):
#'     [lme4::lmer()] (the returned object inherits from class `lmerMod`) for
#'     the [gaussian()] family, [lme4::glmer()] (the returned object inherits
#'     from class `glmerMod`) for all other families. For [lme4::glmer()],
#'     \pkg{projpred} uses `control = lme4::glmerControl(optimizer = "bobyqa")`,
#'     so both optimization stages use BOBYQA. A `control` argument in `...`
#'     replaces this default in the first fit. After some \pkg{lme4} errors,
#'     \pkg{projpred} refits with its own settings and ignores this `control`.
#'     For the [binomial()] family, \pkg{lme4} versions before 2.1-1 round
#'     non-integer responses and weights in the objective of [lme4::glmer()]
#'     with `nAGQ >= 1`. With these \pkg{lme4} versions, \pkg{projpred} uses
#'     `nAGQ = 0` for such submodels unless `nAGQ` is set via `...`. Then the
#'     fixed effects and the group-level effects do not depend on the rounded
#'     objective; the group-level standard deviations still do.
#'     + For the augmented-data projection: [ordinal::clmm()] (the returned
#'     object inherits from class `clmm`) for the [brms::cumulative()] family,
#'     [mclogit::mblogit()] (the returned object inherits from class `mmblogit`)
#'     for the [brms::categorical()] family. For [ordinal::clmm()] with a single
#'     group-level term `(1 | <group>)`, \pkg{projpred} uses
#'     `method = "ucminf"` and `grad = "forward"` (see
#'     [ordinal::clmm.control()]). If this fit fails or does not converge,
#'     \pkg{projpred} refits with the default settings of [ordinal::clmm()]. An
#'     argument `method` or `control` in `...` replaces this default. For the
#'     `cloglog` and `loglog` links, \pkg{projpred} always uses
#'     `useMatrix = TRUE` (see [ordinal::clmm.control()]) and does not use
#'     `ucminf()`. The code of [ordinal::clmm()] for a single group-level term
#'     uses the other of these two links.
#' * Submodel without multilevel but additive terms: [mgcv::gam()] (the returned
#' object inherits from class `gam`).
#' * Submodel with multilevel and additive terms: [gamm4::gamm4()] (within
#' \pkg{projpred}, the returned object inherits from class `gamm4`).
#'
#' # Verbosity, messages, warnings, errors
#'
#' Global option `projpred.verbose` may be used for specifying the value passed
#' to argument `verbose` of [project()], [varsel()], and [cv_varsel()].
#'
#' By default, \pkg{projpred} catches messages and warnings from the draw-wise
#' divergence minimizers and throws their unique collection after performing all
#' draw-wise divergence minimizations (i.e., draw-wise projections). This can be
#' deactivated by setting global option `projpred.warn_proj_drawwise` to
#' `FALSE`.
#'
#' Furthermore, by default, \pkg{projpred} checks the convergence of the
#' draw-wise divergence minimizers and throws a warning if any seem to have not
#' converged. This warning is thrown after the warning message from global
#' option `projpred.warn_proj_drawwise` (see above) and can be deactivated by
#' setting global option `projpred.check_convergence` to `FALSE`.
#'
#' # Parallelization
#'
#' \pkg{projpred} can run two kinds of computations in parallel: the fold-wise
#' computations of [cv_varsel()] (the *CV* parallelization, activated via
#' argument `parallel` of [cv_varsel()], which in turn can be controlled via
#' global option `projpred.parallel_cv`) and the projections of the single
#' draws (or clusters) of the reference model onto a submodel (the
#' *projection* parallelization, see below).
#'
#' Both use the parallel workers of the current plan of the \pkg{future}
#' package, e.g., `future::plan(future::multisession, workers = 4)`. Plans
#' from other packages for \pkg{future} can be used as well, e.g.,
#' `future::plan(future.mirai::mirai_multisession, workers = 4)` from the
#' \pkg{future.mirai} package. In our experiments with 8 workers, the first
#' parallel call with the workers of \pkg{future.mirai} was faster than with
#' those of `future::multisession`, and the workers of \pkg{future.mirai}
#' returned big results faster (see the projection parallelization below).
#' With the sequential plan of \pkg{future} (the default), a parallel run is
#' not possible, so \pkg{projpred} then runs sequentially and
#' throws a warning. The workers must be able to load the same
#' \pkg{projpred} installation as the main R session; \pkg{projpred} checks
#' this. (If the main R session uses \pkg{projpred} via
#' `pkgload::load_all()` or `devtools::load_all()`, the workers need to load
#' the same source directory in the same way.) The parallelization is not
#' tested on Windows.
#'
#' The projection parallelization is controlled by global option
#' `projpred.parallel_proj_trigger`: the number of projected draws below which
#' the projections of the draws are not run in parallel (even if parallel
#' workers are set up). Such a "trigger" threshold exists because of the
#' computational overhead of a parallelization. By default, the projection
#' parallelization is turned off, which can also be achieved by supplying
#' `Inf` (or `NULL`) to option `projpred.parallel_proj_trigger`. The
#' projection parallelization is mainly useful for submodels with multilevel
#' or additive terms: for most other submodels, \pkg{projpred} fits all
#' projected draws in a single call, which does not run in parallel. It is
#' useful for [varsel()], [project()], and [cv_varsel()] without fold-wise
#' computations; with the CV parallelization, the workers do not run the
#' projections of their folds in parallel. Note that the fitted submodel
#' objects of multilevel and additive submodels are quite big,
#' which---when running in parallel---may lead to excessive memory usage
#' which in turn may crash the R session (on Unix systems, setting an
#' appropriate memory limit via [unix::rlimit_as()] may avoid crashing the
#' whole machine). The workers also need time to return these objects to the
#' main R session: in our experiments with 16 workers of
#' `future::multisession`, this took longer than a sequential run, whereas
#' 16 workers of \pkg{future.mirai} gave a clear speedup. So for the
#' projection parallelization, we recommend the workers of
#' \pkg{future.mirai}.
#'
#' In the fold-wise computations of [cv_varsel()], each CV fold has its own
#' random number stream (L'Ecuyer-CMRG), which \pkg{projpred} derives from the
#' random number generator of the main R session, with and without the CV
#' parallelization. So for a given seed, the results depend neither on
#' argument `parallel`, nor on the number of workers, nor on the plan of the
#' \pkg{future} package.
#'
#' A function goes to a parallel worker with its enclosing environments up
#' to, but not including, the global environment. So if a function of a
#' custom reference model is defined in a script and uses other objects of
#' the global environment, the workers need these objects.
#' If the \pkg{globals} package is installed, \pkg{projpred} finds them from
#' the code of such functions and sends them to the workers. This cannot find
#' objects whose names the code creates at run time (e.g., with [get()]).
#' Global option `projpred.export_to_workers` may be set to a character vector
#' of names of such objects, to export them from the global environment to
#' the parallel workers as well. Setting that option to `FALSE` turns off any
#' export of objects of the global environment.
#'
#' During parallelization (either of the projection or the CV), progression
#' updates can be received via the \pkg{progressr} package, e.g., by calling
#' `progressr::handlers(global = TRUE)` before running the projection or the
#' CV in parallel. The \pkg{projpred} package also offers the global option
#' `projpred.use_progressr` for controlling whether to use the \pkg{progressr}
#' package (`TRUE` or `FALSE`), but since that global option defaults to
#' `requireNamespace("progressr", quietly = TRUE) && interactive()`, it usually
#' does not need to be set by the user. Messages and warnings from the workers
#' are shown in the main R session.
#'
#' # Multilevel models: "Integrating out" group-level effects
#'
#' In case of multilevel models, \pkg{projpred} offers two global options for
#' "integrating out" group-level effects: `projpred.mlvl_pred_new` and
#' `projpred.mlvl_proj_ref_new`. When setting `projpred.mlvl_pred_new` to `TRUE`
#' (default is `FALSE`), then at
#' *prediction* time, \pkg{projpred} will treat group levels existing in the
#' training data as *new* group levels, implying that their group-level effects
#' are drawn randomly from a (multivariate) Gaussian distribution. This concerns
#' both, the reference model and the (i.e., any) submodel. Furthermore, setting
#' `projpred.mlvl_pred_new` to `TRUE` causes [as.matrix.projection()] and
#' [as_draws_matrix.projection()] to omit the projected group-level effects (for
#' the group levels from the original dataset). When setting
#' `projpred.mlvl_proj_ref_new` to `TRUE` (default is `FALSE`), then at
#' *projection* time, the reference model's fitted values (that the submodels
#' fit to) will be computed by treating the group levels from the original
#' dataset as *new* group levels, implying that their group-level effects will
#' be drawn randomly from a (multivariate) Gaussian distribution (as long as the
#' reference model is a multilevel model, which---for custom reference
#' models---does not need to be the case). This also affects the latent response
#' values for a latent projection correspondingly. Setting
#' `projpred.mlvl_pred_new` to `TRUE` makes sense, e.g., when the prediction
#' task is such that any group level will be treated as a new one. Typically,
#' setting `projpred.mlvl_proj_ref_new` to `TRUE` only makes sense when
#' `projpred.mlvl_pred_new` is already set to `TRUE`. In that case, the default
#' of `FALSE` for `projpred.mlvl_proj_ref_new` ensures that at projection time,
#' the submodels fit to the best possible fitted values from the reference
#' model, and setting `projpred.mlvl_proj_ref_new` to `TRUE` would make sense if
#' the group-level effects should be integrated out completely.
#'
#' # Memory usage
#'
#' By setting the global option `projpred.run_gc` to `TRUE`, \pkg{projpred} will
#' call [gc()] at some places (e.g., after each size that the forward search
#' passes through) to free up some memory. These [gc()] calls are not always
#' necessary to reduce the peak memory usage, but they add runtime (hence the
#' default of `FALSE` for that global option).
#'
#' During a call of [varsel()], [cv_varsel()] or [project()], \pkg{projpred}
#' caches the submodel formulas and model matrices per set of predictor terms
#' (for submodels without multilevel or additive terms). Global option
#' `projpred.proj_cache_mb` (default `512`) limits the size of this cache, in
#' units of 2^20 bytes as estimated by [utils::object.size()]. The cache holds
#' at most about twice this limit. The limit applies to each \R process, so
#' forked parallel workers can each have a cache of this size. Setting global
#' option `projpred.proj_cache` to `FALSE` switches the cache off. The cache
#' changes only the runtime, not the results.
#'
#' # Other notes
#'
#' Global option `projpred.digits` controls arguments `digits` of
#' [print.vselsummary()] and [print.vsel()].
#'
#' There are several global options to control arguments of [plot.vsel()] and
#' [plot.cv_proportions()] globally, see section "Usage" of the help pages of
#' these two functions.
#'
#' Global option `projpred.warn_L1_interactions` may be set to `FALSE` to
#' deactivate a warning that an L1 search selected an interaction term before
#' all involved lower-order interaction terms (including main-effect terms) were
#' selected (in which case the predictor ranking is automatically modified by
#' \pkg{projpred} so that the lower-order interaction terms come before this
#' interaction term).
#'
#' Most examples are not executed when called via [example()]. To execute them,
#' their code has to be copied and pasted manually to the console.
#'
#' # Functions
#'
#' \describe{
#'   \item{[init_refmodel()], [get_refmodel()]}{For setting up an object
#'   containing information about the reference model, the submodels, and how
#'   the projection should be carried out. Explicit calls to [init_refmodel()]
#'   and [get_refmodel()] are only rarely needed.}
#'   \item{[varsel()], [cv_varsel()]}{For running the *search* part and the
#'   *evaluation* part for a projection predictive variable selection, possibly
#'   with cross-validation (CV).}
#'   \item{[summary.vsel()], [print.vsel()], [plot.vsel()],
#'   [suggest_size.vsel()], [ranking()], [cv_proportions()],
#'   [plot.cv_proportions()], [performances()]}{For post-processing the results
#'   from [varsel()] and [cv_varsel()].}
#'   \item{[project()]}{For projecting the reference model onto submodel(s).
#'   Typically, this follows the variable selection, but it can also be applied
#'   directly (without a variable selection).}
#'   \item{[as.matrix.projection()] and [as_draws_matrix.projection()]}{For
#'   extracting projected parameter draws.}
#'   \item{[proj_linpred()], [proj_predict()]}{For making predictions from a
#'   submodel (after projecting the reference model onto it).}
#' }
#'
"_PACKAGE"

#' Internal global options
#'
#' The following global options are for internal use:
#' * `projpred.mssg_ndraws`, `projpred.mssg_cut_search`, `projpred.mssg_time`,
#' `projpred.warn_wobs_ppd`, `projpred.warn_additive_experimental`,
#' `projpred.warn_allrandom_dis`, `projpred.warn_instable_projections`,
#' `projpred.warn_cvrefbuilder_NULL`, `projpred.warn_kfold_refits`: A single
#' logical value indicating whether to throw certain messages or warnings
#' (depending on the midfix `mssg` or `warn`, respectively). For the exact
#' meaning of these global options, see their occurrences in the codebase. With
#' the exception of `projpred.warn_allrandom_dis`, these global options are
#' currently used in the unit tests to deactivate these messages and warnings.
#' Global option `projpred.warn_instable_projections` is also used (invisibly)
#' in the latent vignette to suppress the corresponding warnings while
#' illustrating the underlying issue (instable projections).
#' * `projpred.additional_checks`: A single logical value indicating whether to
#' run some additional checks that are not necessary to be run when users call
#' the corresponding \pkg{projpred} functions. Currently, these checks are
#' activated during the unit tests.
#' * `projpred.glm_fitter`: A character string naming the function to be used as
#' the submodel fitter for non-multilevel, non-additive projections. Currently,
#' this is an experimental feature and allowed values are
#' `"fit_glm_ridge_callback"` (the default) and `"fit_glm_callback"`.
#' * `projpred.gaussian_not_as_generalized`: A single logical value indicating
#' whether to treat the [gaussian()] family not as a family for a *generalized
#' linear* model (i.e., for which [glm()] would typically be used as a model
#' fitting function outside of \pkg{projpred}), but as the family for an
#' explicit *linear* model (i.e., for which [lm()] would typically be used as a
#' model fitting function outside of \pkg{projpred}). This also holds for models
#' with multilevel terms (because \pkg{lme4} offers both [lme4::glmer()] and
#' [lme4::lmer()]). Currently, this is an experimental feature.
#' * `projpred.PQL`: A single logical value indicating whether to use
#' [MASS::glmmPQL()] as the submodel fitter for multilevel (non-additive)
#' projections (see GitHub issue
#' [#207](https://github.com/stan-dev/projpred/issues/207) and GitHub pull
#' request [#353](https://github.com/stan-dev/projpred/pull/353)). Currently,
#' this is an experimental feature.
#' * `projpred.search_gauss_crossprod`: A single logical value indicating
#' whether to compute the forward search and the evaluation in the LOO CV with
#' validated search for Gaussian submodels with the identity link (including
#' the latent projection) from cross-products of the model matrix (default
#' `TRUE`). With `FALSE`, \pkg{projpred} projects onto each candidate submodel
#' draw by draw. The cross-products are used only with a positive `regul`.
#' * `projpred.categ_mlvl_mblogit`: A single logical value indicating whether to
#' fit the multilevel submodels of the augmented-data projection for the
#' [brms::categorical()] family with [mclogit::mblogit()] (`TRUE`) or with
#' \pkg{projpred}'s own implementation of the PQL algorithm of
#' [mclogit::mblogit()] (`FALSE`, the default).
#'
#' @name internal-global-options
#' @keywords internal
NULL
