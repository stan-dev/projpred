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
#' For the projection of the reference model onto a submodel, \pkg{projpred}
#' currently relies on the following functions (in other words, these are the
#' workhorse functions used by the default divergence minimizers):
#' * Submodel without multilevel or additive terms:
#'     + For the traditional (or latent) projection (or the augmented-data
#'     projection in case of the [binomial()] or [brms::bernoulli()] family): An
#'     internal C++ function which basically serves the same purpose as [lm()]
#'     for the [gaussian()] family and [glm()] for all other families.
#'     + For the augmented-data projection: [MASS::polr()] for the
#'     [brms::cumulative()] family or [rstanarm::stan_polr()] fits,
#'     [nnet::multinom()] for the [brms::categorical()] family.
#' * Submodel with multilevel but no additive terms:
#'     + For the traditional (or latent) projection (or the augmented-data
#'     projection in case of the [binomial()] or [brms::bernoulli()] family):
#'     [lme4::lmer()] for the [gaussian()] family, [lme4::glmer()] for all other
#'     families.
#'     + For the augmented-data projection: [ordinal::clmm()] for the
#'     [brms::cumulative()] family, [mclogit::mblogit()] for the
#'     [brms::categorical()] family.
#' * Submodel without multilevel but additive terms: [mgcv::gam()].
#' * Submodel with multilevel and additive terms: [gamm4::gamm4()].
#'
#' Setting the global option `projpred.extra_verbose` to `TRUE` will print out
#' which submodel \pkg{projpred} is currently projecting onto as well as (if
#' `method = "forward"` and `verbose = TRUE` in `varsel()` or `cv_varsel()`)
#' which submodel has been selected at those steps of the forward search for
#' which a percentage (of the maximum submodel size that the search is run up
#' to) is printed. In general, however, we cannot recommend setting this global
#' option to `TRUE` for `cv_varsel()` with `validate_search = TRUE` (simply due
#' to the amount of information that will be printed, but also due to the
#' progress bar which will not work anymore as intended).
#'
#' The projection of the reference model onto a submodel can be run in parallel
#' (across the projected draws). This is powered by the \pkg{foreach} package.
#' Thus, any parallel (or sequential) backend compatible with \pkg{foreach} can
#' be used, e.g., the backends from packages \pkg{doParallel}, \pkg{doMPI}, or
#' \pkg{doFuture}. Using the global option `projpred.prll_prj_trigger`, the
#' number of projected draws below which no parallelization is applied (even if
#' a parallel backend is registered) can be modified. Such a "trigger" threshold
#' exists because of the computational overhead of a parallelization which makes
#' the projection parallelization only useful for a sufficiently large number of
#' projected draws. By default, the projection parallelization is turned off,
#' which can also be achieved by supplying `Inf` (or `NULL`) to option
#' `projpred.prll_prj_trigger`. Note that we cannot recommend the projection
#' parallelization on Windows because in our experience, the parallelization
#' overhead is larger there, causing a parallel run to take longer than a
#' sequential run. Also note that the projection parallelization works well for
#' submodels which are GLMs (and hence also for the latent projection if the
#' submodel has no multilevel or additive predictor terms), but for all other
#' types of submodels, the fitted submodel objects are quite big, which---when
#' running in parallel---may lead to excessive memory usage which in turn may
#' crash the R session (on Unix systems, setting an appropriate memory limit via
#' [unix::rlimit_as()] may avoid crashing the whole machine). Thus, we currently
#' cannot recommend parallelizing projections onto submodels which are GLMs (in
#' this context, the latent projection onto a submodel without multilevel and
#' without additive terms may be regarded as a projection onto a submodel which
#' is a GLM). However, for [cv_varsel()], there is also a *CV* parallelization
#' (i.e., a parallelization of \pkg{projpred}'s cross-validation) which can be
#' activated via argument `parallel`.
#'
#' In case of multilevel models, \pkg{projpred} offers two global options for
#' "integrating out" group-level effects: `projpred.mlvl_pred_new` and
#' `projpred.mlvl_proj_ref_new`. When setting `projpred.mlvl_pred_new` to `TRUE`
#' (default is `FALSE`), then at
#' *prediction* time, \pkg{projpred} will treat group levels existing in the
#' training data as *new* group levels, implying that their group-level effects
#' are drawn randomly from a (multivariate) Gaussian distribution. This concerns
#' both, the reference model and the (i.e., any) submodel. Furthermore, setting
#' `projpred.mlvl_pred_new` to `TRUE` causes `as.matrix.projection()` to omit
#' the projected group-level effects (for the group levels from the original
#' dataset). When setting `projpred.mlvl_proj_ref_new` to `TRUE` (default is
#' `FALSE`), then at *projection* time, the reference model's fitted values
#' (that the submodels fit to) will be computed by treating the group levels
#' from the original dataset as *new* group levels, implying that their
#' group-level effects will be drawn randomly from a (multivariate) Gaussian
#' distribution (as long as the reference model is a multilevel model,
#' which---for custom reference models---does not need to be the case). This
#' also affects the latent response values for a latent projection
#' correspondingly. Setting `projpred.mlvl_pred_new` to `TRUE` makes sense,
#' e.g., when the prediction task is such that any group level will be treated
#' as a new one. Typically, setting `projpred.mlvl_proj_ref_new` to `TRUE` only
#' makes sense when `projpred.mlvl_pred_new` is already set to `TRUE`. In that
#' case, the default of `FALSE` for `projpred.mlvl_proj_ref_new` ensures that at
#' projection time, the submodels fit to the best possible fitted values from
#' the reference model, and setting `projpred.mlvl_proj_ref_new` to `TRUE` would
#' make sense if the group-level effects should be integrated out completely.
#'
#' @details
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
#'   [plot.cv_proportions()]}{For post-processing the results from [varsel()]
#'   and [cv_varsel()].}
#'   \item{[project()]}{For projecting the reference model onto submodel(s).
#'   Typically, this follows the variable selection, but it can also be applied
#'   directly (without a variable selection).}
#'   \item{[as.matrix.projection()]}{For extracting projected parameter draws.}
#'   \item{[proj_linpred()], [proj_predict()]}{For making predictions from a
#'   submodel (after projecting the reference model onto it).}
#' }
#'
"_PACKAGE"
