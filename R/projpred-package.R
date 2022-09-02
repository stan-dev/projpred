#' Projection predictive feature selection
#'
#' @useDynLib projpred
#' @importFrom Rcpp sourceCpp
#'
#' @import stats
#' @import ggplot2
#' @importFrom rlang .data
#' @importFrom rstantools posterior_linpred
#' @importFrom loo kfold
#'
#' @description
#'
#' \pkg{projpred} is an \R package for performing a projection predictive
#' variable (or "feature") selection for generalized linear models (GLMs),
#' generalized linear multilevel (or "mixed") models (GLMMs), generalized
#' additive models (GAMs), and generalized additive multilevel (or "mixed")
#' models (GAMMs), with the support for additive models still being
#' experimental. Note that the term "generalized" includes the Gaussian family
#' as well. Some regression models whose response family is not an exponential
#' family are supported by the (currently still experimental) augmented-data
#' projection (see below).
#'
#' The package is compatible with \pkg{rstanarm} and \pkg{brms}, but developers
#' of other packages are welcome to add new [get_refmodel()] methods (which
#' enable the compatibility of their packages with \pkg{projpred}). Custom
#' reference models can also be used via [init_refmodel()]. It is via custom
#' reference models that \pkg{projpred} supports the projection onto candidate
#' models whose predictor terms are not a subset of the reference model's
#' predictor terms. However, for \pkg{rstanarm} and \pkg{brms} reference models,
#' \pkg{projpred} only supports the projection onto *submodels* of the reference
#' model. For the sake of simplicity, throughout this package, we use the term
#' "submodel" for all kinds of candidate models onto which the reference model
#' is projected, even though this term is not always appropriate for custom
#' reference models.
#'
#' In the following and throughout the rest of \pkg{projpred}'s documentation,
#' the term "traditional projection" will be used whenever the projection type
#' is not the augmented-data projection (see below for a description of the
#' augmented-data projection).
#'
#' The families supported by the traditional projection are [gaussian()],
#' [binomial()] (and---via [brms::get_refmodel.brmsfit()]---also
#' [brms::bernoulli()]), as well as [poisson()].
#'
#' As soon as possible, a reference for the augmented-data projection will be
#' provided here. For now, [this GitHub
#' issue](https://github.com/stan-dev/projpred/issues/70) provides some basic
#' information. The families supported by the augmented-data projection are
#' [binomial()] (again also [brms::bernoulli()], but note that currently, the
#' support for the [binomial()] family does not include binomial distributions
#' with more than one trial; in such a case, a workaround is to de-aggregate the
#' Bernoulli trials which belong to the same (aggregated) observation, i.e., to
#' use a "long" dataset), [brms::cumulative()], [rstanarm::stan_polr()] fits,
#' and [brms::categorical()]. See [extend_family()] (which is called by
#' [init_refmodel()]) for an explanation how to apply the augmented-data
#' projection to custom reference models. For non-custom reference models (i.e.,
#' those created by [get_refmodel.stanreg()] or [brms::get_refmodel.brmsfit()]),
#' the augmented-data projection is applied automatically if the family is
#' supported by the augmented-data projection and neither [binomial()] nor
#' [brms::bernoulli()]. For applying the augmented-data projection to the
#' [binomial()] (or [brms::bernoulli()]) family, see [extend_family()] as well
#' as [augdat_link_binom()] and [augdat_ilink_binom()]. Note that the
#' augmented-data projection is currently considered as an experimental feature
#' since it has not been tested thoroughly yet and is also subject to some more
#' theoretical investigations.
#'
#' For the projection of the reference model onto a submodel, \pkg{projpred}
#' currently relies on the following functions (in other words, these are the
#' workhorse functions used by the default divergence minimizers):
#' * Submodel without multilevel or additive terms:
#'     + For the traditional projection (or the augmented-data projection in
#'     case of the [binomial()] or [brms::bernoulli()] family): An internal C++
#'     function which basically serves the same purpose as [lm()] for the
#'     [gaussian()] family and [glm()] for all other families.
#'     + For the augmented-data projection: [MASS::polr()] for the
#'     [brms::cumulative()] family or [rstanarm::stan_polr()] fits,
#'     [nnet::multinom()] for the [brms::categorical()] family.
#' * Submodel with multilevel but no additive terms:
#'     + For the traditional projection (or the augmented-data projection in
#'     case of the [binomial()] or [brms::bernoulli()] family): [lme4::lmer()]
#'     for the [gaussian()] family, [lme4::glmer()] for all other families.
#'     + For the augmented-data projection: [ordinal::clmm()] for the
#'     [brms::cumulative()] family, [mclogit::mblogit()] for the
#'     [brms::categorical()] family.
#' * Submodel without multilevel but additive terms: [mgcv::gam()].
#' * Submodel with multilevel and additive terms: [gamm4::gamm4()].
#'
#' The projection of the reference model onto a submodel can be run on multiple
#' CPU cores in parallel (across the projected draws). This is powered by the
#' \pkg{foreach} package. Thus, any parallel (or sequential) backend compatible
#' with \pkg{foreach} can be used, e.g., the backends from packages
#' \pkg{doParallel}, \pkg{doMPI}, or \pkg{doFuture}. Using the global option
#' `projpred.prll_prj_trigger`, the number of projected draws below which no
#' parallelization is applied (even if a parallel backend is registered) can be
#' modified. Such a "trigger" threshold exists because of the computational
#' overhead of a parallelization which makes parallelization only useful for a
#' sufficiently large number of projected draws. By default, parallelization is
#' turned off, which can also be achieved by supplying `Inf` (or `NULL`) to
#' option `projpred.prll_prj_trigger`. Note that we cannot recommend
#' parallelizing the projection on Windows because in our experience, the
#' parallelization overhead is larger there, causing a parallel run to take
#' longer than a sequential run. Also note that the parallelization works well
#' for GLMs, but for all other models, the fitted model objects are quite big,
#' which---when running in parallel---may lead to excessive memory usage which
#' in turn may crash the R session. Thus, we currently cannot recommend the
#' parallelization for models other than GLMs.
#'
#' The [vignettes](https://mc-stan.org/projpred/articles/) (currently, there is
#' only a single one) illustrate how to use the \pkg{projpred} functions in
#' conjunction. Shorter examples are included here in the documentation.
#'
#' Some references relevant for this package are given in section "References"
#' below. See `citation(package = "projpred")` for details on citing
#' \pkg{projpred}.
#'
#' @details
#'
#' # Functions
#'
#' \describe{
#'   \item{[init_refmodel()], [get_refmodel()]}{For setting up a reference model
#'   (only rarely needed explicitly).}
#'   \item{[varsel()], [cv_varsel()]}{For variable selection, possibly with
#'   cross-validation (CV).}
#'   \item{[summary.vsel()], [print.vsel()], [plot.vsel()],
#'   [suggest_size.vsel()], [solution_terms.vsel()]}{For post-processing the
#'   results from the variable selection.}
#'   \item{[project()]}{For projecting the reference model onto submodel(s).
#'   Typically, this follows the variable selection, but it can also be applied
#'   directly (without a variable selection).}
#'   \item{[as.matrix.projection()]}{For extracting projected parameter draws.}
#'   \item{[proj_linpred()], [proj_predict()]}{For making predictions from a
#'   submodel (after projecting the reference model onto it).}
#' }
#'
#' @references
#'
#' Catalina, A., Bürkner, P.-C., and Vehtari, A. (2020). Projection predictive
#' inference for generalized linear and additive multilevel models.
#' *arXiv:2010.06994*. URL: <https://arxiv.org/abs/2010.06994>.
#'
#' Dupuis, J. A. and Robert, C. P. (2003). Variable selection in qualitative
#' models via an entropic explanatory power. *Journal of Statistical Planning
#' and Inference*, **111**(1-2):77–94. \doi{10.1016/S0378-3758(02)00286-0}.
#'
#' Goutis, C. and Robert, C. P. (1998). Model choice in generalised linear
#' models: A Bayesian approach via Kullback–Leibler projections. *Biometrika*,
#' **85**(1):29–37.
#'
#' Piironen, J. and Vehtari, A. (2017). Comparison of Bayesian predictive
#' methods for model selection. *Statistics and Computing*, **27**(3):711-735.
#' \doi{10.1007/s11222-016-9649-y}.
#'
#' Piironen, J., Paasiniemi, M., and Vehtari, A. (2020). Projective inference in
#' high-dimensional problems: Prediction and feature selection. *Electronic
#' Journal of Statistics*, **14**(1):2155-2197. \doi{10.1214/20-EJS1711}.
#'
"_PACKAGE"
