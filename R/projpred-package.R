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
#' as well.
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
#' Currently, the supported families are [gaussian()], [binomial()] (and---via
#' [brms::get_refmodel.brmsfit()]---also [brms::bernoulli()]), as well as
#' [poisson()].
#'
#' The projection of the reference model onto a submodel can be run on multiple
#' CPU cores in parallel (across the projected draws). This is powered by the
#' \pkg{foreach} package. Thus, you can use any parallel (or sequential) backend
#' compatible with \pkg{foreach}, e.g., the backends from packages
#' \pkg{doParallel}, \pkg{doMPI}, or \pkg{doFuture}. Using the global option
#' `projpred.prll_prj_trigger`, you can modify the number of projected draws
#' below which no parallelization is used (even if a parallel backend is
#' registered). Such a "trigger" threshold exists because of the computational
#' overhead of a parallelization which makes parallelization only useful for a
#' sufficiently large number of projected draws. By default, parallelization is
#' turned off, which can also be achieved by supplying `Inf` (or `NULL`) to
#' option `projpred.prll_prj_trigger`. Note that we cannot recommend
#' parallelizing the projection on Windows because in our experience, the
#' parallelization overhead is larger there, causing a parallel run to take
#' longer than a sequential run. Also note that the parallelization works well
#' for GLMs, but for GLMMs, GAMs, and GAMMs, the fitted model objects are quite
#' big, which---when running in parallel---may lead to an excessive memory usage
#' which in turn may crash the R session. Thus, we currently cannot recommend
#' the parallelization for GLMMs, GAMs, and GAMMs.
#'
#' See the [vignette](https://mc-stan.org/projpred/articles/projpred.html) for
#' an example application (which also demonstrates the use of parallelization).
#' Shorter examples are included here in the documentation.
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
#' Goutis, C. and Robert, C. P. (1998). Model choice in generalised linear
#' models: A Bayesian approach via Kullback–Leibler projections. *Biometrika*,
#' **85**(1):29–37.
#'
#' Dupuis, J. A. and Robert, C. P. (2003). Variable selection in qualitative
#' models via an entropic explanatory power. *Journal of Statistical Planning
#' and Inference*, **111**(1-2):77–94. \doi{10.1016/S0378-3758(02)00286-0}.
#'
#' Piironen, J. and Vehtari, A. (2017). Comparison of Bayesian predictive
#' methods for model selection. *Statistics and Computing*, **27**(3):711-735.
#' \doi{10.1007/s11222-016-9649-y}.
#'
#' Piironen, J., Paasiniemi, M., and Vehtari, A. (2020). Projective inference in
#' high-dimensional problems: Prediction and feature selection. *Electronic
#' Journal of Statistics*, **14**(1):2155-2197. \doi{10.1214/20-EJS1711}.
#'
#' Catalina, A., Bürkner, P.-C., and Vehtari, A. (2020). Projection predictive
#' inference for generalized linear and additive multilevel models.
#' *arXiv:2010.06994*. URL: <https://arxiv.org/abs/2010.06994>.
#'
"_PACKAGE"
