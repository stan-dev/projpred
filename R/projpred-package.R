#' Projection predictive feature selection
#'
#' @useDynLib projpred
#' @importFrom Rcpp sourceCpp
#'
#' @import stats
#' @import ggplot2
#' @importFrom loo psis
#' @importFrom rlang .data
#'
#' @description
#'
#' \pkg{projpred} is an \R package for performing a projection predictive
#' variable (or "feature") selection for generalized linear models (GLMs),
#' generalized linear multilevel (or "mixed") models (GLMMs), generalized
#' additive models (GAMs), and generalized additive multilevel (or "mixed")
#' models (GAMMs). Note that the term "generalized" includes the Gaussian family
#' as well.
#'
#' The package is compatible with \pkg{rstanarm} and \pkg{brms}, but developers
#' of other packages are welcome to add new [get_refmodel()] methods (which
#' enable the compatibility of their packages with \pkg{projpred}). Custom
#' reference models can also be used via [init_refmodel()].
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
#' registered). This option exists because of the computational overhead of a
#' parallelization. A value of `Inf` (the default) turns off parallelization.
#'
#' See the vignettes
#' (\href{https://mc-stan.org/projpred/articles/quickstart.html}{quickstart-vignette}
#' and
#' \href{https://mc-stan.org/projpred/articles/quickstart_glmm.html}{quickstart-glmm-vignette})
#' for example applications. Shorter examples are included here in the documentation.
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
#'   \item{[varsel()], [cv_varsel()]}{Perform the variable selection, possibly
#'   with cross-validation (CV).}
#'   \item{[summary.vsel()], [print.vsel()], [plot.vsel()],
#'   [suggest_size.vsel()], [solution_terms.vsel()]}{Post-process the results
#'   from the variable selection.}
#'   \item{[project()]}{Project the reference model onto submodel(s). Typically,
#'   this follows the variable selection, but it can also be applied directly
#'   (without a variable selection).}
#'   \item{[as.matrix.projection()]}{Extract projected parameter draws.}
#'   \item{[proj_linpred()], [proj_predict()]}{Make predictions from a submodel
#'   (after projecting the reference model onto it).}
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
#' and Inference*, **111**(1-2):77–94. DOI:
#' [10.1016/S0378-3758(02)00286-0](https://doi.org/10.1016/S0378-3758(02)00286-0).
#'
#' Piironen, J. and Vehtari, A. (2017). Comparison of Bayesian predictive
#' methods for model selection. *Statistics and Computing*, **27**(3):711-735.
#' DOI: [10.1007/s11222-016-9649-y](https://doi.org/10.1007/s11222-016-9649-y).
#'
#' Piironen, J., Paasiniemi, M., and Vehtari, A. (2020). Projective inference in
#' high-dimensional problems: Prediction and feature selection. *Electronic
#' Journal of Statistics*, **14**(1):2155-2197. DOI:
#' [10.1214/20-EJS1711](https://doi.org/10.1214/20-EJS1711).
#'
#' Catalina, A., Bürkner, P.-C., and Vehtari, A. (2020). Projection predictive
#' inference for generalized linear and additive multilevel models.
#' *arXiv:2010.06994*. URL: <https://arxiv.org/abs/2010.06994>.
#'
"_PACKAGE"
