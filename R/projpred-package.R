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
#' variable (or "feature") selection for (generalized) linear models ((G)LMs),
#' (generalized) linear mixed (or "multilevel") models ((G)LMMs), (generalized)
#' additive models ((G)AMs) and (generalized) additive mixed (or "multilevel")
#' models ((G)AMMs).
#'
#' The package is compatible with \pkg{rstanarm} and \pkg{brms}, but developers
#' of other packages are welcome to add more methods for the [get_refmodel()]
#' generic. Custom reference models can also be used via [init_refmodel()].
#'
#' Currently, the supported families are [gaussian()], [binomial()] (as well as
#' [brms::bernoulli()]), and [poisson()].
#'
#' See the vignettes
#' (\href{https://mc-stan.org/projpred/articles/quickstart.html}{quickstart-vignette}
#' and
#' \href{https://mc-stan.org/projpred/articles/quickstart_glmm.html}{quickstart-glmm-vignette})
#' for example applications. Smaller examples are included here in the
#' documentation.
#'
#' @details # Functions
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
#'   \item{[proj_linpred()], [proj_predict()]}{Make predictions from a submodel
#'   (after projection).}
#' }
#'
#' @details # References
#'
#' See the `CITATION` file for details on citation.
#'
#' * Goutis, C. and Robert, C. P. (1998). Model choice in generalised linear
#' models: a Bayesian approach via Kullback–Leibler projections. *Biometrika*,
#' 85(1):29–37.
#'
#' * Dupuis, J. A. and Robert, C. P. (2003). Variable selection in qualitative
#' models via an entropic explanatory power. *Journal of Statistical Planning
#' and Inference*, 111(1-2):77–94.
#' doi:[10.1016/S0378-3758(02)00286-0](https://doi.org/10.1016/S0378-3758(02)00286-0).
#'
#' * Piironen, J. and Vehtari, A. (2017). Comparison of Bayesian predictive
#' methods for model selection. *Statistics and Computing*, 27(3):711-735.
#' doi:[10.1007/s11222-016-9649-y](https://doi.org/10.1007/s11222-016-9649-y).
#'
#' * Piironen, J., Paasiniemi, M., and Vehtari, A. (2020). Projective inference
#' in high-dimensional problems: Prediction and feature selection. *Electronic
#' Journal of Statistics*, 14(1):2155-2197.
#' doi:[10.1214/20-EJS1711](https://doi.org/10.1214/20-EJS1711).
#'
#' * Catalina, A., Bürkner, P.-C., and Vehtari, A. (2020). Projection predictive
#' inference for generalized linear and additive multilevel models.
#' *arxiv:2010.06994*. URL: <https://arxiv.org/abs/2010.06994>.
#'
"_PACKAGE"
