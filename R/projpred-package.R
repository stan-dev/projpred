#' Projection predictive feature selection
#'
#' @docType package
#' @name projpred
#'
#' @useDynLib projpred
#' @importFrom Rcpp sourceCpp
#'
#' @import stats
#' @import ggplot2
#' @importFrom loo psis
#' 
#' 
#' @description Description
#' 
#' \pkg{projpred} is an R package to perform projection predictive variable
#'   (feature) selection for generalized linear models, generalized linear
#'   multilevel models and generalized additive multilevel models. The package
#'   is aimed to be compatible with \pkg{rstanarm} but also other reference
#'   models can be used (see function \code{\link{init_refmodel}}).
#'
#' Currently, the supported models (family objects in R) include Gaussian,
#'   Binomial and Poisson families, but more will be implemented later. See the
#'   \href{https://mc-stan.org/projpred/articles/quickstart.html}{quickstart-vignette}
#'   and
#'   \href{https://mc-stan.org/projpred/articles/quickstart-glmm.html}{quickstart-glmm-vignette}
#'   for examples.
#' 
#' 
#' @section Functions:
#' 
#' \describe{
#'  \item{\link{varsel}, \link{cv_varsel}, \link{init_refmodel},
#'   \link{suggest_size}}{ Perform and cross-validate the variable selection.
#'   \link{init_refmodel} can be used to initialize a reference model other than
#'   \pkg{rstanarm}-fit.} \item{\link{project}}{ Get the projected posteriors of
#'   the reduced models.} \item{\link{proj_predict}, \link{proj_linpred}}{ Make
#'   predictions with reduced number of features.} \item{\link{plot},
#'   \link{summary}}{ Visualize and get some key statistics about the variable
#'   selection.}
#' }
#' 
#' 
#' 
#' 
#' @section References:
#' 
#' Dupuis, J. A. and Robert, C. P. (2003). Variable selection in qualitative
#'   models via an entropic explanatory power. \emph{Journal of Statistical
#'   Planning and Inference}, 111(1-2):77–94.
#'
#' Goutis, C. and Robert, C. P. (1998). Model choice in generalised linear
#'   models: a Bayesian approach via Kullback–Leibler projections.
#'   \emph{Biometrika}, 85(1):29–37.
#' 
#' Juho Piironen and Aki Vehtari (2017). Comparison of Bayesian predictive
#'   methods for model selection. \emph{Statistics and Computing},
#'   27(3):711-735. doi:10.1007/s11222-016-9649-y.
#'   (\href{https://link.springer.com/article/10.1007/s11222-016-9649-y}{Online}).
NULL
