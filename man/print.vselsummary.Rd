% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/methods.R
\name{print.vselsummary}
\alias{print.vselsummary}
\title{Print summary of a \code{\link[=varsel]{varsel()}} or \code{\link[=cv_varsel]{cv_varsel()}} run}
\usage{
\method{print}{vselsummary}(x, digits = getOption("projpred.digits", 2), ...)
}
\arguments{
\item{x}{An object of class \code{vselsummary}.}

\item{digits}{Passed to \code{\link[=print.data.frame]{print.data.frame()}} (for the table containing the
submodel performance evaluation results) and \code{\link[=print.default]{print.default()}} (for the
vector containing the reference model performance evaluation results).}

\item{...}{Arguments passed to \code{\link[=print.data.frame]{print.data.frame()}} (for the table containing
the submodel performance evaluation results) and \code{\link[=print.default]{print.default()}} (for the
vector containing the reference model performance evaluation results).}
}
\value{
The output of \code{\link[=summary.vsel]{summary.vsel()}} (invisible).
}
\description{
This is the \code{\link[=print]{print()}} method for summary objects created by \code{\link[=summary.vsel]{summary.vsel()}}.
It displays a summary of the results from a \code{\link[=varsel]{varsel()}} or \code{\link[=cv_varsel]{cv_varsel()}} run.
}
\details{
In the submodel predictive performance table printed at (or towards)
the bottom, column \code{ranking_fulldata} contains the full-data predictor
ranking and column \code{cv_proportions_diag} contains the main diagonal of the
matrix returned by \code{\link[=cv_proportions]{cv_proportions()}} (with \code{cumulate} as set in the
\code{\link[=summary.vsel]{summary.vsel()}} call that created \code{x}). To retrieve the fold-wise
predictor rankings, use the \code{\link[=ranking]{ranking()}} function, possibly followed by
\code{\link[=cv_proportions]{cv_proportions()}} for computing the ranking proportions (which can be
visualized by \code{\link[=plot.cv_proportions]{plot.cv_proportions()}}).
}
