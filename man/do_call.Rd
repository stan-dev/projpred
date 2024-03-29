% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/misc.R
\name{do_call}
\alias{do_call}
\title{Execute a function call}
\usage{
do_call(what, args, pkg = NULL)
}
\arguments{
\item{what}{Either a function or a non-empty character string naming the
function to be called.}

\item{args}{A \code{list} of arguments to the function call. The \code{\link{names}}
attribute of \code{args} gives the argument names.}

\item{pkg}{Optional name of the package in which to search for the function
if \code{what} is a character string.}
}
\value{
The result of the (evaluated) function call.
}
\description{
Execute a function call similar to \code{\link[=do.call]{do.call()}}, but without deparsing
function arguments.
}
\keyword{internal}
