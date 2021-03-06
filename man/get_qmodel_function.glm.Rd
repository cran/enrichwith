% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/enrich.glm.R
\name{get_qmodel_function.glm}
\alias{get_qmodel_function.glm}
\title{Function to compute/extract a \code{qmodel} function}
\usage{
\method{get_qmodel_function}{glm}(object, ...)
}
\arguments{
\item{object}{an object of class \code{glm} or\code{enriched_glm}}

\item{...}{currently not used}
}
\description{
Function to compute/extract a \code{qmodel} function
}
\details{
The computed/extracted function has arguments
\describe{

\item{p}{a vector of probabilities with \code{length(p)} equal to
\code{nrow(data)} at which to evaluate quantiles}

\item{data}{a data frame with observations at which to compute the
quantiles. If missing then quantiles are computed at the model
frame extracted from the object (see \code{\link{glm}})}

\item{coefficients}{the regression coefficients at which the
quantiles are computed. If missing then the maximum likelihood
estimates are used}

\item{dispersion}{the dispersion parameter at which the
quantiles are computed. If missing then the maximum
likelihood estimate is used}

\item{log.p}{logical; if \code{TRUE}, the logarithm of the
probabilities is used}

\item{lower.tail}{logical; if \code{TRUE} (default), probabilities
are P[X <= x] otherwise, P[X > x]}

}
}
