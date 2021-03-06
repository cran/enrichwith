% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/enrich.lm.R
\name{enrich.lm}
\alias{enrich.lm}
\title{Enrich objects of class \code{\link{lm}}}
\usage{
\method{enrich}{lm}(object, with = "all", ...)
}
\arguments{
\item{object}{an object of class lm}

\item{with}{a character vector of options for the enrichment of \code{object}}

\item{...}{extra arguments to be passed to the
\code{compute_*} functions}
}
\value{
The object \code{object} of class \code{\link{lm}} with extra
components. \code{get_enrichment_options.lm()} returns the
components and their descriptions.
}
\description{
Enrich objects of class \code{\link{lm}} with any or all of a set
auxiliary functions, the maximum likelihood estimate of the
dispersion parameter, the expected or observed information at the
maximum likelihood estimator, and the first term in the expansion
of the bias of the maximum likelihood estimator.
}
\details{
The auxiliary functions consist of the score functions, the
expected or observed information, the first-order bias of the
maximum likelihood estimator as functions of the model parameters,
and a \code{simulate} function that takes as input the model
parameters (including the dispersion if any). The result from the
\code{simulate} auxiliary function has the same structure to that
of the \code{\link{simulate}} method for \code{\link{lm}} objects.
}
