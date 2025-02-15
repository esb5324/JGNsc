% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/simulation.R
\name{getCountList}
\alias{getCountList}
\title{generate a list of raw count matrices based on the list of covariance matrices}
\usage{
getCountList(sigma.list, nvec = c(500, 500), ngene = NULL, a3 = 2, b3 = 1)
}
\arguments{
\item{sigma.list}{a list of covariance matrices. Please make sure the positive definite property.}

\item{nvec}{a vector of sample size under each condition.}

\item{ngene}{integer. the number of genes}

\item{a3}{the parameter for the non-dropout rate pij~ beta(a3, b3)}

\item{b3}{the parameter for the non-dropout rate pij~ beta(a3, b3)}
}
\value{
a list of count matrices
}
\description{
generate a list of raw count matrices based on the list of covariance matrices
}
