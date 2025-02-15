% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/basic_functions.R
\name{getCondRelated}
\alias{getCondRelated}
\title{get the genes connected to a TARGET GENE under each condition.}
\usage{
getCondRelated(gs, pcor0 = 1e-04, partcorr.list, conditions = NULL)
}
\arguments{
\item{gs}{character, 'gene selected' or 'TARGET GENE': e.g. "GAPDH", "CD4",...}

\item{pcor0}{numeric, the threshold for partial correlations. values greater than pcor0 will be kept.}

\item{partcorr.list}{list of partial correlation matrices}

\item{conditions}{a vector of names for each condition. e.g.: if the list has two matrices, here conditions = c("ConditionA","ConditionB")}
}
\value{
a table of partial correlations for genes connected to TARGET GENE under all conditions
}
\description{
get the genes connected to a TARGET GENE under each condition.
}
