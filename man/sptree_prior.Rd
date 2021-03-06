% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sptree.R
\name{sptree_prior}
\alias{sptree_prior}
\title{Compute the spanning tree prior of a solution over a range of k value}
\usage{
sptree_prior(sol, k_max)
}
\arguments{
\item{sol}{a gtclust bayesian clustering solution}

\item{k_max}{the maximum value of k for which the prior should be calculated}
}
\value{
the log of the prior
}
\description{
This function compute the spanning tree prior of a solution
}
