% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/PlottingFunctions.R
\name{plotInt}
\alias{plotInt}
\title{Plot 2-way interaction matrix}
\usage{
plotInt(NLmod)
}
\arguments{
\item{NLmod}{The fitted model}
}
\description{
This function takes in an object which is a model from the NLint function
and plots the two-way interaction matrix
}
\examples{

n = 200
p = 10
pc = 1

sigma = matrix(0.3, p, p)
diag(sigma) = 1
X = rmvnorm(n, mean=rep(0,p), sigma = sigma)

C = matrix(rnorm(n*pc), nrow=n)

TrueH = function(X) {
  return(0.5*(X[,2]*X[,3]) - 0.6*(X[,4]^2 * X[,5]))
}

Y = 5 + C + TrueH(X) + rnorm(n)

NLmod = NLint(Y=Y, X=X, C=C)

plotInt(NLmod = NLmod)
}
