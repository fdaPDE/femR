% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/operators.R
\name{*,numeric,Function-method}
\alias{*,numeric,Function-method}
\title{product by scalar for FunctionObejct}
\usage{
\S4method{*}{numeric,Function}(e1, e2)
}
\arguments{
\item{e1}{a numeric.}

\item{e2}{a FunctioObject created by \code{Function}.}
}
\value{
A S4 object representing the reaction term of a second order linear differential operator.
}
\description{
product by scalar for FunctionObejct
}
\examples{
\dontrun{
library(femR)
data("unit_square")
mesh <- Mesh(unit_square)
Vh <- FunctionSpace(mesh)
f <- Function(Vh)
reaction <- 2*f
}
}
