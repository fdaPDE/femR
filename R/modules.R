#' @useDynLib femR
#' @import methods Rcpp
#' @import RcppEigen
#' @importFrom Rcpp evalCpp
#' @import plotly
#' @import Matrix
NULL

## load required modules
Rcpp::loadModule("cpp_2d_domain", TRUE)
Rcpp::loadModule("cpp_lagrange_basis_2d_fe1", TRUE)
Rcpp::loadModule("cpp_lagrange_basis_2d_fe2", TRUE)
Rcpp::loadModule("cpp_pde_2d_fe1",     TRUE)
Rcpp::loadModule("cpp_pde_2d_fe2",     TRUE)