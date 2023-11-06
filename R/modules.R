#' @useDynLib femR
#' @import methods Rcpp
#' @import RcppEigen
#' @importFrom Rcpp evalCpp
#' @import plotly
#' @import Matrix
NULL

# #' @exportPattern "^[[:alpha:]]+"


## load required modules
Rcpp::loadModule("PDE_2D_ORDER_1", TRUE)
Rcpp::loadModule("PDE_2D_ORDER_2", TRUE)
Rcpp::loadModule("Mesh_2D", TRUE)

