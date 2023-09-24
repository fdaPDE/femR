#include <RcppEigen.h>
// [[Rcpp::depends(RcppEigen)]]
#include "R_Mesh.h"

// Rcpp module definition
RCPP_MODULE(Mesh_2D_ORDER_1) {
    Rcpp::class_<Mesh_2D_ORDER_1>("Mesh_2D_ORDER_1")
      .constructor<Rcpp::List>()
      .method("nodes", &Mesh_2D_ORDER_1::nodes)
      .method("elements", &Mesh_2D_ORDER_1::elements)
      .method("neighbors", &Mesh_2D_ORDER_1::neighbors)
      .method("boundary", &Mesh_2D_ORDER_1::boundary)
      .method("fe_order", &Mesh_2D_ORDER_1::fe_order)
      .method("eval",     &Mesh_2D_ORDER_1::eval);
}

RCPP_MODULE(Mesh_2D_ORDER_2) {
    Rcpp::class_<Mesh_2D_ORDER_2>("Mesh_2D_ORDER_2")
      .constructor<Rcpp::List>()
      .method("nodes", &Mesh_2D_ORDER_2::nodes)
      .method("elements", &Mesh_2D_ORDER_2::elements)
      .method("neighbors", &Mesh_2D_ORDER_2::neighbors)
      .method("boundary", &Mesh_2D_ORDER_2::boundary)
      .method("fe_order", &Mesh_2D_ORDER_2::fe_order)
      .method("eval",     &Mesh_2D_ORDER_2::eval);
}
