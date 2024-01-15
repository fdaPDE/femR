// This file is part of fdaPDE, a C++ library for physics-informed
// spatial and functional data analysis.
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.

#include <RcppEigen.h>
// [[Rcpp::depends(RcppEigen)]]
#include "include/r_pde.h"

using cpp_pde_2d_fe1 = R_PDE<2,2,1>;
RCPP_MODULE(cpp_pde_2d_fe1) {
    Rcpp::class_<R_PDE<2,2,1>>("cpp_pde_2d_fe1")
      .constructor<Rcpp::Environment, int, Rcpp::Nullable<Rcpp::List>>()
      .method("quadrature_nodes"     , &R_PDE<2,2,1>::quadrature_nodes     )
      .method("dofs_coordinates"     , &R_PDE<2,2,1>::dofs_coordinates     )
      .method("mass"                 , &R_PDE<2,2,1>::mass                 )
      .method("stiff"                , &R_PDE<2,2,1>::stiff                )
      .method("force"                , &R_PDE<2,2,1>::force                )
      .method("set_dirichlet_bc"     , &R_PDE<2,2,1>::set_dirichlet_bc     )
      .method("set_forcing"          , &R_PDE<2,2,1>::set_forcing          )
      .method("set_initial_condition", &R_PDE<2,2,1>::set_initial_condition)
      .method("init"                 , &R_PDE<2,2,1>::init                 )
      .method("solve"                , &R_PDE<2,2,1>::solve                )
      .method("solution"             , &R_PDE<2,2,1>::solution             );
}
using cpp_pde_2d_fe2 = R_PDE<2,2,2>;
RCPP_MODULE(cpp_pde_2d_fe2) {
    Rcpp::class_<R_PDE<2,2,2>>("cpp_pde_2d_fe2")
      .constructor<Rcpp::Environment, int, Rcpp::Nullable<Rcpp::List>>()
      .method("quadrature_nodes"     , &R_PDE<2,2,2>::quadrature_nodes     )
      .method("dofs_coordinates"     , &R_PDE<2,2,2>::dofs_coordinates     )
      .method("mass"                 , &R_PDE<2,2,2>::mass                 )
      .method("stiff"                , &R_PDE<2,2,2>::stiff                )
      .method("force"                , &R_PDE<2,2,2>::force                )
      .method("set_dirichlet_bc"     , &R_PDE<2,2,2>::set_dirichlet_bc     )
      .method("set_forcing"          , &R_PDE<2,2,2>::set_forcing          )
      .method("set_initial_condition", &R_PDE<2,2,2>::set_initial_condition)
      .method("init"                 , &R_PDE<2,2,2>::init                 )
      .method("solve"                , &R_PDE<2,2,2>::solve                )
      .method("solution"             , &R_PDE<2,2,2>::solution             );
}
