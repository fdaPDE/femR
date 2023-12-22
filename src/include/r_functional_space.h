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

#ifndef __R_FUNCTIONAL_SPACE_H__
#define __R_FUNCTIONAL_SPACE_H__

// we only include RcppEigen.h which pulls Rcpp.h in for us
#include <RcppEigen.h>
// [[Rcpp::depends(RcppEigen)]]

#include <fdaPDE/utils.h>
#include <fdaPDE/finite_elements.h>
#include <fdaPDE/mesh.h>
#include "r_mesh.h"
using fdapde::core::Mesh;
using fdapde::core::Integrator;
using fdapde::core::LagrangianBasis;
using fdapde::core::pointwise_evaluation;
using fdapde::core::areal_evaluation;
using fdapde::core::eval;
using fdapde::core::ct_binomial_coefficient;
using fdapde::core::FEM;
using fdapde::core::fem_order;
using fdapde::core::combinations;

// type-erasure wrapper for a function space object
struct I_FunctionalSpace {
    template <typename T>
    using fn_ptrs = fdapde::mem_fn_ptrs<
      &T::size, &T::operator(), &T::template eval<pointwise_evaluation>, &T::template eval<areal_evaluation>>;
    // forwardings
    decltype(auto) size() const { return fdapde::invoke<std::size_t, 0>(*this); }
    decltype(auto) operator()(const DVector<double>& c, const DMatrix<double>& locs) const {
      return fdapde::invoke<DVector<double>, 1>(*this, c, locs);
    }
    decltype(auto) eval(int evaluation_type, const DMatrix<double>& locs) const {
        using RetType = std::pair<SpMatrix<double>, DVector<double>>;
        switch (evaluation_type) {   // run-time switch based on evaluation strategy
        case eval::pointwise: // 0
            return fdapde::invoke<RetType, 2>(*this, locs).first;
        case eval::areal: // 1
            return fdapde::invoke<RetType, 3>(*this, locs).first;
        }
        return SpMatrix<double>{};
    }
};
using FunctionalSpace = fdapde::erase<fdapde::heap_storage, I_FunctionalSpace>;

// supported basis function types
enum space_type { fem_lagrange };

template <int M, int N, int R> class R_FunctionalSpace {
   private:
    using DomainType = Mesh<M, N>;
    using FunctionalSpaceType = FunctionalSpace;
    using QuadratureRule = Integrator<FEM, DomainType::local_dimension, R>; // exact for FEM spaces (TODO: generalize)

    // internal data
    DomainType domain_ {};
    FunctionalSpaceType fun_space_ {}; // wrapped functional space
    QuadratureRule integrator_ {};     // quadrature rule (exact for the provided fem order)
    int space_type_;
    std::size_t n_dofs_ = 0;        // degrees of freedom, i.e. the maximum ID in the dof_table_
    DMatrix<int> dofs_;             // for each element, the degrees of freedom associated to it
    DMatrix<int> boundary_dofs_;    // unknowns on the boundary of the domain, for boundary conditions prescription
    
    void enumerate_dofs(const DomainType& mesh);
   public:
    
    enum {
        fem_order = R,
        n_dof_per_element = ct_nnodes(DomainType::local_dimension, fem_order),
        n_dof_per_edge = fem_order - 1,
        n_dof_internal =
          n_dof_per_element - (DomainType::local_dimension + 1) - DomainType::n_facets_per_element * (fem_order - 1)   // > 0 \iff R > 2
    };
    
    // constructor
    R_FunctionalSpace(Rcpp::Environment mesh, int space_type) : space_type_(space_type) {
        // set domain
        SEXP meshptr = mesh[".pointer"];
        R_Mesh<M, N>* ptr = reinterpret_cast<R_Mesh<M, N>*>(R_ExternalPtrAddr(meshptr));
	domain_ = ptr->domain();
	
	enumerate_dofs(domain_);
	
	// define functional space
	std::size_t size = n_dofs_;
	switch(space_type_) {
	case space_type::fem_lagrange: {
	  fun_space_ = LagrangianBasis<DomainType, R>(domain_, size);
	} break;
        }
    }
    std::size_t size() const { return fun_space_.size(); }   // number of basis functions on phyiscal domain
    // given a coefficient vector c \in \mathbb{R}^size_, evaluates the corresponding basis expansion at locs
    DVector<double> operator()(const DVector<double>& c, const DMatrix<double>& locs) const {
        return fun_space_(c, locs);
    }
    // returns the \Psi is the matrix of basis functions evaluations according
    // to the given policy (see the specific policy for details)
    SpMatrix<double> eval(int evaluation_type, const DMatrix<double>& locs) const {
        return fun_space_.eval(evaluation_type, locs);
    }
    // given a coefficient vector c \in \mathbb{R}^size_, evaluates the integral of the corresponding basis expasion
    // over the whole domain
    double integrate(const DVector<double>& c) const {
      return integrator_.integrate(domain_, fun_space_(c, integrator_.quadrature_nodes(domain_)));
    }
    
    // destructor
    ~R_FunctionalSpace() = default;
};

template<int M, int N, int R>
void R_FunctionalSpace<M,N,R>::enumerate_dofs(const Mesh<M, N>& mesh) {
    if (n_dofs_ != 0) return;   // return early if dofs already computed
    if constexpr (fem_order == 1) {
      n_dofs_ = mesh.n_nodes();
      dofs_ = mesh.elements();
      boundary_dofs_ = mesh.boundary();
    } else {
      dofs_.resize(mesh.n_elements(), n_dof_per_element);
      dofs_.leftCols(DomainType::n_vertices) = mesh.elements();   // copy dofs associated to geometric vertices

      int next = mesh.n_nodes();   // next valid ID to assign
      auto edge_pattern = combinations<DomainType::n_vertices_per_edge, DomainType::n_vertices>();
      std::set<int> boundary_set;

      // cycle over mesh edges
      for (auto edge = mesh.facet_begin(); edge != mesh.facet_end(); ++edge) {
            for (std::size_t i = 0; i < DomainType::n_elements_per_facet; ++i) {
                int element_id = (*edge).adjacent_elements()[i];
                if (element_id >= 0) {
                    // search for dof insertion point
                    std::size_t j = 0;
                    for (; j < edge_pattern.rows(); ++j) {
                        std::array<int, DomainType::n_vertices_per_edge> e {};
                        for (std::size_t k = 0; k < DomainType::n_vertices_per_edge; ++k) {
                            e[k] = mesh.elements()(element_id, edge_pattern(j, k));
                        }
                        std::sort(e.begin(), e.end());   // normalize edge ordering
                        if ((*edge).node_ids() == e) break;
                    }
                    dofs_(element_id, DomainType::n_vertices + j) = next;
                    if ((*edge).on_boundary()) boundary_set.insert(next);

                    // insert any internal dofs, if any (for cubic or higher order) + insert n_dof_per_edge dofs (for
                    // cubic or higher)
                }
            }
            next++;
      }

      n_dofs_ = next;   // store number of unknowns
      // update boundary
      boundary_dofs_ = DMatrix<int>::Zero(n_dofs_, 1);
      boundary_dofs_.topRows(mesh.boundary().rows()) = mesh.boundary();
      for (auto it = boundary_set.begin(); it != boundary_set.end(); ++it) { boundary_dofs_(*it, 0) = 1; }
    }
    return;
}

#endif   // __R_FUNCTIONAL_SPACE_H__
