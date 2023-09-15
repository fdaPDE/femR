// we only include RcppEigen.h which pulls Rcpp.h in for us
#include <RcppEigen.h>
// [[Rcpp::depends(RcppEigen)]]

#include <fdaPDE/utils/symbols.h>
#include <fdaPDE/pde.h>
#include <fdaPDE/mesh.h>
using fdapde::core::PDEBase;
using fdapde::core::PDE;
using fdapde::core::laplacian;
using fdapde::core::advection;
using fdapde::core::diffusion;
using fdapde::core::reaction;
using fdapde::core::FEM;
using fdapde::core::Mesh;
using fdapde::core::Integrator;
using fdapde::core::ADT;
using fdapde::core::LagrangianElement;
using fdapde::core::FiniteElementBasis;

// M : local dimension, N : embedding dimension, R : FEM order
template <unsigned int M, unsigned int N, unsigned int R>
class R_PDE {
private:
  typedef Mesh<M,N,R> DomainType;
  typedef Integrator<DomainType::local_dimension, DomainType::order> QuadratureRule;
  
  // internal data
  Mesh<M,N,R> domain_;
  ADT<M,N,R> point_locator_;
  QuadratureRule integrator_ {};
  PDEBase* pde_ = nullptr;
  int pde_type_;
  bool init_ = false; // asserted true if PDE is init
  Rcpp::Nullable<Rcpp::List> pde_parameters_; // eventual PDE coefficients
  DMatrix<double> u_; // discretized forcing term
   
  Rcpp::Environment function_;

  typedef LagrangianElement<DomainType::local_dimension, DomainType::order> FunctionSpace;
  typedef FiniteElementBasis<FunctionSpace> FunctionBasis;
  FunctionBasis basis_ {};
public:
  // constructor
  R_PDE(const Rcpp::List& mesh_data,
	int pde_type,
	const Rcpp::Nullable<Rcpp::List>& pde_parameters,
	const Rcpp::Environment& function) :
    domain_(Rcpp::as<DMatrix<double>>(mesh_data["nodes"]),
	    Rcpp::as<DMatrix<int>>   (mesh_data["elements"]),
	    Rcpp::as<DMatrix<int>>   (mesh_data["neigh"]),
	    Rcpp::as<DMatrix<int>>   (mesh_data["boundary"])),
    pde_type_(pde_type), pde_parameters_(pde_parameters), function_(function),
    point_locator_(domain_) { }
  
  // setters
  void set_dirichlet_bc(const DMatrix<double>& data) { pde_->set_dirichlet_bc(data); }
  void set_forcing(const DMatrix<double>& data) { u_ = data; }
  // getters
  DMatrix<double> get_quadrature_nodes() const { return integrator_.quadrature_nodes(domain_); };
  DMatrix<double> get_dofs_coordinates() const { return domain_.dof_coords(); };

  SpMatrix<double> R0() const { return pde_->R0(); }
  SpMatrix<double> R1() const { return pde_->R1(); }
  DMatrix<double> solution() const { return pde_->solution(); }
  
  void init() {
    // define at run-time the pde object based on requested type
    switch(pde_type_) {
    case 1:
      { // L = Laplacian
	auto L = -laplacian<FEM>();
	pde_ = new PDE<DomainType,decltype(L),DMatrix<double>,FEM>(domain_, L, u_);
      }
      break;
    case 2:
      { // L = div(K*grad(f)) + b*grad(f) + c*f
	if(pde_parameters_.isNotNull()) {
	  // extract parameters from pack
	  Rcpp::List pde_parameters(pde_parameters_);
	  SMatrix<M> K = Rcpp::as<DMatrix<double>>(pde_parameters["diffusion"]);
	  SVector<M> b = Rcpp::as<DMatrix<double>>(pde_parameters["transport"]);
	  double c = Rcpp::as<double>(pde_parameters["reaction"]);
	  // define bilinear form
	  auto L = -diffusion<FEM>(K) + advection<FEM>(b) + reaction<FEM>(c);
	  pde_ = new PDE<DomainType,decltype(L),DMatrix<double>,FEM>(domain_, L, u_);
	} else {
	  throw std::runtime_error("pde parameters not supplied.");
	}
      }
      break;
    }
    // initialize pointer
    pde_->init();
    return;
  }

  void solve() {
    pde_->solve();
    function_["coeff"] = pde_->solution();
  }

  // eval at point
  DMatrix<double> eval(const DMatrix<double>& coeff, const DMatrix<double>& point) {
    if(point.cols() != M) {
      throw std::runtime_error("evaluation point must be a " + std::to_string(M) + "-dimensional vector");
    }

    DMatrix<double> evals;
    evals.resize(point.rows(), 1);

    for(std::size_t j = 0; j < point.rows(); ++j) { 
      SVector<M> p;
      for(int i = 0; i < M; ++i) p[i] = point(j,i);
      // locate element containing point
      auto e = point_locator_.locate(p);
      // compute value of function at point
      double v = std::numeric_limits<double>::quiet_NaN();
      if(e != nullptr){
	v = 0;
	for(std::size_t j = 0; j < DomainType::n_dof_per_element; ++j) {
	  v += coeff(domain_.elements()(e->ID(), j), 0) * basis_(*e, j)(p);
	}
      }
      evals(j,0) = v;
    }
    return evals;
  }
  
  // destructor
  ~R_PDE() { delete pde_; }
};


// Rcpp modules definition

typedef R_PDE<2,2,1> PDE_2D_ORDER_1;
RCPP_MODULE(PDE_2D_ORDER_1) {
  Rcpp::class_<PDE_2D_ORDER_1>("PDE_2D_ORDER_1")
    .constructor<Rcpp::List, int, Rcpp::Nullable<Rcpp::List>, Rcpp::Environment>()
    // getters
    .method("get_quadrature_nodes", &PDE_2D_ORDER_1::get_quadrature_nodes)
    .method("get_dofs_coordinates", &PDE_2D_ORDER_1::get_dofs_coordinates)
    .method("get_mass",             &PDE_2D_ORDER_1::R0)
    .method("get_stiff",            &PDE_2D_ORDER_1::R1)
    .method("solution",             &PDE_2D_ORDER_1::solution)
    .method("eval",                 &PDE_2D_ORDER_1::eval)
    // setters
    .method("set_dirichlet_bc",     &PDE_2D_ORDER_1::set_dirichlet_bc)
    .method("set_forcing",          &PDE_2D_ORDER_1::set_forcing)
    // init and solve
    .method("init",                 &PDE_2D_ORDER_1::init)
    .method("solve",                &PDE_2D_ORDER_1::solve);
}

typedef R_PDE<2,2,2> PDE_2D_ORDER_2;
RCPP_MODULE(PDE_2D_ORDER_2) {
  Rcpp::class_<PDE_2D_ORDER_2>("PDE_2D_ORDER_2")
    .constructor<Rcpp::List, int, Rcpp::Nullable<Rcpp::List>, Rcpp::Environment>()
    // getters
    .method("get_quadrature_nodes", &PDE_2D_ORDER_2::get_quadrature_nodes)
    .method("get_dofs_coordinates", &PDE_2D_ORDER_2::get_dofs_coordinates)
    .method("get_mass",             &PDE_2D_ORDER_2::R0)
    .method("get_stiff",            &PDE_2D_ORDER_2::R1)
    .method("solution",             &PDE_2D_ORDER_2::solution)
    .method("eval",                 &PDE_2D_ORDER_2::eval)
    // setters
    .method("set_dirichlet_bc",     &PDE_2D_ORDER_2::set_dirichlet_bc)
    .method("set_forcing",          &PDE_2D_ORDER_2::set_forcing)
    // init and solve
    .method("init",                 &PDE_2D_ORDER_2::init)
    .method("solve",                &PDE_2D_ORDER_2::solve);
}
