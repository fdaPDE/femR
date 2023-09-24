#include <fdaPDE/utils/symbols.h>
#include <fdaPDE/mesh.h>
#include<fdaPDE/finite_elements.h>
#include <memory>
using fdapde::core::ADT;
using fdapde::core::Mesh;
using fdapde::core::FiniteElementBasis;
using fdapde::core::LagrangianElement;

// M : local dimension, N : embedding dimension
template <unsigned int M, unsigned int N, unsigned int R> class R_Mesh {
   private:
    typedef Mesh<M, N, R> DomainType;
    typedef ADT<M, N, R> PointLocatorType;
    
    typedef LagrangianElement<DomainType::local_dimension, DomainType::order> FunctionSpace;
    typedef FiniteElementBasis<FunctionSpace> FunctionBasis;
    
    // internal data (use shared_ptr since Rcpp::wrap will make a copy when passed as argument)
    DomainType domain_ {};
    PointLocatorType point_locator_;
    
    FunctionBasis basis_ {};  
    
   public:
    // constructor
    R_Mesh(const Rcpp::List& mesh_data) :
        domain_(
          Rcpp::as<DMatrix<double>>(mesh_data["nodes"]), Rcpp::as<DMatrix<int>>(mesh_data["elements"]),
          Rcpp::as<DMatrix<int>>(mesh_data["neigh"]), Rcpp::as<DMatrix<int>>(mesh_data["boundary"])),
        point_locator_(domain_) {};
    
    static constexpr int fe_order_ = DomainType::order;
    
    // getters
    const DomainType& domain() const { return domain_; }
    const PointLocatorType& point_locator() const { return point_locator_; }
    const DMatrix<double>& nodes() const { return domain_.nodes(); }
    const Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>& elements() const {
        return domain_.elements();
    }
    const DMatrix<int>& neighbors() const { return domain_.neighbors(); }
    const DMatrix<int>& boundary() const { return domain_.boundary(); }
    int fe_order()const {return fe_order_; }
    
    // eval solution at point
    DMatrix<double> eval(const DMatrix<double>& coeff, const DMatrix<double>& point) {
        /*
        SEXP meshptr = mesh[".pointer"];
        R_Mesh<M, N, R>* ptr = (R_Mesh<M, N, R>*)R_ExternalPtrAddr(meshptr);
        if (point.cols() != M) {
            throw std::runtime_error("evaluation point must be a " + std::to_string(M) + "-dimensional vector");
        }
        */
        DMatrix<double> evals;
        evals.resize(point.rows(), 1);
        for (std::size_t j = 0; j < point.rows(); ++j) {
            SVector<N> p;
            for (int i = 0; i < N; ++i) p[i] = point(j, i);
            // locate element containing point
            auto e = point_locator_.locate(p);
            // compute value of function at point
            double v = std::numeric_limits<double>::quiet_NaN();
            if (e != nullptr) {
                v = 0;
                for (std::size_t k = 0; k < DomainType::n_dof_per_element; ++k) {
                    v += coeff(domain_.elements()(e->ID(), k), 0) * basis_(*e, k)(p);
                }
            }
            evals(j, 0) = v;
        }
        return evals;
    }
    // destructor
    ~R_Mesh() = default;
};

typedef R_Mesh<2, 2, 1> Mesh_2D_ORDER_1;
typedef R_Mesh<2, 2, 2> Mesh_2D_ORDER_2;
