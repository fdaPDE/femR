#include <fdaPDE/mesh.h>
#include <fdaPDE/utils/symbols.h>

#include <memory>
using fdapde::core::ADT;
using fdapde::core::Mesh;

// M : local dimension, N : embedding dimension
template <unsigned int M, unsigned int N> class R_Mesh {
   private:
    typedef Mesh<M, N, 1> DomainType;
    typedef ADT<M, N, 1> PointLocatorType;
    // internal data (use shared_ptr since Rcpp::wrap will make a copy when passed as argument)
    DomainType domain_;
    PointLocatorType point_locator_;
   public:
    // constructor
    R_Mesh(const Rcpp::List& mesh_data) :
        domain_(
          Rcpp::as<DMatrix<double>>(mesh_data["nodes"]), Rcpp::as<DMatrix<int>>(mesh_data["elements"]),
          Rcpp::as<DMatrix<int>>(mesh_data["neigh"]), Rcpp::as<DMatrix<int>>(mesh_data["boundary"])),
        point_locator_(domain_) {};

    // getters
    const DomainType& domain() const { return domain_; }
    const PointLocatorType& point_locator() const { return point_locator_; }
    const DMatrix<double>& nodes() const { return domain_.nodes(); }
    const Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>& elements() const {
        return domain_.elements();
    }
    const DMatrix<int>& neighbors() const { return domain_.neighbors(); }
    const DMatrix<int>& boundary() const { return domain_.boundary(); }

    // destructor
    ~R_Mesh() = default;
};

typedef R_Mesh<2, 2> Mesh_2D;
