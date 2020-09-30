#include "lie_groups/SO2.h"

namespace lie_groups {


//----------------------------------------------------------
SO2::SO2(){
 g_ = Eigen::Matrix<double,2,2>::Identity();
}

//----------------------------------------------------------


SO2::SO2(const SO2& G) {
    g_ = G.g_;
}

//----------------------------------------------------------
SO2::SO2(const Eigen::Matrix<double,2,2> & g) {

    // First verify that it is a proper group element.
    if (g.determinant() == 1 && (g.transpose()*g-g.Identity()).norm() == 0) {
        g_ = g;
    }
    else {
        std::cerr << "SO2::Constructor not valide input setting to identity" << std::endl;
        g_ = Eigen::Matrix<double,2,2>::Identity();
    }

}






}