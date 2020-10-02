#include "lie_groups/SO2.h"

namespace lie_groups {


//----------------------------------------------------------

SO2::SO2(const Eigen::Matrix2d & data, bool verify) {

    // First verify that it is a proper group element.
    if (verify ) {
        if (SO2::isElement(data)) {
            data_ = data;
        } else {
            std::cerr << "SO2::Constructor not valid input setting to identity" << std::endl;
            data_.setIdentity();
        }
    } else {
        data_ = data;
    }
}

//----------------------------------------------------------

Eigen::Matrix<double,1,1> SO2::Log() {
    Eigen::Matrix<double,1,1> m;
    m << atan2(data_(1,0),data_(0,0));
    return m;
}

//----------------------------------------------------------

bool SO2::isElement(const Eigen::Matrix2d& data) {
    return (data.transpose()*data - Eigen::Matrix2d::Identity()).norm() < kSO2_threshold_;
}




}