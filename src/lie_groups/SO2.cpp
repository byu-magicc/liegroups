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

bool SO2::isElement(const Eigen::Matrix2d& data) {
    return (data.transpose()*data - Eigen::Matrix2d::Identity()).norm() < kSO2_threshold_;
}




}