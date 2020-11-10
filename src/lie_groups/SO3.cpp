#include "lie_groups/SO3.h"

namespace lie_groups {


//----------------------------------------------------------

SO3::SO3(const Eigen::Matrix3d & data, bool verify) {

    // First verify that it is a proper group element.
    if (verify ) {
        if (SO3::isElement(data)) {
            data_ = data;
        } else {
            std::cerr << "SO3::Constructor not valid input setting to identity" << std::endl;
            data_.setIdentity();
        }
    } else {
        data_ = data;
    }
}


//----------------------------------------------------------

bool SO3::isElement(const Eigen::Matrix3d& data) {
    return (data.transpose()*data - Eigen::Matrix3d::Identity()).norm() < kSO3_threshold_;
}




}