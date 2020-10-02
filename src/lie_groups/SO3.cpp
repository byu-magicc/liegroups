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

Eigen::Matrix<double,3,1> SO3::Log() {
    Eigen::Matrix<double,3,1> u;

    double t = data_.trace();
    if ( fabs(t-3.0) <= kso3_threshold_) { // Rotation matrix is close to identity
    
        u.setZero();

    } else { // Use Rodriguez formula 

        double th = acos( (t-1)/2);
        u = so3::Vee(th*(data_-data_.transpose())/(2*sin(th)));
    }

    return u;
}

//----------------------------------------------------------

bool SO3::isElement(const Eigen::Matrix3d& data) {
    return (data.transpose()*data - Eigen::Matrix3d::Identity()).norm() < kSO3_threshold_;
}




}