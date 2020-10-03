#include "lie_groups/SE3.h"

namespace lie_groups {


//-------------------------------------------------------------------

SE3::SE3(const Eigen::Matrix4d & data, bool verify) : t_(data_.data()+12), R_(data_.block(0,0,3,3)) {

    // First verify that it is a proper group element.
    if (verify ) {
        if (SE3::isElement(data)) {
            data_ = data;
        } else {
            std::cerr << "SE3::Constructor not valid input setting to identity" << std::endl;
            data_.setIdentity();
        }
    } else {
        data_ = data;
    }
}

//-------------------------------------------------------------------

bool SE3::isElement(const Eigen::Matrix4d& data) {
    
    double d = (data.block(0,0,3,3).transpose()*data.block(0,0,3,3)-Eigen::Matrix3d::Identity()).norm();
    
    return d <= kSE3_threshold_ && data(3,0) == 0 && data(3,1)==0 && data(3,2)==0 && data(3,3)==1;
}




}