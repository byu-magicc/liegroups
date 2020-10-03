#include "lie_groups/SE2.h"

namespace lie_groups {


//----------------------------------------------------------

SE2::SE2(const Eigen::Matrix3d & data, bool verify) : t_(data_.data()+6), R_(data_.block(0,0,2,2)) {

    // First verify that it is a proper group element.
    if (verify ) {
        if (SE2::isElement(data)) {
            data_ = data;
        } else {
            std::cerr << "SE2::Constructor not valid input setting to identity" << std::endl;
            data_.setIdentity();
        }
    } else {
        data_ = data;
    }
}

//----------------------------------------------------------

bool SE2::isElement(const Eigen::Matrix3d& data) {
    
    bool d = (data.block(0,0,2,2).transpose()*data.block(0,0,2,2)-Eigen::Matrix2d::Identity()).norm();
    
    return d <= kSE2_threshold_ && data(2,0) == 0 && data(2,1)==0 && data(2,2)==1;
}




}