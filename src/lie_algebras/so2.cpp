#include "lie_algebras/so2.h"


namespace lie_groups {


//---------------------------------------------------------------------

so2::so2(const Eigen::Matrix<double,2,2>& data, bool verify) {

    if(verify)
    {
        if (isElement(data) )
            data_(0) = data(1,0);
        else {
            std::cerr << "so2::Constructor - Input data not valid. Setting to identity element" << std::endl;
            data_ = Eigen::Matrix<double,1,1>::Zero();
        }
    }
    else {
        data_(0) = data(1,0);
    }

}

//---------------------------------------------------------------------

Eigen::Matrix2d so2::Exp(const Eigen::Matrix<double,1,1> &data) {
    Eigen::Matrix2d m;
    m << cos(data(0)), -sin(data(0)), sin(data(0)), cos(data(0));
    return m;
}

//---------------------------------------------------------------------

Eigen::Matrix<double,1,1> so2::Log(const Eigen::Matrix2d& data) {
    Eigen::Matrix<double,1,1> m;
    m(0) = atan2(data(1,0),data(0,0));
    return m;
}

//---------------------------------------------------------------------

bool so2::isElement(const Eigen::Matrix2d& data) {

    if ( (data.transpose() + data).norm() >= kso2_threshold_) {
        return false;
    } else {
        return true;
    }

}
    
}