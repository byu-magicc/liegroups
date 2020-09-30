#include "lie_algebras/so2.h"


namespace lie_groups {

so2::so2() {

    data_ = Eigen::Matrix<double,1,1>::Zero();

}

//---------------------------------------------------------------------

so2::so2(const so2 & u) {
    data_ = u.data_;
}

//---------------------------------------------------------------------

so2::so2(const Eigen::Matrix<double,1,1> data) {
    data_ = data;
}

//---------------------------------------------------------------------

so2::so2(const Eigen::Matrix<double,2,2>& data, bool verify) {

    if(verify)
    {
        if (data.transpose() + data == Eigen::Matrix2d::Zero()  )
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

so2 so2::Bracket(const so2& u) {
    return so2();
}

//---------------------------------------------------------------------

Eigen::Matrix2d so2::Adjoint() {
    return Eigen::Matrix2d::Identity();
}

//---------------------------------------------------------------------

Eigen::Matrix2d so2::Wedge() {
    Eigen::Matrix2d m;
    m << 0, -data_(0), data_(0), 0;
    return m;
}

//---------------------------------------------------------------------

Eigen::Matrix<double,1,1> so2::Vee() {
    return data_;
}

Eigen::Matrix2d so2::Exp() {
    Eigen::Matrix2d m;
    m << cos(data_(0)), -sin(data_(0)), sin(data_(0)), cos(data_(0));
    return m;
}

//---------------------------------------------------------------------

double so2::Norm() {
    return data_.norm();
}

//---------------------------------------------------------------------

Eigen::Matrix2d so2::Jl() {
    return Eigen::Matrix2d::Identity();
}

//---------------------------------------------------------------------

so2 so2::Jl(const so2& u) {
    return u;
}

//---------------------------------------------------------------------

Eigen::Matrix2d so2::JlInv() {
    return Eigen::Matrix2d::Identity();
}

//---------------------------------------------------------------------

so2 so2::JlInv(const so2& u) {
    return u;
}

//---------------------------------------------------------------------

Eigen::Matrix2d so2::Jr() {
    return Eigen::Matrix2d::Identity();
}

//---------------------------------------------------------------------

so2 so2::Jr(const so2& u) {
    return u;
}

//---------------------------------------------------------------------
 
Eigen::Matrix2d so2::JrInv() {
    return Eigen::Matrix2d::Identity();
}

//---------------------------------------------------------------------

so2 so2::JrInv(const so2& u) {
    return u;
}

//---------------------------------------------------------------------

so2 so2::operator + (const so2& u) {
    return so2(data_ + u.data_);
}

//---------------------------------------------------------------------

so2 so2::operator - (const so2& u) {
    return so2(data_ - u.data_);
}

//---------------------------------------------------------------------

void so2::operator = (const so2& u) {
    data_ = u.data_;
}

//---------------------------------------------------------------------

so2 so2::operator * (const double scalar) {
    return so2(scalar*data_);
}

//---------------------------------------------------------------------

void so2::Print() {
    std::cout << "so2: " << data_(0) << std::endl;
}

//---------------------------------------------------------------------

so2  so2::Identity() {
    return so2();
}
    
}