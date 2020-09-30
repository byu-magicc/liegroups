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

so2::so2(const double data) {
    data_ = Eigen::Matrix<double,1,1>(data);
}

//---------------------------------------------------------------------

so2::so2(const Eigen::Matrix<double,2,2>& data, bool verify) {

    if(verify)
    {
        if (data.transpose() + data == Eigen::Matrix2d::Zero()  )
            data_ = data;
        else {
            std::cerr << "so2::Constructor - Input data not valid. Setting to identity element" << std::endl;
            data_ = Eigen::Matrix<double,1,1>::Zero();
        }
    }
    else {
        data_ = data;
    }

}

//---------------------------------------------------------------------

so2* so2::Bracket(const AlgebraBase * u) {

    return new so2;
    
}

// //---------------------------------------------------------------------

// so2::Eigen::MatrixXd Adjoint();

// //---------------------------------------------------------------------

// so2::Eigen::MatrixXd Wedge();

// //---------------------------------------------------------------------

// so2::Eigen::MatrixXd Vee();

// //---------------------------------------------------------------------

// so2::Eigen::MatrixXd Exp();

// //---------------------------------------------------------------------

// so2::double Norm();

// //---------------------------------------------------------------------

// so2::Eigen::MatrixXd Jl();

// //---------------------------------------------------------------------

// so2::AlgebraBase Jl(const AlgebraBase * u);

// //---------------------------------------------------------------------

// so2::Eigen::MatrixXd JlInv();

// //---------------------------------------------------------------------

// so2::AlgebraBase JlInv(const AlgebraBase * u);

// //---------------------------------------------------------------------

// so2::Eigen::MatrixXd Jr();

// //---------------------------------------------------------------------

// so2::AlgebraBase Jr(const AlgebraBase * u);

// //---------------------------------------------------------------------
 
// so2::Eigen::MatrixXd JrInv();

// //---------------------------------------------------------------------
 
// so2::AlgebraBase JrInv(const AlgebraBase * u);

// //---------------------------------------------------------------------

// so2::AlgebraBase operator + (const AlgebraBase * u);

// //---------------------------------------------------------------------
 
// so2::AlgebraBase operator - (const AlgebraBase * u);

// //---------------------------------------------------------------------

// so2::AlgebraBase operator * (const double scalar);

// //---------------------------------------------------------------------
 
// so2::void Print();

// //---------------------------------------------------------------------

// so2::AlgebraBase Identity();







    
}