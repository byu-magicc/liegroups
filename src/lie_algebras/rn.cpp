// #include "lie_algebras/rn.h"


// namespace lie_groups {

// rn<N>::rn() {

//     data_.setZero();

// }

// //---------------------------------------------------------------------

// rn<N>::rn(const rn<N> & u) {
//     data_ = u.data_;
// }

// //---------------------------------------------------------------------

// rn::rn(const Eigen::Matrix<double,N,1> data) {
//     data_ = data;
// }

// //---------------------------------------------------------------------

// rn::rn(const Eigen::Matrix<double,N,N>& data, bool verify) {

//     if(verify)
//     {
//         if (isElement(data) )
//             data_(0) = data(1,0);
//         else {
//             std::cerr << "rn::Constructor - Input data not valid. Setting to identity element" << std::endl;
//             data_ = Eigen::Matrix<double,N,1>::Zero();
//         }
//     }
//     else {
//         data_(0) = data(1,0);
//     }

// }

// //---------------------------------------------------------------------

// rn rn::Bracket(const rn& u) {
//     return rn();
// }

// //---------------------------------------------------------------------

// Eigen::Matrix<double,N,N> rn::Adjoint() {
//     return Eigen::Matrix<double,N,N>::Identity();
// }

// //---------------------------------------------------------------------

// Eigen::Matrix<double,N,N> rn::Wedge() {
//     return Wedge(data_);
// }

// //---------------------------------------------------------------------

// Eigen::Matrix<double,N,N> rn::Wedge(const Eigen::Matrix<double,N,1>& data) {
//     Eigen::Matrix<double,N,N> m;
//     m << 0, -data(0), data(0), 0;
//     return m; 
// }


// //---------------------------------------------------------------------

// Eigen::Matrix<double,N,1> rn::Vee() {
//     return data_;
// }

// //---------------------------------------------------------------------

// Eigen::Matrix<double,N,1> rn::Vee(const Eigen::Matrix<double,N,N>& data) {
//     Eigen::Matrix<double,N,1> m;
//     m << data(1,0);
//     return m;
// }

// //---------------------------------------------------------------------

// Eigen::Matrix<double,N,N> rn::Exp() {
//     Eigen::Matrix<double,N,N> m;
//     m << cos(data_(0)), -sin(data_(0)), sin(data_(0)), cos(data_(0));
//     return m;
// }

// //---------------------------------------------------------------------

// Eigen::Matrix<double,N,1> rn::Log(const Eigen::Matrix<double,N,N>& data) {
//     Eigen::Matrix<double,N,1> m;
//     m(0) = atan2(data(1,0),data(0,0));
//     return m;
// }

// //---------------------------------------------------------------------

// double rn::Norm() {
//     return data_.norm();
// }

// //---------------------------------------------------------------------

// Eigen::Matrix<double,N,N> rn::Jl() {
//     return Eigen::Matrix<double,N,N>::Identity();
// }

// //---------------------------------------------------------------------

// rn rn::Jl(const rn& u) {
//     return u;
// }

// //---------------------------------------------------------------------

// Eigen::Matrix<double,N,N> rn::JlInv() {
//     return Eigen::Matrix<double,N,N>::Identity();
// }

// //---------------------------------------------------------------------

// rn rn::JlInv(const rn& u) {
//     return u;
// }

// //---------------------------------------------------------------------

// Eigen::Matrix<double,N,N> rn::Jr() {
//     return Eigen::Matrix<double,N,N>::Identity();
// }

// //---------------------------------------------------------------------

// rn rn::Jr(const rn& u) {
//     return u;
// }

// //---------------------------------------------------------------------
 
// Eigen::Matrix<double,N,N> rn::JrInv() {
//     return Eigen::Matrix<double,N,N>::Identity();
// }

// //---------------------------------------------------------------------

// rn rn::JrInv(const rn& u) {
//     return u;
// }

// //---------------------------------------------------------------------

// rn rn::operator + (const rn& u) {
//     return rn(data_ + u.data_);
// }

// //---------------------------------------------------------------------

// rn rn::operator - (const rn& u) {
//     return rn(data_ - u.data_);
// }

// //---------------------------------------------------------------------

// void rn::operator = (const rn& u) {
//     data_ = u.data_;
// }

// //---------------------------------------------------------------------

// rn rn::operator * (const double scalar) {
//     return rn(scalar*data_);
// }

// //---------------------------------------------------------------------

// void rn::Print() {
//     std::cout << "rn: " << data_(0) << std::endl;
// }

// //---------------------------------------------------------------------

// rn  rn::Identity() {
//     return rn();
// }

// //---------------------------------------------------------------------

// bool rn::isElement(const Eigen::Matrix<double,N,N>& data) {

//     if ( (data.transpose() + data).norm() >= krn_threshold_) {
//         return false;
//     } else {
//         return true;
//     }

// }
    
// }