#include "lie_algebras/se3.h"


namespace lie_groups {

se3::se3() : p_(data_.data()), th_(data_.data()+3) {

    data_ = Eigen::Matrix<double,6,1>::Zero();

}

//---------------------------------------------------------------------

se3::se3(const se3 & u) : p_(data_.data()), th_(data_.data()+3) {
    data_ = u.data_;
}

//---------------------------------------------------------------------

se3::se3(const Eigen::Matrix<double,6,1> data) : p_(data_.data()), th_(data_.data()+3) {
    data_ = data;
}

//---------------------------------------------------------------------

se3::se3(const Eigen::Matrix<double,6,6>& data, bool verify) : p_(data_.data()), th_(data_.data()+3) {

    if(verify)
    {
        
        if (se3::isElement(data)) {
            data_ = se3::Vee(data);
        }
        else {
            std::cerr << "se3::Constructor - Input data not valid. Setting to identity element" << std::endl;
            data_ = Eigen::Matrix<double,3,1>::Zero();
        }
    }
    else {
        data_ = se3::Vee(data)
    }

}

//---------------------------------------------------------------------

se3 se3::Bracket(const se3& u) {   

    return se3(this->Adjoint()*u.data_);;
}

//---------------------------------------------------------------------

Eigen::Matrix<double,6,6> se3::Adjoint() {
    Eigen::Matrix<double,6,6> m = Eigen::Matrix<double,6,6>::Zero();
    m.block(0,0,3,3) = se3::SSM(th_);
    m.block(3,3,3,3) = se3::SSM(th_);
    m.block(0,3,3,3) = se3::SSM(p_);

    return m;
}

//---------------------------------------------------------------------

Eigen::Matrix4d se3::Wedge() {
    return Wedge(data_);
}

//---------------------------------------------------------------------

Eigen::Matrix4d se3::Wedge(const Eigen::Matrix<double,6,1>& data) {
    Eigen::Matrix4d m = Eigen::Matrix4d::Zero();
    m.block(0,0,3,3) = se3::SSM(data.block(3,0,3,1));
    m.block(0,3,3,1) = data.block(0,0,3,1);
    return m;
}

//---------------------------------------------------------------------

Eigen::Matrix<double,6,1> se3::Vee() {
    return data_;
}

//---------------------------------------------------------------------

Eigen::Matrix<double,6,1> se3::Vee(const Eigen::Matrix4d& data) {
    Eigen::Matrix<double,6,1> m;
    m.block(0,0,3,1) = data.block(0,3,3,1);
    m.block(3,0,3,1) << data(2,1), data(0,2), data(1,0);
    return m;
}

//---------------------------------------------------------------------

Eigen::Matrix3d se3::Exp() {
    
    Eigen::Matrix3d m;
    m.block(0,0,2,2) << cos(th_(0)), - sin(th_(0)), sin(th_(0)), cos(th_(0));
    m.block(0,2,2,1) = this->Wl()*p_;
    m.block(2,0,1,2).setZero();
    m(2,2) = 1;
    return m;
}

//---------------------------------------------------------------------

Eigen::Matrix<double,3,1> se3::Log(const Eigen::Matrix3d data) {
    Eigen::Matrix<double,3,1> u;
    u(2) = atan2(data(1,0),data(0,0)); // Compute the angle
    Eigen::Matrix2d wl = Wl(u(2));
    u.block(0,0,2,1) = wl.inverse()*data.block(0,2,2,1);
    return u;    
}

//---------------------------------------------------------------------

double se3::Norm() {
    return data_.norm();
}

//---------------------------------------------------------------------

Eigen::Matrix3d se3::Jl() {

    Eigen::Matrix3d m;
    m.setIdentity();
    m.block(0,0,2,2) = this->Wl();
    m.block(0,2,2,1) = this->Dl()*p_;

    return m;
}

//---------------------------------------------------------------------

se3 se3::Jl(const se3& u) {
    return se3(this->Jl()*u.data_);
}

//---------------------------------------------------------------------

Eigen::Matrix3d se3::JlInv() {

    Eigen::Matrix2d w_inv = (this->Wl()).inverse();
    Eigen::Matrix3d m;
    m.setIdentity();
    m.block(0,0,2,2) = w_inv;
    m.block(0,2,2,1) = -w_inv*this->Dl()*p_;

    return m;
}

//---------------------------------------------------------------------

se3 se3::JlInv(const se3& u) {
    return se3(this->JlInv()*u.data_);
}

//---------------------------------------------------------------------

Eigen::Matrix3d se3::Jr() {

    Eigen::Matrix3d m;
    m.setIdentity();
    m.block(0,0,2,2) = this->Wr();
    m.block(0,2,2,1) = this->Dr()*p_;

    return m;
}

//---------------------------------------------------------------------

se3 se3::Jr(const se3& u) {
    return se3(this->Jr()*u.data_);
}

//---------------------------------------------------------------------
 
Eigen::Matrix3d se3::JrInv() {
    
    Eigen::Matrix2d w_inv = (this->Wr()).inverse();
    Eigen::Matrix3d m;
    m.setIdentity();
    m.block(0,0,2,2) = w_inv;
    m.block(0,2,2,1) = -w_inv*this->Dr()*p_;

    return m;
}

//---------------------------------------------------------------------

se3 se3::JrInv(const se3& u) {
    return se3(this->JrInv()*u.data_);
}

//---------------------------------------------------------------------

se3 se3::operator + (const se3& u) {
    return se3(data_ + u.data_);
}

//---------------------------------------------------------------------

se3 se3::operator - (const se3& u) {
    return se3(data_ - u.data_);
}

//---------------------------------------------------------------------

void se3::operator = (const se3& u) {
    data_ = u.data_;
}

//---------------------------------------------------------------------

se3 se3::operator * (const double scalar) {
    return se3(scalar*data_);
}

//---------------------------------------------------------------------

void se3::Print() {
    std::cout << "se3: " << std::endl << data_ << std::endl;
}

//---------------------------------------------------------------------

se3  se3::Identity() {
    return se3();
}

//---------------------------------------------------------------------

Eigen::Matrix3d se3::SSM(const Eigen::Matrix<double,3,1>& x) {

    Eigen::Matrix3d m;
    m << 0, -x(2), x(1), x(2), 0, -x(0), -x(1), x(0), 0;
    return m;
}

//---------------------------------------------------------------------

bool se3::isElement(const Eigen::Matrix3d& data) {

    bool is_element = true;
     
    if ( (data.block(0,0,2,2).transpose() + data.block(0,0,2,2)).norm() >= kse3_threshold_) {
        is_element = false;
    }
    else if (data.block(2,0,1,3) != Eigen::Matrix<double,1,3>::Zero()) {
        is_element = false;
    }

    return is_element;
}


//--------------------------------------------------------
//                  Private functions
//--------------------------------------------------------
Eigen::Matrix2d se3::Wl(){return Wl(th_(0));}
Eigen::Matrix2d se3::Wr(){return Wr(th_(0));}
Eigen::Matrix2d se3::Dl(){return Dl(th_(0));}
Eigen::Matrix2d se3::Dr(){return Dr(th_(0));}


Eigen::Matrix2d se3::Wl(const double th) {

    Eigen::Matrix2d m;

    if (fabs(th) > kse3_threshold_) {
        double a = (1-cos(th))/th;
        double b = sin(th)/th;
        m = a*se3::SSM(1) + b*Eigen::Matrix2d::Identity();
    }
    else
    {
        m.setIdentity();
    }
    
return m;
}

//---------------------------------------------------
Eigen::Matrix2d se3::Wr(const double th) {

    Eigen::Matrix2d m;

    if (fabs(th) > kse3_threshold_) {
        double a = (cos(th)-1)/th;
        double b = sin(th)/th;
        m = a*se3::SSM(1) + b*Eigen::Matrix2d::Identity();
    }
    else
    {
        m.setIdentity();
    }
    
return m;
}

//---------------------------------------------------
Eigen::Matrix2d se3::Dl(const double th) {

    Eigen::Matrix2d m;

    if (fabs(th) > kse3_threshold_) {
        double a = (cos(th)-1)/(th*th);
        double b = (th-sin(th))/(th*th);
        m = a*se3::SSM(1) + b*Eigen::Matrix2d::Identity();
    }
    else
    {
        m.setIdentity();
    }

    return m;
}

//---------------------------------------------------
Eigen::Matrix2d se3::Dr(const double th) {

    Eigen::Matrix2d m;

    if (fabs(th) > kse3_threshold_) {
        double a = (1-cos(th))/(th*th);
        double b = (th-sin(th))/(th*th);
        m = a*se3::SSM(1) + b*Eigen::Matrix2d::Identity();
    }
    else
    {
        m.setIdentity();
    }

    return m;
}

//---------------------------------------------------
    
}