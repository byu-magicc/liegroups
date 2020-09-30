#include "lie_algebras/se2.h"


namespace lie_groups {

se2::se2() : p_(data_.data()), th_(data_.data()+2) {

    data_ = Eigen::Matrix<double,3,1>::Zero();

}

//---------------------------------------------------------------------

se2::se2(const se2 & u) : p_(data_.data()), th_(data_.data()+2) {
    data_ = u.data_;
}

//---------------------------------------------------------------------

se2::se2(const Eigen::Matrix<double,3,1> data) : p_(data_.data()), th_(data_.data()+2) {
    data_ = data;
}

//---------------------------------------------------------------------

se2::se2(const Eigen::Matrix<double,3,3>& data, bool verify) : p_(data_.data()), th_(data_.data()+2) {

    if(verify)
    {
        
        if (se2::isElement(data)) {
            data_(0) = data(0,2);
            data_(1) = data(1,2);
            data_(2) = data(1,0);
        }
        else {
            std::cerr << "se2::Constructor - Input data not valid. Setting to identity element" << std::endl;
            data_ = Eigen::Matrix<double,3,1>::Zero();
        }
    }
    else {
        data_(0) = data(0,2);
        data_(1) = data(1,2);
        data_(2) = data(1,0);
    }

}

//---------------------------------------------------------------------

se2 se2::Bracket(const se2& u) {   

    return se2(this->Adjoint()*u.data_);;
}

//---------------------------------------------------------------------

Eigen::Matrix3d se2::Adjoint() {
    Eigen::Matrix3d m = Eigen::Matrix3d::Zero();
    m.block(0,0,2,2) = se2::SSM(th_(0));
    m.block(0,2,2,1) = -se2::SSM(1)*p_;

    return m;
}

//---------------------------------------------------------------------

Eigen::Matrix3d se2::Wedge() {
    Eigen::Matrix3d m = Eigen::Matrix3d::Zero();
    m.block(0,0,2,2) = se2::SSM(th_(0));
    m.block(0,2,2,1) = p_;
    return m;
}

//---------------------------------------------------------------------

Eigen::Matrix<double,3,1> se2::Vee() {
    return data_;
}

Eigen::Matrix3d se2::Exp() {
    
    Eigen::Matrix3d m;
    m.block(0,0,2,2) << cos(th_(0)), - sin(th_(0)), sin(th_(0)), cos(th_(0));
    m.block(0,2,2,1) = this->Wl()*p_;
    m.block(2,0,1,2).setZero();
    m(2,2) = 1;
    return m;
}

//---------------------------------------------------------------------

double se2::Norm() {
    return data_.norm();
}

//---------------------------------------------------------------------

Eigen::Matrix3d se2::Jl() {

    Eigen::Matrix3d m;
    m.setIdentity();
    m.block(0,0,2,2) = this->Wl();
    m.block(0,2,2,1) = this->Dl()*p_;

    return m;
}

//---------------------------------------------------------------------

se2 se2::Jl(const se2& u) {
    return se2(this->Jl()*u.data_);
}

//---------------------------------------------------------------------

Eigen::Matrix3d se2::JlInv() {

    Eigen::Matrix2d w_inv = (this->Wl()).inverse();
    Eigen::Matrix3d m;
    m.setIdentity();
    m.block(0,0,2,2) = w_inv;
    m.block(0,2,2,1) = -w_inv*this->Dl()*p_;

    return m;
}

//---------------------------------------------------------------------

se2 se2::JlInv(const se2& u) {
    return se2(this->JlInv()*u.data_);
}

//---------------------------------------------------------------------

Eigen::Matrix3d se2::Jr() {

    Eigen::Matrix3d m;
    m.setIdentity();
    m.block(0,0,2,2) = this->Wr();
    m.block(0,2,2,1) = this->Dr()*p_;

    return m;
}

//---------------------------------------------------------------------

se2 se2::Jr(const se2& u) {
    return se2(this->Jr()*u.data_);
}

//---------------------------------------------------------------------
 
Eigen::Matrix3d se2::JrInv() {
    
    Eigen::Matrix2d w_inv = (this->Wr()).inverse();
    Eigen::Matrix3d m;
    m.setIdentity();
    m.block(0,0,2,2) = w_inv;
    m.block(0,2,2,1) = -w_inv*this->Dr()*p_;

    return m;
}

//---------------------------------------------------------------------

se2 se2::JrInv(const se2& u) {
    return se2(this->JrInv()*u.data_);
}

//---------------------------------------------------------------------

se2 se2::operator + (const se2& u) {
    return se2(data_ + u.data_);
}

//---------------------------------------------------------------------

se2 se2::operator - (const se2& u) {
    return se2(data_ - u.data_);
}

//---------------------------------------------------------------------

void se2::operator = (const se2& u) {
    data_ = u.data_;
}

//---------------------------------------------------------------------

se2 se2::operator * (const double scalar) {
    return se2(scalar*data_);
}

//---------------------------------------------------------------------

void se2::Print() {
    std::cout << "se2: " << std::endl << data_ << std::endl;
}

//---------------------------------------------------------------------

se2  se2::Identity() {
    return se2();
}

//---------------------------------------------------------------------

Eigen::Matrix2d se2::SSM(double x) {
    Eigen::Matrix2d m;
    m << 0, -x, x, 0;
    return m;
}

//---------------------------------------------------------------------

bool se2::isElement(const Eigen::Matrix3d& data) {

    bool is_element = true;
     
    if (data.block(0,0,2,2).transpose() + data.block(0,0,2,2) != Eigen::Matrix2d::Zero()) {
        is_element = false;
    }
    else if (data.block(2,0,1,3) != Eigen::Matrix<double,1,3>::Zero()) {
        is_element = false;
    }

    return is_element;
}


Eigen::Matrix2d se2::getWl() {return this->Wl();}
Eigen::Matrix2d se2::getWr() {return this->Wr();}
Eigen::Matrix2d se2::getDl() {return this->Dl();}
Eigen::Matrix2d se2::getDr() {return this->Dr();}

//--------------------------------------------------------
//                  Private functions
//--------------------------------------------------------

Eigen::Matrix2d se2::Wl() {

    Eigen::Matrix2d m;

    if (th_.norm() > 1e-6) {
        double th = th_(0);
        double a = (1-cos(th))/th;
        double b = sin(th)/th;
        m = a*se2::SSM(1) + b*Eigen::Matrix2d::Identity();
    }
    else
    {
        m.setIdentity();
    }
    
return m;
}

//---------------------------------------------------
Eigen::Matrix2d se2::Wr() {

    Eigen::Matrix2d m;

    if (th_.norm() > 1e-6) {
        double th = th_(0);
        double a = (cos(th)-1)/th;
        double b = sin(th)/th;
        m = a*se2::SSM(1) + b*Eigen::Matrix2d::Identity();
    }
    else
    {
        m.setIdentity();
    }
    
return m;
}

//---------------------------------------------------
Eigen::Matrix2d se2::Dl() {

    Eigen::Matrix2d m;

    if (th_.norm() > 1e-6) {
        double th = th_(0);
        double a = (cos(th)-1)/(th*th);
        double b = (th-sin(th))/(th*th);
        m = a*se2::SSM(1) + b*Eigen::Matrix2d::Identity();
    }
    else
    {
        m.setIdentity();
    }

    return m;
}

//---------------------------------------------------
Eigen::Matrix2d se2::Dr() {

    Eigen::Matrix2d m;

    if (th_.norm() > 1e-6) {
        double th = th_(0);
        double a = (1-cos(th))/(th*th);
        double b = (th-sin(th))/(th*th);
        m = a*se2::SSM(1) + b*Eigen::Matrix2d::Identity();
    }
    else
    {
        m.setIdentity();
    }

    return m;
}

//---------------------------------------------------
    
}