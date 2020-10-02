#include "lie_algebras/so3.h"


namespace lie_groups {



//---------------------------------------------------------------------

so3::so3(const Eigen::Matrix<double,3,3>& data, bool verify) {

    if(verify)
    {
        
        if (so3::isElement(data)) {
            data_(0) = data(2,1);
            data_(1) = data(0,2);
            data_(2) = data(1,0);
        }
        else {
            std::cerr << "so3::Constructor - Input data not valid. Setting to identity element" << std::endl;
            data_ = Eigen::Matrix<double,3,1>::Zero();
        }
    }
    else {
        data_(0) = data(2,1);
        data_(1) = data(0,2);
        data_(2) = data(1,0);
    }

}

Eigen::Matrix3d so3::Exp(const Eigen::Matrix<double,3,1>& data) {
    Eigen::Matrix3d m;
    double th = data.norm();

    if (th < kso3_threshold_) { // See if the element is close to the identity element.
        m.setIdentity();
    } else {  // Use Rodriguez formula 
        double a = sin(th)/th;
        double b = (1.0-cos(th))/pow(th,2);
        m = Eigen::Matrix3d::Identity() + a*Wedge(data) + b*Wedge(data)*Wedge(data);
    }
    return m;

}

//---------------------------------------------------------------------

Eigen::Matrix<double,3,1> so3::Log(const Eigen::Matrix3d& data) {

    Eigen::Matrix<double,3,1> u;

    double t = data.trace();
    if ( fabs(t-3.0) <= kso3_threshold_) { // Rotation matrix is close to identity
    
        u.setZero();

    } else { // Use Rodriguez formula 

        double th = acos( (t-1)/2);
        u = so3::Vee(th*(data-data.transpose())/(2*sin(th)));
    }

    return u;    
}

//---------------------------------------------------------------------

Eigen::Matrix3d so3::Jl() {

    Eigen::Matrix3d m;
    
    double th = data_.norm();

    if (th < kso3_threshold_) { // See if the element is close to the identity element.
        m.setIdentity();
    } else {   
        double a = (1.0-cos(th))/pow(th,2);
        double b = (th-sin(th))/pow(th,3);
        m = Eigen::Matrix3d::Identity() + a*this->Wedge() + b*this->Wedge()*this->Wedge();
    }

    return m;
}


//---------------------------------------------------------------------

Eigen::Matrix3d so3::JlInv() {


    Eigen::Matrix3d m;

    double th = data_.norm();

    if (th < kso3_threshold_ || sin(th/2) < kso3_threshold_) { // See if the element is close to the identity element.
        m.setIdentity();
    } else {   
        double a = -0.5;
        double cot = cos(th/2.0)/sin(th/2.0);
        double b = -(th*cot-2.0)/(2.0*pow(th,2));
        m = Eigen::Matrix3d::Identity() + a*this->Wedge() + b*this->Wedge()*this->Wedge();
    }

    return m;
}

//---------------------------------------------------------------------

Eigen::Matrix3d so3::Jr() {

    Eigen::Matrix3d m;
    
    double th = data_.norm();

    if (th < kso3_threshold_) { // See if the element is close to the identity element.
        m.setIdentity();
    } else {   
        double a = (cos(th)-1.0)/pow(th,2);
        double b = (th-sin(th))/pow(th,3);
        m = Eigen::Matrix3d::Identity() + a*this->Wedge() + b*this->Wedge()*this->Wedge();
    }

    return m;
}


//---------------------------------------------------------------------
 
Eigen::Matrix3d so3::JrInv() {
    
    Eigen::Matrix3d m;

    double th = data_.norm();

    if (th < kso3_threshold_ || sin(th/2) < kso3_threshold_) { // See if the element is close to the identity element.
        m.setIdentity();
    } else {   
        double a = 0.5;
        double cot = cos(th/2.0)/sin(th/2.0);
        double b = -(th*cot-2.0)/(2.0*pow(th,2));
        m = Eigen::Matrix3d::Identity() + a*this->Wedge() + b*this->Wedge()*this->Wedge();
    }

    return m;
}

//---------------------------------------------------------------------

bool so3::isElement(const Eigen::Matrix3d& data) {

    bool is_element = true;
     
    if ( (data.transpose()+data).norm() >= kso3_threshold_) {
        is_element = false;
    }

    return is_element;
}


//--------------------------------------------------------
//                  Private functions
//--------------------------------------------------------


    
}