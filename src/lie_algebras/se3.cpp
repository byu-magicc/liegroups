#include "lie_algebras/se3.h"


namespace lie_groups {




//---------------------------------------------------------------------

se3::se3(const Eigen::Matrix<double,4,4>& data, bool verify) : p_(data_.data()), th_(data_.data()+3) {

    if(verify)
    {
        
        if (se3::isElement(data)) {
            data_ = se3::Vee(data);
        }
        else {
            std::cerr << "se3::Constructor - Input data not valid. Setting to identity element" << std::endl;
            data_ = Eigen::Matrix<double,6,1>::Zero();
        }
    }
    else {
        data_ = se3::Vee(data);
    }

}



//---------------------------------------------------------------------

Eigen::Matrix<double,4,4> se3::Exp(const Eigen::Matrix<double,6,1>& data) {
    Eigen::Matrix<double,4,4> m;
    so3 omega(data.block(3,0,3,1));
    m.block(0,0,3,3) = omega.Exp();
    m.block(0,3,3,1) = omega.Jl()*data.block(0,0,3,1);
    m.block(3,0,1,4) << 0,0,0,1;
    return m;  
}


//---------------------------------------------------------------------

Eigen::Matrix<double,6,1> se3::Log(const Eigen::Matrix<double,4,4>& data) {
    
    Eigen::Matrix<double,6,1> u;
    so3 omega(so3::Log(data.block(0,0,3,3)));
    u.block(3,0,3,1) = omega.Vee();
    u.block(0,0,3,1) = omega.JlInv()*data.block(0,3,3,1);
    
    return u;    
}

//---------------------------------------------------------------------

Eigen::Matrix<double,6,6> se3::Jl() {

    so3 omega(th_);

    Eigen::Matrix<double,6,6> m;
    m.block(0,0,3,3) = omega.Jl();
    m.block(0,3,3,3) = Bl(data_);
    m.block(3,0,3,3).setZero();
    m.block(3,3,3,3) =  m.block(0,0,3,3);

    return m;
}

//---------------------------------------------------------------------

Eigen::Matrix<double,6,6> se3::JlInv() {

    so3 omega(th_);

    Eigen::Matrix<double,6,6> m;
    m.block(0,0,3,3) = omega.JlInv();
    m.block(0,3,3,3) = -m.block(0,0,3,3)*Bl(data_)*m.block(0,0,3,3);
    m.block(3,0,3,3).setZero();
    m.block(3,3,3,3) =  m.block(0,0,3,3);

    return m;
}

//---------------------------------------------------------------------

Eigen::Matrix<double,6,6> se3::Jr() {

    so3 omega(th_);

    Eigen::Matrix<double,6,6> m;
    m.block(0,0,3,3) = omega.Jr();
    m.block(0,3,3,3) = Br(data_);
    m.block(3,0,3,3).setZero();
    m.block(3,3,3,3) =  m.block(0,0,3,3);

    return m;
}


//---------------------------------------------------------------------
 
Eigen::Matrix<double,6,6> se3::JrInv() {
    
    so3 omega(th_);

    Eigen::Matrix<double,6,6> m;
    m.block(0,0,3,3) = omega.JrInv();
    m.block(0,3,3,3) = -m.block(0,0,3,3)*Br(data_)*m.block(0,0,3,3);
    m.block(3,0,3,3).setZero();
    m.block(3,3,3,3) =  m.block(0,0,3,3);

    return m;

    return m;
}

//---------------------------------------------------------------------

bool se3::isElement(const Eigen::Matrix<double,4,4>& data) {

    bool is_element = true;
     
    if ( !so3::isElement(data.block(0,0,3,3))) {
        is_element = false;
    }
    else if (data.block(3,0,1,4) != Eigen::Matrix<double,1,4>::Zero()) {
        is_element = false;
    }

    return is_element;
}


/////////////////////////////////////////////////
//                  Private Functions
/////////////////////////////////////////////////

Eigen::Matrix3d se3::Bl(const Eigen::Matrix<double, 6,1>& u) {

Eigen::Matrix<double,3,3> m;
Eigen::Map<const Eigen::Matrix<double,3,1>> p(u.data());
Eigen::Map<const Eigen::Matrix<double,3,1>> w(u.data()+3);

double th = w.norm();


if (th <= kse3_threshold_) { // Close to the identity element;
    m.setIdentity();
} else {
    double th2 = pow(th,2);
    double th3 = pow(th,3);
    double th4 = pow(th,4);
    double a = (cos(th)-1)/th2;
    double b = (th - sin(th))/th3;
    double c = -sin(th)/th3 + 2*(1-cos(th))/th4;
    double d = -2/th4 + 3*sin(th)/pow(th,5) - cos(th)/th4;
    Eigen::Matrix3d q;
    q = w.dot(p)*(-c*SSM(w) + d*SSM(w)*SSM(w));

    m = -a*SSM(p) + b*(SSM(w)*SSM(p) + SSM(p)*SSM(w)) + q;
}

return m;

}

//---------------------------------------------------------------------

Eigen::Matrix3d se3::Br(const Eigen::Matrix<double, 6,1>& u) {

Eigen::Matrix<double,3,3> m;
Eigen::Map<const Eigen::Matrix<double,3,1>> p(u.data());
Eigen::Map<const Eigen::Matrix<double,3,1>> w(u.data()+3);

double th = w.norm();


if (th <= kse3_threshold_) { // Close to the identity element;
    m.setIdentity();
} else {
    double th2 = pow(th,2);
    double th3 = pow(th,3);
    double th4 = pow(th,4);
    double a = (cos(th)-1)/th2;
    double b = (th - sin(th))/th3;
    double c = -sin(th)/th3 + 2*(1-cos(th))/th4;
    double d = -2/th4 + 3*sin(th)/pow(th,5) - cos(th)/th4;
    Eigen::Matrix3d q;
    q = w.dot(p)*(c*SSM(w) + d*SSM(w)*SSM(w));
    m = a*SSM(p) + b*(SSM(w)*SSM(p) + SSM(p)*SSM(w)) + q;
}

return m;


}


}