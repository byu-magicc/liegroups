#ifndef _LIEGROUPS_INCLUDE_LIEALGEBRAS_SE3_
#define _LIEGROUPS_INCLUDE_LIEALGEBRAS_SE3_

#include <Eigen/Dense>
#include <iostream>
#include "lie_algebras/so3.h"

namespace lie_groups {


constexpr double kse3_threshold_=1e-7; /** < If two values are within this threshold, they are considered equal.*/

template< typename tDataType = double>
class se3  {

public:

Eigen::Matrix<tDataType,6,1> data_; /** < The vector is translational velocity followed by angular velocity*/
Eigen::Map<Eigen::Matrix<tDataType,3,1>> p_;  /** < The translational velocity */
Eigen::Map<Eigen::Matrix<tDataType,3,1>> th_; /** < The angular velocity */
static constexpr unsigned int dim_ = 6;
static constexpr unsigned int dim_t_vel_=3; /** < The dimension of the translational velocity */
static constexpr unsigned int dim_a_vel_=3; /** < The dimension of the angular velocity */
static constexpr unsigned int size1_ = 6;
static constexpr unsigned int size2_ = 1;
typedef Eigen::Matrix<tDataType,3,1> Vec3d;
typedef Eigen::Matrix<tDataType,6,1> Vec6d;
typedef Eigen::Matrix<tDataType,3,3> Mat3d;
typedef Eigen::Matrix<tDataType,4,4> Mat4d;
typedef Eigen::Matrix<tDataType,6,6> Mat6d;

/**
 * Default constructor. Initializes algebra element to identity.
 */
se3() : p_(data_.data()), th_(data_.data()+3), data_(Vec6d::Zero()){}


/**
 * Copy constructor.
 */ 
se3(const se3 & u) : p_(data_.data()), th_(data_.data()+3), data_(u.data_) {}

/**
 *Copy assignment
 */
void operator = (const se3& u){data_ = u.data_;}

/**
 * Move constructor.
 */ 
se3(const se3 && u) : p_(data_.data()), th_(data_.data()+3), data_(u.data_) {}

/**
 *Move assignment
 */
void operator = (const se3&& u){data_ = u.data_;}


/**
* Initializes algebra element to the one given. 
* @param[in] data The data of an element of Cartesian space of \f$se(3)\f$
*/
se3(const Vec6d data) : p_(data_.data()), th_(data_.data()+3), data_(data){}

/**
* Initializes algebra element to the one given. If verify is set to true,
* it will verify that the element provided is an element of the Lie algebra.
* @param[in] u The data of an element of \f$se(3)\f$
* @param verify If true, the constructor will verify that the element given is an element of the Lie algebra.
*/
se3(const Mat4d & data, bool verify);

//---------------------------------------------------------------------

/**
 * Performs the Lie bracket.. 
 * The input is the right parameter 
 * of the Lie bracket function. \f$ \[v,u\] \f$
 * @param u An element of the Lie algebra.
 * @return The result of the Lie bracket operation.
 */ 
se3 Bracket(const se3& u){return se3(this->Adjoint()*u.data_);}

/**
 * Computes and returns the matrix adjoint representation of the Lie algebra.
 * This is always the Identity element.
 */ 
Mat6d Adjoint() {    
    Mat6d m = Mat6d::Zero();
    m.block(0,0,3,3) = se3<tDataType>::SSM(th_);
    m.block(3,3,3,3) = se3<tDataType>::SSM(th_);
    m.block(0,3,3,3) = se3<tDataType>::SSM(p_);

    return m;}

/**
 * Computes the Wedge operation which maps an element of the Cartesian space to the Lie algebra.
 * @return The result of the Wedge operation.
 */
Mat4d Wedge(){return Wedge(data_);}

/**
 * Computes the Wedge operation which maps an element of the Cartesian space to the Lie algebra.
 * @param data The data of an element  of the Cartesian space isomorphic to the Lie algebra
 * @return The result of the Wedge operation.
 */
static Mat4d Wedge(const Vec6d& data) {   
    Mat4d m = Mat4d::Zero();
    m.block(0,0,3,3) = se3<tDataType>::SSM(data.block(3,0,3,1));
    m.block(0,3,3,1) = data.block(0,0,3,1);
    return m;}

/**
 * Computes the Vee operation which maps an element of the Lie algebra to the Cartesian space.
 * @return The result of the Vee operation.
 */
Vec6d Vee(){return data_;}

/**
 * Computes the Vee operation which maps an element of the Lie algebra to the Cartesian space.
 * @param data The data of an element
 * @return The result of the Vee operation.
 */
static Vec6d Vee(const Mat4d& data){
    Vec6d m;
    m.block(0,0,3,1) = data.block(0,3,3,1);
    m.block(3,0,3,1) << data(2,1), data(0,2), data(1,0);
    return m;
}

/**
 * Computes the exponential of the element of the Lie algebra.
 */
Mat4d Exp(){return se3<tDataType>::Exp(this->data_);}

/**
 * Computes the exponential of the element of the Lie algebra.
 * @return The data associated to the group element.
 */
static Mat4d Exp(const Vec6d& data);

/**
 * Computes the logaritm of the element of the Lie algebra.
 * @param data The data associated with an element of \f$ SE(3) \f$
 * @return The data of an element of the Cartesian space associated with the Lie algebra
 */
static Vec6d Log(const Mat4d& data);

/**
 * Computes and returns the Euclidean norm of the element of the Lie algebra
 */ 
tDataType Norm(){return data_.norm();}

/**
 * Computes and returns the matrix of the Left Jacobian.
 */ 
Mat6d Jl();

/**
 * Computes the left Jacobian using the element of *this and applies it to the
 * parameter provided. \f$ J_l(v)u \f$.
 * Since the left Jacobian is the identity map for \f$se(3)\f$, it will simply return
 * the parameter u.
 * @param u An element of the Lie algebra.
 * @return It will return the parameter u.
 */ 
se3 Jl(const se3& u){return se3(this->Jl()*u.data_);}

/**
 * Computes and returns the matrix of the Left Jacobian inverse.
 */ 
Mat6d JlInv();

/**
 * Computes the left Jacobian inverse using the element of *this and applies it to the
 * parameter provided. \f$ J_l(v)u \f$
 * Since the left Jacobian inverse is the identity map for \f$se(3)\f$, it will simply return
 * the parameter u.
 * @param u An element of the Lie algebra.
 */ 
se3 JlInv(const se3& u){return se3(this->JlInv()*u.data_);}

/**
 * Computes and returns the matrix of the Right Jacobian
 */ 
Mat6d Jr();

/**
 * Computes the right Jacobian using the element of *this and applies it to the
 * parameter provided. \f$ J_l(v)u \f$
 * Since the right Jacobian is the identity map for \f$se(3)\f$, it will simply return
 * the parameter u.
 * @param u An element of the Lie algebra.
 */ 
se3 Jr(const se3& u){return se3(this->Jr()*u.data_);}

/**
 * Computes and returns the matrix of the right Jacobian inverse.
 */ 
Mat6d JrInv();

/**
 * Computes the right Jacobian inverse using the element of *this and applies it to the
 * parameter provided. \f$ J_l(v)u \f$
 * Since the right Jacobian inverse is the identity map for \f$se(3)\f$, it will simply return
 * the parameter u.
 * @param u An element of the Lie algebra.
 */ 
se3 JrInv(const se3& u){ return se3(this->JrInv()*u.data_);}

/**
 * Adds two elements of the Algebra together
 * @param u An element of the Lie algebra.
 */ 
se3 operator + (const se3& u){return se3(data_ + u.data_);}

/**
 * Subtracts two elements of the Algebra together
 * @param u An element of the Lie algebra.
 */ 
se3 operator - (const se3& u){return se3(data_ - u.data_);}

/**
 * Performs Scalar multiplication and returns the result.
 * @param scalar The scalar that will scale the element of the Lie algebra
 */ 
se3 operator * (const tDataType scalar) const {return se3(scalar*data_);}

/**
 * Prints the data of the element.
 */ 
void Print(){std::cout << data_ << std::endl;}

/**
 * Returns the Identity element.
 */
static se3 Identity(){return se3();}

/**
 * Computes and returns the skew symmetric matrix of the
 * parameter.
 * @param x 
 */ 
static Mat3d SSM(const Vec3d& x){
    Mat3d m;
    m << static_cast<tDataType>(0.0), -x(2), x(1), x(2), static_cast<tDataType>(0.0), -x(0), -x(1), x(0), static_cast<tDataType>(0.0);
    return m;
}

/**
 * Verifies that the parameter data belongs to 
 * an element of \f$se(3)\f$
 */ 
static bool isElement(const Mat4d& data);


private:

// The following are used to compute the Jacobians
static Mat3d Bl(const Vec6d& u);
static Mat3d Br(const Vec6d& u);


};

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                    Definitions
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////


//---------------------------------------------------------------------
template< typename tDataType>
se3<tDataType>::se3(const Mat4d& data, bool verify) : p_(data_.data()), th_(data_.data()+3) {

    if(verify)
    {
        
        if (se3<tDataType>::isElement(data)) {
            data_ = se3<tDataType>::Vee(data);
        }
        else {
            std::cerr << "se3<tDataType>::Constructor - Input data not valid. Setting to identity element" << std::endl;
            data_ = Vec6d::Zero();
        }
    }
    else {
        data_ = se3<tDataType>::Vee(data);
    }

}



//---------------------------------------------------------------------
template< typename tDataType>
Eigen::Matrix<tDataType,4,4> se3<tDataType>::Exp(const Vec6d& data) {
    Mat4d m;
    so3<tDataType> omega(data.block(3,0,3,1));
    m.block(0,0,3,3) = omega.Exp();
    m.block(0,3,3,1) = omega.Jl()*data.block(0,0,3,1);
    m.block(3,0,1,4) << static_cast<tDataType>(0.0),static_cast<tDataType>(0.0),static_cast<tDataType>(0.0),static_cast<tDataType>(1.0);
    return m;  
}


//---------------------------------------------------------------------
template< typename tDataType>
Eigen::Matrix<tDataType,6,1> se3<tDataType>::Log(const Mat4d& data) {
    
    Vec6d u;
    so3<tDataType>  omega(so3<tDataType>::Log(data.block(0,0,3,3)));
    u.block(3,0,3,1) = omega.Vee();
    u.block(0,0,3,1) = omega.JlInv()*data.block(0,3,3,1);
    
    return u;    
}

//---------------------------------------------------------------------
template< typename tDataType>
Eigen::Matrix<tDataType,6,6> se3<tDataType>::Jl() {

    so3<tDataType> omega(th_);

    Eigen::Matrix<tDataType,6,6> m;
    m.block(0,0,3,3) = omega.Jl();
    m.block(0,3,3,3) = Bl(data_);
    m.block(3,0,3,3).setZero();
    m.block(3,3,3,3) =  m.block(0,0,3,3);

    return m;
}

//---------------------------------------------------------------------
template< typename tDataType>
Eigen::Matrix<tDataType,6,6> se3<tDataType>::JlInv() {

    so3<tDataType> omega(th_);

    Eigen::Matrix<tDataType,6,6> m;
    m.block(0,0,3,3) = omega.JlInv();
    m.block(0,3,3,3) = -m.block(0,0,3,3)*Bl(data_)*m.block(0,0,3,3);
    m.block(3,0,3,3).setZero();
    m.block(3,3,3,3) =  m.block(0,0,3,3);

    return m;
}

//---------------------------------------------------------------------
template< typename tDataType>
Eigen::Matrix<tDataType,6,6> se3<tDataType>::Jr() {

    so3<tDataType> omega(th_);

    Eigen::Matrix<tDataType,6,6> m;
    m.block(0,0,3,3) = omega.Jr();
    m.block(0,3,3,3) = Br(data_);
    m.block(3,0,3,3).setZero();
    m.block(3,3,3,3) =  m.block(0,0,3,3);

    return m;
}


//---------------------------------------------------------------------
template< typename tDataType>
Eigen::Matrix<tDataType,6,6> se3<tDataType>::JrInv() {
    
    so3<tDataType> omega(th_);

    Eigen::Matrix<tDataType,6,6> m;
    m.block(0,0,3,3) = omega.JrInv();
    m.block(0,3,3,3) = -m.block(0,0,3,3)*Br(data_)*m.block(0,0,3,3);
    m.block(3,0,3,3).setZero();
    m.block(3,3,3,3) =  m.block(0,0,3,3);

    return m;

    return m;
}

//---------------------------------------------------------------------
template< typename tDataType>
bool se3<tDataType>::isElement(const Mat4d& data) {

    bool is_element = true;
     
    if ( !so3<tDataType>::isElement(data.block(0,0,3,3))) {
        is_element = false;
    }
    else if (data.block(3,0,1,4) != Eigen::Matrix<tDataType,1,4>::Zero()) {
        is_element = false;
    }

    return is_element;
}


/////////////////////////////////////////////////
//                  Private Functions
/////////////////////////////////////////////////
template< typename tDataType>
Eigen::Matrix<tDataType,3,3> se3<tDataType>::Bl(const Eigen::Matrix<tDataType, 6,1>& u) {

Eigen::Matrix<tDataType,3,3> m;
Eigen::Map<const Eigen::Matrix<tDataType,3,1>> p(u.data());
Eigen::Map<const Eigen::Matrix<tDataType,3,1>> w(u.data()+3);

tDataType th = w.norm();


if (th <= kse3_threshold_) { // Close to the identity element;
    m = Eigen::Matrix<tDataType,3,3>::Identity() + SSM(p)/static_cast<tDataType>(2.0);
} else {
    tDataType th2 = pow(th,2);
    tDataType th3 = pow(th,3);
    tDataType th4 = pow(th,4);
    tDataType a = (cos(th)-1.0)/th2;
    tDataType b = (th - sin(th))/th3;
    tDataType c = -sin(th)/th3 + 2.0*(1.0-cos(th))/th4;
    tDataType d = -2.0/th4 + 3.0*sin(th)/pow(th,5) - cos(th)/th4;
    Eigen::Matrix<tDataType,3,3> q;
    q = w.dot(p)*(-c*SSM(w) + d*SSM(w)*SSM(w));

    m = -a*SSM(p) + b*(SSM(w)*SSM(p) + SSM(p)*SSM(w)) + q;
}

return m;

}

//---------------------------------------------------------------------
template< typename tDataType>
Eigen::Matrix<tDataType,3,3> se3<tDataType>::Br(const Eigen::Matrix<tDataType, 6,1>& u) {

Eigen::Matrix<tDataType,3,3> m;
Eigen::Map<const Eigen::Matrix<tDataType,3,1>> p(u.data());
Eigen::Map<const Eigen::Matrix<tDataType,3,1>> w(u.data()+3);

tDataType th = w.norm();


if (th <= kse3_threshold_) { // Close to the identity element;
    m = Eigen::Matrix<tDataType,3,3>::Identity() - SSM(p)/static_cast<tDataType>(2.0);
} else {
    tDataType th2 = pow(th,2);
    tDataType th3 = pow(th,3);
    tDataType th4 = pow(th,4);
    tDataType a = (cos(th)-1.0)/th2;
    tDataType b = (th - sin(th))/th3;
    tDataType c = -sin(th)/th3 + 2.0*(1.0-cos(th))/th4;
    tDataType d = -2.0/th4 + 3.0*sin(th)/pow(th,5) - cos(th)/th4;
    Eigen::Matrix<tDataType,3,3> q;
    q = w.dot(p)*(c*SSM(w) + d*SSM(w)*SSM(w));
    m = a*SSM(p) + b*(SSM(w)*SSM(p) + SSM(p)*SSM(w)) + q;
}

return m;


}


}

#endif //_LIEGROUPS_INCLUDE_LIEALGEBRAS_SE3_