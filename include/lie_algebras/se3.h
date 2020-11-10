#ifndef _LIEGROUPS_INCLUDE_LIEALGEBRAS_SE3_
#define _LIEGROUPS_INCLUDE_LIEALGEBRAS_SE3_

#include <Eigen/Dense>
#include <iostream>
#include "lie_algebras/so3.h"

namespace lie_groups {


constexpr double kse3_threshold_=1e-7; /** < If two values are within this threshold, they are considered equal.*/


class se3  {

public:

Eigen::Matrix<double,6,1> data_; /** < The vector is translational velocity followed by angular velocity*/
Eigen::Map<Eigen::Matrix<double,3,1>> p_;  /** < The translational velocity */
Eigen::Map<Eigen::Matrix<double,3,1>> th_; /** < The angular velocity */
static constexpr unsigned int dim_ = 6;
static constexpr unsigned int size1_ = 6;
static constexpr unsigned int size2_ = 1;


/**
 * Default constructor. Initializes algebra element to identity.
 */
se3() : p_(data_.data()), th_(data_.data()+3), data_(Eigen::Matrix<double,6,1>::Zero()){}


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
se3(const Eigen::Matrix<double,6,1> data) : p_(data_.data()), th_(data_.data()+3), data_(data){}

/**
* Initializes algebra element to the one given. If verify is set to true,
* it will verify that the element provided is an element of the Lie algebra.
* @param[in] u The data of an element of \f$se(3)\f$
* @param verify If true, the constructor will verify that the element given is an element of the Lie algebra.
*/
se3(const Eigen::Matrix<double,4,4> & data, bool verify);

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
Eigen::Matrix<double,6,6> Adjoint() {    
    Eigen::Matrix<double,6,6> m = Eigen::Matrix<double,6,6>::Zero();
    m.block(0,0,3,3) = se3::SSM(th_);
    m.block(3,3,3,3) = se3::SSM(th_);
    m.block(0,3,3,3) = se3::SSM(p_);

    return m;}

/**
 * Computes the Wedge operation which maps an element of the Cartesian space to the Lie algebra.
 * @return The result of the Wedge operation.
 */
Eigen::Matrix4d Wedge(){return Wedge(data_);}

/**
 * Computes the Wedge operation which maps an element of the Cartesian space to the Lie algebra.
 * @param data The data of an element  of the Cartesian space isomorphic to the Lie algebra
 * @return The result of the Wedge operation.
 */
static Eigen::Matrix4d Wedge(const Eigen::Matrix<double,6,1>& data) {   
    Eigen::Matrix4d m = Eigen::Matrix4d::Zero();
    m.block(0,0,3,3) = se3::SSM(data.block(3,0,3,1));
    m.block(0,3,3,1) = data.block(0,0,3,1);
    return m;}

/**
 * Computes the Vee operation which maps an element of the Lie algebra to the Cartesian space.
 * @return The result of the Vee operation.
 */
Eigen::Matrix<double,6,1> Vee(){return data_;}

/**
 * Computes the Vee operation which maps an element of the Lie algebra to the Cartesian space.
 * @param data The data of an element
 * @return The result of the Vee operation.
 */
static Eigen::Matrix<double,6,1> Vee(const Eigen::Matrix4d& data){
    Eigen::Matrix<double,6,1> m;
    m.block(0,0,3,1) = data.block(0,3,3,1);
    m.block(3,0,3,1) << data(2,1), data(0,2), data(1,0);
    return m;
}

/**
 * Computes the exponential of the element of the Lie algebra.
 */
Eigen::Matrix<double,4,4> Exp(){return se3::Exp(this->data_);}

/**
 * Computes the exponential of the element of the Lie algebra.
 * @return The data associated to the group element.
 */
static Eigen::Matrix<double,4,4> Exp(const Eigen::Matrix<double,6,1>& data);

/**
 * Computes the logaritm of the element of the Lie algebra.
 * @param data The data associated with an element of \f$ SE(3) \f$
 * @return The data of an element of the Cartesian space associated with the Lie algebra
 */
static Eigen::Matrix<double,6,1> Log(const Eigen::Matrix<double,4,4>& data);

/**
 * Computes and returns the Euclidean norm of the element of the Lie algebra
 */ 
double Norm(){return data_.norm();}

/**
 * Computes and returns the matrix of the Left Jacobian.
 */ 
Eigen::Matrix<double,6,6> Jl();

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
Eigen::Matrix<double,6,6> JlInv();

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
Eigen::Matrix<double,6,6> Jr();

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
Eigen::Matrix<double,6,6> JrInv();

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
se3 operator * (const double scalar){return se3(scalar*data_);}

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
static Eigen::Matrix3d SSM(const Eigen::Matrix<double,3,1>& x){
    Eigen::Matrix3d m;
    m << 0, -x(2), x(1), x(2), 0, -x(0), -x(1), x(0), 0;
    return m;
}

/**
 * Verifies that the parameter data belongs to 
 * an element of \f$se(3)\f$
 */ 
static bool isElement(const Eigen::Matrix<double,4,4>& data);


private:

// The following are used to compute the Jacobians
static Eigen::Matrix3d Bl(const Eigen::Matrix<double, 6,1>& u);
static Eigen::Matrix3d Br(const Eigen::Matrix<double, 6,1>& u);


};

}

#endif //_LIEGROUPS_INCLUDE_LIEALGEBRAS_SE3_