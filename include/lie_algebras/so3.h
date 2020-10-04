#ifndef _LIEGROUPS_INCLUDE_LIEALGEBRAS_SO3_
#define _LIEGROUPS_INCLUDE_LIEALGEBRAS_SO3_

#include <Eigen/Dense>
#include <iostream>

namespace lie_groups {

constexpr double kso3_threshold_=1e-7; /** < If two values are within this threshold, they are considered equal.*/

class so3  {

public:

Eigen::Matrix<double,3,1> data_; /** < The vector is translational velocity followed by angular velocity*/
static const unsigned int dim_ = 3;
static const unsigned int size1_ = 3;
static const unsigned int size2_ = 1;


/**
 * Default constructor. Initializes algebra element to identity.
 */
so3() : data_(Eigen::Matrix<double,3,1>::Zero()){}

/**
 * Copy constructor.
 */ 
so3(const so3 & u) : data_(u.data_) {}

/**
* Initializes algebra element to the one given. 
* @param[in] data The data of an element of Cartesian space of \f$so(3)\f$
*/
so3(const Eigen::Matrix<double,3,1> data) : data_(data) {}

/**
* Initializes algebra element to the one given. If verify is set to true,
* it will verify that the element provided is an element of the Lie algebra.
* @param[in] u The data of an element of \f$so(3)\f$
* @param verify If true, the constructor will verify that the element given is an element of the Lie algebra.
*/
so3(const Eigen::Matrix3d & data, bool verify);

//---------------------------------------------------------------------

/**
 * Performs the Lie bracket.. 
 * The input is the right parameter 
 * of the Lie bracket function. \f$ \[v,u\] \f$
 * @param u An element of the Lie algebra.
 * @return The result of the Lie bracket operation.
 */ 
so3 Bracket(const so3& u) {return so3(this->Adjoint()*u.data_);}

/**
 * Computes and returns the matrix adjoint representation of the Lie algebra.
 * This is always the Identity element.
 */ 
Eigen::Matrix3d Adjoint() {return this->Wedge();}

/**
 * Computes the Wedge operation which maps an element of the Cartesian space to the Lie algebra.
 * @return The result of the Wedge operation.
 */
Eigen::Matrix3d Wedge() {return Wedge(data_);}

/**
 * Computes the Wedge operation which maps an element of the Cartesian space to the Lie algebra.
 * @param data The data of an element  of the Cartesian space isomorphic to the Lie algebra
 * @return The result of the Wedge operation.
 */
static Eigen::Matrix3d Wedge(const Eigen::Matrix<double,3,1>& data) {    Eigen::Matrix3d m;
    m << 0, -data(2), data(1), data(2), 0, -data(0), -data(1), data(0), 0;
    return m;
}

/**
 * Computes the Vee operation which maps an element of the Lie algebra to the Cartesian space.
 * @return The result of the Vee operation.
 */
Eigen::Matrix<double,3,1> Vee() {return data_;}

/**
 * Computes the Vee operation which maps an element of the Lie algebra to the Cartesian space.
 * @param data The data of an element
 * @return The result of the Vee operation.
 */
static Eigen::Matrix<double,3,1> Vee(const Eigen::Matrix3d& data) {    
    Eigen::Matrix<double,3,1> m;
    m << data(2,1), data(0,2), data(1,0);
    return m;}

/**
 * Computes the exponential of the element of the Lie algebra.
 */
Eigen::Matrix3d Exp(){ return so3::Exp(this->data_);}

/**
 * Computes the exponential of the element of the Lie algebra.
 * @return The data associated to the group element.
 */
static Eigen::Matrix3d Exp(const Eigen::Matrix<double,3,1>& data);

/**
 * Computes the logaritm of the element of the Lie algebra.
 * @param data The data associated with an element of \f$ SO(3) \f$
 * @return The data of an element of the Cartesian space associated with the Lie algebra
 */
static Eigen::Matrix<double,3,1> Log(const Eigen::Matrix3d& data);

/**
 * Computes and returns the Euclidean norm of the element of the Lie algebra
 */ 
double Norm(){return data_.norm();}

/**
 * Computes and returns the matrix of the Left Jacobian.
 */ 
Eigen::Matrix3d Jl();

/**
 * Computes the left Jacobian using the element of *this and applies it to the
 * parameter provided. \f$ J_l(v)u \f$.
 * Since the left Jacobian is the identity map for \f$so(3)\f$, it will simply return
 * the parameter u.
 * @param u An element of the Lie algebra.
 * @return It will return the parameter u.
 */ 
so3 Jl(const so3& u){return so3(this->Jl()*u.data_);}

/**
 * Computes and returns the matrix of the Left Jacobian inverse.
 */ 
Eigen::Matrix3d JlInv();

/**
 * Computes the left Jacobian inverse using the element of *this and applies it to the
 * parameter provided. \f$ J_l(v)u \f$
 * Since the left Jacobian inverse is the identity map for \f$so(3)\f$, it will simply return
 * the parameter u.
 * @param u An element of the Lie algebra.
 */ 
so3 JlInv(const so3& u){return so3(this->JlInv()*u.data_);}

/**
 * Computes and returns the matrix of the Right Jacobian
 */ 
Eigen::Matrix3d Jr();

/**
 * Computes the right Jacobian using the element of *this and applies it to the
 * parameter provided. \f$ J_l(v)u \f$
 * Since the right Jacobian is the identity map for \f$so(3)\f$, it will simply return
 * the parameter u.
 * @param u An element of the Lie algebra.
 */ 
so3 Jr(const so3& u){return so3(this->Jr()*u.data_);}

/**
 * Computes and returns the matrix of the right Jacobian inverse.
 */ 
Eigen::Matrix3d JrInv();

/**
 * Computes the right Jacobian inverse using the element of *this and applies it to the
 * parameter provided. \f$ J_l(v)u \f$
 * Since the right Jacobian inverse is the identity map for \f$so(3)\f$, it will simply return
 * the parameter u.
 * @param u An element of the Lie algebra.
 */ 
so3 JrInv(const so3& u){return so3(this->JrInv()*u.data_);}

/**
 * Adds two elements of the Algebra together
 * @param u An element of the Lie algebra.
 */ 
so3 operator + (const so3& u){return so3(data_ + u.data_);}

/**
 * Subtracts two elements of the Algebra together
 * @param u An element of the Lie algebra.
 */ 
so3 operator - (const so3& u){return so3(data_ - u.data_);}

/**
 * Creates a deep copy of the element
 * @param u An element of the Lie algebra.
 */
void operator = (const so3& u){data_ = u.data_;}

/**
 * Performs Scalar multiplication and returns the result.
 * @param scalar The scalar that will scale the element of the Lie algebra
 */ 
so3 operator * (const double scalar){return so3(scalar*data_);}

/**
 * Prints the data of the element.
 */ 
void Print(){std::cout << data_ << std::endl;}

/**
 * Returns the Identity element.
 */
static so3 Identity(){return so3();}

/**
 * Computes and returns the skew symmetric matrix of the
 * parameter.
 * @param x 
 */ 
static Eigen::Matrix2d SSM(double x);

/**
 * Verifies that the parameter data belongs to 
 * an element of \f$so(3)\f$
 */ 
static bool isElement(const Eigen::Matrix3d& data);


};

}

#endif //_LIEGROUPS_INCLUDE_LIEALGEBRAS_SO3_