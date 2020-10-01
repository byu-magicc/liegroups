#ifndef _LIEGROUPS_INCLUDE_LIEALGEBRAS_SE3_
#define _LIEGROUPS_INCLUDE_LIEALGEBRAS_SE3_

#include <Eigen/Dense>
#include <iostream>

namespace lie_groups {


constexpr double kse3_threshold_=1e-7; /** < If two values are within this threshold, they are considered equal.*/


class se3  {

public:

Eigen::Matrix<double,6,1> data_; /** < The vector is translational velocity followed by angular velocity*/
Eigen::Map<Eigen::Matrix<double,3,1>> p_;  /** < The translational velocity */
Eigen::Map<Eigen::Matrix<double,3,1>> th_; /** < The angular velocity */


/**
 * Default constructor. Initializes algebra element to identity.
 */
se3();


/**
 * Copy constructor.
 */ 
se3(const se3 & u);

/**
* Initializes algebra element to the one given. 
* @param[in] data The data of an element of Cartesian space of \f$se(3)\f$
*/
se3(const Eigen::Matrix<double,6,1> data);

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
se3 Bracket(const se3& u);

/**
 * Computes and returns the matrix adjoint representation of the Lie algebra.
 * This is always the Identity element.
 */ 
Eigen::Matrix<double,6,6> Adjoint();

/**
 * Computes the Wedge operation which maps an element of the Cartesian space to the Lie algebra.
 * @return The result of the Wedge operation.
 */
Eigen::Matrix4d Wedge();

/**
 * Computes the Wedge operation which maps an element of the Cartesian space to the Lie algebra.
 * @param data The data of an element  of the Cartesian space isomorphic to the Lie algebra
 * @return The result of the Wedge operation.
 */
static Eigen::Matrix4d Wedge(const Eigen::Matrix<double,6,1>& data);

/**
 * Computes the Vee operation which maps an element of the Lie algebra to the Cartesian space.
 * @return The result of the Vee operation.
 */
Eigen::Matrix<double,6,1> Vee();

/**
 * Computes the Vee operation which maps an element of the Lie algebra to the Cartesian space.
 * @param data The data of an element
 * @return The result of the Vee operation.
 */
static Eigen::Matrix<double,6,1> Vee(const Eigen::Matrix4d& data);

/**
 * Computes the exponential of the element of the Lie algebra.
 */
Eigen::Matrix<double,4,4> Exp();

/**
 * Computes the logaritm of the element of the Lie algebra.
 * @param data The data associated with an element of \f$ SE(3) \f$
 * @return The data of an element of the Cartesian space associated with the Lie algebra
 */
static Eigen::Matrix<double,6,1> Log(const Eigen::Matrix<double,6,6> data);

/**
 * Computes and returns the Euclidean norm of the element of the Lie algebra
 */ 
double Norm();

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
se3 Jl(const se3& u);

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
se3 JlInv(const se3& u);

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
se3 Jr(const se3& u);

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
se3 JrInv(const se3& u);

/**
 * Adds two elements of the Algebra together
 * @param u An element of the Lie algebra.
 */ 
se3 operator + (const se3& u);

/**
 * Subtracts two elements of the Algebra together
 * @param u An element of the Lie algebra.
 */ 
se3 operator - (const se3& u);

/**
 * Creates a deep copy of the element
 * @param u An element of the Lie algebra.
 */
void operator = (const se3& u);

/**
 * Performs Scalar multiplication and returns the result.
 * @param scalar The scalar that will scale the element of the Lie algebra
 */ 
se3 operator * (const double scalar);

/**
 * Prints the data of the element.
 */ 
void Print();

/**
 * Returns the Identity element.
 */
static se3 Identity();

/**
 * Computes and returns the skew symmetric matrix of the
 * parameter.
 * @param x 
 */ 
static Eigen::Matrix3d SSM(const Eigen::Matrix<double,3,1>& x);

/**
 * Verifies that the parameter data belongs to 
 * an element of \f$se(3)\f$
 */ 
static bool isElement(const Eigen::Matrix<double,6,6>& data);


private:

// The following are used to compute the Jacobians
static Eigen::Matrix3d Wl(const double th);
static Eigen::Matrix3d Wr(const double th);
static Eigen::Matrix3d Dl(const double th);
static Eigen::Matrix3d Dr(const double th);
Eigen::Matrix3d Wl();
Eigen::Matrix3d Wr();
Eigen::Matrix3d Dl();
Eigen::Matrix3d Dr();


};

}

#endif //_LIEGROUPS_INCLUDE_LIEALGEBRAS_SE3_