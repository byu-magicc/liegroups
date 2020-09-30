#ifndef _LIEGROUPS_INCLUDE_LIEALGEBRAS_SE2_
#define _LIEGROUPS_INCLUDE_LIEALGEBRAS_SE2_

#include <Eigen/Dense>
#include <iostream>

namespace lie_groups {

class se2  {

public:

Eigen::Matrix<double,3,1> data_; /** < The vector is translational velocity followed by angular velocity*/
Eigen::Map<Eigen::Matrix<double,2,1>> p_;  /** < The translational velocity */
Eigen::Map<Eigen::Matrix<double,1,1>> th_; /** < The angular velocity */


/**
 * Default constructor. Initializes algebra element to identity.
 */
se2();


/**
 * Copy constructor.
 */ 
se2(const se2 & u);

/**
* Initializes algebra element to the one given. 
* @param[in] data The data of an element of Cartesian space of \f$se(2)\f$
*/
se2(const Eigen::Matrix<double,3,1> data);

/**
* Initializes algebra element to the one given. If verify is set to true,
* it will verify that the element provided is an element of the Lie algebra.
* @param[in] u The data of an element of \f$se(2)\f$
* @param verify If true, the constructor will verify that the element given is an element of the Lie algebra.
*/
se2(const Eigen::Matrix3d & data, bool verify);

//---------------------------------------------------------------------

/**
 * Performs the Lie bracket.. 
 * The input is the right parameter 
 * of the Lie bracket function. \f$ \[v,u\] \f$
 * @param u An element of the Lie algebra.
 * @return The result of the Lie bracket operation.
 */ 
se2 Bracket(const se2& u);

/**
 * Computes and returns the matrix adjoint representation of the Lie algebra.
 * This is always the Identity element.
 */ 
Eigen::Matrix3d Adjoint();

/**
 * Computes the Wedge operation which maps an element of the Cartesian space to the Lie algebra.
 * @return The result of the Wedge operation.
 */
Eigen::Matrix3d Wedge();

/**
 * Computes the Vee operation which maps an element of the Lie algebra to the Cartesian space.
 * @return The result of the Vee operation.
 */
Eigen::Matrix<double,3,1> Vee();

/**
 * Computes the exponential of the element of the Lie algebra.
 * @return The data associated to the group element.
 */
Eigen::Matrix3d Exp();

/**
 * Computes and returns the Euclidean norm of the element of the Lie algebra
 */ 
double Norm();

/**
 * Computes and returns the matrix of the Left Jacobian.
 * This will always return the identity map for \f$se(2)\f$.
 */ 
Eigen::Matrix3d Jl();

/**
 * Computes the left Jacobian using the element of *this and applies it to the
 * parameter provided. \f$ J_l(v)u \f$.
 * Since the left Jacobian is the identity map for \f$se(2)\f$, it will simply return
 * the parameter u.
 * @param u An element of the Lie algebra.
 * @return It will return the parameter u.
 */ 
se2 Jl(const se2& u);

/**
 * Computes and returns the matrix of the Left Jacobian inverse.
 * This will always return the identity map for \f$se(2)\f$.
 */ 
Eigen::Matrix3d JlInv();

/**
 * Computes the left Jacobian inverse using the element of *this and applies it to the
 * parameter provided. \f$ J_l(v)u \f$
 * Since the left Jacobian inverse is the identity map for \f$se(2)\f$, it will simply return
 * the parameter u.
 * @param u An element of the Lie algebra.
 */ 
se2 JlInv(const se2& u);

/**
 * Computes and returns the matrix of the Right Jacobian
 * This will always return the identity map for \f$se(2)\f$.
 */ 
Eigen::Matrix3d Jr();

/**
 * Computes the right Jacobian using the element of *this and applies it to the
 * parameter provided. \f$ J_l(v)u \f$
 * Since the right Jacobian is the identity map for \f$se(2)\f$, it will simply return
 * the parameter u.
 * @param u An element of the Lie algebra.
 */ 
se2 Jr(const se2& u);

/**
 * Computes and returns the matrix of the right Jacobian inverse.
 * This will always return the identity map for \f$se(2)\f$.
 */ 
Eigen::Matrix3d JrInv();

/**
 * Computes the right Jacobian inverse using the element of *this and applies it to the
 * parameter provided. \f$ J_l(v)u \f$
 * Since the right Jacobian inverse is the identity map for \f$se(2)\f$, it will simply return
 * the parameter u.
 * @param u An element of the Lie algebra.
 */ 
se2 JrInv(const se2& u);

/**
 * Adds two elements of the Algebra together
 * @param u An element of the Lie algebra.
 */ 
se2 operator + (const se2& u);

/**
 * Subtracts two elements of the Algebra together
 * @param u An element of the Lie algebra.
 */ 
se2 operator - (const se2& u);

/**
 * Creates a deep copy of the element
 * @param u An element of the Lie algebra.
 */
void operator = (const se2& u);

/**
 * Performs Scalar multiplication and returns the result.
 * @param scalar The scalar that will scale the element of the Lie algebra
 */ 
se2 operator * (const double scalar);

/**
 * Prints the data of the element.
 */ 
void Print();

/**
 * Returns the Identity element.
 */
static se2 Identity();

/**
 * Computes and returns the skew symmetric matrix of the
 * parameter.
 * @param x 
 */ 
static Eigen::Matrix2d SSM(double x);

/**
 * Verifies that the parameter data belongs to 
 * an element of \f$se(2)\f$
 */ 
static bool isElement(const Eigen::Matrix3d& data);

// Helper functions for testing
Eigen::Matrix2d getWl();
Eigen::Matrix2d getWr();
Eigen::Matrix2d getDl();
Eigen::Matrix2d getDr();



private:

// The following are used to compute the Jacobians
Eigen::Matrix2d Wl();
Eigen::Matrix2d Wr();
Eigen::Matrix2d Dl();
Eigen::Matrix2d Dr();


};

}

#endif //_LIEGROUPS_INCLUDE_LIEALGEBRAS_SE2_