#ifndef _LIEGROUPS_INCLUDE_LIEALGEBRAS_SO2_
#define _LIEGROUPS_INCLUDE_LIEALGEBRAS_SO2_

#include "lie_algebras/algebra_base.h"

namespace lie_groups {

class so2 : public AlgebraBase {

public:

/**
 * Default constructor. Initializes algebra element to identity.
 */
so2();


/**
 * Copy constructor.
 */ 
so2(const so2 & u);

/**
* Initializes algebra element to the one given. 
* @param[in] data The data of an element of Cartesian space of \f$so(2)\f$
*/
so2(const double data);

/**
* Initializes algebra element to the one given. If verify is set to true,
* it will verify that the element provided is an element of the Lie algebra.
* @param[in] u The data of an element of \f$so(2)\f$
* @param verify If true, the constructor will verify that the element given is an element of the Lie algebra.
*/
so2(const Eigen::Matrix<double,2,2>& data, bool verify);

//---------------------------------------------------------------------

/**
 * Performs the Lie bracket with is always the identity element. 
 * The input is the right parameter 
 * of the Lie bracket function. \f$ \[v,u\] \f$
 * @param u An element of the Lie algebra.
 * @return The result of the Lie bracket operation.
 */ 
so2 * Bracket(const AlgebraBase* u) override;

/**
 * Computes and returns the matrix adjoint representation of the Lie algebra.
 * This is always the Identity element.
 */ 
Eigen::MatrixXd Adjoint() override;

/**
 * Computes the Wedge operation which maps an element of the Cartesian space to the Lie algebra.
 * @return The result of the Wedge operation.
 */
Eigen::MatrixXd Wedge() override;

/**
 * Computes the Vee operation which maps an element of the Lie algebra to the Cartesian space.
 * @return The result of the Vee operation.
 */
Eigen::MatrixXd Vee() override;

/**
 * Computes the exponential of the element of the Lie algebra.
 * @return The data associated to the group element.
 */
Eigen::MatrixXd Exp() override;

/**
 * Computes and returns the Euclidean norm of the element of the Lie algebra
 */ 
double Norm() override;

/**
 * Computes and returns the matrix of the Left Jacobian.
 * This will always return the identity map for \f$so(2)\f$.
 */ 
Eigen::MatrixXd Jl() override;

/**
 * Computes the left Jacobian using the element of *this and applies it to the
 * parameter provided. \f$ J_l(v)u \f$.
 * Since the left Jacobian is the identity map for \f$so(2)\f$, it will simply return
 * the parameter u.
 * @param u An element of the Lie algebra.
 * @return It will return the parameter u.
 */ 
so2 * Jl(const AlgebraBase * u) override;

/**
 * Computes and returns the matrix of the Left Jacobian inverse.
 * This will always return the identity map for \f$so(2)\f$.
 */ 
Eigen::MatrixXd JlInv() override;

/**
 * Computes the left Jacobian inverse using the element of *this and applies it to the
 * parameter provided. \f$ J_l(v)u \f$
 * Since the left Jacobian inverse is the identity map for \f$so(2)\f$, it will simply return
 * the parameter u.
 * @param u An element of the Lie algebra.
 */ 
so2 * JlInv(const AlgebraBase * u) override;

/**
 * Computes and returns the matrix of the Right Jacobian
 * This will always return the identity map for \f$so(2)\f$.
 */ 
Eigen::MatrixXd Jr() override;

/**
 * Computes the right Jacobian using the element of *this and applies it to the
 * parameter provided. \f$ J_l(v)u \f$
 * Since the right Jacobian is the identity map for \f$so(2)\f$, it will simply return
 * the parameter u.
 * @param u An element of the Lie algebra.
 */ 
so2 * Jr(const AlgebraBase * u) override;

/**
 * Computes and returns the matrix of the right Jacobian inverse.
 * This will always return the identity map for \f$so(2)\f$.
 */ 
Eigen::MatrixXd JrInv() override;

/**
 * Computes the right Jacobian inverse using the element of *this and applies it to the
 * parameter provided. \f$ J_l(v)u \f$
 * Since the right Jacobian inverse is the identity map for \f$so(2)\f$, it will simply return
 * the parameter u.
 * @param u An element of the Lie algebra.
 */ 
so2 * JrInv(const AlgebraBase * u) override;

/**
 * Adds two elements of the Algebra together
 * @param u An element of the Lie algebra.
 */ 
so2 * operator + (const AlgebraBase * u);

/**
 * Subtracts two elements of the Algebra together
 * @param u An element of the Lie algebra.
 */ 
so2 * operator - (const AlgebraBase * u);

/**
 * Performs Scalar multiplication and returns the result.
 * @param scalar The scalar that will scale the element of the Lie algebra
 */ 
so2 * operator * (const double scalar) override;

/**
 * Prints the data of the element.
 */ 
void Print();

/**
 * Returns the Identity element.
 */
so2 * Identity();


};

}

#endif //_LIEGROUPS_INCLUDE_LIEALGEBRAS_SO2_