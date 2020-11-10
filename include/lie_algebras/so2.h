#ifndef _LIEGROUPS_INCLUDE_LIEALGEBRAS_SO2_
#define _LIEGROUPS_INCLUDE_LIEALGEBRAS_SO2_

#include <Eigen/Dense>
#include <iostream>

namespace lie_groups {

constexpr double kso2_threshold_=1e-7; /** < If two values are within this threshold, they are considered equal.*/


class so2  {

public:

Eigen::Matrix<double,1,1> data_;

static constexpr unsigned int dim_ = 1;
static constexpr unsigned int size1_ = 1;
static constexpr unsigned int size2_ = 1;
/**
 * Default constructor. Initializes algebra element to identity.
 */
so2() : data_(Eigen::Matrix<double,1,1>::Zero()){}


/**
 * Copy constructor.
 */ 
so2(const so2 & u) : data_(u.data_){}

/**
 * Copy assignment.
 */
void operator = (const so2& u){data_ = u.data_;}

/**
 * Move constructor.
 */ 
so2(const so2 && u) : data_(u.data_){}

/**
 * Move assignment.
 */
void operator = (const so2&& u){data_ = u.data_;}

/**
* Initializes algebra element to the one given. 
* @param[in] data The data of an element of Cartesian space of \f$so(2)\f$
*/
so2(const Eigen::Matrix<double,1,1> data) : data_(data){}

/**
* Initializes algebra element to the one given. If verify is set to true,
* it will verify that the element provided is an element of the Lie algebra.
* @param[in] u The data of an element of \f$so(2)\f$
* @param verify If true, the constructor will verify that the element given is an element of the Lie algebra.
*/
so2(const Eigen::Matrix2d & data, bool verify);



//---------------------------------------------------------------------

/**
 * Performs the Lie bracket with is always the identity element. 
 * The input is the right parameter 
 * of the Lie bracket function. \f$ \[v,u\] \f$
 * @param u An element of the Lie algebra.
 * @return The result of the Lie bracket operation.
 */ 
so2 Bracket(const so2& u) {return so2();}

/**
 * Computes and returns the matrix adjoint representation of the Lie algebra.
 * This is always the Identity element.
 */ 
Eigen::Matrix2d Adjoint() {return Eigen::Matrix2d::Identity();}

/**
 * Computes the Wedge operation which maps an element of the Cartesian space to the Lie algebra.
 * @return The result of the Wedge operation.
 */
Eigen::Matrix2d Wedge() {return Wedge(data_);}

/**
 * Computes the Wedge operation which maps an element of the Cartesian space to the Lie algebra.
 * @param data The data of an element  of the Cartesian space isomorphic to the Lie algebra
 * @return The result of the Wedge operation.
 */
static Eigen::Matrix2d Wedge(const Eigen::Matrix<double,1,1>& data) {
    Eigen::Matrix2d m;
    m << 0, -data(0), data(0), 0;
    return m; 
}


/**
 * Computes the Vee operation which maps an element of the Lie algebra to the Cartesian space.
 * @return The result of the Vee operation.
 */
Eigen::Matrix<double,1,1> Vee() {return data_;}

/**
 * Computes the Vee operation which maps an element of the Lie algebra to the Cartesian space.
 * @param data The data of an element
 * @return The result of the Vee operation.
 */
static Eigen::Matrix<double,1,1> Vee(const Eigen::Matrix2d& data) {
    Eigen::Matrix<double,1,1> m;
    m << data(1,0);
    return m;  
}

/**
 * Computes the exponential of the element of the Lie algebra.
 * @param data The data pertaining to an element of the Cartesian 
 * space isomorphic to the Lie algebra
 * @return The data associated to the group element.
 */
static Eigen::Matrix2d Exp(const Eigen::Matrix<double,1,1> &data);

/**
 * Computes the exponential of the element of the Lie algebra.
 * @return The data associated to the group element.
 */
Eigen::Matrix2d Exp(){ return so2::Exp(this->data_);}

/**
 * Computes the logaritm of the element of the Lie algebra.
 * @param data The data associated with an element of \f$ SO(2) \f$
 * @return The data of an element of the Cartesian space associated with the Lie algebra
 */
static Eigen::Matrix<double,1,1> Log(const Eigen::Matrix2d& data);

/**
 * Computes and returns the Euclidean norm of the element of the Lie algebra
 */ 
double Norm() {return data_.norm();}

/**
 * Computes and returns the matrix of the Left Jacobian.
 * This will always return the identity map for \f$so(2)\f$.
 */ 
Eigen::Matrix2d Jl(){return Eigen::Matrix2d::Identity();}

/**
 * Computes the left Jacobian using the element of *this and applies it to the
 * parameter provided. \f$ J_l(v)u \f$.
 * Since the left Jacobian is the identity map for \f$so(2)\f$, it will simply return
 * the parameter u.
 * @param u An element of the Lie algebra.
 * @return It will return the parameter u.
 */ 
so2 Jl(const so2& u){return u;}

/**
 * Computes and returns the matrix of the Left Jacobian inverse.
 * This will always return the identity map for \f$so(2)\f$.
 */ 
Eigen::Matrix2d JlInv(){return Eigen::Matrix2d::Identity();}

/**
 * Computes the left Jacobian inverse using the element of *this and applies it to the
 * parameter provided. \f$ J_l(v)u \f$
 * Since the left Jacobian inverse is the identity map for \f$so(2)\f$, it will simply return
 * the parameter u.
 * @param u An element of the Lie algebra.
 */ 
so2 JlInv(const so2& u){return u;}

/**
 * Computes and returns the matrix of the Right Jacobian
 * This will always return the identity map for \f$so(2)\f$.
 */ 
Eigen::Matrix2d Jr(){return Eigen::Matrix2d::Identity();}

/**
 * Computes the right Jacobian using the element of *this and applies it to the
 * parameter provided. \f$ J_l(v)u \f$
 * Since the right Jacobian is the identity map for \f$so(2)\f$, it will simply return
 * the parameter u.
 * @param u An element of the Lie algebra.
 */ 
so2 Jr(const so2& u){return u;}

/**
 * Computes and returns the matrix of the right Jacobian inverse.
 * This will always return the identity map for \f$so(2)\f$.
 */ 
Eigen::Matrix2d JrInv(){return Eigen::Matrix2d::Identity();}

/**
 * Computes the right Jacobian inverse using the element of *this and applies it to the
 * parameter provided. \f$ J_l(v)u \f$
 * Since the right Jacobian inverse is the identity map for \f$so(2)\f$, it will simply return
 * the parameter u.
 * @param u An element of the Lie algebra.
 */ 
so2 JrInv(const so2& u) {return u;}

/**
 * Adds two elements of the Algebra together
 * @param u An element of the Lie algebra.
 */ 
so2 operator + (const so2& u){return so2(data_ + u.data_);}

/**
 * Subtracts two elements of the Algebra together
 * @param u An element of the Lie algebra.
 */ 
so2 operator - (const so2& u){return so2(data_ - u.data_);}

/**
 * Performs Scalar multiplication and returns the result.
 * @param scalar The scalar that will scale the element of the Lie algebra
 */ 
so2 operator * (const double scalar){return so2(scalar*data_);}

/**
 * Prints the data of the element.
 */ 
void Print(){std::cout << data_ << std::endl;}

/**
 * Returns the Identity element.
 */
static so2 Identity(){return so2();}

/**
 * Verifies that the parameter data belongs to 
 * an element of \f$so(2)\f$
 */ 
static bool isElement(const Eigen::Matrix2d& data);



};

}

#endif //_LIEGROUPS_INCLUDE_LIEALGEBRAS_SO2_