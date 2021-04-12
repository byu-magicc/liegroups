#ifndef _LIEGROUPS_INCLUDE_LIEALGEBRAS_TRN_
#define _LIEGROUPS_INCLUDE_LIEALGEBRAS_TRN_

#include <Eigen/Dense>
#include <iostream>

namespace lie_groups {

constexpr double krn_threshold_=1e-7; /** < If two values are within this threshold, they are considered equal.*/

template <typename tDataType=double, int N=2, int T=1>
class rn  {


public:

Eigen::Matrix<tDataType,T*N,1> data_;
static constexpr unsigned int dim_ = T*N;              /** < The dimension of the base manifold and tangent space. */ 
static constexpr unsigned int size1_ = T*N;
static constexpr unsigned int size2_ = 1;
static constexpr unsigned int dim_manifold_ = N;       /** < The dimension of the manifold */ 
static constexpr unsigned int deg_tangent_space_ = T;  /** < The degree of the tangent space. For example, we are working with  T^(deg_tangent_space_)M */

// The rule of 5!!

/**
 * Default constructor. Initializes algebra element to identity.
 */
rn() : data_(Eigen::Matrix<tDataType,dim_,1>::Zero()) {}

/**
 * Copy constructor.
 */ 
rn(const rn & u) : data_(u.data_) {}

/**
 * Copy assignment.
 */
void operator = (const rn& u){data_ = u.data_;}

/**
 * Move constructor.
 */ 
rn(const rn && u) : data_(u.data_) {}

/**
 * Move assignment.
 */
void operator = (const rn&& u){data_ = u.data_;}

/**
* Initializes algebra element to the one given. 
* @param[in] data The data of an element of \f$T\mathbb{R}^n\f$
*/
rn(const Eigen::Matrix<tDataType,dim_,1> data) : data_(data) {}

/**
* Initializes algebra element to the one given. If verify is set to true,
* it will verify that the element provided is an element of the Lie algebra.
* @param[in] u The data of an element of \f$ T\mathbb{R}^n\f$
* @param verify If true, the constructor will verify that the element given is an element of the Lie algebra.
*/
rn(const Eigen::Matrix<tDataType,dim_,1> & data, bool verify) : data_(data){}


/**
 * Performs the Lie bracket which is always the identity element. 
 * The input is the right parameter 
 * of the Lie bracket function. \f$ \[v,u\] \f$
 * @param u An element of the Lie algebra.
 * @return The result of the Lie bracket operation.
 */ 
rn Bracket(const rn& u) {return rn();}

/**
 * Computes and returns the matrix adjoint representation of the Lie algebra.
 * This is always the Identity map for \f$ T\mathbb{R}^n\f$
 */ 
Eigen::Matrix<tDataType,dim_,dim_> Adjoint() {return Eigen::Matrix<tDataType,dim_,dim_>::Identity();}

/**
 * Computes the Wedge operation which maps an element of the Cartesian space to the Lie algebra.
 * This will just the data.
 * @return The result of the Wedge operation.
 */
Eigen::Matrix<tDataType,dim_,1> Wedge(){return data_;}

/**
 * Computes the Wedge operation which maps an element of the Cartesian space to the Lie algebra.
 * @param data The data of an element  of the Cartesian space isomorphic to the Lie algebra
 * @return The result of the Wedge operation.
 */
static Eigen::Matrix<tDataType,dim_,1> Wedge(const Eigen::Matrix<tDataType,dim_,1>& data){return data;}


/**
 * Computes the Vee operation which maps an element of the Lie algebra to the Cartesian space.
 * This will just return itself
 * @return The result of the Vee operation.
 */
Eigen::Matrix<tDataType,dim_,1> Vee(){return data_;}

/**
 * Computes the Vee operation which maps an element of the Lie algebra to the Cartesian space.
 * This will just return the data
 * @param data The data of an element
 * @return The result of the Vee operation.
 */
static Eigen::Matrix<tDataType,dim_,1> Vee(const Eigen::Matrix<tDataType,dim_,1>& data){return data;}

/**
 * Computes the exponential of the element of the Lie algebra.
 * The exponential map is the identity map for \f$ T \mathbb{R}^n\f$
 * @return The data associated to the group element.
 */
Eigen::Matrix<tDataType,dim_,1> Exp(){return Exp(data_);}


/**
 * Computes the exponential of the data of an element of the Lie algebra.
 * The exponential map is the identity map for \f$ T \mathbb{R}^n\f$ * 
 * @return The data associated to the group element.
 */
static Eigen::Matrix<tDataType,dim_,1> Exp(const Eigen::Matrix<tDataType,dim_,1>& data ){return data;}

/**
 * Computes the logaritm of the element of the Lie algebra.
 * The logaritm map is the identity map for \f$ T \mathbb{R}^n\f$
 * @param data The data associated with an element of \f$ T \mathbb{R}^n \f$
 * @return The data of an element of the Cartesian space associated with the Lie algebra
 */
static Eigen::Matrix<tDataType,dim_,1> Log(const Eigen::Matrix<tDataType,dim_,1>& data) {return data;}

/**
 * Computes and returns the Euclidean norm of the element of the Lie algebra
 */ 
tDataType Norm(){return data_.norm();}

/**
 * Computes and returns the matrix of the Left Jacobian.
 * This will always return the identity map for \f$ T \mathbb{R}^n\f$.
 */ 
Eigen::Matrix<tDataType,dim_,dim_> Jl(){return Eigen::Matrix<tDataType,dim_,dim_>::Identity();}

/**
 * Computes the left Jacobian using the element of *this and applies it to the
 * parameter provided. \f$ J_l(v)u \f$.
 * Since the left Jacobian is the identity map for \f$ T \mathbb{R}^n\f$, it will simply return
 * the parameter u.
 * @param u An element of the Lie algebra.
 * @return It will return the parameter u.
 */ 
rn Jl(const rn& u){return u;}

/**
 * Computes and returns the matrix of the Left Jacobian inverse.
 * This will always return the identity map for \f$ \mathbb{R}^n\f$$.
 */ 
Eigen::Matrix<tDataType,dim_,dim_> JlInv(){return Eigen::Matrix<tDataType,dim_,dim_>::Identity();}

/**
 * Computes the left Jacobian inverse using the element of *this and applies it to the
 * parameter provided. \f$ J_l(v)u \f$
 * Since the left Jacobian inverse is the identity map for \f$ T \mathbb{R}^n\f$, it will simply return
 * the parameter u.
 * @param u An element of the Lie algebra.
 */ 
rn JlInv(const rn& u){return u;}

/**
 * Computes and returns the matrix of the Right Jacobian
 * This will always return the identity map for \f$ T \mathbb{R}^n\f$.
 */ 
Eigen::Matrix<tDataType,dim_,dim_> Jr(){return Eigen::Matrix<tDataType,dim_,dim_>::Identity();}

/**
 * Computes the right Jacobian using the element of *this and applies it to the
 * parameter provided. \f$ J_l(v)u \f$
 * Since the right Jacobian is the identity map for \f$ T \mathbb{R}^n\f$, it will simply return
 * the parameter u.
 * @param u An element of the Lie algebra.
 */ 
rn Jr(const rn& u){return u;}

/**
 * Computes and returns the matrix of the right Jacobian inverse.
 * This will always return the identity map for \f$ \mathbb{R}^n\f$.
 */ 
Eigen::Matrix<tDataType,dim_,dim_> JrInv(){return Eigen::Matrix<tDataType,dim_,dim_>::Identity();}

/**
 * Computes the right Jacobian inverse using the element of *this and applies it to the
 * parameter provided. \f$ J_l(v)u \f$
 * Since the right Jacobian inverse is the identity map for \f$ T \mathbb{R}^n\f$, it will simply return
 * the parameter u.
 * @param u An element of the Lie algebra.
 */ 
rn JrInv(const rn& u){return u;}

/**
 * Adds two elements of the Algebra together
 * @param u An element of the Lie algebra.
 */ 
rn operator + (const rn& u){return rn(data_ + u.data_);}

/**
 * Subtracts two elements of the Algebra together
 * @param u An element of the Lie algebra.
 */ 
rn operator - (const rn& u){return rn(data_ - u.data_);}

/**
 * Performs Scalar multiplication and returns the result.
 * @param scalar The scalar that will scale the element of the Lie algebra
 */ 
rn operator * (const tDataType scalar) const {return rn(data_*scalar);}

/**
 * Prints the data of the element.
 */ 
void Print(){std::cout << data_ << std::endl;}

/**
 * Returns the Identity element.
 */
static rn Identity(){return rn();}

/**
 * Verifies that the parameter data belongs to 
 * an element of \f$ T \mathbb{R}^n\f$
 */ 
static bool isElement(const Eigen::Matrix<tDataType,dim_,1>& data){return true;}



};

}

#endif //_LIEGROUPS_INCLUDE_LIEALGEBRAS_TRN_