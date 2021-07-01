#ifndef _LIEGROUPS_INCLUDE_LIEALGEBRAS_RN_
#define _LIEGROUPS_INCLUDE_LIEALGEBRAS_RN_

#include <Eigen/Dense>
#include <iostream>

namespace lie_groups {

constexpr double krn_threshold_=1e-7; /** < If two values are within this threshold, they are considered equal.*/

template <typename tDataType=double, int tNumDimensions=2, int tNumTangentSpaces=1>
class rn  {

public:

/*
* If the tangent space is grater than 1, this theory isn't fully developed so be careful on how you use it. 
*/
static_assert(tNumTangentSpaces > 0, "lie_groups::rn the number of tangent spaces must be greater than 0.");


static constexpr unsigned int dim_ = tNumDimensions;
static constexpr unsigned int size1_ = tNumDimensions*tNumTangentSpaces;
static constexpr unsigned int size2_ = 1;
static constexpr unsigned int num_tangent_spaces_ =tNumTangentSpaces;
static constexpr unsigned int total_num_dim_ = tNumTangentSpaces*tNumDimensions;
typedef Eigen::Matrix<tDataType,total_num_dim_,1> MatData;
typedef Eigen::Matrix<tDataType,total_num_dim_,total_num_dim_> MatNN;
MatData data_;

// The rule of 5!!

/**
 * Default constructor. Initializes algebra element to identity.
 */
rn() : data_(MatData::Zero()) {}

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
* @param[in] data The data of an element of \f$\mathbb{R}^2\f$
*/
rn(const MatData data) : data_(data) {}

/**
* Initializes algebra element to the one given. If verify is set to true,
* it will verify that the element provided is an element of the Lie algebra.
* @param[in] u The data of an element of \f$ \mathbb{R}^n\f$
* @param verify If true, the constructor will verify that the element given is an element of the Lie algebra.
*/
rn(const MatData & data, bool verify) : data_(data){}


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
 * This is always the Identity map for \f$ \mathbb{R}^n\f$
 */ 
MatNN Adjoint() {return MatNN::Identity();}

/**
 * Computes the Wedge operation which maps an element of the Cartesian space to the Lie algebra.
 * This will just the data.
 * @return The result of the Wedge operation.
 */
MatData Wedge(){return data_;}

/**
 * Computes the Wedge operation which maps an element of the Cartesian space to the Lie algebra.
 * @param data The data of an element  of the Cartesian space isomorphic to the Lie algebra
 * @return The result of the Wedge operation.
 */
static MatData Wedge(const MatData& data){return data;}


/**
 * Computes the Vee operation which maps an element of the Lie algebra to the Cartesian space.
 * This will just return itself
 * @return The result of the Vee operation.
 */
MatData Vee(){return data_;}

/**
 * Computes the Vee operation which maps an element of the Lie algebra to the Cartesian space.
 * This will just return the data
 * @param data The data of an element
 * @return The result of the Vee operation.
 */
static MatData Vee(const MatData& data){return data;}

/**
 * Computes the exponential of the element of the Lie algebra.
 * The exponential map is the identity map for \f$ \mathbb{R}^n\f$
 * @return The data associated to the group element.
 */
MatData Exp(){return Exp(data_);}


/**
 * Computes the exponential of the data of an element of the Lie algebra.
 * The exponential map is the identity map for \f$ \mathbb{R}^n\f$ * 
 * @return The data associated to the group element.
 */
static MatData Exp(const MatData& data ){return data;}

/**
 * Computes the logaritm of the element of the Lie algebra.
 * The logaritm map is the identity map for \f$ \mathbb{R}^n\f$
 * @param data The data associated with an element of \f$ SO(2) \f$
 * @return The data of an element of the Cartesian space associated with the Lie algebra
 */
static MatData Log(const MatData& data) {return data;}

/**
 * Computes and returns the Euclidean norm of the element of the Lie algebra
 */ 
tDataType Norm(){return data_.norm();}

/**
 * Computes and returns the matrix of the Left Jacobian.
 * This will always return the identity map for \f$ \mathbb{R}^n\f$.
 */ 
MatNN Jl(){return MatNN::Identity();}

/**
 * Computes the left Jacobian using the element of *this and applies it to the
 * parameter provided. \f$ J_l(v)u \f$.
 * Since the left Jacobian is the identity map for \f$ \mathbb{R}^n\f$, it will simply return
 * the parameter u.
 * @param u An element of the Lie algebra.
 * @return It will return the parameter u.
 */ 
rn Jl(const rn& u){return u;}

/**
 * Computes and returns the matrix of the Left Jacobian inverse.
 * This will always return the identity map for \f$ \mathbb{R}^n\f$$.
 */ 
MatNN JlInv(){return MatNN::Identity();}

/**
 * Computes the left Jacobian inverse using the element of *this and applies it to the
 * parameter provided. \f$ J_l(v)u \f$
 * Since the left Jacobian inverse is the identity map for \f$ \mathbb{R}^n\f$, it will simply return
 * the parameter u.
 * @param u An element of the Lie algebra.
 */ 
rn JlInv(const rn& u){return u;}

/**
 * Computes and returns the matrix of the Right Jacobian
 * This will always return the identity map for \f$ \mathbb{R}^n\f$.
 */ 
MatNN Jr(){return MatNN::Identity();}

/**
 * Computes the right Jacobian using the element of *this and applies it to the
 * parameter provided. \f$ J_l(v)u \f$
 * Since the right Jacobian is the identity map for \f$ \mathbb{R}^n\f$, it will simply return
 * the parameter u.
 * @param u An element of the Lie algebra.
 */ 
rn Jr(const rn& u){return u;}

/**
 * Computes and returns the matrix of the right Jacobian inverse.
 * This will always return the identity map for \f$ \mathbb{R}^n\f$.
 */ 
MatNN JrInv(){return MatNN::Identity();}

/**
 * Computes the right Jacobian inverse using the element of *this and applies it to the
 * parameter provided. \f$ J_l(v)u \f$
 * Since the right Jacobian inverse is the identity map for \f$ \mathbb{R}^n\f$, it will simply return
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
 * an element of \f$ \mathbb{R}^n\f$
 */ 
static bool isElement(const MatData& data){return true;}



};

}

#endif //_LIEGROUPS_INCLUDE_LIEALGEBRAS_RN_