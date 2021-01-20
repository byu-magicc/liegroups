#ifndef _LIEGROUPS_INCLUDE_LIEALGEBRAS_SO2_
#define _LIEGROUPS_INCLUDE_LIEALGEBRAS_SO2_

#include <Eigen/Dense>
#include <iostream>

namespace lie_groups {

constexpr double kso2_threshold_=1e-7; /** < If two values are within this threshold, they are considered equal.*/

template< typename tDataType = double>
class so2  {

public:

static constexpr unsigned int dim_ = 1;
static constexpr unsigned int size1_ = 1;
static constexpr unsigned int size2_ = 1;

typedef Eigen::Matrix<tDataType,1,1> Mat1d;
typedef Eigen::Matrix<tDataType,2,2> Mat2d;

Mat1d data_;


/**
 * Default constructor. Initializes algebra element to identity.
 */
so2() : data_(Mat1d::Zero()){}


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
so2(const Mat1d data) : data_(data){}

/**
* Initializes algebra element to the one given. If verify is set to true,
* it will verify that the element provided is an element of the Lie algebra.
* @param[in] u The data of an element of \f$so(2)\f$
* @param verify If true, the constructor will verify that the element given is an element of the Lie algebra.
*/
so2(const Mat2d & data, bool verify);



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
Mat2d Adjoint() {return Mat2d::Identity();}

/**
 * Computes the Wedge operation which maps an element of the Cartesian space to the Lie algebra.
 * @return The result of the Wedge operation.
 */
Mat2d Wedge() {return Wedge(data_);}

/**
 * Computes the Wedge operation which maps an element of the Cartesian space to the Lie algebra.
 * @param data The data of an element  of the Cartesian space isomorphic to the Lie algebra
 * @return The result of the Wedge operation.
 */
static Mat2d Wedge(const Mat1d& data) {
    Mat2d m;
    m << static_cast<tDataType>(0), -data(0), data(0), static_cast<tDataType>(0);
    return m; 
}


/**
 * Computes the Vee operation which maps an element of the Lie algebra to the Cartesian space.
 * @return The result of the Vee operation.
 */
Mat1d Vee() {return data_;}

/**
 * Computes the Vee operation which maps an element of the Lie algebra to the Cartesian space.
 * @param data The data of an element
 * @return The result of the Vee operation.
 */
static Mat1d Vee(const Mat2d& data) {
    Mat1d m;
    m << data(1,0);
    return m;  
}

/**
 * Computes the exponential of the element of the Lie algebra.
 * @param data The data pertaining to an element of the Cartesian 
 * space isomorphic to the Lie algebra
 * @return The data associated to the group element.
 */
static Mat2d Exp(const Mat1d &data);

/**
 * Computes the exponential of the element of the Lie algebra.
 * @return The data associated to the group element.
 */
Mat2d Exp(){ return so2<tDataType>::Exp(this->data_);}

/**
 * Computes the logaritm of the element of the Lie algebra.
 * @param data The data associated with an element of \f$ SO(2) \f$
 * @return The data of an element of the Cartesian space associated with the Lie algebra
 */
static Mat1d Log(const Mat2d& data);

/**
 * Computes and returns the Euclidean norm of the element of the Lie algebra
 */ 
double Norm() {return data_.norm();}

/**
 * Computes and returns the matrix of the Left Jacobian.
 * This will always return the identity map for \f$so(2)\f$.
 */ 
Mat2d Jl(){return Mat2d::Identity();}

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
Mat2d JlInv(){return Mat2d::Identity();}

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
Mat2d Jr(){return Mat2d::Identity();}

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
Mat2d JrInv(){return Mat2d::Identity();}

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
so2 operator * (const double scalar) const {return so2(scalar*data_);}

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
static bool isElement(const Mat2d& data);



};



///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                    Definitions
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//---------------------------------------------------------------------
template< typename tDataType>
so2<tDataType>::so2(const Eigen::Matrix<tDataType,2,2>& data, bool verify) {

    if(verify)
    {
        if (isElement(data) )
            data_(0) = data(1,0);
        else {
            std::cerr << "so2<tDataType>::Constructor - Input data not valid. Setting to identity element" << std::endl;
            data_ = Mat1d::Zero();
        }
    }
    else {
        data_(0) = data(1,0);
    }

}

//---------------------------------------------------------------------
template< typename tDataType>
Eigen::Matrix<tDataType,2,2> so2<tDataType>::Exp(const Mat1d &data) {
    Mat2d m;
    m << cos(data(0)), -sin(data(0)), sin(data(0)), cos(data(0));
    return m;
}

//---------------------------------------------------------------------
template< typename tDataType>
Eigen::Matrix<tDataType,1,1> so2<tDataType>::Log(const Mat2d& data) {
    Mat1d m;
    m(0) = atan2(data(1,0),data(0,0));
    return m;
}

//---------------------------------------------------------------------
template< typename tDataType>
bool so2<tDataType>::isElement(const Mat2d& data) {

    if ( (data.transpose() + data).norm() >= kso2_threshold_) {
        return false;
    } else {
        return true;
    }

}


}

#endif //_LIEGROUPS_INCLUDE_LIEALGEBRAS_SO2_