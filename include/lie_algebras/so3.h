#ifndef _LIEGROUPS_INCLUDE_LIEALGEBRAS_SO3_
#define _LIEGROUPS_INCLUDE_LIEALGEBRAS_SO3_

#include <Eigen/Dense>
#include <iostream>

namespace lie_groups {

constexpr double kso3_threshold_=1e-7; /** < If two values are within this threshold, they are considered equal.*/


template <typename tDataType=double>
class so3  {

public:

typedef Eigen::Matrix<tDataType,3,1> Vec3d;
typedef Eigen::Matrix<tDataType,2,2> Mat2d;
typedef Eigen::Matrix<tDataType,3,3> Mat3d;

static constexpr unsigned int dim_ = 3;
static constexpr unsigned int size1_ = 3;
static constexpr unsigned int size2_ = 1;

Vec3d data_; /** < The vector is translational velocity followed by angular velocity*/

/**
 * Default constructor. Initializes algebra element to identity.
 */
so3() : data_(Vec3d::Zero()){}

/**
 * Copy constructor.
 */ 
so3(const so3 & u) : data_(u.data_) {}

/**
 * Copy assignment.
 */
void operator = (const so3& u){data_ = u.data_;}

/**
 * Move constructor.
 */ 
so3(const so3 && u) : data_(u.data_) {}

/**
 * Move assignment.
 */
void operator = (const so3&& u){data_ = u.data_;}

/**
* Initializes algebra element to the one given. 
* @param[in] data The data of an element of Cartesian space of \f$so(3)\f$
*/
so3(const Vec3d data) : data_(data) {}

/**
* Initializes algebra element to the one given. If verify is set to true,
* it will verify that the element provided is an element of the Lie algebra.
* @param[in] u The data of an element of \f$so(3)\f$
* @param verify If true, the constructor will verify that the element given is an element of the Lie algebra.
*/
so3(const Mat3d & data, bool verify);

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
Mat3d Adjoint() {return this->Wedge();}

/**
 * Computes the Wedge operation which maps an element of the Cartesian space to the Lie algebra.
 * @return The result of the Wedge operation.
 */
Mat3d Wedge() {return Wedge(data_);}

/**
 * Computes the Wedge operation which maps an element of the Cartesian space to the Lie algebra.
 * @param data The data of an element  of the Cartesian space isomorphic to the Lie algebra
 * @return The result of the Wedge operation.
 */
static Mat3d Wedge(const Vec3d& data) {    Mat3d m;
    m << static_cast<tDataType>(0), -data(2), data(1), data(2), static_cast<tDataType>(0), -data(0), -data(1), data(0), static_cast<tDataType>(0);
    return m;
}

/**
 * Computes the Vee operation which maps an element of the Lie algebra to the Cartesian space.
 * @return The result of the Vee operation.
 */
Vec3d Vee() {return data_;}

/**
 * Computes the Vee operation which maps an element of the Lie algebra to the Cartesian space.
 * @param data The data of an element
 * @return The result of the Vee operation.
 */
static Vec3d Vee(const Mat3d& data) {    
    Vec3d m;
    m << data(2,1), data(0,2), data(1,0);
    return m;}

/**
 * Computes the exponential of the element of the Lie algebra.
 */
Mat3d Exp(){ return so3<tDataType>::Exp(this->data_);}

/**
 * Computes the exponential of the element of the Lie algebra.
 * @return The data associated to the group element.
 */
static Mat3d Exp(const Vec3d& data);

/**
 * Computes the logaritm of the element of the Lie algebra.
 * @param data The data associated with an element of \f$ SO(3) \f$
 * @return The data of an element of the Cartesian space associated with the Lie algebra
 */
static Vec3d Log(const Mat3d& data);

/**
 * Computes and returns the Euclidean norm of the element of the Lie algebra
 */ 
double Norm(){return data_.norm();}

/**
 * Computes and returns the matrix of the Left Jacobian.
 */ 
Mat3d Jl();

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
Mat3d JlInv();

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
Mat3d Jr();

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
Mat3d JrInv();

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
 * Performs Scalar multiplication and returns the result.
 * @param scalar The scalar that will scale the element of the Lie algebra
 */ 
so3 operator * (const double scalar) const {return so3(scalar*data_);}

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
static Mat2d SSM(tDataType x);

/**
 * Verifies that the parameter data belongs to 
 * an element of \f$so(3)\f$
 */ 
static bool isElement(const Mat3d& data);


};


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                    Definitions
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////


template <typename tDataType>
so3<tDataType>::so3(const Eigen::Matrix<tDataType,3,3>& data, bool verify) {

    if(verify)
    {
        
        if (so3<tDataType>::isElement(data)) {
            data_(0) = data(2,1);
            data_(1) = data(0,2);
            data_(2) = data(1,0);
        }
        else {
            std::cerr << "so3<tDataType>::Constructor - Input data not valid. Setting to identity element" << std::endl;
            data_ = Vec3d::Zero();
        }
    }
    else {
        data_(0) = data(2,1);
        data_(1) = data(0,2);
        data_(2) = data(1,0);
    }

}

//--------------------------------------------------------------------------------------------------------

template <typename tDataType>
Eigen::Matrix<tDataType,3,3> so3<tDataType>::Exp(const Eigen::Matrix<tDataType,3,1>& data) {
    Mat3d m;
    tDataType th = data.norm();

    if (th < kso3_threshold_) { // See if the element is close to the identity element.
        m.setIdentity();
    } else {  // Use Rodriguez formula 
        tDataType a = sin(th)/th;
        tDataType b = (1.0-cos(th))/pow(th,2);
        m = Mat3d::Identity() + a*Wedge(data) + b*Wedge(data)*Wedge(data);
    }
    return m;

}

//---------------------------------------------------------------------
template <typename tDataType>
Eigen::Matrix<tDataType,3,1> so3<tDataType>::Log(const Eigen::Matrix<tDataType,3,3>& data) {

    Vec3d u;

    tDataType t = data.trace();
    if ( fabs(t-3.0) <= kso3_threshold_) { // Rotation matrix is close to identity
    
        u.setZero();

    } else { // Use Rodriguez formula 

        tDataType th = acos( (t-1)/2);
        u = so3<tDataType>::Vee(th*(data-data.transpose())/(2*sin(th)));
    }

    return u;    
}

//---------------------------------------------------------------------
template <typename tDataType>
Eigen::Matrix<tDataType,3,3> so3<tDataType>::Jl() {

    Mat3d m;
    
    tDataType th = data_.norm();

    if (th < kso3_threshold_) { // See if the element is close to the identity element.
        m.setIdentity();
    } else {   
        tDataType a = (static_cast<tDataType>(1.0)-cos(th))/pow(th,2);
        tDataType b = (th-sin(th))/pow(th,3);
        m = Mat3d::Identity() + a*this->Wedge() + b*this->Wedge()*this->Wedge();
    }

    return m;
}


//---------------------------------------------------------------------
template <typename tDataType>
Eigen::Matrix<tDataType,3,3> so3<tDataType>::JlInv() {


    Mat3d m;

    tDataType th = data_.norm();

    if (th < kso3_threshold_ || sin(th/2) < kso3_threshold_) { // See if the element is close to the identity element.
        m.setIdentity();
    } else {   
        tDataType a = -0.5;
        tDataType cot = cos(th/2.0)/sin(th/2.0);
        tDataType b = -(th*cot-2.0)/(2.0*pow(th,2));
        m = Mat3d::Identity() + a*this->Wedge() + b*this->Wedge()*this->Wedge();
    }

    return m;
}

//---------------------------------------------------------------------
template <typename tDataType>
Eigen::Matrix<tDataType,3,3> so3<tDataType>::Jr() {

    Mat3d m;
    
    tDataType th = data_.norm();

    if (th < kso3_threshold_) { // See if the element is close to the identity element.
        m.setIdentity();
    } else {   
        tDataType a = (cos(th)-1.0)/pow(th,2);
        tDataType b = (th-sin(th))/pow(th,3);
        m = Mat3d::Identity() + a*this->Wedge() + b*this->Wedge()*this->Wedge();
    }

    return m;
}


//---------------------------------------------------------------------
template <typename tDataType>
Eigen::Matrix<tDataType,3,3> so3<tDataType>::JrInv() {
    
    Mat3d m;

    tDataType th = data_.norm();

    if (th < kso3_threshold_ || sin(th/2) < kso3_threshold_) { // See if the element is close to the identity element.
        m.setIdentity();
    } else {   
        tDataType a = 0.5;
        tDataType cot = cos(th/2.0)/sin(th/2.0);
        tDataType b = -(th*cot-2.0)/(2.0*pow(th,2));
        m = Mat3d::Identity() + a*this->Wedge() + b*this->Wedge()*this->Wedge();
    }

    return m;
}

//---------------------------------------------------------------------
template <typename tDataType>
bool so3<tDataType>::isElement(const Eigen::Matrix<tDataType,3,3>& data) {

    bool is_element = true;
     
    if ( (data.transpose()+data).norm()/static_cast<tDataType>(2.0) >= kso3_threshold_) {
        is_element = false;
    }

    return is_element;
}


} // lie_groups

#endif //_LIEGROUPS_INCLUDE_LIEALGEBRAS_SO3_