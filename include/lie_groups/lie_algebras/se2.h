#ifndef _LIEGROUPS_INCLUDE_LIEALGEBRAS_SE2_
#define _LIEGROUPS_INCLUDE_LIEALGEBRAS_SE2_

#include <Eigen/Dense>
#include <iostream>

namespace lie_groups {

constexpr double kse2_threshold_=1e-7; /** < If two values are within this threshold, they are considered equal.*/

template <typename tDataType=double, int tNumDimensions=3, int tNumTangentSpaces=1>
class se2  {

static_assert(tNumTangentSpaces == 1, "lie_groups::se2 the number of tangent spaces must be 1.");
static_assert(tNumDimensions == 3, "lie_groups::se2 the number of dimensions must be 3.");


public:

Eigen::Matrix<tDataType,3,1> data_; /** < The vector is translational velocity followed by angular velocity*/
Eigen::Map<Eigen::Matrix<tDataType,2,1>> p_;  /** < The translational velocity */
Eigen::Map<Eigen::Matrix<tDataType,1,1>> th_; /** < The angular velocity */

typedef Eigen::Matrix<tDataType,3,1> Vec3d;
typedef Eigen::Matrix<tDataType,3,3> Mat3d;
typedef Eigen::Matrix<tDataType,2,2> Mat2d;

static constexpr unsigned int dim_ = tNumDimensions;
static constexpr unsigned int dim_t_vel_=2; /** < The dimension of the translational velocity */
static constexpr unsigned int dim_a_vel_=1; /** < The dimension of the angular velocity */
static constexpr unsigned int size1_ = tNumDimensions;
static constexpr unsigned int size2_ = 1;
static constexpr unsigned int total_num_dim_ = tNumTangentSpaces*tNumDimensions;



/**
 * Default constructor. Initializes algebra element to identity.
 */
se2(): p_(data_.data()), th_(data_.data()+2), data_(Vec3d::Zero()){}


/**
 * Copy constructor.
 */ 
se2(const se2 & u) : p_(data_.data()), th_(data_.data()+2), data_(u.data_){}

/**
 * Copy assignment
 */
void operator = (const se2& u){data_ = u.data_;}

/**
 * Move constructor.
 */ 
se2(const se2 && u) : p_(data_.data()), th_(data_.data()+2), data_(u.data_){}

/**
 * Move assignment.
 */
void operator = (const se2&& u){data_ = u.data_;}



/**
* Initializes algebra element to the one given. 
* @param[in] data The data of an element of Cartesian space of \f$se(2)\f$
*/
se2(const Vec3d data) : p_(data_.data()), th_(data_.data()+2), data_(data){}

/**
* Initializes algebra element to the one given. If verify is set to true,
* it will verify that the element provided is an element of the Lie algebra.
* @param[in] u The data of an element of \f$se(2)\f$
* @param verify If true, the constructor will verify that the element given is an element of the Lie algebra.
*/
se2(const Mat3d & data, bool verify);

//---------------------------------------------------------------------

/**
 * Performs the Lie bracket.. 
 * The input is the right parameter 
 * of the Lie bracket function. \f$ \[v,u\] \f$
 * @param u An element of the Lie algebra.
 * @return The result of the Lie bracket operation.
 */ 
se2 Bracket(const se2& u){return se2(this->Adjoint()*u.data_);}

/**
 * Computes and returns the matrix adjoint representation of the Lie algebra.
 * This is always the Identity element.
 */ 
Mat3d Adjoint(){
    Mat3d m = Mat3d::Zero();
    m.block(0,0,2,2) = se2::SSM(th_(0));
    m.block(0,2,2,1) = -se2::SSM(1)*p_;
    return m;
}

/**
 * Computes the Wedge operation which maps an element of the Cartesian space to the Lie algebra.
 * @return The result of the Wedge operation.
 */
Mat3d Wedge(){return Wedge(data_);}

/**
 * Computes the Wedge operation which maps an element of the Cartesian space to the Lie algebra.
 * @param data The data of an element  of the Cartesian space isomorphic to the Lie algebra
 * @return The result of the Wedge operation.
 */
static Mat3d Wedge(const Vec3d& data){   
    Mat3d m = Mat3d::Zero();
    m.block(0,0,2,2) = se2::SSM(data(2));
    m.block(0,2,2,1) = data.block(0,0,2,1);
    return m;}

/**
 * Computes the Vee operation which maps an element of the Lie algebra to the Cartesian space.
 * @return The result of the Vee operation.
 */
Vec3d Vee(){return data_;}

/**
 * Computes the Vee operation which maps an element of the Lie algebra to the Cartesian space.
 * @param data The data of an element
 * @return The result of the Vee operation.
 */
static Vec3d Vee(const Mat3d& data) {
    Vec3d m;
    m << data(0,2), data(1,2), data(1,0);
    return m;
}

/**
 * Computes the exponential of the element of the Lie algebra.
 */
Mat3d Exp(){return se2::Exp(this->data_);}

/**
 * Computes the exponential of the element of the Lie algebra.
 * @return The data associated to the group element.
 */
static Mat3d Exp(const Vec3d& data);

/**
 * Computes the logaritm of the element of the Lie algebra.
 * @param data The data associated with an element of \f$ SE(2) \f$
 * @return The data of an element of the Cartesian space associated with the Lie algebra
 */
static Vec3d Log(const Mat3d& data);

/**
 * Computes and returns the Euclidean norm of the element of the Lie algebra
 */ 
tDataType Norm(){return data_.norm();}

/**
 * Computes and returns the matrix of the Left Jacobian.
 */ 
Mat3d Jl();

/**
 * Computes the left Jacobian using the element of *this and applies it to the
 * parameter provided. \f$ J_l(v)u \f$.
 * Since the left Jacobian is the identity map for \f$se(2)\f$, it will simply return
 * the parameter u.
 * @param u An element of the Lie algebra.
 * @return It will return the parameter u.
 */ 
se2 Jl(const se2& u){return se2(this->Jl()*u.data_);}

/**
 * Computes and returns the matrix of the Left Jacobian inverse.
 */ 
Mat3d JlInv();

/**
 * Computes the left Jacobian inverse using the element of *this and applies it to the
 * parameter provided. \f$ J_l(v)u \f$
 * Since the left Jacobian inverse is the identity map for \f$se(2)\f$, it will simply return
 * the parameter u.
 * @param u An element of the Lie algebra.
 */ 
se2 JlInv(const se2& u){return se2(this->JlInv()*u.data_);}

/**
 * Computes and returns the matrix of the Right Jacobian
 */ 
Mat3d Jr();

/**
 * Computes the right Jacobian using the element of *this and applies it to the
 * parameter provided. \f$ J_l(v)u \f$
 * Since the right Jacobian is the identity map for \f$se(2)\f$, it will simply return
 * the parameter u.
 * @param u An element of the Lie algebra.
 */ 
se2 Jr(const se2& u){return se2(this->Jr()*u.data_);}

/**
 * Computes and returns the matrix of the right Jacobian inverse.
 */ 
Mat3d JrInv();

/**
 * Computes the right Jacobian inverse using the element of *this and applies it to the
 * parameter provided. \f$ J_l(v)u \f$
 * Since the right Jacobian inverse is the identity map for \f$se(2)\f$, it will simply return
 * the parameter u.
 * @param u An element of the Lie algebra.
 */ 
se2 JrInv(const se2& u){return se2(this->JrInv()*u.data_);}

/**
 * Adds two elements of the Algebra together
 * @param u An element of the Lie algebra.
 */ 
se2 operator + (const se2& u){return se2(data_ + u.data_);}

/**
 * Subtracts two elements of the Algebra together
 * @param u An element of the Lie algebra.
 */ 
se2 operator - (const se2& u){return se2(data_ - u.data_);}

/**
 * Performs Scalar multiplication and returns the result.
 * @param scalar The scalar that will scale the element of the Lie algebra
 */ 
se2 operator * (const tDataType scalar) const {return se2(scalar*data_);}

/**
 * Prints the data of the element.
 */ 
void Print(){std::cout << data_ << std::endl;}

/**
 * Returns the Identity element.
 */
static se2 Identity(){return se2();}

/**
 * Computes and returns the skew symmetric matrix of the
 * parameter.
 * @param x 
 */ 
static Mat2d SSM(tDataType x) {
    Mat2d m;
    m << static_cast<tDataType>(0), -x, x, static_cast<tDataType>(0);
    return m;
}

/**
 * Verifies that the parameter data belongs to 
 * an element of \f$se(2)\f$
 */ 
static bool isElement(const Mat3d& data);


private:

// The following are used to compute the Jacobians
static Mat2d Wl(const tDataType th);
static Mat2d Wr(const tDataType th);
static Mat2d Dl(const tDataType th);
static Mat2d Dr(const tDataType th);
Mat2d Wl(){return Wl(th_(0));}
Mat2d Wr(){return Wr(th_(0));}
Mat2d Dl(){return Dl(th_(0));}
Mat2d Dr(){return Dr(th_(0));}


};

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                    Definitions
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//---------------------------------------------------------------------
template <typename tDataType, int tNumDimensions, int tNumTangentSpaces>
se2<tDataType,tNumDimensions,tNumTangentSpaces>::se2(const Eigen::Matrix<tDataType,3,3>& data, bool verify) : p_(data_.data()), th_(data_.data()+2) {

    if(verify)
    {
        
        if (se2<tDataType,tNumDimensions,tNumTangentSpaces>::isElement(data)) {
            data_(0) = data(0,2);
            data_(1) = data(1,2);
            data_(2) = data(1,0);
        }
        else {
            std::cerr << "se2::Constructor - Input data not valid. Setting to identity element" << std::endl;
            data_ = Eigen::Matrix<tDataType,3,1>::Zero();
        }
    }
    else {
        data_(0) = data(0,2);
        data_(1) = data(1,2);
        data_(2) = data(1,0);
    }

}

//---------------------------------------------------------------------
template <typename tDataType, int tNumDimensions, int tNumTangentSpaces>
Eigen::Matrix<tDataType,3,3> se2<tDataType,tNumDimensions,tNumTangentSpaces>::Exp(const Eigen::Matrix<tDataType,3,1>& data) {
    
    Eigen::Matrix<tDataType,3,3> m;
    m.block(0,0,2,2) << cos(data(2)), - sin(data(2)), sin(data(2)), cos(data(2));
    m.block(0,2,2,1) = Wl(data(2))*data.block(0,0,2,1);
    m.block(2,0,1,2).setZero();
    m(2,2) = static_cast<tDataType>(1.0);
    return m;
}

//---------------------------------------------------------------------
template <typename tDataType, int tNumDimensions, int tNumTangentSpaces>
Eigen::Matrix<tDataType,3,1> se2<tDataType,tNumDimensions,tNumTangentSpaces>::Log(const Eigen::Matrix<tDataType,3,3>& data) {
    Eigen::Matrix<tDataType,3,1> u;
    u(2) = atan2(data(1,0),data(0,0)); // Compute the angle
    Eigen::Matrix<tDataType,2,2> wl = Wl(u(2));
    u.block(0,0,2,1) = wl.inverse()*data.block(0,2,2,1);
    return u;    
}

//---------------------------------------------------------------------
template <typename tDataType, int tNumDimensions, int tNumTangentSpaces>
Eigen::Matrix<tDataType,3,3> se2<tDataType,tNumDimensions,tNumTangentSpaces>::Jl() {

    Eigen::Matrix<tDataType,3,3> m;
    m.setIdentity();
    m.block(0,0,2,2) = this->Wl();
    m.block(0,2,2,1) = this->Dl()*p_;

    return m;
}

//---------------------------------------------------------------------
template <typename tDataType, int tNumDimensions, int tNumTangentSpaces>
Eigen::Matrix<tDataType,3,3> se2<tDataType,tNumDimensions,tNumTangentSpaces>::JlInv() {

    Eigen::Matrix<tDataType,2,2> w_inv = (this->Wl()).inverse();
    Eigen::Matrix<tDataType,3,3> m;
    m.setIdentity();
    m.block(0,0,2,2) = w_inv;
    m.block(0,2,2,1) = -w_inv*this->Dl()*p_;

    return m;
}

//---------------------------------------------------------------------
template <typename tDataType, int tNumDimensions, int tNumTangentSpaces>
Eigen::Matrix<tDataType,3,3> se2<tDataType,tNumDimensions,tNumTangentSpaces>::Jr() {

    Eigen::Matrix<tDataType,3,3> m;
    m.setIdentity();
    m.block(0,0,2,2) = this->Wr();
    m.block(0,2,2,1) = this->Dr()*p_;

    return m;
}

//---------------------------------------------------------------------
template <typename tDataType, int tNumDimensions, int tNumTangentSpaces>
Eigen::Matrix<tDataType,3,3> se2<tDataType,tNumDimensions,tNumTangentSpaces>::JrInv() {
    
    Eigen::Matrix<tDataType,2,2> w_inv = (this->Wr()).inverse();
    Eigen::Matrix<tDataType,3,3> m;
    m.setIdentity();
    m.block(0,0,2,2) = w_inv;
    m.block(0,2,2,1) = -w_inv*this->Dr()*p_;

    return m;
}

//---------------------------------------------------------------------
template <typename tDataType, int tNumDimensions, int tNumTangentSpaces>
bool se2<tDataType,tNumDimensions,tNumTangentSpaces>::isElement(const Eigen::Matrix<tDataType,3,3>& data) {

    bool is_element = true;
     
    if ( (data.block(0,0,2,2).transpose() + data.block(0,0,2,2)).norm() >= kse2_threshold_) {
        is_element = false;
    }
    else if (data.block(2,0,1,3) != Eigen::Matrix<tDataType,1,3>::Zero()) {
        is_element = false;
    }

    return is_element;
}


//--------------------------------------------------------
//                  Private functions
//--------------------------------------------------------
template <typename tDataType, int tNumDimensions, int tNumTangentSpaces>
Eigen::Matrix<tDataType,2,2> se2<tDataType,tNumDimensions,tNumTangentSpaces>::Wl(const tDataType th) {

    Eigen::Matrix<tDataType,2,2> m;

    if (th > kse2_threshold_ || th < -kse2_threshold_) {
        tDataType a = (1.0-cos(th))/th;
        tDataType b = sin(th)/th;
        m = a*se2<tDataType,tNumDimensions,tNumTangentSpaces>::SSM(static_cast<tDataType>(1.0)) + b*Eigen::Matrix<tDataType,2,2>::Identity();
    }
    else
    {
        m = Eigen::Matrix<tDataType,2,2>::Identity() + se2<tDataType,tNumDimensions,tNumTangentSpaces>::SSM(static_cast<tDataType>(1.0))*th/2.0;
    }
    
return m;
}

//---------------------------------------------------
template <typename tDataType, int tNumDimensions, int tNumTangentSpaces>
Eigen::Matrix<tDataType,2,2> se2<tDataType,tNumDimensions,tNumTangentSpaces>::Wr(const tDataType th) {

    Eigen::Matrix<tDataType,2,2> m;

    if (th > kse2_threshold_ || th < -kse2_threshold_) {
        tDataType a = (cos(th)-1.0)/th;
        tDataType b = sin(th)/th;
        m = a*se2<tDataType,tNumDimensions,tNumTangentSpaces>::SSM(static_cast<tDataType>(1.0)) + b*Eigen::Matrix<tDataType,2,2>::Identity();
    }
    else
    {
        m = Eigen::Matrix<tDataType,2,2>::Identity() - se2<tDataType,tNumDimensions,tNumTangentSpaces>::SSM(static_cast<tDataType>(1.0))*th/2.0;
    }
    
return m;
}

//---------------------------------------------------
template <typename tDataType, int tNumDimensions, int tNumTangentSpaces>
Eigen::Matrix<tDataType,2,2> se2<tDataType,tNumDimensions,tNumTangentSpaces>::Dl(const tDataType th) {

    Eigen::Matrix<tDataType,2,2> m;

    if (fabs(th) > kse2_threshold_) {
        tDataType a = (cos(th)-1.0)/(th*th);
        tDataType b = (th-sin(th))/(th*th);
        m = a*se2<tDataType,tNumDimensions,tNumTangentSpaces>::SSM(1) + b*Eigen::Matrix<tDataType,2,2>::Identity();
    }
    else
    {
        m = Eigen::Matrix<tDataType,2,2>::Identity()*th/6.0 - se2<tDataType,tNumDimensions,tNumTangentSpaces>::SSM(static_cast<tDataType>(1.0))*th/2.0;
    }

    return m;
}

//---------------------------------------------------
template <typename tDataType, int tNumDimensions, int tNumTangentSpaces>
Eigen::Matrix<tDataType,2,2> se2<tDataType,tNumDimensions,tNumTangentSpaces>::Dr(const tDataType th) {

    Eigen::Matrix<tDataType,2,2> m;

    if (th > kse2_threshold_ || th < -kse2_threshold_) {
        tDataType a = (1.0-cos(th))/(th*th);
        tDataType b = (th-sin(th))/(th*th);
        m = a*se2<tDataType,tNumDimensions,tNumTangentSpaces>::SSM(static_cast<tDataType>(1.0)) + b*Eigen::Matrix<tDataType,2,2>::Identity();
    }
    else
    {
        m = Eigen::Matrix<tDataType,2,2>::Identity()*th/6.0 + se2<tDataType,tNumDimensions,tNumTangentSpaces>::SSM(static_cast<tDataType>(1.0))*th/2.0;
    }

    return m;
}



} // namespace

#endif //_LIEGROUPS_INCLUDE_LIEALGEBRAS_SE2_