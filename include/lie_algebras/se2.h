#ifndef _LIEGROUPS_INCLUDE_LIEALGEBRAS_SE2_
#define _LIEGROUPS_INCLUDE_LIEALGEBRAS_SE2_

#include <Eigen/Dense>
#include <iostream>

namespace lie_groups {

constexpr double kse2_threshold_=1e-7; /** < If two values are within this threshold, they are considered equal.*/


class se2  {

public:

Eigen::Matrix<double,3,1> data_; /** < The vector is translational velocity followed by angular velocity*/
Eigen::Map<Eigen::Matrix<double,2,1>> p_;  /** < The translational velocity */
Eigen::Map<Eigen::Matrix<double,1,1>> th_; /** < The angular velocity */


/**
 * Default constructor. Initializes algebra element to identity.
 */
se2(): p_(data_.data()), th_(data_.data()+2), data_(Eigen::Matrix<double,3,1>::Zero()){}


/**
 * Copy constructor.
 */ 
se2(const se2 & u) : p_(data_.data()), th_(data_.data()+2), data_(u.data_){}

/**
* Initializes algebra element to the one given. 
* @param[in] data The data of an element of Cartesian space of \f$se(2)\f$
*/
se2(const Eigen::Matrix<double,3,1> data) : p_(data_.data()), th_(data_.data()+2), data_(data){}

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
se2 Bracket(const se2& u){return se2(this->Adjoint()*u.data_);}

/**
 * Computes and returns the matrix adjoint representation of the Lie algebra.
 * This is always the Identity element.
 */ 
Eigen::Matrix3d Adjoint(){
    Eigen::Matrix3d m = Eigen::Matrix3d::Zero();
    m.block(0,0,2,2) = se2::SSM(th_(0));
    m.block(0,2,2,1) = -se2::SSM(1)*p_;
    return m;
}

/**
 * Computes the Wedge operation which maps an element of the Cartesian space to the Lie algebra.
 * @return The result of the Wedge operation.
 */
Eigen::Matrix3d Wedge(){return Wedge(data_);}

/**
 * Computes the Wedge operation which maps an element of the Cartesian space to the Lie algebra.
 * @param data The data of an element  of the Cartesian space isomorphic to the Lie algebra
 * @return The result of the Wedge operation.
 */
static Eigen::Matrix3d Wedge(const Eigen::Matrix<double,3,1>& data){   
    Eigen::Matrix3d m = Eigen::Matrix3d::Zero();
    m.block(0,0,2,2) = se2::SSM(data(2));
    m.block(0,2,2,1) = data.block(0,0,2,1);
    return m;}

/**
 * Computes the Vee operation which maps an element of the Lie algebra to the Cartesian space.
 * @return The result of the Vee operation.
 */
Eigen::Matrix<double,3,1> Vee(){return data_;}

/**
 * Computes the Vee operation which maps an element of the Lie algebra to the Cartesian space.
 * @param data The data of an element
 * @return The result of the Vee operation.
 */
static Eigen::Matrix<double,3,1> Vee(const Eigen::Matrix3d& data) {
    Eigen::Matrix<double,3,1> m;
    m << data(0,2), data(1,2), data(1,0);
    return m;
}

/**
 * Computes the exponential of the element of the Lie algebra.
 */
Eigen::Matrix3d Exp(){return se2::Exp(this->data_);}

/**
 * Computes the exponential of the element of the Lie algebra.
 * @return The data associated to the group element.
 */
static Eigen::Matrix3d Exp(const Eigen::Matrix<double,3,1>& data);

/**
 * Computes the logaritm of the element of the Lie algebra.
 * @param data The data associated with an element of \f$ SE(2) \f$
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
 * Since the left Jacobian is the identity map for \f$se(2)\f$, it will simply return
 * the parameter u.
 * @param u An element of the Lie algebra.
 * @return It will return the parameter u.
 */ 
se2 Jl(const se2& u){return se2(this->Jl()*u.data_);}

/**
 * Computes and returns the matrix of the Left Jacobian inverse.
 */ 
Eigen::Matrix3d JlInv();

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
Eigen::Matrix3d Jr();

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
Eigen::Matrix3d JrInv();

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
 * Creates a deep copy of the element
 * @param u An element of the Lie algebra.
 */
void operator = (const se2& u){data_ = u.data_;}

/**
 * Performs Scalar multiplication and returns the result.
 * @param scalar The scalar that will scale the element of the Lie algebra
 */ 
se2 operator * (const double scalar){return se2(scalar*data_);}

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
static Eigen::Matrix2d SSM(double x) {
    Eigen::Matrix2d m;
    m << 0, -x, x, 0;
    return m;
}

/**
 * Verifies that the parameter data belongs to 
 * an element of \f$se(2)\f$
 */ 
static bool isElement(const Eigen::Matrix3d& data);


private:

// The following are used to compute the Jacobians
static Eigen::Matrix2d Wl(const double th);
static Eigen::Matrix2d Wr(const double th);
static Eigen::Matrix2d Dl(const double th);
static Eigen::Matrix2d Dr(const double th);
Eigen::Matrix2d Wl(){return Wl(th_(0));}
Eigen::Matrix2d Wr(){return Wr(th_(0));}
Eigen::Matrix2d Dl(){return Dl(th_(0));}
Eigen::Matrix2d Dr(){return Dr(th_(0));}


};

}

#endif //_LIEGROUPS_INCLUDE_LIEALGEBRAS_SE2_