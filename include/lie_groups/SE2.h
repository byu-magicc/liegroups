#ifndef _LIEGROUPS_INCLUDE_LIEGROUPS_SE2_
#define _LIEGROUPS_INCLUDE_LIEGROUPS_SE2_

#include <Eigen/Dense>
#include <iostream>
#include <lie_algebras/se2.h>

namespace lie_groups {

constexpr double kSE2_threshold_ = 1e-7;

class SE2 {

public:

Eigen::Matrix3d data_;
Eigen::Map<Eigen::Matrix<double,2,1>> t_;  /** < The position */
Eigen::Ref<Eigen::Matrix2d> R_; /** < The rotation */

/**
 * Default constructor. Initializes group element to identity.
 */
SE2() : t_(data_.data()+6), R_(data_.block(0,0,2,2)), data_(Eigen::Matrix3d::Identity() ) {}


/**
 * Copy constructor.
 */ 
SE2(const SE2 & g) : t_(data_.data()+6), R_(data_.block(0,0,2,2)), data_(g.data_) {}

/**
* Initializes group element to the one given. If verify is true
* it will check that the input is an element of \f$SE(2)\f$
* @param[in] data The data pertaining to an element of \f$SE(2)\f$
* @param verify If true, the constructor will verify that the provided 
* element is a member of \f$SE(2)\f$
*/
SE2(const Eigen::Matrix3d & data, bool verify);

/**
* Initializes group element to the data of the one given. 
* @param[in] data  The data pertaining to an element of \f$SE(2)\f$
*/
SE2(const Eigen::Matrix3d & data) :data_(data): t_(data_.data()+6), R_(data_.block(0,0,2,2)),{}

/**
 * Returns the inverse of the element
 */ 
SE2 Inverse(){  
    Eigen::Matrix3d m;
    m.block(0,0,2,2) = R_.transpose();
    m.block(0,2,2,1) = -m.block(0,0,2,2)*t_;
    m.block(2,0,3,1) << 0,0,1;
    return SE2(data_.transpose());}

/**
 * Returns the identity element
 */ 
static SE2 Identity(){return SE2();}

/**
 * Returns the matrix adjoint map.
 * For this Lie group it is the identity map.
 */ 
Eigen::Matrix3d SE2::Adjoint(){
    Eigen::Matrix3d m = data_;
    m.block(0,2,2,1) << t_(1), -t_(0);
    return m;
}

/**
 * Computes the log of the element.
 */ 
Eigen::Matrix<double,3,1> Log() {return se2::Log(this->data_);}

/**
 * Performs the left group action on itself. i.e. this is on the left of
 * the bilinear operation.
 */ 
SE2 operator * (const SE2& g){ return SE2(data_ *g.data_);}


/**
 * Assignment Operator. Deep copy of the input parameter
 * @param g The element to be copied.
 */ 
void operator = (const SE2& g){ this->data_ = g.data_; }

/**
 * Performs the OPlus operation 
 * @param g_data The data belonging to the group element.
 * @param u_data The data belonging to the Cartesian space that is isomorphic to the Lie algebra.
 * @return The result of the BoxPlus operation.
 */ 
static Eigen::Matrix3d OPlus(const Eigen::Matrix3d& g_data, const Eigen::Matrix<double,3,1>& u_data)
{return g_data*se2::Exp(u_data);}

/**
 * Performs the OPlus operation 
 * @param u_data The data belonging to the Cartesian space that is isomorphic to the Lie algebra.
 * @return The result of the BoxPlus operation.
 */ 
Eigen::Matrix3d OPlus(const Eigen::Matrix<double,3,1>& u_data)
{return OPlus(this->data_,u_data);}

/**
 * Performs the OPlus operation and assigns the result to the group element
 * @param u_data The data belonging to the Cartesian space that is isomorphic to the Lie algebra.
 * @return The result of the BoxPlus operation.
 */ 
void OPlusEq(const Eigen::Matrix<double,3,1>& u_data)
{data_ = this->OPlus(u_data);}

/**
 * Performs the BoxPlus operation 
 * @param g_data The data belonging to the group element.
 * @param u_data The data belonging to an element of the Lie algebra.
 * @return The result of the BoxPlus operation.
 */ 
static Eigen::Matrix3d BoxPlus(const Eigen::Matrix3d& g_data, const Eigen::Matrix3d& u_data)
{return OPlus(g_data,se2::Vee(u_data));}

/**
 * Performs the BoxPlus operation 
 * @param u_data The data belonging to an element of the Lie algebra.
 * @return The result of the BoxPlus operation.
 */ 
Eigen::Matrix3d BoxPlus(const Eigen::Matrix3d& u_data)
{return this->OPlus(se2::Vee(u_data));}

/**
 * Performs the BoxPlus operation and assigns the result to the group element
 * @param u_data The data belonging to an element of the Lie algebra.
 * @return The result of the BoxPlus operation.
 */ 
void BoxPlusEq(const Eigen::Matrix3d& u_data)
{data_ = this->BoxPlus(u_data);}

/**
 * Performs the BoxPlus operation 
 * @param g An element of the group
 * @param u An element of the Lie algebra.
 * @return The result of the BoxPlus operation.
 */ 
static SE2 BoxPlus(const SE2& g, const se2& u)
{return SE2(OPlus(g.data_,u.data_));}

/**
 * Performs the BoxPlus operation 
 * @param u An element of the Lie algebra.
 * @return The result of the BoxPlus operation.
 */ 
SE2 BoxPlus(const se2& u)
{return SE2::BoxPlus(*this,u);}

/**
 * Performs the BoxPlus operation and assigns the result to the group element
 * @param u An element of the Lie algebra.
 * @return The result of the BoxPlus operation.
 */ 
void BoxPlusEq(const se2& u)
{data_ = (this->BoxPlus(u)).data_;}


/**
 * Performs the O-minus operation \f$ \log(g_1^-1*g_2) \f$
 * @param g_data1 The data of  \f$ g_1 \f$
 * @param g_data2 The data of \f$ g_2 \f$
 * @return The data of an element of the Cartesian space isomorphic to the Lie algebra
 */ 
static Eigen::Matrix<double,3,1> OMinus(const Eigen::Matrix3d& g1_data,const Eigen::Matrix3d& g2_data)
{return se2::Log(g1_data.transpose()*g2_data);}


/**
 * Performs the O-minus operation with this being \f$ g_1 \f$ in the equation \f$ \log(g_1^-1*g_2) \f$
 * @param g_data The data of  \f$ g_2 \f$
 * @return The data of an element of the Cartesian space isomorphic to the Lie algebra
 */ 
Eigen::Matrix<double,3,1> OMinus(const Eigen::Matrix3d& g_data)
{return SE2::OMinus(data_,g_data);}

/**
 * Performs the O-minus operation \f$ \log(g_1^-1*g_2) \f$
 * @param g_data1 The data of  \f$ g_1 \f$
 * @param g_data2 The data of \f$ g_2 \f$
 * @return The data of an element of the Lie algebra
 */ 
static Eigen::Matrix3d BoxMinus(const Eigen::Matrix3d& g1_data,const Eigen::Matrix3d& g2_data)
{return se2::Wedge(SE2::OMinus(g1_data, g2_data));}


/**
 * Performs the O-minus operation with this being \f$ g_1 \f$ in the equation \f$ \log(g_1^-1*g_2) \f$
 * @param g_data The data of  \f$ g_2 \f$
 * @return The data of an element of the Lie algebra
 */ 
Eigen::Matrix3d BoxMinus(const Eigen::Matrix3d& g_data)
{return BoxMinus(this->data_,g_data);}

/**
 * Performs the O-minus operation with this being \f$ g_1 \f$ in the equation \f$ \log(g_1^-1*g_2) \f$
 * @param g An element of the group
 * @return An element of the Lie algebra
 */ 
se2 BoxMinus(const SE2& g){return se2(se2::Vee(this->BoxMinus(g.data_)));}

/**
 * Prints the content of the data
 */ 
void Print() {
    std::cout << data_ << std::endl; 
}

/**
 * Verifies that the data of an element properly corresponds to the set. 
 */ 
static bool isElement(const Eigen::Matrix3d& data);

};

}


#endif // _LIEGROUPS_INCLUDE_LIEGROUPS_SE2_