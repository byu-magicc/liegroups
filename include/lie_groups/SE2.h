#ifndef _LIEGROUPS_INCLUDE_LIEGROUPS_SE2_
#define _LIEGROUPS_INCLUDE_LIEGROUPS_SE2_

#include <Eigen/Dense>
#include <iostream>
#include <lie_algebras/se2.h>
#include <lie_algebras/so2.h>
#include "lie_groups/group_base.h"

namespace lie_groups {

constexpr double kSE2_threshold_ = 1e-7;

// The GroupBase template needs the Group, Algebra, Group Data, and Cartesian Data
class SE2 : public GroupBase<SE2,se2, Eigen::Matrix3d, Eigen::Matrix<double,3,1>>{

public:

Eigen::Matrix3d data_;
Eigen::Map<Eigen::Matrix<double,2,1>> t_;  /** < The position */
Eigen::Ref<Eigen::Matrix2d> R_; /** < The rotation */

// These are evaluated at compile time
static constexpr unsigned int dim_ = 3;
static constexpr unsigned int dim_pos_=2;
static constexpr unsigned int dim_rot_=1; /** < The dimension of the rotation */
static constexpr unsigned int size1_ = 3;
static constexpr unsigned int size2_ = 3;
typedef GroupBase<SE2,se2, Eigen::Matrix3d, Eigen::Matrix<double,3,1>> Base; 
typedef se2 algebra;
typedef so2 rot_algebra;
using Base::BoxPlus;
using Base::BoxMinus;


/**
 * Default constructor. Initializes group element to identity.
 */
SE2() : t_(data_.data()+6), R_(data_.block(0,0,2,2)), data_(Eigen::Matrix3d::Identity() ) {}


/**
 * Copy constructor.
 */ 
SE2(const SE2 & g) : t_(data_.data()+6), R_(data_.block(0,0,2,2)), data_(g.data_) {}

/**
 * Copy assignment
 */ 
void operator = (const SE2& g){ this->data_ = g.data_; }

/**
 * Move constructor.
 */ 
SE2(const SE2 && g) : t_(data_.data()+6), R_(data_.block(0,0,2,2)), data_(g.data_) {}

/**
 * Move assignment
 */ 
void operator = (const SE2&& g){ this->data_ = g.data_; }

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
SE2(const Eigen::Matrix3d & data) :data_(data), t_(data_.data()+6), R_(data_.block(0,0,2,2)){}

/*
 * Returns the inverse of the element
 */ 
SE2 Inverse(){ return SE2::Inverse(this->data_);}


/*
 * Returns the inverse of the data of an element
 */ 
static Eigen::Matrix3d Inverse(const Eigen::Matrix3d& data){  
    return data.inverse();}

/**
 * Returns the identity element
 */ 
static SE2 Identity(){return SE2();}

/**
 * Returns the matrix adjoint map.
 * For this Lie group it is the identity map.
 */ 
Eigen::Matrix3d Adjoint(){
    Eigen::Matrix3d m = data_;
    m.block(0,2,2,1) << t_(1), -t_(0);
    return m;
}


/**
 * Performs the left group action 
 */
SE2 operator * (const SE2& g){ return SE2(data_ *g.data_) ;}

/**
 * Performs the group operation between the data of two elements
 */ 
static Eigen::Matrix3d Mult(const Eigen::Matrix3d& data1, const Eigen::Matrix3d& data2 ){
    return data1*data2;
}

/**
 * Performs the BoxPlus operation 
 * @param g An element of the group
 * @param u An element of the Lie algebra.
 * @return The result of the BoxPlus operation.typename
 */ 
static SE2 BoxPlus(const SE2& g, const se2& u)
{return SE2(OPlus(g.data_,u.data_));}


/**
 * Performs the BoxPlus operation 
 * @param u An element of the Lie algebra.
 * @return The result of the BoxPlus operation.
 */ 
SE2 BoxPlus(const se2& u)
{return BoxPlus(*this,u);}

/**
 * Performs the O-minus operation with this being \f$ g_1 \f$ in the equation \f$ \log(g_2^-1*g_1) \f$
 * @param g An element of the group
 * @return An element of the Lie algebra
 */ 
se2 BoxMinus(const SE2& g){ return se2( se2::Vee(BoxMinus(g.data_)));}


/**
 * Verifies that the data of an element properly corresponds to the set. 
 */ 
static bool isElement(const Eigen::Matrix3d& data);

};

}


#endif // _LIEGROUPS_INCLUDE_LIEGROUPS_SE2_