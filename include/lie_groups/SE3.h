#ifndef _LIEGROUPS_INCLUDE_LIEGROUPS_SE3_
#define _LIEGROUPS_INCLUDE_LIEGROUPS_SE3_

#include <Eigen/Dense>
#include <iostream>
#include <lie_algebras/se3.h>
#include "lie_groups/group_base.h"

namespace lie_groups {

constexpr double kSE3_threshold_ = 1e-7;

class SE3 : public GroupBase<SE3,se3, Eigen::Matrix4d, Eigen::Matrix<double,6,1>> {

public:

Eigen::Matrix<double,4,4> data_;
Eigen::Map<Eigen::Matrix<double,3,1>> t_;  /** < The position */
Eigen::Ref<Eigen::Matrix3d> R_; /** < The rotation */
static constexpr unsigned int dim_ = 6;
static constexpr unsigned int size1_ = 4;
static constexpr unsigned int size2_ = 4;
typedef se3 algebra;
typedef GroupBase<SE3,se3, Eigen::Matrix4d, Eigen::Matrix<double,6,1>> Base; 
using Base::BoxPlus;
using Base::BoxMinus;


/**
 * Default constructor. Initializes group element to identity.
 */
SE3() : t_(data_.data()+12), R_(data_.block(0,0,3,3)), data_(Eigen::Matrix<double,4,4>::Identity() ) {}


/**
 * Copy constructor.
 */ 
SE3(const SE3 & g) : t_(data_.data()+12), R_(data_.block(0,0,3,3)), data_(g.data_) {}

/**
 * Copy assignment
 */ 
void operator = (const SE3& g){ this->data_ = g.data_; }

/**
 * Move constructor.
 */ 
SE3(const SE3 && g) : t_(data_.data()+12), R_(data_.block(0,0,3,3)), data_(g.data_) {}

/**
 * Move assignment
 */ 
void operator = (const SE3&& g){ this->data_ = g.data_; }

/**
* Initializes group element to the one given. If verify is true
* it will check that the input is an element of \f$SE(3)\f$
* @param[in] data The data pertaining to an element of \f$SE(3)\f$
* @param verify If true, the constructor will verify that the provided 
* element is a member of \f$SE(3)\f$
*/
SE3(const Eigen::Matrix<double,4,4> & data, bool verify);

/**
* Initializes group element to the data of the one given. 
* @param[in] data  The data pertaining to an element of \f$SE(3)\f$
*/
SE3(const Eigen::Matrix<double,4,4> & data) :data_(data), t_(data_.data()+12), R_(data_.block(0,0,3,3)){}

/*
 * Returns the inverse of the element
 */ 
SE3 Inverse(){ return SE3::Inverse(this->data_);}


/*
 * Returns the inverse of the data of an element
 */ 

static Eigen::Matrix<double,4,4> Inverse(const Eigen::Matrix<double,4,4>& data){  
    return data.inverse();}

/**
 * Returns the identity element
 */ 
static SE3 Identity(){return SE3();}

/**
 * Returns the matrix adjoint map.
 * For this Lie group it is the identity map.
 */ 
Eigen::Matrix<double,6,6> Adjoint(){
    Eigen::Matrix<double,6,6> m;
    m.block(0,0,3,3) = R_;
    m.block(3,3,3,3) = R_;
    m.block(0,3,3,3) = se3::SSM(t_)*R_;
    m.block(3,0,3,3).setZero(); 
    return m;
}

/**
 * Performs the left group action 
 */
SE3 operator * (const SE3& g){ return SE3(data_ *g.data_) ;}

/**
 * Performs the group operation between the data of two elements
 */ 
static Eigen::Matrix4d Mult(const Eigen::Matrix4d& data1, const Eigen::Matrix4d& data2 ){
    return data1*data2;
}

/**
 * Performs the BoxPlus operation 
 * @param g An element of the group
 * @param u An element of the Lie algebra.
 * @return The result of the BoxPlus operation.typename
 */ 
static SE3 BoxPlus(const SE3& g, const se3& u)
{return SE3(OPlus(g.data_,u.data_));}


/**
 * Performs the BoxPlus operation 
 * @param u An element of the Lie algebra.
 * @return The result of the BoxPlus operation.
 */ 
SE3 BoxPlus(const se3& u)
{return BoxPlus(*this,u);}

/**
 * Performs the O-minus operation with this being \f$ g_1 \f$ in the equation \f$ \log(g_2^-1*g_1) \f$
 * @param g An element of the group
 * @return An element of the Lie algebra
 */ 
se3 BoxMinus(const SE3& g){ return se3( se3::Vee(BoxMinus(g.data_)));}



/**
 * Prints the content of the data
 */ 
void Print() {
    std::cout << data_ << std::endl; 
}

/**
 * Verifies that the data of an element properly corresponds to the set. 
 */ 
static bool isElement(const Eigen::Matrix<double,4,4>& data);

};

}


#endif // _LIEGROUPS_INCLUDE_LIEGROUPS_SE3_