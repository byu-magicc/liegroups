#ifndef _LIEGROUPS_INCLUDE_LIEGROUPS_SO3_
#define _LIEGROUPS_INCLUDE_LIEGROUPS_SO3_

#include <Eigen/Dense>
#include <iostream>
#include <lie_algebras/so3.h>
#include "lie_groups/group_base.h"

namespace lie_groups {

constexpr double kSO3_threshold_ = 1e-7;

class SO3 : public GroupBase<SO3,so3, Eigen::Matrix3d, Eigen::Matrix<double,3,1>> {

public:

Eigen::Matrix3d data_;
static constexpr unsigned int dim_ = 3;
static constexpr unsigned int size1_ = 3;
static constexpr unsigned int size2_ = 3;
typedef GroupBase<SO3,so3, Eigen::Matrix3d, Eigen::Matrix<double,3,1>> Base;
typedef so3 algebra;
using Base::BoxPlus;
using Base::BoxMinus;

/**
 * Default constructor. Initializes group element to identity.
 */
SO3() : data_(Eigen::Matrix3d::Identity()){}


/**
 * Copy constructor.
 */ 
SO3(const SO3 & g) : data_(g.data_) {}

/**
 * copy assignment
 */ 
void operator = (const SO3& g){ this->data_ = g.data_; }

/**
 * Move constructor.
 */ 
SO3(const SO3 && g) : data_(g.data_) {}

/**
 * Move assignment
 */ 
void operator = (const SO3&& g){ this->data_ = g.data_; }

/**
* Initializes group element to the one given. If verify is true
* it will check that the input is an element of \f$SO(3)\f$
* @param[in] data The data pertaining to an element of \f$SO(3)\f$
* @param verify If true, the constructor will verify that the provided 
* element is a member of \f$SO(3)\f$
*/
SO3(const Eigen::Matrix3d & data, bool verify);

/**
* Initializes group element to the data of the one given. 
* @param[in] data  The data pertaining to an element of \f$SO(3)\f$
*/
SO3(const Eigen::Matrix3d & data) :data_(data){}

/*
 * Returns the inverse of the element
 */ 
SO3 Inverse(){ return SO3::Inverse(this->data_);}

/*
 * Returns the inverse of the data of an element
 */ 
static Eigen::Matrix3d Inverse(const Eigen::Matrix3d& data){  return data.transpose();}

/**
 * Returns the identity element
 */ 
static SO3 Identity(){return SO3();}

/**
 * Returns the matrix adjoint map.
 * For this Lie group it is the identity map.
 */ 
Eigen::Matrix3d Adjoint(){return data_;}

/**
 * Performs the left group action on itself. i.e. this is on the left of
 * the bilinear operation.
 */ 
SO3 operator * (const SO3& g){ return SO3(data_ *g.data_);}

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
static SO3 BoxPlus(const SO3& g, const so3& u)
{return SO3(OPlus(g.data_,u.data_));}


/**
 * Performs the BoxPlus operation 
 * @param u An element of the Lie algebra.
 * @return The result of the BoxPlus operation.
 */ 
SO3 BoxPlus(const so3& u)
{return BoxPlus(*this,u);}

/**
 * Performs the O-minus operation with this being \f$ g_1 \f$ in the equation \f$ \log(g_2^-1*g_1) \f$
 * @param g An element of the group
 * @return An element of the Lie algebra
 */ 
so3 BoxMinus(const SO3& g){ return so3( so3::Vee(BoxMinus(g.data_)));}


/**
 * Verifies that the data of an element properly corresponds to the set. 
 */ 
static bool isElement(const Eigen::Matrix3d& data);

};

}


#endif // _LIEGROUPS_INCLUDE_LIEGROUPS_SO3_