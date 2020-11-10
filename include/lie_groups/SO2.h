#ifndef _LIEGROUPS_INCLUDE_LIEGROUPS_SO2_
#define _LIEGROUPS_INCLUDE_LIEGROUPS_SO2_

#include <Eigen/Dense>
#include <utility>
#include <string>
#include <iostream>
#include <lie_algebras/so2.h>
#include "lie_groups/group_base.h"

namespace lie_groups {

constexpr double kSO2_threshold_ = 1e-7;

class SO2 : public GroupBase<SO2,so2, Eigen::Matrix2d, Eigen::Matrix<double,1,1>>{

public:

Eigen::Matrix2d data_;

static constexpr unsigned int dim_ = 1;
static constexpr unsigned int size1_ = 2;
static constexpr unsigned int size2_ = 2;
typedef so2 algebra;
typedef GroupBase<SO2,so2, Eigen::Matrix2d, Eigen::Matrix<double,1,1>> Base;
using Base::BoxPlus;
using Base::BoxMinus;

/**
 * Default constructor. Initializes group element to identity.
 */
SO2() : data_(Eigen::Matrix<double,2,2>::Identity()){}


/**
 * Copy constructor.
 */ 
SO2(const SO2 & g) : data_(g.data_) {}

/**
 * Copy assignment
 */ 
void operator = (const SO2& g){ this->data_ = g.data_; }

/**
 * Move constructor.
 */ 
SO2(const SO2 && g) : data_(g.data_) {}

/**
 * Move assignment
 */ 
void operator = (const SO2&& g){ this->data_ = g.data_; }

/**
* Initializes group element to the one given. If verify is true
* it will check that the input is an element of \f$SO(2)\f$
* @param[in] data The data pertaining to an element of \f$SO(2)\f$
* @param verify If true, the constructor will verify that the provided 
* element is a member of \f$SO(2)\f$
*/
SO2(const Eigen::Matrix2d & data, bool verify);

/**
* Initializes group element to the data of the one given. 
* @param[in] data  The data pertaining to an element of \f$SO(2)\f$
*/
SO2(const Eigen::Matrix2d & data) :data_(data){}

/*
 * Returns the inverse of the element
 */ 
SO2 Inverse(){ return SO2::Inverse(this->data_);}

/*
 * Returns the inverse of the data of an element
 */ 
static Eigen::Matrix2d Inverse(const Eigen::Matrix2d& data){  return data.transpose();}


/**
 * Returns the identity element
 */ 
static SO2 Identity(){return SO2();}

/**
 * Returns the matrix adjoint map.
 * For this Lie group it is the identity map.
 */ 
Eigen::Matrix<double,1,1> Adjoint(){return Eigen::Matrix<double,1,1>::Identity();}

/**
 * Performs the left group action on itself. i.e. this is on the left of
 * the bilinear operation.
 */ 
SO2 operator * (const SO2& g){ return SO2(data_ *g.data_);}

/**
 * Performs the group operation between the data of two elements
 */ 
static Eigen::Matrix2d Mult(const Eigen::Matrix2d& data1, const Eigen::Matrix2d& data2 ){
    return data1*data2;
}

/**
 * Performs the BoxPlus operation 
 * @param g An element of the group
 * @param u An element of the Lie algebra.
 * @return The result of the BoxPlus operation.typename
 */ 
static SO2 BoxPlus(const SO2& g, const so2& u)
{return SO2(OPlus(g.data_,u.data_));}


/**
 * Performs the BoxPlus operation 
 * @param u An element of the Lie algebra.
 * @return The result of the BoxPlus operation.
 */ 
SO2 BoxPlus(const so2& u)
{return BoxPlus(*this,u);}

/**
 * Performs the O-minus operation with this being \f$ g_1 \f$ in the equation \f$ \log(g_2^{-1}*g_1) \f$
 * @param g An element of the group
 * @return An element of the Lie algebra
 */ 
so2 BoxMinus(const SO2& g){ return so2( so2::Vee(BoxMinus(g.data_)));}

/**
 * Prints the content of the data
 */ 
void Print() {
    std::cout << data_ << std::endl; 
}

/**
 * Verifies that the data of an element properly corresponds to the set. 
 */ 
static bool isElement(const Eigen::Matrix2d& data);

};

}


#endif // _LIEGROUPS_INCLUDE_LIEGROUPS_SO2_