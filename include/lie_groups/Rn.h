#ifndef _LIEGROUPS_INCLUDE_LIEGROUPS_RN_
#define _LIEGROUPS_INCLUDE_LIEGROUPS_RN_

#include <Eigen/Dense>
#include <iostream>
#include <lie_algebras/rn.h>
#include "lie_groups/group_base.h"

namespace lie_groups {

constexpr double kRn_threshold_ = 1e-7;

template <int N>
class Rn : public GroupBase<Rn<N>, rn<N>, Eigen::Matrix<double,N,1>,Eigen::Matrix<double,N,1>>{

public:

Eigen::Matrix<double,N,1> data_;

static constexpr unsigned int dim_ = N;
static constexpr unsigned int size1_ = N;
static constexpr unsigned int size2_ = 1;
typedef GroupBase<Rn<N>, rn<N>, Eigen::Matrix<double,N,1>,Eigen::Matrix<double,N,1>> Base; 
typedef rn<N> algebra;
typedef Eigen::Matrix<double,N,1> Mat_G;
typedef Eigen::Matrix<double,N,1> Mat_A;
typedef Eigen::Matrix<double,N,1> Mat_C;
using Base::BoxPlus;
using Base::BoxMinus;



/**
 * Default constructor. Initializes group element to identity.
 */
Rn() : data_(Mat_G::Zero()) {}


/**
 * Copy constructor.
 */ 
Rn(const Rn & g) : data_(g.data_) {}

/**
 * Copy assignment
 */ 
void operator = (const Rn& g){ this->data_ = g.data_; }

/**
 * Move constructor.
 */ 
Rn(const Rn && g) : data_(g.data_) {}

/**
 * Move assignment
 */ 
void operator = (const Rn&& g){ this->data_ = g.data_; }

/**
* Initializes group element to the one given. If verify is true
* it will check that the input is an element of \f$\mathbb{R}^n\f$
* @param[in] data The data pertaining to an element of \f$\mathbb{R}^n\f$
* @param verify If true, the constructor will verify that the provided 
* element is a member of \f$\mathbb{R}^n\f$
*/
Rn(const Mat_G& data, bool verify) : data_(data) {}

/**
* Initializes group element to the data of the one given. 
* @param[in] data  The data pertaining to an element of \f$\mathbb{R}^n\f$
*/
Rn(const Mat_G& data) :data_(data){}

/*
 * Returns the inverse of the element
 */ 
Rn Inverse(){ return Rn(-this->data_);}

/*
 * Returns the inverse of the data of an element
 */ 
static Mat_G Inverse(const Mat_G& data){  return -data;}


/**
 * Returns the identity element
 */ 
static Rn Identity(){return Rn();}

/**
 * Returns the matrix adjoint map.
 * For this Lie group it is the identity map.
 */ 
Eigen::Matrix<double,N,N> Adjoint(){return Eigen::Matrix<double,N,N>::Identity();}

/**
 * Computes the log of the element.
 */ 
Mat_C Log(){return data_;}

/**
 * Performs the left group action on itself. i.e. this is on the left of
 * the bilinear operation.
 */ 
Rn operator * (const Rn& g){ return Rn(data_  + g.data_);}

/**
 * Performs the group operation between the data of two elements
 */ 
static Mat_G Mult(const Mat_G& data1, const Mat_G& data2 ){
    return data1 + data2;
}

/**
 * Performs the BoxPlus operation 
 * @param g An element of the group
 * @param u An element of the Lie algebra.
 * @return The result of the BoxPlus operation.typename
 */ 
static Rn<N>  BoxPlus(const Rn<N> & g, const rn<N>& u)
{return Rn<N> (OPlus(g.data_,u.data_));}


/**
 * Performs the BoxPlus operation 
 * @param u An element of the Lie algebra.
 * @return The result of the BoxPlus operation.
 */ 
Rn<N> BoxPlus(const rn<N>& u)
{return BoxPlus(*this,u);}

/**
 * Performs the O-minus operation with this being \f$ g_1 \f$ in the equation \f$ \log(g_2^-1*g_1) \f$
 * @param g An element of the group
 * @return An element of the Lie algebra
 */ 
rn<N> BoxMinus(const Rn<N>& g){ return rn<N>( rn<N>::Vee(BoxMinus(g.data_)));}

/**
 * Verifies that the data of an element properly corresponds to the set. 
 */ 
static bool isElement(const Mat_G& data) {return true;}

};

}


#endif // _LIEGROUPS_INCLUDE_LIEGROUPS_RN_