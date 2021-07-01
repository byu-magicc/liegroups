#ifndef _LIEGROUPS_INCLUDE_LIEGROUPS_RN_
#define _LIEGROUPS_INCLUDE_LIEGROUPS_RN_

#include <Eigen/Dense>
#include <iostream>

#include "lie_groups/lie_algebras/rn.h"
#include "lie_groups/lie_groups/group_base.h"

namespace lie_groups {

constexpr double kRn_threshold_ = 1e-7;

template <typename tDataType=double, int tNumDimensions=2, int tNumTangentSpaces=1>
class Rn : public GroupBase<Rn<tDataType, tNumDimensions,tNumTangentSpaces>, rn<tDataType,tNumDimensions,tNumTangentSpaces>, Eigen::Matrix<tDataType,tNumDimensions,1>,Eigen::Matrix<tDataType,tNumDimensions*tNumTangentSpaces,1>, tDataType>{

static_assert(tNumTangentSpaces > 0, "lie_groups::Rn the number of tangent spaces must be greater than 0.");


public:


static constexpr unsigned int num_tangent_spaces_ = tNumTangentSpaces;
static constexpr unsigned int dim_ = tNumDimensions;
static constexpr unsigned int size1_ = tNumDimensions;
static constexpr unsigned int size2_ = 1;
typedef GroupBase<Rn<tDataType, tNumDimensions,tNumTangentSpaces>, rn<tDataType,tNumDimensions,tNumTangentSpaces>, Eigen::Matrix<tDataType,tNumDimensions,1>,Eigen::Matrix<tDataType,tNumDimensions*tNumTangentSpaces,1>, tDataType> Base; 
typedef rn<tDataType,tNumDimensions,tNumTangentSpaces> Algebra;
typedef Abelian GroupType;
typedef Eigen::Matrix<tDataType,tNumDimensions,1> MatNd;
using Base::BoxPlus;
using Base::BoxMinus;

Eigen::Matrix<tDataType,tNumDimensions,1> data_;

/**
 * Default constructor. Initializes group element to identity.
 */
Rn() : data_(MatNd::Zero()) {}


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
Rn(const MatNd& data, bool verify) : data_(data) {}

/**
* Initializes group element to the data of the one given. 
* @param[in] data  The data pertaining to an element of \f$\mathbb{R}^n\f$
*/
Rn(const MatNd& data) :data_(data){}

/*
 * Returns the inverse of the element
 */ 
Rn Inverse() const { return Rn(-this->data_);}

/*
 * Returns the inverse of the data of an element
 */ 
static MatNd Inverse(const MatNd& data){  return -data;}


/**
 * Returns the identity element
 */ 
static Rn Identity(){return Rn();}

/**
 * Returns the matrix adjoint map.
 * For this Lie group it is the identity map.
 */ 
Eigen::Matrix<tDataType,dim_,dim_> Adjoint(){return Eigen::Matrix<tDataType,dim_,dim_>::Identity();}

/**
 * Computes the log of the element.
 */ 
MatNd Log() const {return data_;}

/**
 * Performs the left group action on itself. i.e. this is on the left of
 * the bilinear operation.
 */ 
Rn operator * (const Rn& g){ return Rn(data_  + g.data_);}

/**
 * Performs the group operation between the data of two elements
 */ 
static MatNd Mult(const MatNd& data1, const MatNd& data2 ){
    return data1 + data2;
}

/**
 * Performs the BoxPlus operation 
 * @param g An element of the group
 * @param u An element of the Lie algebra.
 * @return The result of the BoxPlus operation.typename
 */ 
static Rn  BoxPlus(const Rn & g, const Algebra& u)
{return Rn (OPlus(g.data_,u.data_));}


/**
 * Performs the BoxPlus operation 
 * @param u An element of the Lie algebra.
 * @return The result of the BoxPlus operation.
 */ 
Rn BoxPlus(const Algebra& u) const
{return BoxPlus(*this,u);}

/**
 * Performs the O-minus operation with this being \f$ g_1 \f$ in the equation \f$ \log(g_2^-1*g_1) \f$
 * @param g An element of the group
 * @return An element of the Lie algebra
 */ 
Algebra BoxMinus(const Rn& g) const { return Algebra( Algebra::Vee(BoxMinus(g.data_)));}

/**
 * Verifies that the data of an element properly corresponds to the set. 
 */ 
static bool isElement(const MatNd& data) {return true;}

};

} // namespace lie_groups


#endif // _LIEGROUPS_INCLUDE_LIEGROUPS_RN_