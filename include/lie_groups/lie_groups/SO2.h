#ifndef _LIEGROUPS_INCLUDE_LIEGROUPS_SO2_
#define _LIEGROUPS_INCLUDE_LIEGROUPS_SO2_

#include <Eigen/Dense>
#include <utility>
#include <string>
#include <iostream>

#include "lie_groups/lie_algebras/so2.h"
#include "lie_groups/lie_groups/group_base.h"

namespace lie_groups {

constexpr double kSO2_threshold_ = 1e-7;

template <typename tDataType=double, int tNumDimensions=1, int tNumTangentSpaces=1>
class SO2 : public GroupBase<SO2<tDataType,tNumDimensions,tNumTangentSpaces>,so2<tDataType,tNumDimensions,tNumTangentSpaces>, Eigen::Matrix<tDataType,2,2>, Eigen::Matrix<tDataType,2,2>, Eigen::Matrix<tDataType,1,1>>{

static_assert(tNumTangentSpaces == 1, "lie_groups::SO2 the number of tangent spaces must be greater than 0.");


public:


typedef Eigen::Matrix<tDataType,1,1> Mat1d;
typedef Eigen::Matrix<tDataType,2,2> Mat2d;

static constexpr unsigned int dim_ = tNumDimensions;
static constexpr unsigned int size1_ = 2;
static constexpr unsigned int size2_ = 2;
typedef so2<tDataType,tNumDimensions,tNumTangentSpaces> Algebra;
typedef Abelian GroupType;
typedef GroupBase<SO2<tDataType,tNumDimensions,tNumTangentSpaces>,so2<tDataType,tNumDimensions,tNumTangentSpaces>, Eigen::Matrix<tDataType,2,2>, Eigen::Matrix<tDataType,2,2>, Eigen::Matrix<tDataType,1,1>> Base;
using Base::BoxPlus;
using Base::BoxMinus;

Mat2d data_;

/**
 * Default constructor. Initializes group element to identity.
 */
SO2() : data_(Mat2d::Identity()){}


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
SO2(const Mat2d & data, bool verify);

/**
* Initializes group element to the data of the one given. 
* @param[in] data  The data pertaining to an element of \f$SO(2)\f$
*/
SO2(const Mat2d & data) :data_(data){}

/*
 * Returns the inverse of the element
 */ 
SO2 Inverse() const { return SO2<tDataType>::Inverse(this->data_);}

/*
 * Returns the inverse of the data of an element
 */ 
static Mat2d Inverse(const Mat2d& data){  return data.transpose();}


/**
 * Returns the identity element
 */ 
static SO2 Identity(){return SO2();}

/**
 * Returns the matrix adjoint map.
 * For this Lie group it is the identity map.
 */ 
Mat1d Adjoint(){return Mat1d::Identity();}

/**
 * Performs the left group action on itself. i.e. this is on the left of
 * the bilinear operation.
 */ 
SO2 operator * (const SO2& g){ return SO2(data_ *g.data_);}

/**
 * Performs the group operation between the data of two elements
 */ 
static Mat2d Mult(const Mat2d& data1, const Mat2d& data2 ){
    return data1*data2;
}

/**
 * Performs the BoxPlus operation 
 * @param g An element of the group
 * @param u An element of the Lie algebra.
 * @return The result of the BoxPlus operation.typename
 */ 
static SO2 BoxPlus(const SO2& g, const Algebra& u)
{return SO2(OPlus(g.data_,u.data_));}


/**
 * Performs the BoxPlus operation 
 * @param u An element of the Lie algebra.
 * @return The result of the BoxPlus operation.
 */ 
SO2 BoxPlus(const Algebra& u) const
{return BoxPlus(*this,u);}

/**
 * Performs the O-minus operation with this being \f$ g_1 \f$ in the equation \f$ \log(g_2^{-1}*g_1) \f$
 * @param g An element of the group
 * @return An element of the Lie algebra
 */ 
so2<tDataType> BoxMinus(const SO2& g) const { return Algebra( Algebra::Vee(BoxMinus(g.data_)));}

/**
 * Prints the content of the data
 */ 
void Print() {
    std::cout << data_ << std::endl; 
}

/**
 * Verifies that the data of an element properly corresponds to the set. 
 */ 
static bool isElement(const Mat2d& data);

};

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                    Definitions
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template <typename tDataType, int tNumDimensions, int tNumTangentSpaces>
SO2<tDataType,tNumDimensions,tNumTangentSpaces>::SO2(const Eigen::Matrix<tDataType,2,2> & data, bool verify) {

    // First verify that it is a proper group element.
    if (verify ) {
        if (SO2<tDataType>::isElement(data)) {
            data_ = data;
        } else {
            std::cerr << "SO2<tDataType>::Constructor not valid input setting to identity" << std::endl;
            data_.setIdentity();
        }
    } else {
        data_ = data;
    }
}

//----------------------------------------------------------
template <typename tDataType, int tNumDimensions, int tNumTangentSpaces>
bool SO2<tDataType,tNumDimensions,tNumTangentSpaces>::isElement(const Eigen::Matrix<tDataType,2,2>& data) {
    return (data.transpose()*data - Mat2d::Identity()).norm() < kSO2_threshold_;
}


} // namespace lie_groups


#endif // _LIEGROUPS_INCLUDE_LIEGROUPS_SO2_