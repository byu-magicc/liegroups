#ifndef _LIEGROUPS_INCLUDE_LIEGROUPS_SE2_
#define _LIEGROUPS_INCLUDE_LIEGROUPS_SE2_

#include <Eigen/Dense>
#include <iostream>

#include "lie_groups/lie_algebras/se2.h"
#include "lie_groups/lie_algebras/so2.h"
#include "lie_groups/lie_groups/group_base.h"

namespace lie_groups {

constexpr double kSE2_threshold_ = 1e-6;

// The GroupBase template needs the Group, Algebra, Group Data, and Cartesian Data
template <typename tDataType=double, int tNumDimensions=3, int tNumTangentSpaces=1>
class SE2 : public GroupBase<SE2<tDataType,tNumDimensions,tNumTangentSpaces>,se2<tDataType,tNumDimensions,tNumTangentSpaces>, Eigen::Matrix<tDataType,3,3>, Eigen::Matrix<tDataType,3,3>, Eigen::Matrix<tDataType,3,1>, tDataType>{

static_assert(tNumTangentSpaces == 1, "lie_groups::SE2 the number of tangent spaces must be 1.");


public:

typedef Eigen::Matrix<tDataType,2,1> Vec2d;
typedef Eigen::Matrix<tDataType,3,1> Vec3d;
typedef Eigen::Matrix<tDataType,2,2> Mat2d;
typedef Eigen::Matrix<tDataType,3,3> Mat3d;

// These are evaluated at compile time
static constexpr unsigned int dim_ = 3;
static constexpr unsigned int dim_pos_=2;
static constexpr unsigned int dim_rot_=1; /** < The dimension of the rotation */
static constexpr unsigned int size1_ = 3;
static constexpr unsigned int size2_ = 3;
typedef GroupBase<SE2<tDataType,tNumDimensions,tNumTangentSpaces>,se2<tDataType,tNumDimensions,tNumTangentSpaces>, Eigen::Matrix<tDataType,3,3>, Eigen::Matrix<tDataType,3,3>, Eigen::Matrix<tDataType,3,1>, tDataType> Base; 
typedef se2<tDataType,tNumDimensions,tNumTangentSpaces> Algebra;
typedef NonAbelian GroupType;

typedef so2<tDataType> RotAlgebra;
using Base::BoxPlus;
using Base::BoxMinus;


Mat3d data_;
Eigen::Map<Vec2d> t_;  /** < The position */
Eigen::Ref<Mat2d> R_; /** < The rotation */

/**
 * Default constructor. Initializes group element to identity.
 */
SE2() : t_(data_.data()+6), R_(data_.block(0,0,2,2)), data_(Mat3d::Identity() ) {}


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
SE2(const Mat3d & data, bool verify);

/**
* Initializes group element to the data of the one given. 
* @param[in] data  The data pertaining to an element of \f$SE(2)\f$
*/
SE2(const Mat3d & data) :data_(data), t_(data_.data()+6), R_(data_.block(0,0,2,2)){}

/*
 * Returns the inverse of the element
 */ 
SE2 Inverse() const { return SE2::Inverse(this->data_);}


/*
 * Returns the inverse of the data of an element
 */ 
static Mat3d Inverse(const Mat3d& data){  
    return data.inverse();}

/**
 * Returns the identity element
 */ 
static SE2 Identity(){return SE2();}

/**
 * Returns the matrix adjoint map.
 * For this Lie group it is the identity map.
 */ 
Mat3d Adjoint(){
    Mat3d m = data_;
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
static Mat3d Mult(const Mat3d& data1, const Mat3d& data2 ){
    return data1*data2;
}

/**
 * Performs the BoxPlus operation 
 * @param g An element of the group
 * @param u An element of the Lie algebra.
 * @return The result of the BoxPlus operation.typename
 */ 
static SE2 BoxPlus(const SE2& g, const Algebra& u)
{return SE2(OPlus(g.data_,u.data_));}


/**
 * Performs the BoxPlus operation 
 * @param u An element of the Lie algebra.
 * @return The result of the BoxPlus operation.
 */ 
SE2 BoxPlus(const Algebra& u) const
{return BoxPlus(*this,u);}

/**
 * Performs the O-minus operation with this being \f$ g_1 \f$ in the equation \f$ \log(g_2^-1*g_1) \f$
 * @param g An element of the group
 * @return An element of the Lie algebra
 */ 
Algebra BoxMinus(const SE2& g) const { return Algebra( Algebra::Vee(BoxMinus(g.data_)));}


/**
 * Verifies that the data of an element properly corresponds to the set. 
 */ 
static bool isElement(const Mat3d& data);

};



///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                    Definitions
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template <typename tDataType, int tNumDimensions, int tNumTangentSpaces>
SE2<tDataType,tNumDimensions,tNumTangentSpaces>::SE2(const Eigen::Matrix<tDataType,3,3> & data, bool verify) : t_(data_.data()+6), R_(data_.block(0,0,2,2)) {

    // First verify that it is a proper group element.
    if (verify ) {
        if (SE2::isElement(data)) {
            data_ = data;
        } else {
            std::cerr << "SE2::Constructor not valid input setting to identity" << std::endl;
            data_.setIdentity();
        }
    } else {
        data_ = data;
    }
}

//----------------------------------------------------------
template <typename tDataType, int tNumDimensions, int tNumTangentSpaces>
bool SE2<tDataType,tNumDimensions,tNumTangentSpaces>::isElement(const Eigen::Matrix<tDataType,3,3>& data) {
    
    tDataType d = (data.block(0,0,2,2).transpose()*data.block(0,0,2,2)-Eigen::Matrix<tDataType,2,2>::Identity()).norm();
    
    return d <= kSE2_threshold_ && data(2,0) == 0 && data(2,1)==0 && data(2,2)==1;
}

} // namespace lie_groups


#endif // _LIEGROUPS_INCLUDE_LIEGROUPS_SE2_