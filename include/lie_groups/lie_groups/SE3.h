#ifndef _LIEGROUPS_INCLUDE_LIEGROUPS_SE3_
#define _LIEGROUPS_INCLUDE_LIEGROUPS_SE3_

#include <Eigen/Dense>
#include <iostream>

#include "lie_groups/lie_algebras/se3.h"
#include "lie_groups/lie_algebras/so3.h"
#include "lie_groups/lie_groups/group_base.h"

namespace lie_groups {

constexpr double kSE3_threshold_ = 1e-6;
template<typename tDataType=double, int tN=6>
class SE3 : public GroupBase<SE3<tDataType>,se3<tDataType>, Eigen::Matrix<tDataType,4,4>, Eigen::Matrix<tDataType,6,1>,tDataType> {

public:

typedef Eigen::Matrix<tDataType,3,1> Vec3d;
typedef Eigen::Matrix<tDataType,6,1> Vec6d;
typedef Eigen::Matrix<tDataType,3,3> Mat3d;
typedef Eigen::Matrix<tDataType,4,4> Mat4d;
typedef Eigen::Matrix<tDataType,6,6> Mat6d;

static constexpr unsigned int dim_ = 6;
static constexpr unsigned int dim_pos_ = 3; /** < The dimension of the position */
static constexpr unsigned int dim_rot_ = 3; /** < The dimension of the rotation */
static constexpr unsigned int size1_ = 4;
static constexpr unsigned int size2_ = 4;
typedef se3<tDataType> Algebra;
typedef NonAbelian GroupType;
typedef so3<tDataType> RotAlgebra;
typedef GroupBase<SE3<tDataType>,se3<tDataType>, Eigen::Matrix<tDataType,4,4>, Eigen::Matrix<tDataType,6,1>,tDataType> Base; 
using Base::BoxPlus;
using Base::BoxMinus;


Mat4d data_;
Eigen::Map<Vec3d> t_;  /** < The position */
Eigen::Ref<Mat3d> R_; /** < The rotation */


/**
 * Default constructor. Initializes group element to identity.
 */
SE3() : t_(data_.data()+12), R_(data_.block(0,0,3,3)), data_(Mat4d::Identity() ) {}


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
SE3(const Mat4d & data, bool verify);

/**
* Initializes group element to the data of the one given. 
* @param[in] data  The data pertaining to an element of \f$SE(3)\f$
*/
SE3(const Mat4d & data) :data_(data), t_(data_.data()+12), R_(data_.block(0,0,3,3)){}

/*
 * Returns the inverse of the element
 */ 
SE3 Inverse() const { return SE3::Inverse(this->data_);}


/*
 * Returns the inverse of the data of an element
 */ 

static Mat4d Inverse(const Mat4d& data){  
    return data.inverse();}

/**
 * Returns the identity element
 */ 
static SE3 Identity(){return SE3();}

/**
 * Returns the matrix adjoint map.
 * For this Lie group it is the identity map.
 */ 
Mat6d Adjoint(){
    Mat6d m;
    m.block(0,0,3,3) = R_;
    m.block(3,3,3,3) = R_;
    m.block(0,3,3,3) = se3<tDataType>::SSM(t_)*R_;
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
static Mat4d Mult(const Mat4d& data1, const Mat4d& data2 ){
    return data1*data2;
}

/**
 * Performs the BoxPlus operation 
 * @param g An element of the group
 * @param u An element of the Lie algebra.
 * @return The result of the BoxPlus operation.typename
 */ 
static SE3 BoxPlus(const SE3& g, const Algebra& u)
{return SE3(OPlus(g.data_,u.data_));}


/**
 * Performs the BoxPlus operation 
 * @param u An element of the Lie algebra.
 * @return The result of the BoxPlus operation.
 */ 
SE3 BoxPlus(const Algebra& u) const
{return BoxPlus(*this,u);}

/**
 * Performs the O-minus operation with this being \f$ g_1 \f$ in the equation \f$ \log(g_2^-1*g_1) \f$
 * @param g An element of the group
 * @return An element of the Lie algebra
 */ 
Algebra BoxMinus(const SE3& g) const { return Algebra( Algebra::Vee(BoxMinus(g.data_)));}



/**
 * Prints the content of the data
 */ 
void Print() {
    std::cout << data_ << std::endl; 
}

/**
 * Verifies that the data of an element properly corresponds to the set. 
 */ 
static bool isElement(const Mat4d& data);

};

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                    Definitions
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename tDataType, int tN>
SE3<tDataType, tN>::SE3(const Eigen::Matrix<tDataType,4,4> & data, bool verify) : t_(data_.data()+12), R_(data_.block(0,0,3,3)) {

    // First verify that it is a proper group element.
    if (verify ) {
        if (SE3::isElement(data)) {
            data_ = data;
        } else {
            std::cerr << "SE3::Constructor not valid input setting to identity" << std::endl;
            data_.setIdentity();
        }
    } else {
        data_ = data;
    }
}

//-------------------------------------------------------------------
template<typename tDataType, int tN>
bool SE3<tDataType, tN>::isElement(const Eigen::Matrix<tDataType,4,4>& data) {
    
    double d = (data.block(0,0,3,3).transpose()*data.block(0,0,3,3)-Mat3d::Identity()).norm();
    
    return d <= kSE3_threshold_ && data(3,0) == 0 && data(3,1)==0 && data(3,2)==0 && data(3,3)==1;
}


} // namespace lie_groups


#endif // _LIEGROUPS_INCLUDE_LIEGROUPS_SE3_