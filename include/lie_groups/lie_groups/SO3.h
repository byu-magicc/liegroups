#ifndef _LIEGROUPS_INCLUDE_LIEGROUPS_SO3_
#define _LIEGROUPS_INCLUDE_LIEGROUPS_SO3_

#include <Eigen/Dense>
#include <iostream>

#include "lie_groups/lie_algebras/so3.h"
#include "lie_groups/lie_groups/group_base.h"

namespace lie_groups {

constexpr double kSO3_threshold_ = 1e-6;

template<typename tDataType=double, int tN=3>
class SO3 : public GroupBase<SO3<tDataType>,so3<tDataType>, Eigen::Matrix<tDataType,3,3>, Eigen::Matrix<tDataType,3,1>,tDataType> {

public:

typedef Eigen::Matrix<tDataType,3,1> Vec3d;
typedef Eigen::Matrix<tDataType,3,3> Mat3d;

static constexpr unsigned int dim_ = 3;
static constexpr unsigned int size1_ = 3;
static constexpr unsigned int size2_ = 3;
typedef GroupBase<SO3<tDataType>,so3<tDataType>, Eigen::Matrix<tDataType,3,3>, Eigen::Matrix<tDataType,3,1>,tDataType> Base;
typedef so3<tDataType> Algebra;
typedef NonAbelian GroupType;
using Base::BoxPlus;
using Base::BoxMinus;

 Mat3d data_;

/**
 * Default constructor. Initializes group element to identity.
 */
SO3() : data_( Mat3d::Identity()){}


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
SO3(const  Mat3d & data, bool verify);

/**
* Initializes group element to the data of the one given. 
* @param[in] data  The data pertaining to an element of \f$SO(3)\f$
*/
SO3(const  Mat3d & data) :data_(data){}

/*
 * Returns the inverse of the element
 */ 
SO3 Inverse() const { return SO3::Inverse(this->data_);}

/*
 * Returns the inverse of the data of an element
 */ 
static  Mat3d Inverse(const  Mat3d& data){  return data.transpose();}

/**
 * Returns the identity element
 */ 
static SO3 Identity(){return SO3();}

/**
 * Returns the matrix adjoint map.
 * For this Lie group it is the identity map.
 */ 
 Mat3d Adjoint(){return data_;}

/**
 * Performs the left group action on itself. i.e. this is on the left of
 * the bilinear operation.
 */ 
SO3 operator * (const SO3& g){ return SO3(data_ *g.data_);}

/**
 * Performs the group operation between the data of two elements
 */ 
static  Mat3d Mult(const  Mat3d& data1, const  Mat3d& data2 ){
    return data1*data2;
}

/**
 * Performs the BoxPlus operation 
 * @param g An element of the group
 * @param u An element of the Lie algebra.
 * @return The result of the BoxPlus operation.typename
 */ 
static SO3 BoxPlus(const SO3& g, const Algebra& u)
{return SO3(OPlus(g.data_,u.data_));}


/**
 * Performs the BoxPlus operation 
 * @param u An element of the Lie algebra.
 * @return The result of the BoxPlus operation.
 */ 
SO3 BoxPlus(const Algebra& u) const
{return BoxPlus(*this,u);}

/**
 * Performs the O-minus operation with this being \f$ g_1 \f$ in the equation \f$ \log(g_2^-1*g_1) \f$
 * @param g An element of the group
 * @return An element of the Lie algebra
 */ 
Algebra BoxMinus(const SO3& g) const { return Algebra( Algebra::Vee(BoxMinus(g.data_)));}


/**
 * Verifies that the data of an element properly corresponds to the set. 
 */ 
static bool isElement(const  Mat3d& data);

};

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                    Definitions
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename tDataType, int tN>
SO3<tDataType,tN>::SO3(const  Eigen::Matrix<tDataType,3,3> & data, bool verify) {

    // First verify that it is a proper group element.
    if (verify ) {
        if (SO3::isElement(data)) {
            data_ = data;
        } else {
            std::cerr << "SO3::Constructor not valid input setting to identity" << std::endl;
            data_.setIdentity();
        }
    } else {
        data_ = data;
    }
}


//----------------------------------------------------------
template<typename tDataType, int tN>
bool SO3<tDataType,tN>::isElement(const  Eigen::Matrix<tDataType,3,3>& data) {
    return (data.transpose()*data -  Mat3d::Identity()).norm() < kSO3_threshold_;
}


} // namespace lie_groups


#endif // _LIEGROUPS_INCLUDE_LIEGROUPS_SO3_