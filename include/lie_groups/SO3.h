#ifndef _LIEGROUPS_INCLUDE_LIEGROUPS_SO3_
#define _LIEGROUPS_INCLUDE_LIEGROUPS_SO3_

#include <Eigen/Dense>
#include <iostream>
#include <lie_algebras/so3.h>

namespace lie_groups {

constexpr double kSO3_threshold_ = 1e-7;

class SO3 {

public:

Eigen::Matrix3d data_;
static const unsigned int dim_ = 3;
static const unsigned int size1_ = 3;
static const unsigned int size2_ = 3;


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
 * Computes the log of the element.
 */ 
Eigen::Matrix<double,3,1> Log();

/**
 * Performs the left group action on itself. i.e. this is on the left of
 * the bilinear operation.
 */ 
SO3 operator * (const SO3& g){ return SO3(data_ *g.data_);}

/**
 * Performs the OPlus operation 
 * @param g_data The data belonging to the group element.
 * @param u_data The data belonging to the Cartesian space that is isomorphic to the Lie algebra.
 * @return The result of the BoxPlus operation.
 */ 
static Eigen::Matrix3d OPlus(const Eigen::Matrix3d& g_data, const Eigen::Matrix<double,3,1>& u_data)
{return g_data*so3::Exp(u_data);}

/**
 * Performs the OPlus operation 
 * @param u_data The data belonging to the Cartesian space that is isomorphic to the Lie algebra.
 * @return The result of the BoxPlus operation.
 */ 
Eigen::Matrix3d OPlus(const Eigen::Matrix<double,3,1>& u_data)
{return OPlus(this->data_,u_data);}

/**
 * Performs the OPlus operation and assigns the result to the group element
 * @param u_data The data belonging to the Cartesian space that is isomorphic to the Lie algebra.
 * @return The result of the BoxPlus operation.
 */ 
void OPlusEq(const Eigen::Matrix<double,3,1>& u_data)
{data_ = this->OPlus(u_data);}

/**
 * Performs the BoxPlus operation 
 * @param g_data The data belonging to the group element.
 * @param u_data The data belonging to an element of the Lie algebra.
 * @return The result of the BoxPlus operation.
 */ 
static Eigen::Matrix3d BoxPlus(const Eigen::Matrix3d& g_data, const Eigen::Matrix3d& u_data)
{return OPlus(g_data,so3::Vee(u_data));}

/**
 * Performs the BoxPlus operation 
 * @param u_data The data belonging to an element of the Lie algebra.
 * @return The result of the BoxPlus operation.
 */ 
Eigen::Matrix3d BoxPlus(const Eigen::Matrix3d& u_data)
{return this->OPlus(so3::Vee(u_data));}

/**
 * Performs the BoxPlus operation and assigns the result to the group element
 * @param u_data The data belonging to an element of the Lie algebra.
 * @return The result of the BoxPlus operation.
 */ 
void BoxPlusEq(const Eigen::Matrix3d& u_data)
{data_ = this->BoxPlus(u_data);}

/**
 * Performs the BoxPlus operation 
 * @param g An element of the group
 * @param u An element of the Lie algebra.
 * @return The result of the BoxPlus operation.
 */ 
static SO3 BoxPlus(const SO3& g, const so3& u)
{return SO3(OPlus(g.data_,u.data_));}

/**
 * Performs the BoxPlus operation 
 * @param u An element of the Lie algebra.
 * @return The result of the BoxPlus operation.
 */ 
SO3 BoxPlus(const so3& u)
{return SO3::BoxPlus(*this,u);}

/**
 * Performs the BoxPlus operation and assigns the result to the group element
 * @param u An element of the Lie algebra.
 * @return The result of the BoxPlus operation.
 */ 
void BoxPlusEq(const so3& u)
{data_ = (this->BoxPlus(u)).data_;}


/**
 * Performs the O-minus operation \f$ \log(g_1^-1*g_2) \f$
 * @param g_data1 The data of  \f$ g_1 \f$
 * @param g_data2 The data of \f$ g_2 \f$
 * @return The data of an element of the Cartesian space isomorphic to the Lie algebra
 */ 
static Eigen::Matrix<double,3,1> OMinus(const Eigen::Matrix3d& g1_data,const Eigen::Matrix3d& g2_data)
{return so3::Log(SO3::Inverse(g1_data)*g2_data);}


/**
 * Performs the O-minus operation with this being \f$ g_1 \f$ in the equation \f$ \log(g_1^-1*g_2) \f$
 * @param g_data The data of  \f$ g_2 \f$
 * @return The data of an element of the Cartesian space isomorphic to the Lie algebra
 */ 
Eigen::Matrix<double,3,1> OMinus(const Eigen::Matrix3d& g_data)
{return SO3::OMinus(data_,g_data);}

/**
 * Performs the O-minus operation \f$ \log(g_1^-1*g_2) \f$
 * @param g_data1 The data of  \f$ g_1 \f$
 * @param g_data2 The data of \f$ g_2 \f$
 * @return The data of an element of the Lie algebra
 */ 
static Eigen::Matrix3d BoxMinus(const Eigen::Matrix3d& g1_data,const Eigen::Matrix3d& g2_data)
{return so3::Wedge(SO3::OMinus(g1_data, g2_data));}


/**
 * Performs the O-minus operation with this being \f$ g_1 \f$ in the equation \f$ \log(g_1^-1*g_2) \f$
 * @param g_data The data of  \f$ g_2 \f$
 * @return The data of an element of the Lie algebra
 */ 
Eigen::Matrix3d BoxMinus(const Eigen::Matrix3d& g_data)
{return BoxMinus(this->data_,g_data);}

/**
 * Performs the O-minus operation with this being \f$ g_1 \f$ in the equation \f$ \log(g_1^-1*g_2) \f$
 * @param g An element of the group
 * @return An element of the Lie algebra
 */ 
so3 BoxMinus(const SO3& g){return so3(so3::Vee(this->BoxMinus(g.data_)));}

/**
 * Prints the content of the data
 */ 
void Print() {
    std::cout << data_ << std::endl; 
}

/**
 * Verifies that the data of an element properly corresponds to the set. 
 */ 
static bool isElement(const Eigen::Matrix3d& data);

};

}


#endif // _LIEGROUPS_INCLUDE_LIEGROUPS_SO3_