#ifndef _LIEGROUPS_INCLUDE_LIEGROUPS_SO2_
#define _LIEGROUPS_INCLUDE_LIEGROUPS_SO2_

#include <Eigen/Dense>
#include <utility>
#include <string>
#include <iostream>
#include <lie_algebras/so2.h>

namespace lie_groups {

constexpr double kSO2_threshold_ = 1e-7;

class SO2 {

public:

Eigen::Matrix2d data_;

static const unsigned int dim_ = 1;
static const unsigned int size1_ = 2;
static const unsigned int size2_ = 2;
// std::pair <std::string,double> product2;

/**
 * Default constructor. Initializes group element to identity.
 */
SO2() : data_(Eigen::Matrix<double,2,2>::Identity()){}


/**
 * Copy constructor.
 */ 
SO2(const SO2 & g) : data_(g.data_) {}

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
 * Computes the log of the element.
 */ 
Eigen::Matrix<double,1,1> Log();

/**
 * Performs the left group action on itself. i.e. this is on the left of
 * the bilinear operation.
 */ 
SO2 operator * (const SO2& g){ return SO2(data_ *g.data_);}


/**
 * Assignment Operator. Deep copy of the input parameter
 * @param g The element to be copied.
 */ 
void operator = (const SO2& g){ this->data_ = g.data_; }

/**
 * Performs the OPlus operation 
 * @param g_data The data belonging to the group element.
 * @param u_data The data belonging to the Cartesian space that is isomorphic to the Lie algebra.
 * @return The result of the BoxPlus operation.
 */ 
static Eigen::Matrix2d OPlus(const Eigen::Matrix2d& g_data, const Eigen::Matrix<double,1,1>& u_data)
{return g_data*so2::Exp(u_data);}

/**
 * Performs the OPlus operation 
 * @param u_data The data belonging to the Cartesian space that is isomorphic to the Lie algebra.
 * @return The result of the BoxPlus operation.
 */ 
Eigen::Matrix2d OPlus(const Eigen::Matrix<double,1,1>& u_data)
{return OPlus(this->data_,u_data);}

/**
 * Performs the OPlus operation and assigns the result to the group element
 * @param u_data The data belonging to the Cartesian space that is isomorphic to the Lie algebra.
 * @return The result of the BoxPlus operation.
 */ 
void OPlusEq(const Eigen::Matrix<double,1,1>& u_data)
{data_ = this->OPlus(u_data);}

/**
 * Performs the BoxPlus operation 
 * @param g_data The data belonging to the group element.
 * @param u_data The data belonging to an element of the Lie algebra.
 * @return The result of the BoxPlus operation.
 */ 
static Eigen::Matrix2d BoxPlus(const Eigen::Matrix2d& g_data, const Eigen::Matrix2d& u_data)
{return OPlus(g_data,so2::Vee(u_data));}

/**
 * Performs the BoxPlus operation 
 * @param u_data The data belonging to an element of the Lie algebra.
 * @return The result of the BoxPlus operation.
 */ 
Eigen::Matrix2d BoxPlus(const Eigen::Matrix2d& u_data)
{return this->OPlus(so2::Vee(u_data));}

/**
 * Performs the BoxPlus operation and assigns the result to the group element
 * @param u_data The data belonging to an element of the Lie algebra.
 * @return The result of the BoxPlus operation.
 */ 
void BoxPlusEq(const Eigen::Matrix2d& u_data)
{data_ = this->BoxPlus(u_data);}

/**
 * Performs the BoxPlus operation 
 * @param g An element of the group
 * @param u An element of the Lie algebra.
 * @return The result of the BoxPlus operation.
 */ 
static SO2 BoxPlus(const SO2& g, const so2& u)
{return SO2(OPlus(g.data_,u.data_));}

/**
 * Performs the BoxPlus operation 
 * @param u An element of the Lie algebra.
 * @return The result of the BoxPlus operation.
 */ 
SO2 BoxPlus(const so2& u)
{return SO2::BoxPlus(*this,u);}

/**
 * Performs the BoxPlus operation and assigns the result to the group element
 * @param u An element of the Lie algebra.
 * @return The result of the BoxPlus operation.
 */ 
void BoxPlusEq(const so2& u)
{data_ = (this->BoxPlus(u)).data_;}


/**
 * Performs the O-minus operation \f$ \log(g_1^-1*g_2) \f$
 * @param g_data1 The data of  \f$ g_1 \f$
 * @param g_data2 The data of \f$ g_2 \f$
 * @return The data of an element of the Cartesian space isomorphic to the Lie algebra
 */ 
static Eigen::Matrix<double,1,1> OMinus(const Eigen::Matrix2d& g1_data,const Eigen::Matrix2d& g2_data)
{return so2::Log(SO2::Inverse(g1_data)*g2_data);}


/**
 * Performs the O-minus operation with this being \f$ g_1 \f$ in the equation \f$ \log(g_1^-1*g_2) \f$
 * @param g_data The data of  \f$ g_2 \f$
 * @return The data of an element of the Cartesian space isomorphic to the Lie algebra
 */ 
Eigen::Matrix<double,1,1> OMinus(const Eigen::Matrix2d& g_data)
{return SO2::OMinus(data_,g_data);}

/**
 * Performs the O-minus operation \f$ \log(g_1^-1*g_2) \f$
 * @param g_data1 The data of  \f$ g_1 \f$
 * @param g_data2 The data of \f$ g_2 \f$
 * @return The data of an element of the Lie algebra
 */ 
static Eigen::Matrix2d BoxMinus(const Eigen::Matrix2d& g1_data,const Eigen::Matrix2d& g2_data)
{return so2::Wedge(SO2::OMinus(g1_data, g2_data));}


/**
 * Performs the O-minus operation with this being \f$ g_1 \f$ in the equation \f$ \log(g_1^-1*g_2) \f$
 * @param g_data The data of  \f$ g_2 \f$
 * @return The data of an element of the Lie algebra
 */ 
Eigen::Matrix2d BoxMinus(const Eigen::Matrix2d& g_data)
{return BoxMinus(this->data_,g_data);}

/**
 * Performs the O-minus operation with this being \f$ g_1 \f$ in the equation \f$ \log(g_1^-1*g_2) \f$
 * @param g An element of the group
 * @return An element of the Lie algebra
 */ 
so2 BoxMinus(const SO2& g){return so2(so2::Vee(this->BoxMinus(g.data_)));}

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