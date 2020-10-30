#ifndef _LIEGROUPS_INCLUDE_LIEGROUPS_RN_
#define _LIEGROUPS_INCLUDE_LIEGROUPS_RN_

#include <Eigen/Dense>
#include <iostream>
#include <lie_algebras/rn.h>

namespace lie_groups {

constexpr double kRn_threshold_ = 1e-7;

template <int N>
class Rn {

public:

Eigen::Matrix<double,N,1> data_;

static const unsigned int dim_ = N;
static const unsigned int size1_ = N;
static const unsigned int size2_ = 1;


/**
 * Default constructor. Initializes group element to identity.
 */
Rn() : data_(Eigen::Matrix<double,N,1>::Zero()) {}


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
Rn(const Eigen::Matrix<double,N,1> & data, bool verify) : data_(data) {}

/**
* Initializes group element to the data of the one given. 
* @param[in] data  The data pertaining to an element of \f$\mathbb{R}^n\f$
*/
Rn(const Eigen::Matrix<double,N,1> & data) :data_(data){}

/*
 * Returns the inverse of the element
 */ 
Rn Inverse(){ return Rn(-this->data_);}

/*
 * Returns the inverse of the data of an element
 */ 
static Eigen::Matrix<double,N,1> Inverse(const Eigen::Matrix<double,N,1>& data){  return -data;}


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
Eigen::Matrix<double,N,1> Log(){return data_;}

/**
 * Performs the left group action on itself. i.e. this is on the left of
 * the bilinear operation.
 */ 
Rn operator * (const Rn& g){ return Rn(data_  + g.data_);}

/**
 * Performs the OPlus operation 
 * @param g_data The data belonging to the group element.
 * @param u_data The data belonging to the Cartesian space that is isomorphic to the Lie algebra.
 * @return The result of the BoxPlus operation.
 */ 
static Eigen::Matrix<double,N,1> OPlus(const Eigen::Matrix<double,N,1>& g_data, const Eigen::Matrix<double,N,1>& u_data)
{return g_data + u_data;}

/**
 * Performs the OPlus operation 
 * @param u_data The data belonging to the Cartesian space that is isomorphic to the Lie algebra.
 * @return The result of the BoxPlus operation.
 */ 
Eigen::Matrix<double,N,1> OPlus(const Eigen::Matrix<double,N,1>& u_data)
{return OPlus(this->data_,u_data);}

/**
 * Performs the OPlus operation and assigns the result to the group element
 * @param u_data The data belonging to the Cartesian space that is isomorphic to the Lie algebra.
 * @return The result of the BoxPlus operation.
 */ 
void OPlusEq(const Eigen::Matrix<double,N,1>& u_data)
{data_ = this->OPlus(u_data);}

/**
 * Performs the BoxPlus operation 
 * @param g_data The data belonging to the group element.
 * @param u_data The data belonging to an element of the Lie algebra.
 * @return The result of the BoxPlus operation.
 */ 
static Eigen::Matrix<double,N,1> BoxPlus(const Eigen::Matrix<double,N,1>& g_data, const Eigen::Matrix<double,N,1>& u_data)
{return OPlus(g_data,u_data);}

/**
 * Performs the BoxPlus operation 
 * @param u_data The data belonging to an element of the Lie algebra.
 * @return The result of the BoxPlus operation.
 */ 
Eigen::Matrix<double,N,1> BoxPlus(const Eigen::Matrix<double,N,1>& u_data)
{return this->OPlus(u_data);}

/**
 * Performs the BoxPlus operation and assigns the result to the group element
 * @param u_data The data belonging to an element of the Lie algebra.
 * @return The result of the BoxPlus operation.
 */ 
void BoxPlusEq(const Eigen::Matrix<double,N,1>& u_data)
{data_ = this->BoxPlus(u_data);}

/**
 * Performs the BoxPlus operation 
 * @param g An element of the group
 * @param u An element of the Lie algebra.
 * @return The result of the BoxPlus operation.
 */ 
static Rn BoxPlus(const Rn& g, const rn<N>& u)
{return Rn(OPlus(g.data_,u.data_));}

/**
 * Performs the BoxPlus operation 
 * @param u An element of the Lie algebra.
 * @return The result of the BoxPlus operation.
 */ 
Rn BoxPlus(const rn<N>& u)
{return Rn::BoxPlus(*this,u);}

/**
 * Performs the BoxPlus operation and assigns the result to the group element
 * @param u An element of the Lie algebra.
 * @return The result of the BoxPlus operation.
 */ 
void BoxPlusEq(const rn<N>& u)
{data_ = (this->BoxPlus(u)).data_;}


/**
 * Performs the O-minus operation \f$ \log(g_1^-1*g_2) \f$
 * @param g_data1 The data of  \f$ g_1 \f$
 * @param g_data2 The data of \f$ g_2 \f$
 * @return The data of an element of the Cartesian space isomorphic to the Lie algebra
 */ 
static Eigen::Matrix<double,N,1> OMinus(const Eigen::Matrix<double,N,1>& g1_data,const Eigen::Matrix<double,N,1>& g2_data)
{return g2_data - g1_data;}


/**
 * Performs the O-minus operation with this being \f$ g_1 \f$ in the equation \f$ \log(g_1^-1*g_2) \f$
 * @param g_data The data of  \f$ g_2 \f$
 * @return The data of an element of the Cartesian space isomorphic to the Lie algebra
 */ 
Eigen::Matrix<double,N,1> OMinus(const Eigen::Matrix<double,N,1>& g_data)
{return Rn::OMinus(data_,g_data);}

/**
 * Performs the O-minus operation \f$ \log(g_1^-1*g_2) \f$
 * @param g_data1 The data of  \f$ g_1 \f$
 * @param g_data2 The data of \f$ g_2 \f$
 * @return The data of an element of the Lie algebra
 */ 
static Eigen::Matrix<double,N,1> BoxMinus(const Eigen::Matrix<double,N,1>& g1_data,const Eigen::Matrix<double,N,1>& g2_data)
{return Rn::OMinus(g1_data, g2_data);}


/**
 * Performs the O-minus operation with this being \f$ g_1 \f$ in the equation \f$ \log(g_1^-1*g_2) \f$
 * @param g_data The data of  \f$ g_2 \f$
 * @return The data of an element of the Lie algebra
 */ 
Eigen::Matrix<double,N,1> BoxMinus(const Eigen::Matrix<double,N,1>& g_data)
{return BoxMinus(this->data_,g_data);}

/**
 * Performs the O-minus operation with this being \f$ g_1 \f$ in the equation \f$ \log(g_1^-1*g_2) \f$
 * @param g An element of the group
 * @return An element of the Lie algebra
 */ 
rn<N> BoxMinus(const Rn& g){return rn<N>(this->BoxMinus(g.data_));}

/**
 * Prints the content of the data
 */ 
void Print() {
    std::cout << data_ << std::endl; 
}

/**
 * Verifies that the data of an element properly corresponds to the set. 
 */ 
static bool isElement(const Eigen::Matrix<double,N,1>& data) {return true;}

};

}


#endif // _LIEGROUPS_INCLUDE_LIEGROUPS_RN_