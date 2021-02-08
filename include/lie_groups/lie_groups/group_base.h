#ifndef _LIEGROUPS_INCLUDE_LIEGROUPS_GROUPBASE_
#define _LIEGROUPS_INCLUDE_LIEGROUPS_GROUPBASE_


#include <Eigen/Dense>
#include <iostream>
#include "lie_algebras/se2.h"





namespace lie_groups {

// These classes are used to give addition information about a state
struct Abelian {};
struct NonAbelian{};



template<typename Group, typename Algebra, typename Mat_G, typename Mat_C, typename tDataType = double>
class GroupBase{

private:
GroupBase()=default;
~GroupBase()=default;
friend Group;

public:

typedef Mat_G Mat_A;                                               // Lie algebra matrix

/**
 * Computes the log of the element.
 */ 
Mat_C Log() const {return Algebra::Log(static_cast<const Group*>(this)->data_);}
// Mat_C Log() {return se2::Log(static_cast<Group*>(this)->data_);}


// static Mat_G Mult(const Mat_G& data1, const Mat_G& data2 ){
//     return Group::Mult(data1,data2);
// }

/**
 * Creates a random element
 */ 
static Mat_G Random() {return Algebra::Exp(Mat_C::Random());}

/**
 * Performs the OPlus operation 
 * @param g_data The data belonging to the group element.
 * @param u_data The data belonging to the Cartesian space that is isomorphic to the Lie algebra.
 * @return The result of the BoxPlus operation.
 */ 
static Mat_G OPlus(const Mat_G& g_data, const Mat_C& u_data)
{return Group::Mult(g_data,Algebra::Exp(u_data));}


/**
 * Performs the OPlus operation 
 * @param u_data The data belonging to the Cartesian space that is isomorphic to the Lie algebra.
 * @return The result of the BoxPlus operation.
 */ 
Mat_G OPlus(const Mat_C& u_data) const
{return OPlus(static_cast<const Group*>(this)->data_ ,u_data);}

/**
 * Performs the OPlus operation and assigns the result to the group element
 * @param u_data The data belonging to the Cartesian space that is isomorphic to the Lie algebra.
 * @return The result of the BoxPlus operation.
 */ 
void OPlusEq(const Mat_C& u_data)
{static_cast<Group*>(this)->data_ = this->OPlus(u_data);}

/**
 * Performs the BoxPlus operation 
 * @param g_data The data belonging to the group element.
 * @param u_data The data belonging to an element of the Lie algebra.
 * @return The result of the BoxPlus operation.
 */ 
static Mat_G BoxPlus(const Mat_G& g_data, const Mat_A& u_data)
{return OPlus(g_data,Algebra::Vee(u_data));}

/**
 * Performs the BoxPlus operation 
 * @param u_data The data belonging to an element of the Lie algebra.
 * @return The result of the BoxPlus operation.
 */ 
Mat_G BoxPlus(const Mat_A& u_data) const
{return this->OPlus(Algebra::Vee(u_data));}

/**
 * Performs the BoxPlus operation and assigns the result to the group element
 * @param u_data The data belonging to an element of the Lie algebra.
 * @return The result of the BoxPlus operation.
 */ 
void BoxPlusEq(const Mat_A& u_data)
{static_cast<Group* >(this)->data_ = this->BoxPlus(u_data);}

/**
 * Performs the BoxPlus operation and assigns the result to the group element
 * @param u An element of the Lie algebra.
 * @return The result of the BoxPlus operation.
 */ 
void BoxPlusEq(const Algebra& u)
{static_cast<Group* >(this)->data_ = (this->OPlus(u.data_));}


/**
 * Performs the O-minus operation \f$ \log(g_2^-1*g_1) \f$
 * @param g_data1 The data of  \f$ g_1 \f$
 * @param g_data2 The data of \f$ g_2 \f$
 * @return The data of an element of the Cartesian space isomorphic to the Lie algebra
 */ 
static Mat_C OMinus(const Mat_G& g1_data,const Mat_G& g2_data) 
{return Algebra::Log(Group::Mult(Group::Inverse(g2_data),g1_data));}


/**
 * Performs the O-minus operation with this being \f$ g_1 \f$ in the equation \f$ \log(g_2^-1*g_1) \f$
 * @param g_data The data of  \f$ g_1 \f$
 * @return The data of an element of the Cartesian space isomorphic to the Lie algebra
 */ 
Mat_C OMinus(const Mat_G& g_data) const
{return OMinus(static_cast<const Group*>(this)->data_,g_data);}

/**
 * Performs the O-minus operation \f$ \log(g_2^-1*g_1) \f$
 * @param g_data1 The data of  \f$ g_1 \f$
 * @param g_data2 The data of \f$ g_2 \f$
 * @return The data of an element of the Lie algebra
 */ 
static Mat_A BoxMinus(const Mat_G& g1_data,const Mat_G& g2_data)
{return Algebra::Wedge(OMinus(g1_data, g2_data));}

/**
 * Performs the O-minus operation with this being \f$ g_1 \f$ in the equation \f$ \log(g_2^-1*g_1) \f$
 * @param g_data The data of  \f$ g_1 \f$
 * @return The data of an element of the Lie algebra
 */ 
Mat_A BoxMinus(const Mat_G& g_data) const
{return BoxMinus(static_cast<const Group*>(this)->data_,g_data);}


/**
 * Prints the content of the data
 */ 
void Print() {
    std::cout << static_cast<Group*>(this)->data_ << std::endl; 
}


};


}


#endif // _LIEGROUPS_INCLUDE_LIEGROUPS_GROUPBASE_