#ifndef _LIEGROUPS_INCLUDE_LIEGROUPS_STATE_
#define _LIEGROUPS_INCLUDE_LIEGROUPS_STATE_

#include <Eigen/Dense>
#include <iostream>

// Lie algebras
#include "lie_algebras/rn.h"
#include "lie_algebras/so2.h"
#include "lie_algebras/so3.h"
#include "lie_algebras/se2.h"
#include "lie_algebras/se3.h"

// Lie groups
#include "lie_groups/Rn.h"
#include "lie_groups/SO2.h"
#include "lie_groups/SO3.h"
#include "lie_groups/SE2.h"
#include "lie_groups/SE3.h"

namespace lie_groups {

template <class G> 
class State {

public:
typedef typename G::algebra U;
G g_;   /** < The pose of the object.*/
U u_;   /** < The twist (velocity) of the object.*/

static constexpr unsigned int dim_ = G::dim_ + U::dim_;

typedef G g_type_; /** < The group type .*/
typedef U u_type_; /** < The algebra type .*/
typedef Eigen::Matrix<double,G::size1_, G::size2_> Mat_G;     /**< The group data type. */
typedef Eigen::Matrix<double,U::size1_, U::size2_> Mat_C;     /**< The Cartesian space data type. */
typedef Mat_G Mat_A;                                          /**< The Lie algebra data type. */
typedef typename G::GroupType StateType;
typedef Eigen::Matrix<double,2*G::dim_, G::dim_> Mat_Adj; /**< The Adjoint data type. */
typedef Eigen::Matrix<double,2*U::size1_,1> Mat_SC;           /**< The State Cartesian space data type. */

/**
 * Default constructor. Initializes group element to identity.
 */
State()=default;

/**
 * Copy constructor.
 */ 
State(const State& s) : g_(s.g_), u_(s.u_) {};

/**
 * Copy assignment.
 */ 
void operator = (const State& s){g_ = s.g_; u_ = s.u_;};

/**
 * Move constructor.
 */ 
State(const State&& s) : g_(s.g_), u_(s.u_) {};

/**
 * Move assignment.
 */ 
void operator = (const State&& s){g_ = s.g_; u_ = s.u_;};

/**
 * Copy constructor using group and algebra elements.
 */ 
State(const G& g, const U& u) : g_(g), u_(u) {};

/**
 * Move constructor using group and algebra elements.
 */ 
State(const G&& g, const U&& u) : g_(g), u_(u) {};

/**
* Initializes state using the data provied. 
* @param g_data The data pertaining to the group.
* @param u_data The data pertaining to a element of the Cartesian space isomorphic Lie algebra
* are elements of the group and Lie algebra.
*/
State(bool verify, const Mat_G & g_data, const Mat_C & u_data) : g_(g_data,verify), u_(u_data) {}


/**
* Initializes state using the data provied. If verify is true
* it will check that the inputs are elements of the group and Lie algebra.
* @param g_data The data pertaining to the group.
* @param u_data The data pertaining to the Lie algebra
* @param verify If true, the constructor will verify that the provided data 
* are elements of the group and Lie algebra.
*/
State(const Mat_G & g_data, const Mat_A & u_data, bool verify) : g_(g_data,verify), u_(u_data,verify) {}

/*
 * Returns the inverse of the element
 */ 
State Inverse(){ return State(g_.Inverse(),U(-u_.data_));}
/**
 * Returns the identity element
 */ 
static State Identity(){return State();}

/**
 * Left group action multiplication. (this*s)
 */ 
State operator * (const State& s){
  return State(this->g_*s.g_,this->u_+s.u_);
}

/**
 * Return a random state element
 */
static State Random() {
  return State(false,G::Random(),Mat_C::Random());
} 

/**
 * Returns the state adjoint
 */ 
Mat_Adj Adjoint(){
  Mat_Adj tmp;
  tmp.block(0,0,G::dim_,G::dim_) = g_.Adjoint();
  tmp.block(G::dim_,0,G::dim_,G::dim_) = u_.Adjoint();
  return tmp;}

/**
 * Performs the O-minus operation \f$ \log(S_2^{-1}*S_1) i.e. S_1-S_2\f$
 * @param g1_data The group data of  \f$ g_1 \f$
 * @param g2_data The group data of \f$ g_2 \f$
 * @param u1_data The Cartesian data of \f$ u_1 \f$
 * @param u2_data The Cartesian data of \f$ u_2 \f$
 * @return The data of an element of the Cartesian space isomorphic to the Lie algebra
 */ 
static Mat_SC OMinus(const Mat_G& g1_data,const Mat_G& g2_data,const Mat_C & u1_data,const Mat_C & u2_data)
{ Mat_SC tmp;
  tmp.block(0,0,G::dim_,1) = G::OMinus(g1_data,g2_data);
  tmp.block(G::dim_,0,U::dim_,1) = u1_data - u2_data;
    return tmp;}

/**
 * Performs the O-minus operation \f$ \log(S_2^{-1}*S_1) i.e. S_1-S_2\f$
 * @param s1 The state  \f$ s_1 \f$
 * @param s2 The state  \f$ s_2 \f$
 * @return The data of an element of the Cartesian space isomorphic to the Lie algebra
 */ 
static Mat_SC OMinus(const State<G>& s1, const State<G>&s2)
{ return OMinus(s1.g_.data_,s2.g_.data_, s1.u_.data_,s2.u_.data_);}


/**
 * Performs the O-minus operation \f$ \log(S_2^{-1}*S_1) i.e. S_2-S_1\f$ with this being S_1
 * @param s2 The state  \f$ s_2 \f$
 * @return The data of an element of the Cartesian space isomorphic to the Lie algebra
 */ 
Mat_SC OMinus(const State<G>& s2) const {
  return OMinus(*this,s2);
}


/**
 * Performs the O-Plus operation \f$ \text{state} \exp{\text{cartesian}}\f$ 
 * @param state The state  
 * @param cartesian An in the Cartesian space
 * @return A state that is the result of the O-Plus operation
 */ 
static State OPlus(const State& state, Mat_SC cartesian) {
  State tmp;
  tmp.g_.data_ = state.g_.OPlus(cartesian.block(0,0,G::dim_,1));
  tmp.u_.data_ = state.u_.data_ + cartesian.block(G::dim_,0,G::dim_,1);
  return tmp;
}

/**
 * Performs the O-Plus operation \f$ \text{this} \exp{\text{cartesian}}\f$ 
 * @param cartesian An in the Cartesian space
 * @return A state that is the result of the O-Plus operation
 */ 
State OPlus( Mat_SC cartesian) const{
  return OPlus(*this,cartesian);
}

/**
 * Performs the O-Plus operation \f$ \text{this} \exp{\text{cartesian}}\f$ and sets the state to the result.
 * @param cartesian An in the Cartesian space
 */ 
void OPlusEQ( Mat_SC cartesian) {
  *this = OPlus(*this,cartesian);
}

// /**
//  * Performs the O-minus operation with this being \f$ g_1 \f$ in the equation \f$ \log(g_1^-1*g_2) \f$
//  * @param g_data The data of  \f$ g_2 \f$
//  * @return The data of an element of the Cartesian space isomorphic to the Lie algebra
//  */ 
// Eigen::Matrix<double,3,1> OMinus(const Eigen::Matrix3d& g_data)
// {return SE2::OMinus(data_,g_data);}

// /**
//  * Performs the O-minus operation \f$ \log(g_1^-1*g_2) \f$
//  * @param g_data1 The data of  \f$ g_1 \f$
//  * @param g_data2 The data of \f$ g_2 \f$
//  * @return The data of an element of the Lie algebra
//  */ 
// static Eigen::Matrix3d BoxMinus(const Eigen::Matrix3d& g1_data,const Eigen::Matrix3d& g2_data)
// {return se2::Wedge(SE2::OMinus(g1_data, g2_data));}


// /**
//  * Performs the O-minus operation with this being \f$ g_1 \f$ in the equation \f$ \log(g_1^-1*g_2) \f$
//  * @param g_data The data of  \f$ g_2 \f$
//  * @return The data of an element of the Lie algebra
//  */ 
// Eigen::Matrix3d BoxMinus(const Eigen::Matrix3d& g_data)
// {return BoxMinus(this->data_,g_data);}

// /**
//  * Performs the O-minus operation with this being \f$ g_1 \f$ in the equation \f$ \log(g_1^-1*g_2) \f$
//  * @param g An element of the group
//  * @return An element of the Lie algebra
//  */ 
// se2 BoxMinus(const SE2& g){return se2(se2::Vee(this->BoxMinus(g.data_)));}


};

typedef State<Rn<2>> R2_r2;
typedef State<Rn<3>> R3_r3;
typedef State<SO2> SO2_so2;
typedef State<SO3> SO3_so3;
typedef State<SE2> SE2_se2;
typedef State<SE3> SE3_se3;

}

#endif //_LIEGROUPS_INCLUDE_LIEGROUPS_STATE