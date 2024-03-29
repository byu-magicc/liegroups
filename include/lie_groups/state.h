#ifndef _LIEGROUPS_INCLUDE_LIEGROUPS_STATE_
#define _LIEGROUPS_INCLUDE_LIEGROUPS_STATE_

#include <Eigen/Dense>
#include <iostream>

// Lie algebras
#include "lie_groups/lie_algebras/rn.h"
#include "lie_groups/lie_algebras/so2.h"
#include "lie_groups/lie_algebras/so3.h"
#include "lie_groups/lie_algebras/se2.h"
#include "lie_groups/lie_algebras/se3.h"

// Lie groups
#include "lie_groups/lie_groups/Rn.h"
#include "lie_groups/lie_groups/SO2.h"
#include "lie_groups/lie_groups/SO3.h"
#include "lie_groups/lie_groups/SE2.h"
#include "lie_groups/lie_groups/SE3.h"

namespace lie_groups {

template <template<typename , int, int > class tG, typename tDataType = double,int tGroupDim =2, int tNumTangentSpaces = 1> 
class State {

public:
static constexpr int N = tGroupDim;                              /**< The dimensions of the group. */
static constexpr int NumTangentSpaces = tNumTangentSpaces;       /**< The number of tangent spaces to consider. For example, if the number is 1, then only velocity is considered. If the number is 2, then velocity and acceleration are considered. 
                                                                      Not every group can support more than one tangent space; currently only RN can.*/

typedef tDataType DataType;
typedef tG<tDataType,tGroupDim,tNumTangentSpaces> G;
typedef typename G::Algebra U;
typedef G Group;
typedef U Algebra;
typedef G g_type_; /** < The group type .*/
typedef U u_type_; /** < The algebra type .*/
typedef Eigen::Matrix<tDataType,G::size1_, G::size2_> Mat_G;      /**< The group data type. */
typedef Eigen::Matrix<tDataType,U::size1_, U::size2_> Mat_C;      /**< The Cartesian space data type. */
typedef typename G::Base::Mat_A Mat_A;                            /**< The Lie algebra data type. */
typedef typename G::GroupType StateType;
// typedef Eigen::Matrix<tDataType,G::dim_ + U::total_num_dim_, G::dim_> Mat_Adj;    /**< The Adjoint data type. */
static constexpr unsigned int dim_ = G::dim_ + U::total_num_dim_;
typedef Eigen::Matrix<tDataType,dim_,1> Vec_SC;                   /**< The State Cartesian space data type. */
typedef Eigen::Matrix<tDataType,dim_,dim_> Mat_SC;                /**< The State Cartesian space matrix data type. */

template<typename T>
using StateTemplate = State<tG, T, tGroupDim, tNumTangentSpaces>;

template< typename T1, int T2, int T3>
using GroupTemplate = tG<T1,T2,T3>;



G g_;   /** < The pose of the object.*/
U u_;   /** < The twist (velocity) of the object.*/

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
static State Random(const DataType scalar = static_cast<DataType>(1.0)) {
  return State(false,G::Random(scalar),Mat_C::Random()*scalar);
} 

// /**
//  * Returns the state adjoint
//  */ 
// Mat_Adj Adjoint(){
//   Mat_Adj tmp;
//   tmp.block(0,0,G::dim_,G::dim_) = g_.Adjoint();
//   tmp.block(G::dim_,0,U::total_num_dim_,G::dim_) = u_.Adjoint();
//   return tmp;}

/**
 * Performs the O-minus operation \f$ \log(S_2^{-1}*S_1) i.e. S_1-S_2\f$
 * @param g1_data The group data of  \f$ g_1 \f$
 * @param g2_data The group data of \f$ g_2 \f$
 * @param u1_data The Cartesian data of \f$ u_1 \f$
 * @param u2_data The Cartesian data of \f$ u_2 \f$
 * @return The data of an element of the Cartesian space isomorphic to the Lie algebra
 */ 
static Vec_SC OMinus(const Mat_G& g1_data,const Mat_G& g2_data,const Mat_C & u1_data,const Mat_C & u2_data)
{ Vec_SC tmp;
  tmp.block(0,0,G::dim_,1) = G::OMinus(g1_data,g2_data).block(0,0,G::dim_,1);
  tmp.block(G::dim_,0,U::total_num_dim_,1) = u1_data - u2_data;
    return tmp;}

/**
 * Performs the O-minus operation \f$ \log(S_2^{-1}*S_1) i.e. S_1-S_2\f$
 * @param s1 The state  \f$ s_1 \f$
 * @param s2 The state  \f$ s_2 \f$
 * @return The data of an element of the Cartesian space isomorphic to the Lie algebra
 */ 
static Vec_SC OMinus(const State& s1, const State& s2)
{ return OMinus(s1.g_.data_,s2.g_.data_, s1.u_.data_,s2.u_.data_);}


/**
 * Performs the O-minus operation \f$ \log(S_2^{-1}*S_1) i.e. S_2-S_1\f$ with this being S_1
 * @param s2 The state  \f$ s_2 \f$
 * @return The data of an element of the Cartesian space isomorphic to the Lie algebra
 */ 
Vec_SC OMinus(const State& s2) const {
  return OMinus(*this,s2);
}


/**
 * Performs the O-Plus operation \f$ \text{state} \exp{\text{cartesian}}\f$ 
 * @param state The state  
 * @param cartesian An in the Cartesian space
 * @return A state that is the result of the O-Plus operation
 */ 
static State OPlus(const State& state, Vec_SC cartesian) {
  State tmp;
  tmp.g_.data_ = state.g_.OPlus(cartesian.block(0,0,U::total_num_dim_,1));
  tmp.u_.data_ = state.u_.data_ +  cartesian.block(G::dim_,0,U::total_num_dim_,1);
  return tmp;
}

/**
 * Performs the O-Plus operation \f$ \text{this} \exp{\text{cartesian}}\f$ 
 * @param cartesian An in the Cartesian space
 * @return A state that is the result of the O-Plus operation
 */ 
State OPlus( Vec_SC cartesian) const{
  return OPlus(*this,cartesian);
}

/**
 * Performs the O-Plus operation \f$ \text{this} \exp{\text{cartesian}}\f$ and sets the state to the result.
 * @param cartesian An in the Cartesian space
 */ 
void OPlusEQ( Vec_SC cartesian) {
  *this = OPlus(*this,cartesian);
}


/**
 * Computes the right Jacobian of the states Lie algebra
 * @param cartesian An element in the state's Cartesian space
 */ 
static Mat_SC Jr(const Vec_SC& cartesian) {
  
  Eigen::Matrix<DataType, Algebra::total_num_dim_,1> vec_algebra;
  vec_algebra.setZero();
  vec_algebra.block(0,0,Group::dim_,1) = cartesian.block(0,0,Group::dim_,1);
  Algebra group(vec_algebra);
  Mat_SC jacobian = Mat_SC::Identity();
  jacobian.block(0,0,Group::dim_,Group::dim_) = group.Jr().block(0,0,Group::dim_,Group::dim_);
  return jacobian;

}

/**
 * Computes the left Jacobian of the states Lie algebra
 * @param cartesian An element in the state's Cartesian space
 */ 
static Mat_SC Jl(const Vec_SC& cartesian) {
  Eigen::Matrix<DataType, Algebra::total_num_dim_,1> vec_algebra;
  vec_algebra.setZero();
  vec_algebra.block(0,0,Group::dim_,1) = cartesian.block(0,0,Group::dim_,1);
  Algebra group(vec_algebra);
  Mat_SC jacobian = Mat_SC::Identity();
  jacobian.block(0,0,Group::dim_,Group::dim_) = group.Jl().block(0,0,Group::dim_,Group::dim_);
  return jacobian;
}

/**
 * Computes the inverse of the right Jacobian of the states Lie algebra
 * @param cartesian An element in the state's Cartesian space
 */ 
static Mat_SC JrInv(const Vec_SC& cartesian) {
  Eigen::Matrix<DataType, Algebra::total_num_dim_,1> vec_algebra;
  vec_algebra.setZero();
  vec_algebra.block(0,0,Group::dim_,1) = cartesian.block(0,0,Group::dim_,1);
  Algebra group(vec_algebra);
  Mat_SC jacobian = Mat_SC::Identity();
  jacobian.block(0,0,Group::dim_,Group::dim_) = group.JrInv().block(0,0,Group::dim_,Group::dim_);
  return jacobian;
}

/**
 * Computes the inverse of the left Jacobian of the states Lie algebra
 * @param cartesian An element in the state's Cartesian space
 */ 
static Mat_SC JlInv(const Vec_SC& cartesian) {
  Eigen::Matrix<DataType, Algebra::total_num_dim_,1> vec_algebra;
  vec_algebra.setZero();
  vec_algebra.block(0,0,Group::dim_,1) = cartesian.block(0,0,Group::dim_,1);
  Algebra group(vec_algebra);
  Mat_SC jacobian = Mat_SC::Identity();
  jacobian.block(0,0,Group::dim_,Group::dim_) = group.JlInv().block(0,0,Group::dim_,Group::dim_);
  return jacobian;
}


/**
 * Computes the exponential of an element in the cartesian space
 */ 
static State Exp(const Vec_SC& cartesian) {

  State state;
  state.g_.data_ = State::Algebra::Exp(cartesian.block(0,0,State::Algebra::total_num_dim_,1));
  state.u_.data_ = cartesian.block(State::Group::dim_,0,State::Algebra::total_num_dim_,1);
  return state;
}

/**
 * Computes the Log of the state
 */
static Vec_SC Log(const State& state) {
  Vec_SC cartesian;
  cartesian.block(0,0,State::Group::dim_,1) = State::Algebra::Log(state.g_.data_);
  cartesian.block(State::Group::dim_,0, State::Algebra::total_num_dim_,1) = state.u_.data_;
  return cartesian;
}




};

typedef State<Rn,  double,2,1> R2_r2;
typedef State<Rn,  double,3,1> R3_r3;
typedef State<SO2, double,1,1> SO2_so2;
typedef State<SO3, double,3,1> SO3_so3;
typedef State<SE2, double,3,1> SE2_se2;
typedef State<SE3, double,6,1> SE3_se3;

}

#endif //_LIEGROUPS_INCLUDE_LIEGROUPS_STATE