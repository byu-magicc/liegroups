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

template <class G, class U> 
class State {

public:

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
 * Copy constructor.
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
State(const Eigen::Matrix<double,G::size1_,G::size2_> & g_data, const Eigen::Matrix<double,U::size1_,U::size2_> & u_data, bool verify) : g_(g_data,verify), u_(u_data,verify) {}


/**
* Initializes state using the data provied. If verify is true
* it will check that the inputs are elements of the group and Lie algebra.
* @param g_data The data pertaining to the group.
* @param u_data The data pertaining to the Lie algebra
* @param verify If true, the constructor will verify that the provided data 
* are elements of the group and Lie algebra.
*/
State(const Eigen::Matrix<double,G::size1_,G::size2_> & g_data, const Eigen::Matrix<double,U::size1_,U::size1_> & u_data, bool verify) : g_(g_data,verify), u_(u_data,verify) {}

/*
 * Returns the inverse of the element
 */ 
State Inverse(){ return State(g_.Inverse(),u_.Inverse());}

/**
 * Returns the identity element
 */ 
static State Identity(){return State();}

/**
 * Returns the matrix group adjoint map.
 * For this Lie group it is the identity map.
 */ 
Eigen::Matrix<double,G::dim_,G::dim_> Adjoint(){return g_.Adjoint();}


};

typedef State<Rn<2>,rn<2>> R2_r2;
typedef State<Rn<3>,rn<3>> R3_r3;
typedef State<SE2,se2> SE2_se2;
typedef State<SE3,se3> SE3_se3;

}

#endif //_LIEGROUPS_INCLUDE_LIEGROUPS_STATE