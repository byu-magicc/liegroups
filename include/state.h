#ifndef _LIEGROUPS_INCLUDE_LIEGROUPS_STATE_
#define _LIEGROUPS_INCLUDE_LIEGROUPS_STATE

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

G g_;
U u_;

/**
 * Default constructor. Initializes group element to identity.
 */
State()=default;

/**
 * Copy constructor.
 */ 
State(const G& g, const U& u) : g_(g), u_(u) {};

/**
* Initializes state using the data provied. 
* @param g_data The data pertaining to the group.
* @param u_data The data pertaining to a element of the Cartesian space isomorphic Lie algebra
* are elements of the group and Lie algebra.
*/
State(const Eigen::Matrix<double,G::size1_,G::size2_> & g_data, const Eigen::Matrix<double,U::size1_,U::size2_> & u_data, bool verify) : g(g_data,verify), u(u_data,verify) {}


/**
* Initializes state using the data provied. If verify is true
* it will check that the inputs are elements of the group and Lie algebra.
* @param g_data The data pertaining to the group.
* @param u_data The data pertaining to the Lie algebra
* @param verify If true, the constructor will verify that the provided data 
* are elements of the group and Lie algebra.
*/
State(const Eigen::Matrix<double,G::size1_,G::size2_> & g_data, const Eigen::Matrix<double,U::size1_,U::size1_> & u_data, bool verify) : g(g_data,verify), u(u_data,verify) {}

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

}

#endif //_LIEGROUPS_INCLUDE_LIEGROUPS_STATE