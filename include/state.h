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

G g;
U u;

State()=default;



};

}

#endif //_LIEGROUPS_INCLUDE_LIEGROUPS_STATE