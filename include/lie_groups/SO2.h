#ifndef _LIEGROUPS_INCLUDE_LIEGROUPS_SO2_
#define _LIEGROUPS_INCLUDE_LIEGROUPS_SO2_

#include "lie_groups/group_base.h"

namespace lie_groups {


class SO2 : public GroupBase {

public:

/**
 * Default constructor. Initializes group element to identity.
 */
SO2();


/**
 * Copy constructor.
 */ 
SO2(const SO2 & g);

/**
* Initializes group element to the one given.
*/
SO2(const Eigen::Matrix<double,2,2> & g);


// SO2 operator * (const SO2* g){
//     return g_*g.g_;
// }

void Print() {

std::cout << "SO(2): " << g_(0,0) << " , " << g_(0,1) << std::endl << "\t" << g_(1,0) << " , " << g_(1,1) << std::endl; 

}

private:

};

}


#endif // _LIEGROUPS_INCLUDE_LIEGROUPS_SO2_