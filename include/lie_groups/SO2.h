#ifndef _LIEGROUPS_INCLUDE_LIEGROUPS_SO2_
#define _LIEGROUPS_INCLUDE_LIEGROUPS_SO2_

#include "lie_groups/group_base.h"

namespace lie_groups {


class SO2 : public GroupBase {

public:

SO2(){

g_ = Eigen::Matrix<double,2,2>::Identity();


}

void Print() {

std::cout << "SO(2): " << g_(0,0) << " , " << g_(0,1) << std::endl << "\t" << g_(1,0) << " , " << g_(1,1) << std::endl; 

}

private:

};

}


#endif // _LIEGROUPS_INCLUDE_LIEGROUPS_SO2_