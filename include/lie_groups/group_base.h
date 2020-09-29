#ifndef _LIEGROUPS_INCLUDE_LIEGROUPS_GROUPBASE_
#define _LIEGROUPS_INCLUDE_LIEGROUPS_GROUPBASE_


#include <Eigen/Dense>
#include <iostream>

namespace lie_groups {

class GroupBase {

public:

    Eigen::MatrixXd g_; 

    // virtual GroupBase operator*(const GroupBase& b)=0;

    virtual void Print()=0;

private:


};


}


#endif // _LIEGROUPS_INCLUDE_LIEGROUPS_GROUPBASE_