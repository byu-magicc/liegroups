#ifndef _LIEGROUPS_INCLUDE_LIEGROUPS_GROUPBASE_
#define _LIEGROUPS_INCLUDE_LIEGROUPS_GROUPBASE_


#include <Eigen/Dense>
#include <iostream>

namespace lie_groups {

class GroupBase {

public:


    Eigen::MatrixXd g_; 

    /**
    * This operator performs the left group action on the element passed in.
    * @param g The group element that the left action is being applied to.
    * @return The result of the left group action.
    */
    virtual GroupBase operator * (const GroupBase& g)=0;

    /**
     * Computes the matrix adjoint representation.
     * @return The matrix adjoint representation.
     */
    virtual Eigen::MatrixXd Adjoint()=0;

    /**
    * Computes and returns the identity element.
    */
    virtual Eigen::MatrixXd Identity()=0;

    /**
    * Computes and returns the inverse element.
    */
    virtual Eigen::MatrixXd Inverse()=0; 

    /**
     * Applies the log map to the element.
     * @return The result of the log function.
     */ 
    virtual Eigen::MatrixXd Log()=0;

    // /**
    //  * Performs the left box-plus operation 
    //  * \f$ g\exp{u} \f$.
    //  * @param u The data of an element of the Lie Algebra
    //  * @return The rusult of the operation.
    //  */ 
    // virtual Eigen::MatrixXd LBoxPlus(const Eigen::MatrixXd& u);

    // /**
    //  * Performs the right box-plus operation 
    //  * \f$ \exp{u}g \f$.
    //  * @param u The data of an element of the Lie Algebra
    //  * @return The rusult of the operation.
    //  */ 
    // virtual Eigen::MatrixXd RBoxPlus(const Eigen::MatrixXd& u);

    // /**
    //  * Performs the left O-plus operation 
    //  * \f$ g\exp{u^{\wedge}} \f$.
    //  * @param u The data of an element of the Cartesian space
    //  * @return The rusult of the operation.
    //  */ 
    // virtual Eigen::MatrixXd LOPlus(const Eigen::MatrixXd& u);

    // /**
    //  * Performs the right O-plus operation 
    //  * \f$ \exp{u^{\wedge}}g \f$.
    //  * @param u A element of the Cartesian space
    //  * @return The rusult of the operation.
    //  */ 
    // virtual Eigen::MatrixXd ROPlus(const Eigen::MatrixXd& u);

    

    virtual void Print()=0;

private:


};


}


#endif // _LIEGROUPS_INCLUDE_LIEGROUPS_GROUPBASE_