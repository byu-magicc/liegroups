#ifndef _LIEGROUPS_INCLUDE_LIEALGEBRAS_ALGEBRABASE_
#define _LIEGROUPS_INCLUDE_LIEALGEBRAS_ALGEBRABASE_


#include <Eigen/Dense>
#include <iostream>

namespace lie_groups {

/**
 * \class The AlgebraBase class defines the required functions for a Lie algebra.
 */ 

class AlgebraBase {

public:

    Eigen::Matrix<double,Eigen::Dynamic,1> data_; /** < The element of the Lie algebra represented in Cartesian space */

    /**
     * Performs the Lie bracket. The input is the right parameter 
     * of the Lie bracket function. \f$ \[v,u\] \f$
     * @param u An element of the Lie algebra.
     * @return The result of the Lie bracket operation.
     */ 
    virtual AlgebraBase * Bracket(const AlgebraBase* u)=0;

    /**
     * Computes and returns the matrix adjoint representation of the Lie algebra.
     */ 
    virtual Eigen::MatrixXd Adjoint()=0;

    /**
     * Computes the Wedge operation which maps an element of the Cartesian space to the Lie algebra.
     * @return The result of the Wedge operation.
     */
    virtual Eigen::MatrixXd Wedge()=0;

    /**
     * Computes the Vee operation which maps an element of the Lie algebra to the Cartesian space.
     * @return The result of the Vee operation.
     */
    virtual Eigen::MatrixXd Vee()=0;

    /**
     * Computes the exponential of the element of the Lie algebra.
     * @return The data associated to the group element.
     */
    virtual Eigen::MatrixXd Exp()=0;

    /**
     * Computes and returns the Euclidean norm of the element of the Lie algebra
     */ 
    virtual double Norm()=0;

    /**
     * Computes and returns the matrix of the Left Jacobian
     */ 
    virtual Eigen::MatrixXd Jl()=0;

    /**
     * Computes the left Jacobian using the element of *this and applies it to the
     * parameter provided. \f$ J_l(v)u \f$
     * @param u An element of the Lie algebra.
     */ 
    virtual AlgebraBase * Jl(const AlgebraBase * u)=0;

    /**
     * Computes and returns the matrix of the Left Jacobian inverse
     */ 
    virtual Eigen::MatrixXd JlInv()=0;

    /**
     * Computes the left Jacobian inverse using the element of *this and applies it to the
     * parameter provided. \f$ J_l(v)u \f$
     * @param u An element of the Lie algebra.
     */ 
    virtual AlgebraBase * JlInv(const AlgebraBase * u)=0;

    /**
     * Computes and returns the matrix of the Right Jacobian
     */ 
    virtual Eigen::MatrixXd Jr()=0;

    /**
     * Computes the right Jacobian using the element of *this and applies it to the
     * parameter provided. \f$ J_l(v)u \f$
     * @param u An element of the Lie algebra.
     */ 
    virtual AlgebraBase * Jr(const AlgebraBase * u)=0;

    /**
     * Computes and returns the matrix of the right Jacobian inverse
     */ 
    virtual Eigen::MatrixXd JrInv()=0;

    /**
     * Computes the right Jacobian inverse using the element of *this and applies it to the
     * parameter provided. \f$ J_l(v)u \f$
     * @param u An element of the Lie algebra.
     */ 
    virtual AlgebraBase * JrInv(const AlgebraBase * u)=0;

    /**
     * Adds two elements of the Algebra together
     * @param u An element of the Lie algebra.
     */ 
    virtual AlgebraBase * operator + (const AlgebraBase * u)=0;

    /**
     * Subtracts two elements of the Algebra together
     * @param u An element of the Lie algebra.
     */ 
    virtual AlgebraBase * operator - (const AlgebraBase * u)=0;

    /**
     * Performs Scalar multiplication and returns the result.
     * @param scalar The scalar that will scale the element of the Lie algebra
     */ 
    virtual AlgebraBase * operator * (const double scalar)=0;

    /**
     * Prints the data of the element.
     */ 
    virtual void Print()=0;

    /**
     * Returns the Identity element.
     */
    virtual AlgebraBase Identity()=0; 
};


}


#endif //_LIEGROUPS_INCLUDE_LIEALGEBRAS_ALGEBRABASE_