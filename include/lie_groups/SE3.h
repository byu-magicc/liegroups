#ifndef _LIEGROUPS_INCLUDE_LIEGROUPS_SE3_
#define _LIEGROUPS_INCLUDE_LIEGROUPS_SE3_

#include <Eigen/Dense>
#include <iostream>
#include <lie_algebras/se3.h>

namespace lie_groups {

constexpr double kSE3_threshold_ = 1e-7;

class SE3 {

public:

Eigen::Matrix<double,4,4> data_;
Eigen::Map<Eigen::Matrix<double,3,1>> t_;  /** < The position */
Eigen::Ref<Eigen::Matrix3d> R_; /** < The rotation */
static const unsigned int dim_ = 6;
static const unsigned int size1_ = 4;
static const unsigned int size2_ = 4;


/**
 * Default constructor. Initializes group element to identity.
 */
SE3() : t_(data_.data()+12), R_(data_.block(0,0,3,3)), data_(Eigen::Matrix<double,4,4>::Identity() ) {}


/**
 * Copy constructor.
 */ 
SE3(const SE3 & g) : t_(data_.data()+12), R_(data_.block(0,0,3,3)), data_(g.data_) {}

/**
* Initializes group element to the one given. If verify is true
* it will check that the input is an element of \f$SE(3)\f$
* @param[in] data The data pertaining to an element of \f$SE(3)\f$
* @param verify If true, the constructor will verify that the provided 
* element is a member of \f$SE(3)\f$
*/
SE3(const Eigen::Matrix<double,4,4> & data, bool verify);

/**
* Initializes group element to the data of the one given. 
* @param[in] data  The data pertaining to an element of \f$SE(3)\f$
*/
SE3(const Eigen::Matrix<double,4,4> & data) :data_(data), t_(data_.data()+12), R_(data_.block(0,0,3,3)){}

/*
 * Returns the inverse of the element
 */ 
SE3 Inverse(){ return SE3::Inverse(this->data_);}


/*
 * Returns the inverse of the data of an element
 */ 

static Eigen::Matrix<double,4,4> Inverse(const Eigen::Matrix<double,4,4>& data){  
    return data.inverse();}

/**
 * Returns the identity element
 */ 
static SE3 Identity(){return SE3();}

/**
 * Returns the matrix adjoint map.
 * For this Lie group it is the identity map.
 */ 
Eigen::Matrix<double,6,6> Adjoint(){
    Eigen::Matrix<double,6,6> m;
    m.block(0,0,3,3) = R_;
    m.block(3,3,3,3) = R_;
    m.block(0,3,3,3) = se3::SSM(t_)*R_;
    m.block(3,0,3,3).setZero(); 
    return m;
}

/**
 * Computes the log of the element.
 */ 
Eigen::Matrix<double,6,1> Log() {return se3::Log(this->data_);}

/**
 * Performs the left group action on itself. i.e. this is on the left of
 * the bilinear operation.
 */ 
SE3 operator * (const SE3& g){ return SE3(data_ *g.data_);}

/**
 * Assignment Operator. Deep copy of the input parameter
 * @param g The element to be copied.
 */ 
void operator = (const SE3& g){ this->data_ = g.data_; }

/**
 * Performs the OPlus operation 
 * @param g_data The data belonging to the group element.
 * @param u_data The data belonging to the Cartesian space that is isomorphic to the Lie algebra.
 * @return The result of the BoxPlus operation.
 */ 
static Eigen::Matrix<double,4,4> OPlus(const Eigen::Matrix<double,4,4>& g_data, const Eigen::Matrix<double,6,1>& u_data)
{return g_data*se3::Exp(u_data);}

/**
 * Performs the OPlus operation 
 * @param u_data The data belonging to the Cartesian space that is isomorphic to the Lie algebra.
 * @return The result of the BoxPlus operation.
 */ 
Eigen::Matrix<double,4,4> OPlus(const Eigen::Matrix<double,6,1>& u_data)
{return OPlus(this->data_,u_data);}

/**
 * Performs the OPlus operation and assigns the result to the group element
 * @param u_data The data belonging to the Cartesian space that is isomorphic to the Lie algebra.
 * @return The result of the BoxPlus operation.
 */ 
void OPlusEq(const Eigen::Matrix<double,6,1>& u_data)
{data_ = this->OPlus(u_data);}

/**
 * Performs the BoxPlus operation 
 * @param g_data The data belonging to the group element.
 * @param u_data The data belonging to an element of the Lie algebra.
 * @return The result of the BoxPlus operation.
 */ 
static Eigen::Matrix<double,4,4> BoxPlus(const Eigen::Matrix<double,4,4>& g_data, const Eigen::Matrix<double,4,4>& u_data)
{return OPlus(g_data,se3::Vee(u_data));}

/**
 * Performs the BoxPlus operation 
 * @param u_data The data belonging to an element of the Lie algebra.
 * @return The result of the BoxPlus operation.
 */ 
Eigen::Matrix<double,4,4> BoxPlus(const Eigen::Matrix<double,4,4>& u_data)
{return this->OPlus(se3::Vee(u_data));}

/**
 * Performs the BoxPlus operation and assigns the result to the group element
 * @param u_data The data belonging to an element of the Lie algebra.
 * @return The result of the BoxPlus operation.
 */ 
void BoxPlusEq(const Eigen::Matrix<double,4,4>& u_data)
{data_ = this->BoxPlus(u_data);}

/**
 * Performs the BoxPlus operation 
 * @param g An element of the group
 * @param u An element of the Lie algebra.
 * @return The result of the BoxPlus operation.
 */ 
static SE3 BoxPlus(const SE3& g, const se3& u)
{return SE3(OPlus(g.data_,u.data_));}

/**
 * Performs the BoxPlus operation 
 * @param u An element of the Lie algebra.
 * @return The result of the BoxPlus operation.
 */ 
SE3 BoxPlus(const se3& u)
{return SE3::BoxPlus(*this,u);}

/**
 * Performs the BoxPlus operation and assigns the result to the group element
 * @param u An element of the Lie algebra.
 * @return The result of the BoxPlus operation.
 */ 
void BoxPlusEq(const se3& u)
{data_ = (this->BoxPlus(u)).data_;}


/**
 * Performs the O-minus operation \f$ \log(g_1^-1*g_2) \f$
 * @param g_data1 The data of  \f$ g_1 \f$
 * @param g_data2 The data of \f$ g_2 \f$
 * @return The data of an element of the Cartesian space isomorphic to the Lie algebra
 */ 
static Eigen::Matrix<double,6,1> OMinus(const Eigen::Matrix<double,4,4>& g1_data,const Eigen::Matrix<double,4,4>& g2_data)
{return se3::Log(SE3::Inverse(g1_data)*g2_data);}


/**
 * Performs the O-minus operation with this being \f$ g_1 \f$ in the equation \f$ \log(g_1^-1*g_2) \f$
 * @param g_data The data of  \f$ g_2 \f$
 * @return The data of an element of the Cartesian space isomorphic to the Lie algebra
 */ 
Eigen::Matrix<double,6,1> OMinus(const Eigen::Matrix<double,4,4>& g_data)
{return SE3::OMinus(data_,g_data);}

/**
 * Performs the O-minus operation \f$ \log(g_1^-1*g_2) \f$
 * @param g_data1 The data of  \f$ g_1 \f$
 * @param g_data2 The data of \f$ g_2 \f$
 * @return The data of an element of the Lie algebra
 */ 
static Eigen::Matrix<double,4,4> BoxMinus(const Eigen::Matrix<double,4,4>& g1_data,const Eigen::Matrix<double,4,4>& g2_data)
{return se3::Wedge(SE3::OMinus(g1_data, g2_data));}


/**
 * Performs the O-minus operation with this being \f$ g_1 \f$ in the equation \f$ \log(g_1^-1*g_2) \f$
 * @param g_data The data of  \f$ g_2 \f$
 * @return The data of an element of the Lie algebra
 */ 
Eigen::Matrix<double,4,4> BoxMinus(const Eigen::Matrix<double,4,4>& g_data)
{return BoxMinus(this->data_,g_data);}

/**
 * Performs the O-minus operation with this being \f$ g_1 \f$ in the equation \f$ \log(g_1^-1*g_2) \f$
 * @param g An element of the group
 * @return An element of the Lie algebra
 */ 
se3 BoxMinus(const SE3& g){return se3(se3::Vee(this->BoxMinus(g.data_)));}

/**
 * Prints the content of the data
 */ 
void Print() {
    std::cout << data_ << std::endl; 
}

/**
 * Verifies that the data of an element properly corresponds to the set. 
 */ 
static bool isElement(const Eigen::Matrix<double,4,4>& data);

};

}


#endif // _LIEGROUPS_INCLUDE_LIEGROUPS_SE3_