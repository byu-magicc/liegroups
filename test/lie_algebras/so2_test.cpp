#include <Eigen/Dense>

#include "lie_algebras/so2.h"
#include "gtest/gtest.h"

namespace lie_groups {
//#########################################################################
TEST(so2TEST, Constructors) {

// default constructor
so2 u1;
ASSERT_EQ(u1.data_(0), 0 ) << "Element not initialized to identity. i.e. zero.";

// assignment constructor Cartesian space
Eigen::Matrix<double,1,1> v = Eigen::Matrix<double,1,1>::Random();
so2 u2(v);
ASSERT_EQ(u2.data_,v) << "Error with the assignment constructor";

// assignment constructor invaid element
Eigen::Matrix2d V = Eigen::Matrix2d::Random();
V = (V+Eigen::Matrix2d(V.transpose())); // Symmetric matrix
if (V.norm() == 0) // If the zero matrix 
    V(0,0) = 1;
so2 u3(V,true);
ASSERT_NE(u3.data_(0),V(1,0)) << "Error with the assignment constructor";

// assignment constructor valid element
V = V-Eigen::Matrix2d(V.transpose());
so2 u4(V,true);
ASSERT_EQ(u4.data_(0),V(1,0)) << "Error with the assignment constructor";

// copy constructor
so2 u5(u2);
ASSERT_EQ(u5.data_,v) << "Error with the copy constructor";

// identity function test
so2 u6 = so2::Identity();
ASSERT_EQ(u6.data_(0),0) << "Error with the identity function";

}

//#########################################################################


// Tests the Bracket, Adjoint, Wedge, Vee, Exp, and Norm functions
TEST(so2TEST, BAWVEN) {

// Bracket test
so2 u1;
so2 u2;
so2 u3 = u1.Bracket(u2);
ASSERT_EQ(u3.data_(0), 0) << "Error with bracket function";

// Adjoint test
ASSERT_EQ(u1.Adjoint(), Eigen::Matrix2d::Identity()) << "Error with Adjoint function";

// Wedge test
Eigen::Matrix<double,1,1> v = Eigen::Matrix<double,1,1>::Random();
so2 u4(v);
Eigen::Matrix2d V;
V << 0, -v, v, 0;
ASSERT_EQ(u4.Wedge(),V) << "Error with wedge function";

// Vee test
ASSERT_EQ(u4.Vee(),u4.data_) << "Error with vee function";
ASSERT_EQ(so2::Vee(u4.Wedge()),u4.data_) << "Error with vee function";

// Exp test
Eigen::Matrix<double,1,1> v1 = Eigen::Matrix<double,1,1>::Random();
so2 u5(v1);
Eigen::Matrix2d g=u5.Exp();

ASSERT_EQ(g.transpose()*g, Eigen::Matrix2d::Identity()) << "Error with exp function";
ASSERT_EQ(g.determinant(), 1) << "Error with exp function";

// Norm test
Eigen::Matrix<double,1,1> v2 = Eigen::Matrix<double,1,1>::Random();
so2 u6(v2);
ASSERT_EQ(u6.Norm(),v2.norm()) << "Error with norm function";


}

//#########################################################################


TEST(so2Test, OperatorTest) {

Eigen::Matrix<double,1,1> v1 = Eigen::Matrix<double,1,1>::Random();
Eigen::Matrix<double,1,1> v2 = Eigen::Matrix<double,1,1>::Random();
Eigen::Matrix<double,1,1> v3 = Eigen::Matrix<double,1,1>::Random();
Eigen::Matrix<double,1,1> v4 = Eigen::Matrix<double,1,1>::Random();

so2 u1(v1);
so2 u2(v2);
so2 u3(v3);
so2 u4(v4);

// Addition test
so2 u5 = u1+u2;
ASSERT_EQ(u5.data_, u1.data_ + u2.data_) << "Error with addition operator";

// Subtraction test
so2 u6 = u3-u4;
ASSERT_EQ(u6.data_, u3.data_-u4.data_) << "Error with subtraction operator";

// Scalar multiplication test
so2 u7 = u1*6;
ASSERT_EQ(u7.data_,u1.data_*6) << "Error with multiplication test";

// Assignment test
so2 u8 = u4;
ASSERT_EQ(u8.data_,u4.data_) << "Error with assignment test";

}

TEST(so2Test, JacobianTest) {

Eigen::Matrix<double,1,1> v1 = Eigen::Matrix<double,1,1>::Random();
Eigen::Matrix<double,1,1> v2 = Eigen::Matrix<double,1,1>::Random();
Eigen::Matrix<double,1,1> v3 = Eigen::Matrix<double,1,1>::Random();
Eigen::Matrix<double,1,1> v4 = Eigen::Matrix<double,1,1>::Random();

so2 u1(v1);
so2 u2(v2);
so2 u3(v3);
so2 u4(v4);

ASSERT_EQ(u1.Jl(),Eigen::Matrix2d::Identity()) << "Error with left Jacobian";
ASSERT_EQ(u1.Jr(),Eigen::Matrix2d::Identity()) << "Error with right Jacobian";
ASSERT_EQ(u1.JlInv(),Eigen::Matrix2d::Identity()) << "Error with left Jacobian inverse";
ASSERT_EQ(u1.JrInv(),Eigen::Matrix2d::Identity()) << "Error with right Jacobian inverse";

ASSERT_EQ((u1.Jl(u2)).data_,u2.data_) << "Error with left Jacobian";
ASSERT_EQ((u1.Jr(u2)).data_,u2.data_) << "Error with right Jacobian";
ASSERT_EQ((u1.JlInv(u2)).data_,u2.data_) << "Error with left Jacobian inverse";
ASSERT_EQ((u1.JrInv(u2)).data_,u2.data_) << "Error with right Jacobian inverse";


}


}
