#include <Eigen/Dense>

#include "lie_algebras/se2.h"
#include "gtest/gtest.h"

namespace lie_groups {

TEST(se2TEST, Constructors) {

// default constructor
se2 u1;
Eigen::Matrix<double,3,1> Identity = Eigen::Matrix<double,3,1>::Zero();
ASSERT_EQ(u1.data_, Identity ) << "Element not initialized to identity. i.e. zero.";
u1.data_.setRandom();
ASSERT_EQ(u1.data_.block(0,0,2,1),u1.p_) << "Translational velocity map not set properly.";
ASSERT_EQ(u1.data_.block(2,0,1,1),u1.th_) << "Angular velocity map not set properly.";



// assignment constructor Cartesian space
Eigen::Matrix<double,3,1> v = Eigen::Matrix<double,3,1>::Random();
se2 u2(v);
ASSERT_EQ(u2.data_,v) << "Error with the assignment constructor";
ASSERT_EQ(u2.data_.block(0,0,2,1),u2.p_) << "Translational velocity map not set properly.";
ASSERT_EQ(u2.data_.block(2,0,1,1),u2.th_) << "Angular velocity map not set properly.";

// assignment constructor invaid element
Eigen::Matrix3d V1 = Eigen::Matrix3d::Random();
V1 = (V1+Eigen::Matrix3d(V1.transpose())); // Symmetric matrix
if (V1.norm() == 0) // If the zero matrix 
    V1(0,0) = 1;
se2 u3(V1,true);
ASSERT_EQ(u3.data_, Identity) << "Error with the assignment constructor";
ASSERT_EQ(u3.data_.block(0,0,2,1),u3.p_) << "Translational velocity map not set properly.";
ASSERT_EQ(u3.data_.block(2,0,1,1),u3.th_) << "Angular velocity map not set properly.";


// assignment constructor valid element
Eigen::Matrix3d V2;
Eigen::Matrix<double,3,1> v2; 
V2 << 0, -4, 5, 4, 0, 6, 0, 0, 0;
v2 << 5, 6, 4;
se2 u4(V2,true);
ASSERT_EQ(u4.data_,v2) << "Error with the assignment constructor";
ASSERT_EQ(u4.data_.block(0,0,2,1),u4.p_) << "Translational velocity map not set properly.";
ASSERT_EQ(u4.data_.block(2,0,1,1),u4.th_) << "Angular velocity map not set properly.";


// copy constructor
se2 u5(u4);
ASSERT_EQ(u5.data_,u4.data_) << "Error with the copy constructor";
ASSERT_EQ(u5.data_.block(0,0,2,1),u5.p_) << "Translational velocity map not set properly.";
ASSERT_EQ(u5.data_.block(2,0,1,1),u5.th_) << "Angular velocity map not set properly.";


// identity function test
se2 u6 = se2::Identity();
ASSERT_EQ(u6.data_,Identity) << "Error with the identity function";

}

// Tests the Bracket, Adjoint, Wedge, Vee, Exp, and Norm functions
TEST(se2TEST, BAWVEN) {

se2 u1(Eigen::Matrix<double, 3,1>::Random());
se2 u2(Eigen::Matrix<double, 3,1>::Random());
se2 u3(Eigen::Matrix<double, 3,1>::Random());
se2 u4(Eigen::Matrix<double, 3,1>::Random());


// Wedge test
Eigen::Matrix3d m2;
m2 << 0, -u3.th_(0), u3.p_(0), u3.th_(0), 0, u3.p_(1), 0,0,0;
ASSERT_EQ(u3.Wedge(),m2) << "Error with wedge function";

// Adjoint test
Eigen::Matrix3d m3;
m3  << 0, -u3.th_(0), u3.p_(1), u3.th_(0), 0, -u3.p_(0), 0,0,0;
ASSERT_EQ(u3.Adjoint(),m3) << "Error with adjoint test"; 

// Bracket and adjoint test
Eigen::Matrix3d m1;
m1 = u1.Wedge()*u2.Wedge() - u2.Wedge()*u1.Wedge();
ASSERT_EQ( (u1.Bracket(u2)).Wedge(), m1) << "Error with bracket or Adjoint function";


// Vee test
ASSERT_EQ(u4.Vee(),u4.data_) << "Error with vee function";

// Exp test
Eigen::Matrix3d E;
Eigen::Matrix3d A;
double factorial=1;
E.setZero();
A.setIdentity();

for (long int i = 1; i < 1000; ++i) {

    E = E+A/factorial;
    A = A*u4.Wedge();
    factorial = factorial*i;
}
ASSERT_LE((u4.Exp()-E).norm(), 1e-10 ) << "Error with exp function";

// Norm test
ASSERT_EQ(u1.Norm(),u1.data_.norm()) << "Error with norm function";


}

TEST(se2Test, OperatorTest) {

Eigen::Matrix<double,3,1> v1 = Eigen::Matrix<double,3,1>::Random();
Eigen::Matrix<double,3,1> v2 = Eigen::Matrix<double,3,1>::Random();
Eigen::Matrix<double,3,1> v3 = Eigen::Matrix<double,3,1>::Random();
Eigen::Matrix<double,3,1> v4 = Eigen::Matrix<double,3,1>::Random();

se2 u1(v1);
se2 u2(v2);
se2 u3(v3);
se2 u4(v4);

// Addition test
se2 u5 = u1+u2;
ASSERT_EQ(u5.data_, u1.data_ + u2.data_) << "Error with addition operator";

// Subtraction test
se2 u6 = u3-u4;
ASSERT_EQ(u6.data_, u3.data_-u4.data_) << "Error with subtraction operator";

// Scalar multiplication test
se2 u7 = u1*6;
ASSERT_EQ(u7.data_,u1.data_*6) << "Error with multiplication test";

// Assignment test
se2 u8 = u4;
ASSERT_EQ(u8.data_,u4.data_) << "Error with assignment test";

}

TEST(se2Test, JacobianTest) {

Eigen::Matrix<double,3,1> v1 = Eigen::Matrix<double,3,1>::Random();
Eigen::Matrix<double,3,1> v2 = Eigen::Matrix<double,3,1>::Random();
Eigen::Matrix<double,3,1> v3 = Eigen::Matrix<double,3,1>::Random();
Eigen::Matrix<double,3,1> v4 = Eigen::Matrix<double,3,1>::Random();
Eigen::Matrix<double,3,1> v5;
v5 << 1,2,3;



se2 u1(v1);
se2 u2(v2);
se2 u3(v3);
se2 u4(v4);
se2 u5(v5);

std::cerr << u5.Jl() << std::endl;



// ASSERT_EQ(u1.Jl(),Eigen::Matrix2d::Identity()) << "Error with left Jacobian";
// ASSERT_EQ(u1.Jr(),Eigen::Matrix2d::Identity()) << "Error with right Jacobian";
// ASSERT_EQ(u1.JlInv(),Eigen::Matrix2d::Identity()) << "Error with left Jacobian inverse";
// ASSERT_EQ(u1.JrInv(),Eigen::Matrix2d::Identity()) << "Error with right Jacobian inverse";

// ASSERT_EQ((u1.Jl(u2)).data_,u2.data_) << "Error with left Jacobian";
// ASSERT_EQ((u1.Jr(u2)).data_,u2.data_) << "Error with right Jacobian";
// ASSERT_EQ((u1.JlInv(u2)).data_,u2.data_) << "Error with left Jacobian inverse";
// ASSERT_EQ((u1.JrInv(u2)).data_,u2.data_) << "Error with right Jacobian inverse";


}


}
