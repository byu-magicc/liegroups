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

//#########################################################################


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
ASSERT_EQ(u4.Vee(u4.Wedge()),u4.data_) << "Error with vee function";

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

// Log test
ASSERT_LT((se2::Log(u4.Exp())-u4.data_).norm(),1e-8);

// Norm test
ASSERT_EQ(u1.Norm(),u1.data_.norm()) << "Error with norm function";


}

//#########################################################################


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


// These tests are performed in a very specific order. They build
// upon each other. An error in an earlier test can cause failures in
// following tests even if there isn't an error. 
TEST(se2Test, JacobianTest) {

Eigen::Matrix<double,3,1> v1 = Eigen::Matrix<double,3,1>::Random();
Eigen::Matrix<double,3,1> v2 = Eigen::Matrix<double,3,1>::Random();
Eigen::Matrix<double,3,1> e1;
Eigen::Matrix<double,3,1> e2;
Eigen::Matrix<double,3,1> e3;
Eigen::Matrix3d estimated_jr_inv;
e1 << 1,0,0;
e2 << 0,1,0;
e3 << 0,0,1;
double dt = 1e-7;


se2 u1(v1);
se2 u2(-v1);
se2 u3(v2);



estimated_jr_inv.block(0,0,3,1) = (se2(se2::Log(u1.Exp()*se2(e1*dt).Exp())).Vee() - v1)/dt;
estimated_jr_inv.block(0,1,3,1) = (se2(se2::Log(u1.Exp()*se2(e2*dt).Exp())).Vee() - v1)/dt;
estimated_jr_inv.block(0,2,3,1) = (se2(se2::Log(u1.Exp()*se2(e3*dt).Exp())).Vee() - v1)/dt;

// Test JrInverse
ASSERT_LE( (estimated_jr_inv-u1.JrInv()).norm(),1e-8) << "Error with the right Jacobian";

// Test JlInverse: JlInv(-v) = JrInv(v)
ASSERT_LE( (u2.JlInv()-u1.JrInv()).norm(),1e-8) << "Error with the right Jacobian";

// Test Jr
ASSERT_LE( (u1.Jr()*u1.JrInv()-Eigen::Matrix3d::Identity()).norm(),1e-8) << "Error with the right Jacobian";

// Test Jl
ASSERT_LE( (u1.Jl()*u1.JlInv()-Eigen::Matrix3d::Identity()).norm(),1e-8) << "Error with the right Jacobian";


// Test Remaining


ASSERT_LE( ((u1.Jl(u3)).data_-u1.Jl()*u3.data_).norm(), 1e-8) << "Error with left Jacobian";
ASSERT_LE( ((u1.Jr(u3)).data_-u1.Jr()*u3.data_).norm(), 1e-8) << "Error with right Jacobian";
ASSERT_LE( ((u1.JlInv(u3)).data_-u1.JlInv()*u3.data_).norm(), 1e-8) << "Error with left Jacobian inverse";
ASSERT_LE( ((u1.JrInv(u3)).data_-u1.JrInv()*u3.data_).norm(), 1e-8) << "Error with right Jacobian inverse";


}


}
