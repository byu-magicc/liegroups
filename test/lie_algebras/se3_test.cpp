#include <Eigen/Dense>

#include "lie_algebras/se3.h"
#include "gtest/gtest.h"

namespace lie_groups {

TEST(se3TEST, Constructors) {

typedef float DataType;
typedef Eigen::Matrix<DataType,6,1> Vec6d;
typedef Eigen::Matrix<DataType,4,4> Mat4d;


// default constructor
se3<DataType> u1;
Vec6d Identity = Vec6d::Zero();
ASSERT_EQ(u1.data_, Identity ) << "Element not initialized to identity. i.e. zero.";
u1.data_.setRandom();
ASSERT_EQ(u1.data_.block(0,0,3,1),u1.p_) << "Translational velocity map not set properly.";
ASSERT_EQ(u1.data_.block(3,0,3,1),u1.th_) << "Angular velocity map not set properly.";



// assignment constructor Cartesian space
Vec6d v = Vec6d::Random();
se3<DataType> u2(v);
ASSERT_EQ(u2.data_,v) << "Error with the assignment constructor";
ASSERT_EQ(u2.data_.block(0,0,3,1),u2.p_) << "Translational velocity map not set properly.";
ASSERT_EQ(u2.data_.block(3,0,3,1),u2.th_) << "Angular velocity map not set properly.";

// assignment constructor invaid element
Mat4d V1;
V1.setRandom();
V1 = (V1+Mat4d(V1.transpose())); // Symmetric matrix
if (V1.norm() == 0) // If the zero matrix 
    V1(0,0) = 1;
se3<DataType> u3(V1,true);
ASSERT_EQ(u3.data_, Identity) << "Error with the assignment constructor";
ASSERT_EQ(u3.data_.block(0,0,3,1),u3.p_) << "Translational velocity map not set properly.";
ASSERT_EQ(u3.data_.block(3,0,3,1),u3.th_) << "Angular velocity map not set properly.";


// assignment constructor valid element
Mat4d V2;
Vec6d v2; 
V2 << 0, -1, 2, 4, 1, 0, -3, 5, -2, 3, 0, 6,0,0,0,0;
v2 << 4,5,6,3,2,1;
se3<DataType> u4(V2,true);
ASSERT_EQ(u4.data_,v2) << "Error with the assignment constructor";
ASSERT_EQ(u4.data_.block(0,0,3,1),u4.p_) << "Translational velocity map not set properly.";
ASSERT_EQ(u4.data_.block(3,0,3,1),u4.th_) << "Angular velocity map not set properly.";


// copy constructor
se3<DataType> u5(u4);
ASSERT_EQ(u5.data_,u4.data_) << "Error with the copy constructor";
ASSERT_EQ(u5.data_.block(0,0,3,1),u5.p_) << "Translational velocity map not set properly.";
ASSERT_EQ(u5.data_.block(3,0,3,1),u5.th_) << "Angular velocity map not set properly.";


// identity function test
se3<DataType> u6 = se3<DataType>::Identity();
ASSERT_EQ(u6.data_,Identity) << "Error with the identity function";

}

//#########################################################################


// Tests the Bracket, Adjoint, Wedge, Vee, Exp, and Norm functions
TEST(se3TEST, BAWVEN) {

typedef double DataType;
typedef Eigen::Matrix<DataType,3,1> Vec3d;
typedef Eigen::Matrix<DataType,6,1> Vec6d;
typedef Eigen::Matrix<DataType,3,3> Mat3d;
typedef Eigen::Matrix<DataType,4,4> Mat4d;
typedef Eigen::Matrix<DataType,6,6> Mat6d;

se3<DataType> u1(Vec6d::Random());
se3<DataType> u2(Vec6d::Random());
se3<DataType> u3(Vec6d::Random());
se3<DataType> u4(Vec6d::Random());


// Wedge test
Mat4d m1;
m1 << 0, -u1.data_(5), u1.data_(4), u1.data_(0), u1.data_(5), 0, -u1.data_(3), u1.data_(1), -u1.data_(4), u1.data_(3),0,u1.data_(2),0,0,0,0;

ASSERT_EQ(u1.Wedge(),m1) << "Error with wedge function";

// SSM test
Mat3d m2;
Vec3d x;
x.setRandom();
m2 << 0, -x(2), x(1), x(2), 0, -x(0), -x(1), x(0), 0;
ASSERT_EQ(m2,se3<DataType>::SSM(x));


// Adjoint test

Mat6d m3 = Mat6d::Zero();
m3.block(0,0,3,3) = se3<DataType>::SSM(u2.th_);
m3.block(3,3,3,3) = se3<DataType>::SSM(u2.th_);
m3.block(0,3,3,3) = se3<DataType>::SSM(u2.p_);
ASSERT_EQ(u2.Adjoint(),m3) << "Error with adjoint test"; 

// Bracket and adjoint test
Mat4d m4;
m4 = u1.Wedge()*u2.Wedge() - u2.Wedge()*u1.Wedge();
ASSERT_LE( ((u1.Bracket(u2)).Wedge()- m4).norm(), kse3_threshold_) << "Error with bracket or Adjoint function";


// Vee test
ASSERT_EQ(u4.Vee(),u4.data_) << "Error with vee function";
ASSERT_EQ(u4.Vee(u4.Wedge()),u4.data_) << "Error with vee function";

// Exp test
Mat4d E;
Mat4d A;
double factorial=1;
E.setZero();
A.setIdentity();

for (double i = 1; i < 1000; i+=1.0) {

    E +=A/factorial;
    A *=u4.Wedge();
    factorial *=i;
}



ASSERT_LE((u4.Exp()-E).norm(), kse3_threshold_ ) << "Error with exp function";

// Log test
ASSERT_LT((se3<DataType>::Log(u4.Exp())-u4.data_).norm(),kse3_threshold_);

// Norm test
ASSERT_EQ(u1.Norm(),u1.data_.norm()) << "Error with norm function";


}

//#########################################################################


TEST(se3Test, OperatorTest) {

typedef float DataType;
typedef Eigen::Matrix<DataType,3,1> Vec3d;
typedef Eigen::Matrix<DataType,6,1> Vec6d;
typedef Eigen::Matrix<DataType,3,3> Mat3d;
typedef Eigen::Matrix<DataType,4,4> Mat4d;
typedef Eigen::Matrix<DataType,6,6> Mat6d;


Vec6d v1 = Vec6d::Random();
Vec6d v2 = Vec6d::Random();
Vec6d v3 = Vec6d::Random();
Vec6d v4 = Vec6d::Random();

se3<DataType> u1(v1);
se3<DataType> u2(v2);
se3<DataType> u3(v3);
se3<DataType> u4(v4);

// Addition test
se3<DataType> u5 = u1+u2;
ASSERT_EQ(u5.data_, u1.data_ + u2.data_) << "Error with addition operator";

// Subtraction test
se3<DataType> u6 = u3-u4;
ASSERT_EQ(u6.data_, u3.data_-u4.data_) << "Error with subtraction operator";

// Scalar multiplication test
se3<DataType> u7 = u1*6;
ASSERT_EQ(u7.data_,u1.data_*6) << "Error with multiplication test";

// Assignment test
se3<DataType> u8 = u4;
ASSERT_EQ(u8.data_,u4.data_) << "Error with assignment test";

}


// These tests are performed in a very specific order. They build
// upon each other. An error in an earlier test can cause failures in
// following tests even if there isn't an error. 
TEST(se3Test, JacobianTest) {

typedef double DataType;
typedef Eigen::Matrix<DataType,3,1> Vec3d;
typedef Eigen::Matrix<DataType,6,1> Vec6d;
typedef Eigen::Matrix<DataType,3,3> Mat3d;
typedef Eigen::Matrix<DataType,4,4> Mat4d;
typedef Eigen::Matrix<DataType,6,6> Mat6d;

Vec6d v1 = Vec6d::Random();
Vec6d v2 = Vec6d::Random();
Vec6d e1;
Vec6d e2;
Vec6d e3;
Vec6d e4;
Vec6d e5;
Vec6d e6;
Mat6d estimated_jr_inv;
e1 << 1,0,0,0,0,0;
e2 << 0,1,0,0,0,0;
e3 << 0,0,1,0,0,0;
e4 << 0,0,0,1,0,0;
e5 << 0,0,0,0,1,0;
e6 << 0,0,0,0,0,1;
double dt = 1e-7;


se3<DataType> u1(v1);
se3<DataType> u2(-v1);
se3<DataType> u3(v2);



estimated_jr_inv.block(0,0,6,1) = (se3<DataType>(se3<DataType>::Log(u1.Exp()*se3<DataType>(e1*dt).Exp())).Vee() - v1)/dt;
estimated_jr_inv.block(0,1,6,1) = (se3<DataType>(se3<DataType>::Log(u1.Exp()*se3<DataType>(e2*dt).Exp())).Vee() - v1)/dt;
estimated_jr_inv.block(0,2,6,1) = (se3<DataType>(se3<DataType>::Log(u1.Exp()*se3<DataType>(e3*dt).Exp())).Vee() - v1)/dt;
estimated_jr_inv.block(0,3,6,1) = (se3<DataType>(se3<DataType>::Log(u1.Exp()*se3<DataType>(e4*dt).Exp())).Vee() - v1)/dt;
estimated_jr_inv.block(0,4,6,1) = (se3<DataType>(se3<DataType>::Log(u1.Exp()*se3<DataType>(e5*dt).Exp())).Vee() - v1)/dt;
estimated_jr_inv.block(0,5,6,1) = (se3<DataType>(se3<DataType>::Log(u1.Exp()*se3<DataType>(e6*dt).Exp())).Vee() - v1)/dt;

// Test JrInverse
ASSERT_LE( (estimated_jr_inv-u1.JrInv()).norm(),kse3_threshold_) << "Error with the right Jacobian";

// Test JlInverse: JlInv(-v) = JrInv(v)
ASSERT_LE( (u2.JlInv()-u1.JrInv()).norm(),kse3_threshold_) << "Error with the right Jacobian";

// Test Jr
ASSERT_LE( (u1.Jr()*u1.JrInv()-Mat6d::Identity()).norm(),kse3_threshold_) << "Error with the right Jacobian";

// Test Jl
ASSERT_LE( (u1.Jl()*u1.JlInv()-Mat6d::Identity()).norm(),kse3_threshold_) << "Error with the right Jacobian";


// Test Remaining


ASSERT_LE( ((u1.Jl(u3)).data_-u1.Jl()*u3.data_).norm(), kse3_threshold_) << "Error with left Jacobian";
ASSERT_LE( ((u1.Jr(u3)).data_-u1.Jr()*u3.data_).norm(), kse3_threshold_) << "Error with right Jacobian";
ASSERT_LE( ((u1.JlInv(u3)).data_-u1.JlInv()*u3.data_).norm(), kse3_threshold_) << "Error with left Jacobian inverse";
ASSERT_LE( ((u1.JrInv(u3)).data_-u1.JrInv()*u3.data_).norm(), kse3_threshold_) << "Error with right Jacobian inverse";


}


}
