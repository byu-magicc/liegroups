#include <Eigen/Dense>

#include "lie_algebras/so3.h"
#include "gtest/gtest.h"

namespace lie_groups {

TEST(so3TEST, Constructors) {

typedef float DataType;
typedef Eigen::Matrix<DataType,3,1> Vec3d;
typedef Eigen::Matrix<DataType,3,3> Mat3d;

// default constructor
so3<DataType> u1;
Vec3d Identity = Vec3d::Zero();
ASSERT_EQ(u1.data_, Identity ) << "Element not initialized to identity. i.e. zero.";




// assignment constructor Cartesian space
Vec3d v = Vec3d::Random();
so3<DataType> u2(v);
ASSERT_EQ(u2.data_,v) << "Error with the assignment constructor";

// assignment constructor invaid element
Mat3d V1 = Mat3d::Random();
V1 = (V1+Mat3d(V1.transpose())); // Symmetric matrix
if (V1.norm() == 0) // If the zero matrix 
    V1(0,0) = 1;
so3<DataType> u3(V1,true);
ASSERT_EQ(u3.data_, Identity) << "Error with the assignment constructor";


// assignment constructor valid element
Mat3d V2 = Mat3d::Random();
Mat3d tmp = V2.transpose() - V2; // Create a skew symmetric matrix
V2 = tmp;
Vec3d v2;
v2 << V2(2,1), V2(0,2), V2(1,0);
so3<DataType> u4(V2,true);
ASSERT_EQ(u4.data_,v2) << "Error with the assignment constructor";


// copy constructor
so3<DataType> u5(u4);
ASSERT_EQ(u5.data_,u4.data_) << "Error with the copy constructor";


// identity function test
so3<DataType> u6 = so3<DataType>::Identity();
ASSERT_EQ(u6.data_,Identity) << "Error with the identity function";

}

// Tests the Bracket, Adjoint, Wedge, Vee, Exp, and Norm functions
TEST(so3TEST, BAWVEN) {

typedef double DataType;
typedef Eigen::Matrix<DataType,3,1> Vec3d;
typedef Eigen::Matrix<DataType,3,3> Mat3d;

so3<DataType> u1(Vec3d::Random());
so3<DataType> u2(Vec3d::Random());
so3<DataType> u3(Vec3d::Random());
so3<DataType> u4(Vec3d::Random());


// Wedge test
Mat3d m1;
m1 << 0, -u1.data_(2), u1.data_(1), u1.data_(2), 0, -u1.data_(0), -u1.data_(1), u1.data_(0),0;
ASSERT_EQ(u1.Wedge(),m1) << "Error with wedge function";


// Adjoint test
Mat3d m2;
m2 << 0, -u2.data_(2), u2.data_(1), u2.data_(2), 0, -u2.data_(0), -u2.data_(1), u2.data_(0),0;
ASSERT_EQ(u2.Adjoint(),m2) << "Error with adjoint test"; 

// Bracket and adjoint test
Mat3d m3;
m3 = u3.Wedge()*u2.Wedge() - u2.Wedge()*u3.Wedge();
ASSERT_EQ( (u3.Bracket(u2)).Wedge(), m3) << "Error with bracket or Adjoint function";


// Vee test
ASSERT_EQ(u4.Vee(),u4.data_) << "Error with vee function";
ASSERT_EQ(so3<DataType>::Vee(u4.Wedge()),u4.data_) << "Error with vee function";


// Exp test
Mat3d E;
Mat3d A;
DataType factorial=1;
E.setZero();
A.setIdentity();
Mat3d W = u4.Wedge();

for (long int i = 1; i < 1000; ++i) {

    E += A/factorial;
    A *= W;
    factorial *= i;
}
ASSERT_LE((u4.Exp()-E).norm(), 1e-10 ) << "Error with exp function";

// Log test
ASSERT_LT((so3<DataType>::Log(u4.Exp())-u4.data_).norm(),1e-10);

// Norm test
ASSERT_EQ(u1.Norm(),u1.data_.norm()) << "Error with norm function";


}

TEST(so3Test, OperatorTest) {

typedef float DataType;
typedef Eigen::Matrix<DataType,3,1> Vec3d;
typedef Eigen::Matrix<DataType,3,3> Mat3d;

Vec3d v1 = Vec3d::Random();
Vec3d v2 = Vec3d::Random();
Vec3d v3 = Vec3d::Random();
Vec3d v4 = Vec3d::Random();

so3<DataType> u1(v1);
so3<DataType> u2(v2);
so3<DataType> u3(v3);
so3<DataType> u4(v4);

// Addition test
so3<DataType> u5 = u1+u2;
ASSERT_EQ(u5.data_, u1.data_ + u2.data_) << "Error with addition operator";

// Subtraction test
so3<DataType> u6 = u3-u4;
ASSERT_EQ(u6.data_, u3.data_-u4.data_) << "Error with subtraction operator";

// Scalar multiplication test
so3<DataType> u7 = u1*6;
ASSERT_EQ(u7.data_,u1.data_*6) << "Error with multiplication test";

// Assignment test
so3<DataType> u8 = u4;
ASSERT_EQ(u8.data_,u4.data_) << "Error with assignment test";

}


// These tests are performed in a very specific order. They build
// upon each other. An error in an earlier test can cause failures in
// following tests even if there isn't an error. 
TEST(so3Test, JacobianTest) {

typedef double DataType;
typedef Eigen::Matrix<DataType,3,1> Vec3d;
typedef Eigen::Matrix<DataType,3,3> Mat3d;

Vec3d v1 = Vec3d::Random();
Vec3d v2 = Vec3d::Random();
Vec3d e1;
Vec3d e2;
Vec3d e3;
Mat3d estimated_jr_inv;
e1 << 1,0,0;
e2 << 0,1,0;
e3 << 0,0,1;
double dt = 1e-7;


so3<DataType> u1(v1);
so3<DataType> u2(-v1);
so3<DataType> u3(v2);



estimated_jr_inv.block(0,0,3,1) = (so3<DataType>(so3<DataType>::Log(u1.Exp()*so3<DataType>(e1*dt).Exp())).Vee() - v1)/dt;
estimated_jr_inv.block(0,1,3,1) = (so3<DataType>(so3<DataType>::Log(u1.Exp()*so3<DataType>(e2*dt).Exp())).Vee() - v1)/dt;
estimated_jr_inv.block(0,2,3,1) = (so3<DataType>(so3<DataType>::Log(u1.Exp()*so3<DataType>(e3*dt).Exp())).Vee() - v1)/dt;

// Test JrInverse
ASSERT_LE( (estimated_jr_inv-u1.JrInv()).norm(),1e-7) << "Error with the right Jacobian";

// Test JlInverse: JlInv(-v) = JrInv(v)
ASSERT_LE( (u2.JlInv()-u1.JrInv()).norm(),1e-8) << "Error with the right Jacobian";

// Test Jr
ASSERT_LE( (u1.Jr()*u1.JrInv()-Mat3d::Identity()).norm(),1e-8) << "Error with the right Jacobian";

// Test Jl
ASSERT_LE( (u1.Jl()*u1.JlInv()-Mat3d::Identity()).norm(),1e-8) << "Error with the right Jacobian";


// Test Remaining


ASSERT_LE( ((u1.Jl(u3)).data_-u1.Jl()*u3.data_).norm(), 1e-8) << "Error with left Jacobian";
ASSERT_LE( ((u1.Jr(u3)).data_-u1.Jr()*u3.data_).norm(), 1e-8) << "Error with right Jacobian";
ASSERT_LE( ((u1.JlInv(u3)).data_-u1.JlInv()*u3.data_).norm(), 1e-8) << "Error with left Jacobian inverse";
ASSERT_LE( ((u1.JrInv(u3)).data_-u1.JrInv()*u3.data_).norm(), 1e-8) << "Error with right Jacobian inverse";


}


}
