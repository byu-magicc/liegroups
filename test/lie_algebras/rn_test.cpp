#include <Eigen/Dense>

#include "gtest/gtest.h"

#include "lie_groups/lie_algebras/rn.h"

namespace lie_groups {

constexpr int kN = 5;
constexpr int NumTangentSpaces = 2;
constexpr int TotalDimensions = kN*NumTangentSpaces;
typedef float tDataType;

//#########################################################################
TEST(rnTEST, Constructors) {
// default constructor
rn<tDataType,kN,NumTangentSpaces> u1;
Eigen::Matrix<tDataType,TotalDimensions,1> Identity;
Identity.setZero();
ASSERT_EQ(u1.data_, Identity ) << "Element not initialized to identity. i.e. zero.";

// assignment constructor Cartesian space
Eigen::Matrix<tDataType,TotalDimensions,1> v = Eigen::Matrix<tDataType,TotalDimensions,1>::Random();
rn<tDataType,kN,NumTangentSpaces> u2(v);
ASSERT_EQ(u2.data_,v) << "Error with the assignment constructor";

Eigen::Matrix<int,1,1> n;
n.setRandom();
int nn = std::min(abs(n(0)),1000);
for ( unsigned long int i = 0; i< nn; ++i) {
    rn<tDataType,kN,NumTangentSpaces> u;
}

// assignment constructor valid element
Eigen::Matrix<tDataType,TotalDimensions,1> V;
V.setRandom();
rn<tDataType,kN,NumTangentSpaces> u4(V,true);
ASSERT_EQ(u4.data_,V) << "Error with the assignment constructor";
// copy constructor
rn<tDataType,kN,NumTangentSpaces> u5(u2);
ASSERT_EQ(u5.data_,v) << "Error with the copy constructor";

// identity function test
rn<tDataType,kN,NumTangentSpaces> u6 = rn<tDataType,kN,NumTangentSpaces>::Identity();
ASSERT_EQ(u6.data_,Identity) << "Error with the identity function";

}

//#########################################################################


// Tests the Bracket, Adjoint, Wedge, Vee, Exp, and Norm functions
TEST(rnTEST, BAWVEN) {

Eigen::Matrix<tDataType,TotalDimensions,1> Identity;
Identity.setZero();

Eigen::Matrix<tDataType,TotalDimensions,TotalDimensions> I = Eigen::Matrix<tDataType,TotalDimensions,TotalDimensions>::Identity();

// Bracket test
rn<tDataType,TotalDimensions> u1;
rn<tDataType,TotalDimensions> u2;
rn<tDataType,TotalDimensions> u3 = u1.Bracket(u2);
ASSERT_EQ(u3.data_, Identity) << "Error with bracket function";

// Adjoint test
ASSERT_EQ(u1.Adjoint(), I ) << "Error with Adjoint function";

// Wedge test
Eigen::Matrix<tDataType,TotalDimensions,1> v = Eigen::Matrix<tDataType,TotalDimensions,1>::Random();
rn<tDataType,kN,NumTangentSpaces> u4(v);
ASSERT_EQ(u4.Wedge(),v) << "Error with wedge function";

typedef rn<tDataType,kN,NumTangentSpaces> group;

// Vee test
ASSERT_EQ(u4.Vee(),u4.data_) << "Error with vee function";
ASSERT_EQ(group::Vee(u4.Wedge()),u4.data_) << "Error with vee function";

// Exp test
Eigen::Matrix<tDataType,TotalDimensions,1> v1 = Eigen::Matrix<tDataType,TotalDimensions,1>::Random();
rn<tDataType,kN,NumTangentSpaces> u5(v1);

ASSERT_EQ(u5.Exp(), v1) << "Error with exp function";

// Log test
ASSERT_LT((rn<tDataType,kN,NumTangentSpaces>::Log(u5.Exp())-u5.data_).norm(),krn_threshold_) << "Error with log function";

// Norm test
Eigen::Matrix<tDataType,TotalDimensions,1> v2 = Eigen::Matrix<tDataType,TotalDimensions,1>::Random();
rn<tDataType,kN,NumTangentSpaces> u6(v2);
ASSERT_EQ(u6.Norm(),v2.norm()) << "Error with norm function";


}

//#########################################################################


TEST(rnTest, OperatorTest) {

Eigen::Matrix<tDataType,TotalDimensions,1> v1 = Eigen::Matrix<tDataType,TotalDimensions,1>::Random();
Eigen::Matrix<tDataType,TotalDimensions,1> v2 = Eigen::Matrix<tDataType,TotalDimensions,1>::Random();
Eigen::Matrix<tDataType,TotalDimensions,1> v3 = Eigen::Matrix<tDataType,TotalDimensions,1>::Random();
Eigen::Matrix<tDataType,TotalDimensions,1> v4 = Eigen::Matrix<tDataType,TotalDimensions,1>::Random();

rn<tDataType,kN,NumTangentSpaces> u1(v1);
rn<tDataType,kN,NumTangentSpaces> u2(v2);
rn<tDataType,kN,NumTangentSpaces> u3(v3);
rn<tDataType,kN,NumTangentSpaces> u4(v4);

// Addition test
rn<tDataType,kN,NumTangentSpaces> u5 = u1+u2;
ASSERT_EQ(u5.data_, u1.data_ + u2.data_) << "Error with addition operator";

// Subtraction test
rn<tDataType,kN,NumTangentSpaces> u6 = u3-u4;
ASSERT_EQ(u6.data_, u3.data_-u4.data_) << "Error with subtraction operator";

// Scalar multiplication test
rn<tDataType,kN,NumTangentSpaces> u7 = u1*6;
ASSERT_EQ(u7.data_,u1.data_*6) << "Error with multiplication test";

// Assignment test
rn<tDataType,kN,NumTangentSpaces> u8 = u4;
ASSERT_EQ(u8.data_,u4.data_) << "Error with assignment test";

}


TEST(rnTest, JacobianTest) {

Eigen::Matrix<tDataType,TotalDimensions,1> v1 = Eigen::Matrix<tDataType,TotalDimensions,1>::Random();
Eigen::Matrix<tDataType,TotalDimensions,1> v2 = Eigen::Matrix<tDataType,TotalDimensions,1>::Random();
Eigen::Matrix<tDataType,TotalDimensions,1> v3 = Eigen::Matrix<tDataType,TotalDimensions,1>::Random();
Eigen::Matrix<tDataType,TotalDimensions,1> v4 = Eigen::Matrix<tDataType,TotalDimensions,1>::Random();

rn<tDataType,kN,NumTangentSpaces> u1(v1);
rn<tDataType,kN,NumTangentSpaces> u2(v2);
rn<tDataType,kN,NumTangentSpaces> u3(v3);
rn<tDataType,kN,NumTangentSpaces> u4(v4);

Eigen::Matrix<tDataType,TotalDimensions,TotalDimensions> I = Eigen::Matrix<tDataType,TotalDimensions,TotalDimensions>::Identity();

ASSERT_EQ(u1.Jl(),I) << "Error with left Jacobian";
ASSERT_EQ(u1.Jr(),I) << "Error with right Jacobian";
ASSERT_EQ(u1.JlInv(),I) << "Error with left Jacobian inverse";
ASSERT_EQ(u1.JrInv(),I) << "Error with right Jacobian inverse";

ASSERT_EQ((u1.Jl(u2)).data_,u2.data_) << "Error with left Jacobian";
ASSERT_EQ((u1.Jr(u2)).data_,u2.data_) << "Error with right Jacobian";
ASSERT_EQ((u1.JlInv(u2)).data_,u2.data_) << "Error with left Jacobian inverse";
ASSERT_EQ((u1.JrInv(u2)).data_,u2.data_) << "Error with right Jacobian inverse";


}


}
