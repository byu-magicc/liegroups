#include <Eigen/Dense>

#include "lie_algebras/rn.h"
#include "gtest/gtest.h"

namespace lie_groups {

constexpr int kN = 5;
typedef float tDataType;

//#########################################################################
TEST(rnTEST, Constructors) {
// default constructor
rn<tDataType,kN> u1;
Eigen::Matrix<tDataType,kN,1> Identity;
Identity.setZero();
ASSERT_EQ(u1.data_, Identity ) << "Element not initialized to identity. i.e. zero.";

// assignment constructor Cartesian space
Eigen::Matrix<tDataType,kN,1> v = Eigen::Matrix<tDataType,kN,1>::Random();
rn<tDataType,kN> u2(v);
ASSERT_EQ(u2.data_,v) << "Error with the assignment constructor";

Eigen::Matrix<int,1,1> n;
n.setRandom();
int nn = std::min(abs(n(0)),1000);
for ( unsigned long int i = 0; i< nn; ++i) {
    rn<tDataType,kN> u;
}

// assignment constructor valid element
Eigen::Matrix<tDataType,kN,1> V;
V.setRandom();
rn<tDataType,kN> u4(V,true);
ASSERT_EQ(u4.data_,V) << "Error with the assignment constructor";
// copy constructor
rn<tDataType,kN> u5(u2);
ASSERT_EQ(u5.data_,v) << "Error with the copy constructor";

// identity function test
rn<tDataType,kN> u6 = rn<tDataType,kN>::Identity();
ASSERT_EQ(u6.data_,Identity) << "Error with the identity function";

}

//#########################################################################


// Tests the Bracket, Adjoint, Wedge, Vee, Exp, and Norm functions
TEST(rnTEST, BAWVEN) {

Eigen::Matrix<tDataType,kN,1> Identity;
Identity.setZero();

Eigen::Matrix<tDataType,kN,kN> I = Eigen::Matrix<tDataType,kN,kN>::Identity();

// Bracket test
rn<tDataType,kN> u1;
rn<tDataType,kN> u2;
rn<tDataType,kN> u3 = u1.Bracket(u2);
ASSERT_EQ(u3.data_, Identity) << "Error with bracket function";

// Adjoint test
ASSERT_EQ(u1.Adjoint(), I ) << "Error with Adjoint function";

// Wedge test
Eigen::Matrix<tDataType,kN,1> v = Eigen::Matrix<tDataType,kN,1>::Random();
rn<tDataType,kN> u4(v);
ASSERT_EQ(u4.Wedge(),v) << "Error with wedge function";

typedef rn<tDataType,kN> group;

// Vee test
ASSERT_EQ(u4.Vee(),u4.data_) << "Error with vee function";
ASSERT_EQ(group::Vee(u4.Wedge()),u4.data_) << "Error with vee function";

// Exp test
Eigen::Matrix<tDataType,kN,1> v1 = Eigen::Matrix<tDataType,kN,1>::Random();
rn<tDataType,kN> u5(v1);

ASSERT_EQ(u5.Exp(), v1) << "Error with exp function";

// Log test
ASSERT_LT((rn<tDataType,kN>::Log(u5.Exp())-u5.data_).norm(),krn_threshold_) << "Error with log function";

// Norm test
Eigen::Matrix<tDataType,kN,1> v2 = Eigen::Matrix<tDataType,kN,1>::Random();
rn<tDataType,kN> u6(v2);
ASSERT_EQ(u6.Norm(),v2.norm()) << "Error with norm function";


}

//#########################################################################


TEST(rnTest, OperatorTest) {

Eigen::Matrix<tDataType,kN,1> v1 = Eigen::Matrix<tDataType,kN,1>::Random();
Eigen::Matrix<tDataType,kN,1> v2 = Eigen::Matrix<tDataType,kN,1>::Random();
Eigen::Matrix<tDataType,kN,1> v3 = Eigen::Matrix<tDataType,kN,1>::Random();
Eigen::Matrix<tDataType,kN,1> v4 = Eigen::Matrix<tDataType,kN,1>::Random();

rn<tDataType,kN> u1(v1);
rn<tDataType,kN> u2(v2);
rn<tDataType,kN> u3(v3);
rn<tDataType,kN> u4(v4);

// Addition test
rn<tDataType,kN> u5 = u1+u2;
ASSERT_EQ(u5.data_, u1.data_ + u2.data_) << "Error with addition operator";

// Subtraction test
rn<tDataType,kN> u6 = u3-u4;
ASSERT_EQ(u6.data_, u3.data_-u4.data_) << "Error with subtraction operator";

// Scalar multiplication test
rn<tDataType,kN> u7 = u1*6;
ASSERT_EQ(u7.data_,u1.data_*6) << "Error with multiplication test";

// Assignment test
rn<tDataType,kN> u8 = u4;
ASSERT_EQ(u8.data_,u4.data_) << "Error with assignment test";

}


TEST(rnTest, JacobianTest) {

Eigen::Matrix<tDataType,kN,1> v1 = Eigen::Matrix<tDataType,kN,1>::Random();
Eigen::Matrix<tDataType,kN,1> v2 = Eigen::Matrix<tDataType,kN,1>::Random();
Eigen::Matrix<tDataType,kN,1> v3 = Eigen::Matrix<tDataType,kN,1>::Random();
Eigen::Matrix<tDataType,kN,1> v4 = Eigen::Matrix<tDataType,kN,1>::Random();

rn<tDataType,kN> u1(v1);
rn<tDataType,kN> u2(v2);
rn<tDataType,kN> u3(v3);
rn<tDataType,kN> u4(v4);

Eigen::Matrix<tDataType,kN,kN> I = Eigen::Matrix<tDataType,kN,kN>::Identity();

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
