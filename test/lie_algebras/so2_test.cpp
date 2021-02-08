#include <Eigen/Dense>

#include "lie_groups/lie_algebras/so2.h"
#include "gtest/gtest.h"

namespace lie_groups {
//#########################################################################
TEST(so2TEST, Constructors) {

typedef float DataType;
typedef Eigen::Matrix<DataType, 1,1> Mat1d;
typedef Eigen::Matrix<DataType, 2,2> Mat2d;

// default constructor
so2<DataType> u1;
ASSERT_EQ(u1.data_(0), 0 ) << "Element not initialized to identity. i.e. zero.";

// assignment constructor Cartesian space
Mat1d v = Mat1d::Random();
so2<DataType> u2(v);
ASSERT_EQ(u2.data_,v) << "Error with the assignment constructor";

// assignment constructor invaid element
Mat2d V = Mat2d::Random();
V = (V+Mat2d(V.transpose())); // Symmetric matrix
if (V.norm() == 0) // If the zero matrix 
    V(0,0) = 1;
so2<DataType> u3(V,true);
ASSERT_NE(u3.data_(0),V(1,0)) << "Error with the assignment constructor";

// assignment constructor valid element
V = V-Mat2d(V.transpose());
so2<DataType> u4(V,true);
ASSERT_EQ(u4.data_(0),V(1,0)) << "Error with the assignment constructor";

// copy constructor
so2<DataType> u5(u2);
ASSERT_EQ(u5.data_,v) << "Error with the copy constructor";

// identity function test
so2<DataType> u6 = so2<DataType>::Identity();
ASSERT_EQ(u6.data_(0),0) << "Error with the identity function";

}

//#########################################################################


// Tests the Bracket, Adjoint, Wedge, Vee, Exp, and Norm functions
TEST(so2TEST, BAWVEN) {

typedef float DataType;
typedef Eigen::Matrix<DataType, 1,1> Mat1d;
typedef Eigen::Matrix<DataType, 2,2> Mat2d;

// Bracket test
so2<DataType> u1;
so2<DataType> u2;
so2<DataType> u3 = u1.Bracket(u2);
ASSERT_EQ(u3.data_(0), 0) << "Error with bracket function";

// Adjoint test
ASSERT_EQ(u1.Adjoint(), Mat2d::Identity()) << "Error with Adjoint function";

// Wedge test
Mat1d v = Mat1d::Random();
so2<DataType> u4(v);
Mat2d V;
V << 0, -v, v, 0;
ASSERT_EQ(u4.Wedge(),V) << "Error with wedge function";

// Vee test
ASSERT_EQ(u4.Vee(),u4.data_) << "Error with vee function";
ASSERT_EQ(so2<DataType>::Vee(u4.Wedge()),u4.data_) << "Error with vee function";

// Exp test
Mat1d v1 = Mat1d::Random();
so2<DataType> u5(v1);
Mat2d g=u5.Exp();

ASSERT_EQ(g.transpose()*g, Mat2d::Identity()) << "Error with exp function";
ASSERT_EQ(g.determinant(), 1) << "Error with exp function";

// Log test
ASSERT_LT((so2<DataType>::Log(u5.Exp())-u5.data_).norm(),kso2_threshold_) << "Error with log function";

// Norm test
Mat1d v2 = Mat1d::Random();
so2<DataType> u6(v2);
ASSERT_EQ(u6.Norm(),v2.norm()) << "Error with norm function";


}

//#########################################################################


TEST(so2Test, OperatorTest) {

typedef float DataType;
typedef Eigen::Matrix<DataType, 1,1> Mat1d;
typedef Eigen::Matrix<DataType, 2,2> Mat2d;

Mat1d v1 = Mat1d::Random();
Mat1d v2 = Mat1d::Random();
Mat1d v3 = Mat1d::Random();
Mat1d v4 = Mat1d::Random();

so2<DataType> u1(v1);
so2<DataType> u2(v2);
so2<DataType> u3(v3);
so2<DataType> u4(v4);

// Addition test
so2<DataType> u5 = u1+u2;
ASSERT_EQ(u5.data_, u1.data_ + u2.data_) << "Error with addition operator";

// Subtraction test
so2<DataType> u6 = u3-u4;
ASSERT_EQ(u6.data_, u3.data_-u4.data_) << "Error with subtraction operator";

// Scalar multiplication test
so2<DataType> u7 = u1*6;
ASSERT_EQ(u7.data_,u1.data_*6) << "Error with multiplication test";

// Assignment test
so2<DataType> u8 = u4;
ASSERT_EQ(u8.data_,u4.data_) << "Error with assignment test";

}

TEST(so2Test, JacobianTest) {

typedef float DataType;
typedef Eigen::Matrix<DataType, 1,1> Mat1d;
typedef Eigen::Matrix<DataType, 2,2> Mat2d;

Mat1d v1 = Mat1d::Random();
Mat1d v2 = Mat1d::Random();
Mat1d v3 = Mat1d::Random();
Mat1d v4 = Mat1d::Random();

so2<DataType> u1(v1);
so2<DataType> u2(v2);
so2<DataType> u3(v3);
so2<DataType> u4(v4);

ASSERT_EQ(u1.Jl(),Mat2d::Identity()) << "Error with left Jacobian";
ASSERT_EQ(u1.Jr(),Mat2d::Identity()) << "Error with right Jacobian";
ASSERT_EQ(u1.JlInv(),Mat2d::Identity()) << "Error with left Jacobian inverse";
ASSERT_EQ(u1.JrInv(),Mat2d::Identity()) << "Error with right Jacobian inverse";

ASSERT_EQ((u1.Jl(u2)).data_,u2.data_) << "Error with left Jacobian";
ASSERT_EQ((u1.Jr(u2)).data_,u2.data_) << "Error with right Jacobian";
ASSERT_EQ((u1.JlInv(u2)).data_,u2.data_) << "Error with left Jacobian inverse";
ASSERT_EQ((u1.JrInv(u2)).data_,u2.data_) << "Error with right Jacobian inverse";


}


}
