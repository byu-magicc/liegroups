#include "lie_groups/Rn.h"
#include "gtest/gtest.h"
#include <Eigen/Dense>


namespace lie_groups {

constexpr int N = 5;

Eigen::Matrix<double,N,1> GenRandElem() {
    Eigen::Matrix<double,N,1> m;
    m.setRandom();
    return m;
}

Eigen::Matrix<double,N,1> GenElem(Eigen::Matrix<double,N,1> th) {
    return th;
}

////////////////////////////////////////////////////////////////////////

// Test the constructors
TEST(RnTEST, Constructors) {

Eigen::Matrix<double,N,1> Identity;
Identity.setZero();

// Valid element
Eigen::Matrix<double,N,1> data1 =GenRandElem();

Eigen::Matrix<double,N,1> data2;
data2.setRandom();

Rn<N> g1;
Rn<N> g2(data1,true);
Rn<N> g3(data2,true);
Rn<N> g4(g2);
Rn<N> g5(data2);

ASSERT_EQ(g1.data_,Identity) << "Default constructor not set to identity";
ASSERT_EQ(g2.data_,data1) << "Assignment constructor error";
ASSERT_EQ(g3.data_,data2) << "Copy constructor constructor error";
ASSERT_EQ(g4.data_,data1) << "Assignment constructor error: Invalid element not excepted. Should set element to identity.";
ASSERT_EQ(g5.data_,data2) << "Assignment constructor error";


}


////////////////////////////////////////////////////////////////////////

// Test the Inverse, Adjoint, Identity, Log, and operators
TEST(RnTEST, IAILO) {

Eigen::Matrix<double,N,1> Identity;
Identity.setZero();
Eigen::Matrix<double,N,N> I;
I.setIdentity();
Eigen::Matrix<double,N,1> th;
th.setRandom();
Eigen::Matrix<double,N,1> g = GenElem(th);

Rn<N> g1(GenRandElem());
Rn<N> g2(GenRandElem());
Rn<N> g3(g);
Rn<N> g4 = g1;
Rn<N> g5 = g1*g2;

ASSERT_LE( ((g1.Inverse()).data_ + g1.data_).norm(), kRn_threshold_) << " Error with the inverse operation ";
ASSERT_EQ( Rn<N>::Identity().data_, Identity  ) << "Error with identity function ";
ASSERT_EQ( g2.Adjoint(),  I) << "Error with the Adjoint operation";

ASSERT_EQ(th, g3.Log()) << "Error with the log function";

ASSERT_EQ(g4.data_, g1.data_) << "Error with assignmet operator";

ASSERT_EQ(g5.data_, g1.data_ + g2.data_) << "Error with the group operator";

}

////////////////////////////////////////////////////////////////////////

// Tests the box plus and box minus functions
TEST(RnTest, BoxPlus) {

Eigen::Matrix<double,N,1> th1;
Eigen::Matrix<double,N,1> th2;
th1.setRandom();
th2.setRandom();



Eigen::Matrix<double,N,1> data1 = GenElem(th1);
Eigen::Matrix<double,N,1> data2 = GenElem(th2);
Eigen::Matrix<double,N,1> data3 = data1+data2;

Rn<N> g1(data1);
g1.OPlusEq(th2);
Rn<N> g2(data1);
Rn<N> g3(data1);
g3.BoxPlusEq(rn<N>(th2));
Rn<N> g4(data2);
Rn<N> g5(data3);

ASSERT_EQ(g1.data_, data1+data2) << "Error with an OPlus function";

// std::cerr << Rn<N>::BoxPlus(data1,Th2)<< std::endl;
// std::cerr << data1*data2<< std::endl;

ASSERT_EQ(Rn<N>::BoxPlus(data1,th2), data1 + data2) << "Error with the BoxPlus function";
ASSERT_EQ(g2.BoxPlus(th2), data1 + data2) << "Error with the BoxPlus function";

g2.BoxPlusEq(th2);
ASSERT_EQ(g2.data_, data1 + data2) << "Error with the BoxPlus function";

ASSERT_EQ(g3.data_, data1+data2) << "Error with the BoxPlus function";

// std::cerr << g4.OMinus(data3)<< std::endl;
// std::cerr << th1<< std::endl;
ASSERT_LE( (g4.OMinus(data3)- th1).norm(), kRn_threshold_) << "Error with the OPlus function";

ASSERT_LE( (g4.BoxMinus(g5).data_- th1).norm(), kRn_threshold_) << "Error with the OPlus function";

}


}