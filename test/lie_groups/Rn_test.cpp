#include "lie_groups/lie_groups/Rn.h"
#include "gtest/gtest.h"
#include <Eigen/Dense>


namespace lie_groups {

constexpr int N = 5;
constexpr int kNumTangentSpaces = 2;
constexpr int kTotalNumDimensions = N*kNumTangentSpaces;

template<typename tDataType>
Eigen::Matrix<tDataType,N,1> GenRandElem() {
    Eigen::Matrix<tDataType,N,1> m;
    m.setRandom();
    return m;
}

template<typename tDataType>
Eigen::Matrix<tDataType,N,1> GenElem(Eigen::Matrix<tDataType,N,1> th) {
    return th;
}

////////////////////////////////////////////////////////////////////////

// Test the constructors
TEST(RnTEST, Constructors) {

typedef float DataType;
typedef Eigen::Matrix<DataType,N,1> VecNd;
typedef Eigen::Matrix<DataType,N,N> MatNd;


VecNd Identity;
Identity.setZero();

// Valid element
VecNd data1 =GenRandElem<DataType>();

VecNd data2;
data2.setRandom();

Rn<DataType,N,kNumTangentSpaces> g1;
Rn<DataType,N,kNumTangentSpaces> g2(data1,true);
Rn<DataType,N,kNumTangentSpaces> g3(data2,true);
Rn<DataType,N,kNumTangentSpaces> g4(g2);
Rn<DataType,N,kNumTangentSpaces> g5(data2);

ASSERT_EQ(g1.data_,Identity) << "Default constructor not set to identity";
ASSERT_EQ(g2.data_,data1) << "Assignment constructor error";
ASSERT_EQ(g3.data_,data2) << "Copy constructor constructor error";
ASSERT_EQ(g4.data_,data1) << "Assignment constructor error: Invalid element not excepted. Should set element to identity.";
ASSERT_EQ(g5.data_,data2) << "Assignment constructor error";


}


////////////////////////////////////////////////////////////////////////

// Test the Inverse, Adjoint, Identity, Log, and operators
TEST(RnTEST, IAILO) {

typedef float DataType;
typedef Eigen::Matrix<DataType,N,1> VecNd;
typedef Eigen::Matrix<DataType,N,N> MatNd;


VecNd Identity;
Identity.setZero();
MatNd I;
I.setIdentity();
VecNd th;
th.setRandom();
VecNd g = GenElem(th);

Rn<DataType,N,kNumTangentSpaces> g1(GenRandElem<DataType>());
Rn<DataType,N,kNumTangentSpaces> g2(GenRandElem<DataType>());
Rn<DataType,N,kNumTangentSpaces> g3(g);
Rn<DataType,N,kNumTangentSpaces> g4 = g1;
Rn<DataType,N,kNumTangentSpaces> g5 = g1*g2;

ASSERT_LE( ((g1.Inverse()).data_ + g1.data_).norm(), kRn_threshold_) << " Error with the inverse operation ";
ASSERT_EQ( (Rn<DataType,N>::Identity()).data_, Identity  ) << "Error with identity function ";
ASSERT_EQ( g2.Adjoint(),  I) << "Error with the Adjoint operation";

ASSERT_EQ(th, g3.Log()) << "Error with the log function";

ASSERT_EQ(g4.data_, g1.data_) << "Error with assignmet operator";

ASSERT_EQ(g5.data_, g1.data_ + g2.data_) << "Error with the group operator";

}

////////////////////////////////////////////////////////////////////////

// Tests the box plus and box minus functions
TEST(RnTest, BoxPlus) {

typedef float DataType;
typedef Eigen::Matrix<DataType,N,1> VecNd;
typedef Eigen::Matrix<DataType,N*kNumTangentSpaces,1> VecTNd;
typedef Eigen::Matrix<DataType,N,N> MatNd;

VecNd th1;
VecNd th2;
th1.setRandom();
th2.setRandom();

VecTNd Tth1;
VecTNd Tth2;
Tth1.setZero();
Tth1.block(0,0,N,1) = th1;
Tth2.setZero();
Tth2.block(0,0,N,1) = th2;


VecNd data1 = GenElem<DataType>(th1);
VecNd data2 = GenElem<DataType>(th2);
VecNd data3 = data1+data2;

Rn<DataType,N,kNumTangentSpaces> g1(data1);
g1.OPlusEq(Tth2);
Rn<DataType,N,kNumTangentSpaces> g2(data1);
Rn<DataType,N,kNumTangentSpaces> g3(data1);
g3.BoxPlusEq(rn<DataType,N,kNumTangentSpaces>(Tth2));
Rn<DataType,N,kNumTangentSpaces> g4(data2);
Rn<DataType,N,kNumTangentSpaces> g5(data3);

ASSERT_EQ(g1.data_, data1+data2) << "Error with an OPlus function";

// std::cerr << Rn<DataType,N>::BoxPlus(data1,Th2)<< std::endl;
// std::cerr << data1*data2<< std::endl;
typedef Rn<DataType,N,kNumTangentSpaces> Group;

ASSERT_EQ(Group::BoxPlus(data1,Tth2), data3) << "Error with the BoxPlus function";
ASSERT_EQ(g2.BoxPlus(Tth2), data1 + data2) << "Error with the BoxPlus function";

g2.BoxPlusEq(Tth2);
ASSERT_EQ(g2.data_, data1 + data2) << "Error with the BoxPlus function";

ASSERT_EQ(g3.data_, data1+data2) << "Error with the BoxPlus function";

// std::cerr << g4.OMinus(data3)<< std::endl;
// std::cerr << th1<< std::endl;
ASSERT_LE( (g5.OMinus(data2)- Tth1).norm(), kRn_threshold_) << "Error with the OPlus function";

ASSERT_LE( (g5.BoxMinus(g4).data_- Tth1).norm(), kRn_threshold_) << "Error with the OPlus function";

}


}