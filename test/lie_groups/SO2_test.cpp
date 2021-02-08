#include "lie_groups/lie_groups/SO2.h"
#include "gtest/gtest.h"
#include <Eigen/Dense>


namespace lie_groups {

template<typename tDataType>
Eigen::Matrix<tDataType,2,2> GenRandElem() {
    Eigen::Matrix<tDataType,1,1> u;
    u.setRandom();
    Eigen::Matrix<tDataType,2,2> m;
    m = so2<tDataType>::Exp(u);
    return m;
}

template<typename tDataType>
Eigen::Matrix<tDataType,2,2> GenElem(Eigen::Matrix<tDataType,1,1> th) {
    Eigen::Matrix<tDataType,2,2> m;
    m = so2<tDataType>::Exp(th);
    return m;
}

////////////////////////////////////////////////////////////////////////

// Test the constructors
TEST(SO2TEST, Constructors) {

typedef float tDataType;

Eigen::Matrix<tDataType,2,2> Identity;
Identity.setIdentity();

// Valid element
Eigen::Matrix<tDataType,2,2> data1 =GenRandElem<tDataType>();

// Invalid element
Eigen::Matrix<tDataType,2,2> data2;
data2.setRandom();
while  ((data2.transpose()*data2 - Identity).norm() <= kSO2_threshold_ ) {
    data2.setRandom();
}



SO2<tDataType> g1;
SO2<tDataType> g2(data1,true);
SO2<tDataType> g3(data2,true);
SO2<tDataType> g4(g2);
SO2<tDataType> g5(data2);

ASSERT_EQ(g1.data_,Identity) << "Default constructor not set to identity";
ASSERT_EQ(g2.data_,data1) << "Assignment constructor error";
ASSERT_EQ(g3.data_,Identity) << "Copy constructor constructor error";
ASSERT_EQ(g4.data_,data1) << "Assignment constructor error: Invalid element not excepted. Should set element to identity.";
ASSERT_EQ(g5.data_,data2) << "Assignment constructor error";


}


////////////////////////////////////////////////////////////////////////

// Test the Inverse, Adjoint, Identity, Log, and operators
TEST(SO2TEST, IAILO) {

typedef float tDataType;

Eigen::Matrix<tDataType,2,2> Identity;
Identity.setIdentity();
Eigen::Matrix<tDataType,1,1> I;
I << 1;
Eigen::Matrix<tDataType,1,1> th;
th << 0.1;
Eigen::Matrix<tDataType,2,2> g = GenElem(th);

SO2<tDataType> g1(GenRandElem<tDataType>());
SO2<tDataType> g2(GenRandElem<tDataType>());
SO2<tDataType> g3(g);
SO2<tDataType> g4 = g1;
SO2<tDataType> g5 = g1*g2;

ASSERT_LE( ((g1.Inverse()).data_ - g1.data_.inverse()).norm(), kSO2_threshold_) << " Error with the inverse operation ";
ASSERT_EQ( SO2<tDataType>::Identity().data_, Identity  ) << "Error with identity function ";
ASSERT_EQ( g2.Adjoint(),  I) << "Error with the Adjoint operation";

ASSERT_DOUBLE_EQ(th(0), g3.Log()(0)) << "Error with the log function";

ASSERT_EQ(g4.data_, g1.data_) << "Error with assignmet operator";

ASSERT_EQ(g5.data_, g1.data_*g2.data_) << "Error with the group operator";

}

////////////////////////////////////////////////////////////////////////

// Tests the box plus and box minus functions
TEST(SO2Test, BoxPlus) {

typedef float tDataType;

Eigen::Matrix<tDataType,1,1> th1;
Eigen::Matrix<tDataType,1,1> th2;
th1 << 0.1;
th2 << 0.2;
Eigen::Matrix<tDataType,2,2> Th1 = so2<tDataType>::Wedge(th1);
Eigen::Matrix<tDataType,2,2> Th2 = so2<tDataType>::Wedge(th2);



Eigen::Matrix<tDataType,2,2> data1 = GenElem<tDataType>(th1);
Eigen::Matrix<tDataType,2,2> data2 = GenElem<tDataType>(th2);
Eigen::Matrix<tDataType,2,2> data3 = data1*data2;

SO2<tDataType> g1(data1);
g1.OPlusEq(th2);
SO2<tDataType> g2(data1);
SO2<tDataType> g3(data1);
g3.BoxPlusEq(so2<tDataType>(th2));
SO2<tDataType> g4(data1);
SO2<tDataType> g5(data3);

ASSERT_EQ(g1.data_, data1*data2) << "Error with an OPlus function";

// std::cerr << SO2::BoxPlus(data1,Th2)<< std::endl;
// std::cerr << data1*data2<< std::endl;

ASSERT_EQ(SO2<tDataType>::BoxPlus(data1,Th2), data1*data2) << "Error with the BoxPlus function";
ASSERT_EQ(g2.BoxPlus(Th2), data1*data2) << "Error with the BoxPlus function";

g2.BoxPlusEq(Th2);
ASSERT_EQ(g2.data_, data1*data2) << "Error with the BoxPlus function";

ASSERT_EQ(g3.data_, data1*data2) << "Error with the BoxPlus function";

// std::cerr << g4.OMinus(data3)<< std::endl;
// std::cerr << th1<< std::endl;
ASSERT_LE( (g5.OMinus(data1)- th2).norm(), kSO2_threshold_) << "Error with the OPlus function";

ASSERT_LE( (g5.BoxMinus(g4).data_- th2).norm(), kSO2_threshold_) << "Error with the OPlus function";

}

}