#include "lie_groups/SO3.h"
#include "gtest/gtest.h"
#include <Eigen/Dense>


namespace lie_groups {

template<typename tDataType>
Eigen::Matrix<tDataType,3,3> GenRandElem() {
    Eigen::Matrix<tDataType,3,1> u;
    u.setRandom();
    Eigen::Matrix<tDataType,3,3> m;
    m = so3<tDataType>::Exp(u);
    return m;
}

template<typename tDataType>
Eigen::Matrix<tDataType,3,3> GenElem(Eigen::Matrix<tDataType,3,1> th) {
    Eigen::Matrix<tDataType,3,3> m;
    m = so3<tDataType>::Exp(th);
    return m;
}

////////////////////////////////////////////////////////////////////////

// Test the constructors
TEST(SO3TEST, Constructors) {

typedef float tDataType;

Eigen::Matrix<tDataType,3,3> Identity;
Identity.setIdentity();

// Valid element
Eigen::Matrix<tDataType,3,3> data1 =GenRandElem<tDataType>();

// Invalid element
Eigen::Matrix<tDataType,3,3> data2;
data2.setRandom();
while  ((data2.transpose()*data2 - Identity).norm() <= kSO3_threshold_ ) {
    data2.setRandom();
}



SO3<tDataType> g1;
SO3<tDataType> g2(data1,true);
SO3<tDataType> g3(data2,true);
SO3<tDataType> g4(g2);
SO3<tDataType> g5(data2);

ASSERT_EQ(g1.data_,Identity) << "Default constructor not set to identity";
ASSERT_EQ(g2.data_,data1) << "Assignment constructor error";
ASSERT_EQ(g3.data_,Identity) << "Copy constructor constructor error";
ASSERT_EQ(g4.data_,data1) << "Assignment constructor error: Invalid element not excepted. Should set element to identity.";
ASSERT_EQ(g5.data_,data2) << "Assignment constructor error";


}


////////////////////////////////////////////////////////////////////////

// Test the Inverse, Adjoint, Identity, Log, and operators
TEST(SO3TEST, IAILO) {

typedef float tDataType;
typedef Eigen::Matrix<tDataType,3,1> Vec3d;
typedef Eigen::Matrix<tDataType,3,3> Mat3d;

Mat3d Identity;
Identity.setIdentity();
Vec3d th;
th << 0.1,0.2,0.3;
Mat3d g = GenElem(th);

SO3<tDataType> g1(GenRandElem<tDataType>());
SO3<tDataType> g2(GenRandElem<tDataType>());
SO3<tDataType> g3(g);
SO3<tDataType> g4 = g1;
SO3<tDataType> g5 = g1*g2;

ASSERT_LE( ((g1.Inverse()).data_ - g1.data_.inverse()).norm(), kSO3_threshold_) << " Error with the inverse operation ";
ASSERT_EQ( SO3<tDataType>::Identity().data_, Identity  ) << "Error with identity function ";
ASSERT_EQ( g2.Adjoint(),  g2.data_) << "Error with the Adjoint operation";

ASSERT_DOUBLE_EQ(th(0), g3.Log()(0)) << "Error with the log function";

ASSERT_EQ(g4.data_, g1.data_) << "Error with assignmet operator";

ASSERT_EQ(g5.data_, g1.data_*g2.data_) << "Error with the group operator";

}

////////////////////////////////////////////////////////////////////////

// Tests the box plus and box minus functions
TEST(SO3Test, BoxPlus) {

typedef float tDataType;
typedef Eigen::Matrix<tDataType,3,1> Vec3d;
typedef Eigen::Matrix<tDataType,3,3> Mat3d;

Vec3d th1;
Vec3d th2;
th1 << 0.1,0.2,0.3;
th2 << 0.2,0.2,0.1;
Mat3d Th1 = so3<tDataType>::Wedge(th1);
Mat3d Th2 = so3<tDataType>::Wedge(th2);



Mat3d data1 = GenElem<tDataType>(th1);
Mat3d data2 = GenElem<tDataType>(th2);
Mat3d data3 = data1*data2;

SO3<tDataType> g1(data1,true);
g1.OPlusEq(th2);
SO3<tDataType> g2(data1,true);
SO3<tDataType> g3(data1,true);
g3.BoxPlusEq(so3<tDataType>(th2));
SO3<tDataType> g4(data1,true);
SO3<tDataType> g5(data3,true);

ASSERT_EQ(g1.data_, data1*data2) << "Error with an OPlus function";

// std::cerr << SO3::BoxPlus(data1,Th2)<< std::endl;
// std::cerr << data1*data2<< std::endl;

ASSERT_EQ(SO3<tDataType>::BoxPlus(data1,Th2), data1*data2) << "Error with the BoxPlus function";
ASSERT_EQ(g2.BoxPlus(Th2), data1*data2) << "Error with the BoxPlus function";

g2.BoxPlusEq(Th2);
ASSERT_EQ(g2.data_, data1*data2) << "Error with the BoxPlus function";

ASSERT_EQ(g3.data_, data1*data2) << "Error with the BoxPlus function";


// std::cout << "th2: " << std::endl << th2 << std::endl;
// std::cout << "data1: " << std::endl << data1 << std::endl;
// std::cout << "data1 inv: " << std::endl << data1.inverse() << std::endl;
// std::cout << "data1 inv data: " << std::endl << data1.inverse()*data1 << std::endl;
// std::cout << "g5 : " << std::endl << g5.data_ << std::endl;
// std::cout << "data3 : " << std::endl << data3 << std::endl;
// std::cout << "mult: " << std::endl << data1.inverse()*g5.data_ << std::endl;

ASSERT_LE( (g5.OMinus(data1)- th2).norm(), kSO3_threshold_) << "Error with the OPlus function";



ASSERT_LE( (g5.BoxMinus(g4).data_- th2).norm(), kSO3_threshold_) << "Error with the OPlus function";

}

}