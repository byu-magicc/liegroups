#include "lie_groups/SE3.h"
#include "gtest/gtest.h"
#include <Eigen/Dense>


namespace lie_groups {

template<typename tDataType>
Eigen::Matrix<tDataType,4,4> GenRandElem() {
    Eigen::Matrix<tDataType,6,1> u;
    u.setRandom();
    Eigen::Matrix<tDataType,4,4> m;
    m = se3<tDataType>::Exp(u);
    return m;
}

template<typename tDataType>
Eigen::Matrix<tDataType,4,4> GenElem(Eigen::Matrix<tDataType,6,1> th) {
    Eigen::Matrix<tDataType,4,4> m;
    m = se3<tDataType>::Exp(th);
    return m;
}

////////////////////////////////////////////////////////////////////////

// Test the constructors
TEST(SE3TEST, Constructors) {

typedef float tDataType;
typedef Eigen::Matrix<tDataType,3,3> Mat3d;
typedef Eigen::Matrix<tDataType,4,4> Mat4d;

Mat4d Identity;
Identity.setIdentity();

// Valid element
Mat4d data1 =GenRandElem<tDataType>();

// Invalid element
Mat4d data2;
Eigen::Ref<Mat3d> R(data2.block(0,0,3,3));
data2.setRandom();
while  ( ((R.transpose()*R- Identity.block(0,0,3,3)).norm() <= kSE3_threshold_ ) && data2(3,0)==0 && data2(3,1)==0 && data2(3,2)==0 && data2(3,3)==1) {
    data2.setRandom();
}



SE3<tDataType>g1;
SE3<tDataType>g2(data1,true);
SE3<tDataType>g3(data2,true);
SE3<tDataType>g4(g2);
SE3<tDataType>g5(data2);

ASSERT_EQ(g1.data_,Identity) << "Default constructor not set to identity";
ASSERT_EQ(g1.data_.block(0,3,3,1),g1.t_) << "Translational velocity map not set properly.";
ASSERT_EQ(g1.data_.block(0,0,3,3),g1.R_) << "Angular velocity map not set properly.";



ASSERT_EQ(g2.data_,data1) << "Assignment constructor error";
ASSERT_EQ(g2.data_.block(0,3,3,1),g2.t_) << "Translational velocity map not set properly.";
ASSERT_EQ(g2.data_.block(0,0,3,3),g2.R_) << "Angular velocity map not set properly.";



ASSERT_EQ(g3.data_,Identity) << "Copy constructor constructor error";
ASSERT_EQ(g3.data_.block(0,3,3,1),g3.t_) << "Translational velocity map not set properly.";
ASSERT_EQ(g3.data_.block(0,0,3,3),g3.R_) << "Angular velocity map not set properly.";

ASSERT_EQ(g4.data_,data1) << "Assignment constructor error: Invalid element not excepted. Should set element to identity.";
ASSERT_EQ(g4.data_.block(0,3,3,1),g4.t_) << "Translational velocity map not set properly.";
ASSERT_EQ(g4.data_.block(0,0,3,3),g4.R_) << "Angular velocity map not set properly.";


ASSERT_EQ(g5.data_,data2) << "Assignment constructor error";
ASSERT_EQ(g5.data_.block(0,3,3,1),g5.t_) << "Translational velocity map not set properly.";
ASSERT_EQ(g5.data_.block(0,0,3,3),g5.R_) << "Angular velocity map not set properly.";



}


////////////////////////////////////////////////////////////////////////

// Test the Inverse, Adjoint, Identity, Log, and operators
TEST(SE3TEST, IAILO) {

typedef float tDataType;
typedef Eigen::Matrix<tDataType,6,1> Vec6d;
typedef Eigen::Matrix<tDataType,3,3> Mat3d;
typedef Eigen::Matrix<tDataType,4,4> Mat4d;

Mat4d Identity;
Identity.setIdentity();
Vec6d th;
th << -1, 2, 3, 0.1, 0.2, 0.3;
Mat4d g = GenElem(th);

se3<tDataType>u1(Vec6d::Random());

SE3<tDataType>g1(GenRandElem<tDataType>());
SE3<tDataType>g2(GenRandElem<tDataType>());
SE3<tDataType>g3(g);
SE3<tDataType>g4 = g1;
SE3<tDataType>g5 = g1*g2;

ASSERT_LE( ((g1.Inverse()).data_ - (g1.data_).inverse() ).norm(), kSE3_threshold_) << " Error with the inverse operation ";

ASSERT_EQ( SE3<tDataType>::Identity().data_, Identity  ) << "Error with identity function ";

se3<tDataType>u2(g2.Adjoint()*u1.data_);
se3<tDataType>u3(g2.data_*u1.Wedge()*(g2.Inverse()).data_,true);

ASSERT_LE( (u2.data_-u3.data_).norm(), kSE3_threshold_) << "Error with the Adjoint operation";

ASSERT_LE( (se3<tDataType>::Exp(g3.Log())-g3.data_).norm(), kSE3_threshold_ ) << "Error with the log function";

ASSERT_EQ(g4.data_, g1.data_) << "Error with assignmet operator";

ASSERT_EQ(g5.data_, g1.data_*g2.data_) << "Error with the group operator";

}


////////////////////////////////////////////////////////////////////////

// Tests the box plus and box minus functions
TEST(SE3Test, BoxPlus) {

typedef double tDataType;
typedef Eigen::Matrix<tDataType,6,1> Vec6d;
typedef Eigen::Matrix<tDataType,3,3> Mat3d;
typedef Eigen::Matrix<tDataType,4,4> Mat4d;

Vec6d th1;
Vec6d th2;
th1 << 1,2,3, 0.1,0.2,0.3;
th2 << 4,5,6, 0.2,0.2,0.1;
Mat4d Th1 = se3<tDataType>::Wedge(th1);
Mat4d Th2 = se3<tDataType>::Wedge(th2);



Mat4d data1 = GenElem(th1);
Mat4d data2 = GenElem(th2);
Mat4d data3 = data1*data2;

SE3<tDataType>g1(data1,true);
g1.OPlusEq(th2);
SE3<tDataType>g2(data1,true);
SE3<tDataType>g3(data1,true);
g3.BoxPlusEq(se3<tDataType>(th2));
SE3<tDataType>g4(data1,true);
SE3<tDataType>g5(data3,true);

ASSERT_EQ(g1.data_, data1*data2) << "Error with an OPlus function";

// std::cerr << se3<tDataType>::BoxPlus(data1,Th2)<< std::endl;
// std::cerr << data1*data2<< std::endl;

ASSERT_EQ(SE3<tDataType>::BoxPlus(data1,Th2), data1*data2) << "Error with the BoxPlus function";
ASSERT_EQ(g2.BoxPlus(Th2), data1*data2) << "Error with the BoxPlus function";

g2.BoxPlusEq(Th2);
ASSERT_EQ(g2.data_, data1*data2) << "Error with the BoxPlus function";

ASSERT_EQ(g3.data_, data1*data2) << "Error with the BoxPlus function";

ASSERT_LE( (g5.OMinus(data1)- th2).norm(), kSE3_threshold_) << "Error with the OPlus function";

ASSERT_LE( (g5.BoxMinus(g4).data_- th2).norm(), kSE3_threshold_) << "Error with the OPlus function";

}

}