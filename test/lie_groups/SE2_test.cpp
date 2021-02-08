#include "lie_groups/lie_groups/SE2.h"
#include "gtest/gtest.h"
#include <Eigen/Dense>


namespace lie_groups {

template<typename tDataType>
Eigen::Matrix<tDataType,3,3> GenRandElem() {
    Eigen::Matrix<tDataType,3,1> u;
    u.setRandom();
    Eigen::Matrix<tDataType,3,3> m;
    m = se2<tDataType>::Exp(u);
    return m;
}
template<typename tDataType>
Eigen::Matrix<tDataType,3,3> GenElem(Eigen::Matrix<tDataType,3,1> th) {
    Eigen::Matrix<tDataType,3,3> m;
    m = se2<tDataType>::Exp(th);
    return m;
}

////////////////////////////////////////////////////////////////////////

// Test the constructors
TEST(SE2TEST, Constructors) {

typedef float tDataType;
typedef Eigen::Matrix<tDataType,2,2> Mat2d;
typedef Eigen::Matrix<tDataType,3,3> Mat3d;

Mat3d Identity;
Identity.setIdentity();

// Valid element
Mat3d data1 =GenRandElem<tDataType>();

// Invalid element
Mat3d data2;
Eigen::Ref<Mat2d> R(data2.block(0,0,2,2));
data2.setRandom();
while  ( ((R.transpose()*R- Identity.block(0,0,2,2)).norm() <= kSE2_threshold_ ) && data2(2,0)==0 && data2(2,1)==0 && data2(2,2)==1) {
    data2.setRandom();
}



SE2<tDataType> g1;
SE2<tDataType> g2(data1,true);
SE2<tDataType> g3(data2,true);
SE2<tDataType> g4(g2);
SE2<tDataType> g5(data2);

ASSERT_EQ(g1.data_,Identity) << "Default constructor not set to identity";
ASSERT_EQ(g1.data_.block(0,2,2,1),g1.t_) << "Translational velocity map not set properly.";
ASSERT_EQ(g1.data_.block(0,0,2,2),g1.R_) << "Angular velocity map not set properly.";



ASSERT_EQ(g2.data_,data1) << "Assignment constructor error";
ASSERT_EQ(g2.data_.block(0,2,2,1),g2.t_) << "Translational velocity map not set properly.";
ASSERT_EQ(g2.data_.block(0,0,2,2),g2.R_) << "Angular velocity map not set properly.";



ASSERT_EQ(g3.data_,Identity) << "Copy constructor constructor error";
ASSERT_EQ(g3.data_.block(0,2,2,1),g3.t_) << "Translational velocity map not set properly.";
ASSERT_EQ(g3.data_.block(0,0,2,2),g3.R_) << "Angular velocity map not set properly.";

ASSERT_EQ(g4.data_,data1) << "Assignment constructor error: Invalid element not excepted. Should set element to identity.";
ASSERT_EQ(g4.data_.block(0,2,2,1),g4.t_) << "Translational velocity map not set properly.";
ASSERT_EQ(g4.data_.block(0,0,2,2),g4.R_) << "Angular velocity map not set properly.";


ASSERT_EQ(g5.data_,data2) << "Assignment constructor error";
ASSERT_EQ(g5.data_.block(0,2,2,1),g5.t_) << "Translational velocity map not set properly.";
ASSERT_EQ(g5.data_.block(0,0,2,2),g5.R_) << "Angular velocity map not set properly.";


}


////////////////////////////////////////////////////////////////////////

// Test the Inverse, Adjoint, Identity, Log, and operators
TEST(SE2TEST, IAILO) {

typedef float tDataType;
typedef Eigen::Matrix<tDataType,2,2> Mat2d;
typedef Eigen::Matrix<tDataType,3,3> Mat3d;

Mat3d Identity;
Identity.setIdentity();
Eigen::Matrix<tDataType,3,1> th;
th << 0.1,0.2,0.3;
Mat3d g = GenElem<tDataType>(th);

se2<tDataType> u1(Eigen::Matrix<tDataType,3,1>::Random());

SE2<tDataType> g1(GenRandElem<tDataType>());
SE2<tDataType> g2(GenRandElem<tDataType>());
SE2<tDataType> g3(g);
SE2<tDataType> g4 = g1;
SE2<tDataType> g5 = g1*g2;

ASSERT_LE( ((g1.Inverse()).data_ - (g1.data_).inverse() ).norm(), kSE2_threshold_) << " Error with the inverse operation ";

ASSERT_EQ( SE2<tDataType>::Identity().data_, Identity  ) << "Error with identity function ";

se2<tDataType> u2(g2.Adjoint()*u1.data_);
se2<tDataType> u3(g2.data_*u1.Wedge()*(g2.Inverse()).data_,true);

ASSERT_LE( (u2.data_-u3.data_).norm(), kSE2_threshold_) << "Error with the Adjoint operation";

ASSERT_LE( (se2<tDataType>::Exp(g3.Log())-g3.data_).norm(), kSE2_threshold_ ) << "Error with the log function";

ASSERT_EQ(g4.data_, g1.data_) << "Error with assignmet operator";

ASSERT_EQ(g5.data_, g1.data_*g2.data_) << "Error with the group operator";

}

////////////////////////////////////////////////////////////////////////

// Tests the box plus and box minus functions
TEST(SE2Test, BoxPlus) {

typedef float tDataType;
typedef Eigen::Matrix<tDataType,2,2> Mat2d;
typedef Eigen::Matrix<tDataType,3,3> Mat3d;

Eigen::Matrix<tDataType,3,1> th1;
Eigen::Matrix<tDataType,3,1> th2;
th1 << 0.1,0.2,0.3;
th2 << 0.2,0.2,0.1;
Mat3d Th1 = se2<tDataType>::Wedge(th1);
Mat3d Th2 = se2<tDataType>::Wedge(th2);



Mat3d data1 = GenElem<tDataType>(th1);
Mat3d data2 = GenElem<tDataType>(th2);
Mat3d data3 = data1*data2;

SE2<tDataType> g1(data1,true);



g1.OPlusEq(th2);
SE2<tDataType> g2(data1,true);
SE2<tDataType> g3(data1,true);
g3.BoxPlusEq(se2<tDataType>(th2));
SE2<tDataType> g4(data1,true);
SE2<tDataType> g5(data3,true);

ASSERT_EQ(g1.data_, data1*data2) << "Error with an OPlus function";

// std::cerr << SE2::BoxPlus(data1,Th2)<< std::endl;
// std::cerr << data1*data2<< std::endl;

ASSERT_EQ(SE2<tDataType>::BoxPlus(data1,Th2), data1*data2) << "Error with the BoxPlus function";
ASSERT_EQ(g2.BoxPlus(Th2), data1*data2) << "Error with the BoxPlus function";

g2.BoxPlusEq(Th2);
ASSERT_EQ(g2.data_, data1*data2) << "Error with the BoxPlus function";

ASSERT_EQ(g3.data_, data1*data2) << "Error with the BoxPlus function";

ASSERT_LE( (g5.OMinus(data1)- th2).norm(), kSE2_threshold_) << "Error with the OPlus function";

ASSERT_LE( (g5.BoxMinus(g4).data_- th2).norm(), kSE2_threshold_) << "Error with the OPlus function";

}

} // namespace lie_groups