#include "lie_groups/SE2.h"
#include "gtest/gtest.h"
#include <Eigen/Dense>


namespace lie_groups {


Eigen::Matrix3d GenRandElem() {
    Eigen::Matrix<double,3,1> u;
    u.setRandom();
    Eigen::Matrix3d m;
    m = se2::Exp(u);
    return m;
}

Eigen::Matrix3d GenElem(Eigen::Matrix<double,3,1> th) {
    Eigen::Matrix3d m;
    m = se2::Exp(th);
    return m;
}

////////////////////////////////////////////////////////////////////////

// Test the constructors
TEST(SE2TEST, Constructors) {

Eigen::Matrix3d Identity;
Identity.setIdentity();

// Valid element
Eigen::Matrix3d data1 =GenRandElem();

// Invalid element
Eigen::Matrix3d data2;
data2.setRandom();
while  ((data2.transpose()*data2 - Identity).norm() <= kSE2_threshold_ ) {
    data2.setRandom();
}



SE2 g1;
g1.data_ << 1,2,3,4,5,6,7,8,9;

std::cerr << g1.data_ << std::endl;
std::cerr << g1.R_ << std::endl;
std::cerr << g1.t_ << std::endl;

// SE2 g2(data1,true);
// SE2 g3(data2,true);
// SE2 g4(g2);
// SE2 g5(data2);

// ASSERT_EQ(g1.data_,Identity) << "Default constructor not set to identity";
// ASSERT_EQ(g2.data_,data1) << "Assignment constructor error";
// ASSERT_EQ(g3.data_,Identity) << "Copy constructor constructor error";
// ASSERT_EQ(g4.data_,data1) << "Assignment constructor error: Invalid element not excepted. Should set element to identity.";
// ASSERT_EQ(g5.data_,data2) << "Assignment constructor error";


}

/*
////////////////////////////////////////////////////////////////////////

// Test the Inverse, Adjoint, Identity, Log, and operators
TEST(SE2TEST, IAILO) {

Eigen::Matrix3d Identity;
Identity.setIdentity();
Eigen::Matrix<double,3,1> th;
th << 0.1,0.2,0.3;
Eigen::Matrix3d g = GenElem(th);

SE2 g1(GenRandElem());
SE2 g2(GenRandElem());
SE2 g3(g);
SE2 g4 = g1;
// SE2 g5 = g1*g2;

ASSERT_LE( ((g1.Inverse()).data_ - g1.data_.inverse()).norm(), kSE2_threshold_) << " Error with the inverse operation ";
ASSERT_EQ( SE2::Identity().data_, Identity  ) << "Error with identity function ";
ASSERT_EQ( g2.Adjoint(),  g2.data_) << "Error with the Adjoint operation";

ASSERT_DOUBLE_EQ(th(0), g3.Log()(0)) << "Error with the log function";

ASSERT_EQ(g4.data_, g1.data_) << "Error with assignmet operator";

// ASSERT_EQ(g5.data_, g1.data_*g2.data_) << "Error with the group operator";

}

////////////////////////////////////////////////////////////////////////

// Tests the box plus and box minus functions
TEST(SE2Test, BoxPlus) {

Eigen::Matrix<double,3,1> th1;
Eigen::Matrix<double,3,1> th2;
th1 << 0.1,0.2,0.3;
th2 << 0.2,0.2,0.1;
Eigen::Matrix3d Th1 = se2::Wedge(th1);
Eigen::Matrix3d Th2 = se2::Wedge(th2);



Eigen::Matrix3d data1 = GenElem(th1);
Eigen::Matrix3d data2 = GenElem(th2);
Eigen::Matrix3d data3 = data1*data2;

SE2 g1(data1,true);
g1.OPlusEq(th2);
SE2 g2(data1,true);
SE2 g3(data1,true);
g3.BoxPlusEq(se2(th2));
SE2 g4(data1,true);
SE2 g5(data3,true);

ASSERT_EQ(g1.data_, data1*data2) << "Error with an OPlus function";

// std::cerr << SE2::BoxPlus(data1,Th2)<< std::endl;
// std::cerr << data1*data2<< std::endl;

ASSERT_EQ(SE2::BoxPlus(data1,Th2), data1*data2) << "Error with the BoxPlus function";
ASSERT_EQ(g2.BoxPlus(Th2), data1*data2) << "Error with the BoxPlus function";

g2.BoxPlusEq(Th2);
ASSERT_EQ(g2.data_, data1*data2) << "Error with the BoxPlus function";

ASSERT_EQ(g3.data_, data1*data2) << "Error with the BoxPlus function";

ASSERT_LE( (g4.OMinus(data3)- th2).norm(), kSE2_threshold_) << "Error with the OPlus function";

ASSERT_LE( (g4.BoxMinus(g5).data_- th2).norm(), kSE2_threshold_) << "Error with the OPlus function";

}
*/
}