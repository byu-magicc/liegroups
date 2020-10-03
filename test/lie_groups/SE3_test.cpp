#include "lie_groups/SE3.h"
#include "gtest/gtest.h"
#include <Eigen/Dense>


namespace lie_groups {


Eigen::Matrix4d GenRandElem() {
    Eigen::Matrix<double,6,1> u;
    u.setRandom();
    Eigen::Matrix4d m;
    m = se3::Exp(u);
    return m;
}

Eigen::Matrix4d GenElem(Eigen::Matrix<double,6,1> th) {
    Eigen::Matrix4d m;
    m = se3::Exp(th);
    return m;
}

////////////////////////////////////////////////////////////////////////

// Test the constructors
TEST(SE3TEST, Constructors) {

Eigen::Matrix4d Identity;
Identity.setIdentity();

// Valid element
Eigen::Matrix4d data1 =GenRandElem();

// // Invalid element
// Eigen::Matrix4d data2;
// data2.setRandom();
// while  ((data2.transpose()*data2 - Identity).norm() <= kSE3_threshold_ ) {
//     data2.setRandom();
// }


SE3 g1;
SE3 g2;
SE3 g3;
for (unsigned long int i=0; i<10000;i++)
{
 g3 = g1.Inverse();   
}


// SE3 g1;
// SE3 g2(data1,true);
// SE3 g3(data2,true);
// SE3 g4(g2);
// SE3 g5(data2);

// ASSERT_EQ(g1.data_,Identity) << "Default constructor not set to identity";
// ASSERT_EQ(g2.data_,data1) << "Assignment constructor error";
// ASSERT_EQ(g3.data_,Identity) << "Copy constructor constructor error";
// ASSERT_EQ(g4.data_,data1) << "Assignment constructor error: Invalid element not excepted. Should set element to identity.";
// ASSERT_EQ(g5.data_,data2) << "Assignment constructor error";


}

/*
////////////////////////////////////////////////////////////////////////

// Test the Inverse, Adjoint, Identity, Log, and operators
TEST(SE3TEST, IAILO) {

Eigen::Matrix4d Identity;
Identity.setIdentity();
Eigen::Matrix<double,6,1> th;
th << 0.1,0.2,0.3;
Eigen::Matrix4d g = GenElem(th);

se3 u1(Eigen::Matrix<double,6,1>::Random());

SE3 g1(GenRandElem());
SE3 g2(GenRandElem());
SE3 g3(g);
SE3 g4 = g1;
SE3 g5 = g1*g2;

ASSERT_LE( ((g1.Inverse()).data_ - (g1.data_).inverse() ).norm(), kSE3_threshold_) << " Error with the inverse operation ";

ASSERT_EQ( SE3::Identity().data_, Identity  ) << "Error with identity function ";

se3 u2(g2.Adjoint()*u1.data_);
se3 u3(g2.data_*u1.Wedge()*(g2.Inverse()).data_,true);

ASSERT_LE( (u2.data_-u3.data_).norm(), kSE3_threshold_) << "Error with the Adjoint operation";

ASSERT_LE( (se3::Exp(g3.Log())-g3.data_).norm(), kSE3_threshold_ ) << "Error with the log function";

ASSERT_EQ(g4.data_, g1.data_) << "Error with assignmet operator";

ASSERT_EQ(g5.data_, g1.data_*g2.data_) << "Error with the group operator";

}

////////////////////////////////////////////////////////////////////////

// Tests the box plus and box minus functions
TEST(SE3Test, BoxPlus) {

Eigen::Matrix<double,6,1> th1;
Eigen::Matrix<double,6,1> th2;
th1 << 0.1,0.2,0.3;
th2 << 0.2,0.2,0.1;
Eigen::Matrix4d Th1 = se3::Wedge(th1);
Eigen::Matrix4d Th2 = se3::Wedge(th2);



Eigen::Matrix4d data1 = GenElem(th1);
Eigen::Matrix4d data2 = GenElem(th2);
Eigen::Matrix4d data3 = data1*data2;

SE3 g1(data1,true);
g1.OPlusEq(th2);
SE3 g2(data1,true);
SE3 g3(data1,true);
g3.BoxPlusEq(se3(th2));
SE3 g4(data1,true);
SE3 g5(data3,true);

ASSERT_EQ(g1.data_, data1*data2) << "Error with an OPlus function";

// std::cerr << SE3::BoxPlus(data1,Th2)<< std::endl;
// std::cerr << data1*data2<< std::endl;

ASSERT_EQ(SE3::BoxPlus(data1,Th2), data1*data2) << "Error with the BoxPlus function";
ASSERT_EQ(g2.BoxPlus(Th2), data1*data2) << "Error with the BoxPlus function";

g2.BoxPlusEq(Th2);
ASSERT_EQ(g2.data_, data1*data2) << "Error with the BoxPlus function";

ASSERT_EQ(g3.data_, data1*data2) << "Error with the BoxPlus function";

ASSERT_LE( (g4.OMinus(data3)- th2).norm(), kSE3_threshold_) << "Error with the OPlus function";

ASSERT_LE( (g4.BoxMinus(g5).data_- th2).norm(), kSE3_threshold_) << "Error with the OPlus function";

}
*/
}