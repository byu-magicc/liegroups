#include "lie_groups/SO2.h"
#include "gtest/gtest.h"
#include <Eigen/Dense>


namespace lie_groups {


// Test the constructors
TEST(SO2TEST, Constructors) {

Eigen::Matrix2d g;
double th = 0.3;
// Valid element
g <<  cos(th), -sin(th), sin(th), cos(th);

SO2 g1;
SO2 g2(g);
SO2 g3(g2);

ASSERT_EQ(g1.g_,Eigen::Matrix2d::Identity()) << "Default constructor not set to identity";
ASSERT_EQ(g2.g_,g) << "Assignment constructor error";
ASSERT_EQ(g3.g_,g) << "Copy constructor constructor error";

// Not valid element. Determinant is 1 but is not unitary
g <<  2, 0, 0, 0.5;
SO2 g4(g);


ASSERT_NE(g4.g_,g) << "Assignment constructor error: Excepted invalid element";
ASSERT_EQ(g4.g_,Eigen::Matrix2d::Identity()) << "Assignment constructor error: Invalid element not excepted. Should set element to identity.";


}


}