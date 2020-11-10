#include "state.h"
#include "gtest/gtest.h"

#include <Eigen/Dense>
#include <ctime>
#include <chrono>


namespace lie_groups {


using MyTypes = ::testing::Types<R2_r2,R3_r3,SO2_so2,SO3_so3,SE2_se2,SE3_se3>;

template <typename T>
class ContructorTest : public testing::Test {
public:
typedef T type;
};

TYPED_TEST_SUITE(ContructorTest, MyTypes);
TYPED_TEST(ContructorTest, Constructors) {


// Test default constructor
TypeParam state;
std::cout << state.g_.data_ << std::endl << std::endl;
// ASSERT_EQ(state.g_.data_, TypeParam::g_type_::Identity());



}


}