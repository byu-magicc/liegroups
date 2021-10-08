#include "gtest/gtest.h"

#include <Eigen/Dense>
#include <ctime>
#include <chrono>

#include "lie_groups/state.h"

namespace lie_groups {


using MyTypes = ::testing::Types<State<Rn,double,1,3>,State<Rn,double,2,2>,State<Rn,double,2,3>,State<Rn,double,3,2>,State<Rn,double,3,3>,R2_r2,R3_r3,SO2_so2,SO3_so3,SE2_se2,SE3_se3>;
// using MyTypes = ::testing::Types<State<Rn,double,3,3>>;

template <typename T>
class ContructorTest : public testing::Test {
public:
typedef T type;
};

TYPED_TEST_SUITE(ContructorTest, MyTypes);

// Used to test the functions inverse, identity and adjoint
template <typename T>
class OtherFunctionTest : public testing::Test {
    public:
    typedef T type;
};

TYPED_TEST_SUITE(OtherFunctionTest, MyTypes);

// Used to test the boxplus boxminus oplus and ominus operators
template <typename T>
class BoxOTest : public testing::Test {
    public:
    typedef T type;
};

TYPED_TEST_SUITE(BoxOTest, MyTypes);

////////////////////////////////////////////////////////////
//                        Constructor test
////////////////////////////////////////////////////////////

TYPED_TEST(ContructorTest, Constructors) {


// Test default constructor
TypeParam state;
if(TypeParam::g_type_::size2_ == 1) {
    ASSERT_EQ(state.g_.data_, TypeParam::Mat_G::Zero());
    ASSERT_EQ(state.u_.data_, TypeParam::Mat_C::Zero());
} else {
    ASSERT_EQ(state.g_.data_, TypeParam::Mat_G::Identity());
    ASSERT_EQ(state.u_.data_, TypeParam::Mat_C::Zero());

}

// Copy constructor and assignment
state.g_.data_.setRandom();
state.u_.data_.setRandom();
TypeParam state_copy(state);
ASSERT_EQ(state_copy.g_.data_,state.g_.data_);
ASSERT_EQ(state_copy.u_.data_,state.u_.data_);

TypeParam state_copy2 = state;
ASSERT_EQ(state_copy2.g_.data_,state.g_.data_);
ASSERT_EQ(state_copy2.u_.data_,state.u_.data_);

//Move constructor and assignment
TypeParam&& state_tmp = TypeParam();
state_tmp.g_.data_.setRandom();
state_tmp.u_.data_.setRandom();
TypeParam state_move(state_tmp);
TypeParam state_move2 = state_tmp;
ASSERT_EQ(state_move.g_.data_, state_tmp.g_.data_);
ASSERT_EQ(state_move.u_.data_, state_tmp.u_.data_);
ASSERT_EQ(state_move2.g_.data_,state_tmp.g_.data_);
ASSERT_EQ(state_move2.u_.data_,state_tmp.u_.data_);

// Copy and move constructor using group and algebra elements
TypeParam state_copy3(state.g_,state.u_);
TypeParam state_move3(state_tmp.g_, state_tmp.u_);
ASSERT_EQ(state_copy3.g_.data_, state.g_.data_);
ASSERT_EQ(state_copy3.u_.data_, state.u_.data_);
ASSERT_EQ(state_move3.g_.data_,state_tmp.g_.data_);
ASSERT_EQ(state_move3.u_.data_,state_tmp.u_.data_);

// Initialize state using data of elements
typename TypeParam::Mat_G g;
typename TypeParam::Mat_C c;
typename TypeParam::Mat_A a;
bool set = true;
if(TypeParam::g_type_::size2_ != 1) { // Test invalid entries

    g.setOnes();
    a.setOnes();
    c.setRandom();
    TypeParam state_invalid1(set,g,c);
    ASSERT_EQ(state_invalid1.g_.data_,TypeParam::Mat_G::Identity());
    ASSERT_EQ(state_invalid1.u_.data_,c);
    TypeParam state_invalid2(g,a,set);
    ASSERT_EQ(state_invalid2.g_.data_,TypeParam::Mat_G::Identity());
    ASSERT_EQ(state_invalid2.u_.data_,TypeParam::Mat_C::Zero());
    // ASSERT_FALSE(set);
    
} 


g = TypeParam::g_type_::Random();
c.setRandom();
a = TypeParam::u_type_::Wedge(c);

TypeParam state_valid1(set,g,c);
ASSERT_EQ(state_valid1.g_.data_,g);
ASSERT_EQ(state_valid1.u_.data_,c);

TypeParam state_valid2(g,a,set);
ASSERT_EQ(state_valid2.g_.data_,g);
ASSERT_EQ(state_valid2.u_.data_,c);



}


////////////////////////////////////////////////////////////
//               Other Function test
////////////////////////////////////////////////////////////
TYPED_TEST(OtherFunctionTest, OtherFunctions) {


TypeParam state;

state.g_.data_ = TypeParam::g_type_::Random();
state.u_.data_.setRandom();

// Test inverse
TypeParam state_inverse = state.Inverse();
ASSERT_EQ(state_inverse.g_.data_, (state.g_.Inverse()).data_);
ASSERT_EQ(state_inverse.u_.data_, -state.u_.data_);

// // Test state Adjoint
// typename TypeParam::Mat_Adj Adjoint = state.Adjoint();
// typename TypeParam::Mat_Adj Adjoint_verify;
// Adjoint_verify.block(0,0,TypeParam::g_type_::dim_,TypeParam::g_type_::dim_) = state.g_.Adjoint();
// Adjoint_verify.block(TypeParam::g_type_::dim_,0, TypeParam::g_type_::dim_,TypeParam::g_type_::dim_) = state.u_.Adjoint();
// ASSERT_EQ(Adjoint,Adjoint_verify);

// Test identity
state = TypeParam::Identity();
if(TypeParam::g_type_::size2_ == 1) {
    ASSERT_EQ(state.g_.data_, TypeParam::Mat_G::Zero());
    ASSERT_EQ(state.u_.data_, TypeParam::Mat_C::Zero());
} else {
    ASSERT_EQ(state.g_.data_, TypeParam::Mat_G::Identity());
    ASSERT_EQ(state.u_.data_, TypeParam::Mat_C::Zero());

}

// Random state test
TypeParam state_random = TypeParam::Random();
ASSERT_TRUE(TypeParam::g_type_::isElement(state_random.g_.data_));

// Group multiplication
TypeParam state1 = TypeParam::Random();
TypeParam state2 = TypeParam::Random();
TypeParam state3 = state1*state2;
typename TypeParam::g_type_ g = state1.g_*state2.g_;
typename TypeParam::u_type_ u = state1.u_+state2.u_;

ASSERT_EQ(state3.g_.data_, g.data_);
ASSERT_EQ(state3.u_.data_, u.data_);


}

////////////////////////////////////////////////////////////
//               BoxPlus/Minus OPlus/Minus test
////////////////////////////////////////////////////////////
TYPED_TEST(BoxOTest, OtherFunctions) {


TypeParam state1 = TypeParam::Random();
TypeParam state2 = TypeParam::Random();
TypeParam state3 = state1*state2;

typename TypeParam::Vec_SC data2 = state3.OMinus(state1);
typename TypeParam::Vec_SC data2_correct;
data2_correct.block(0,0,TypeParam::g_type_::dim_,1) = state2.g_.Log().block(0,0,TypeParam::g_type_::dim_,1);
data2_correct.block(TypeParam::g_type_::dim_,0,TypeParam::u_type_::total_num_dim_,1) = state2.u_.data_;

ASSERT_LE( (data2-data2_correct).norm(), 1e-10 );

TypeParam state4 = state1;
state4.OPlusEQ(data2);

ASSERT_LE( (state4.g_.data_ - state3.g_.data_).norm(), 1e-10);
ASSERT_LE( (state4.u_.data_ - state3.u_.data_).norm(), 1e-10);

state4 = state1.OPlus(data2);
ASSERT_LE( (state4.g_.data_ - state3.g_.data_).norm(), 1e-10);
ASSERT_LE( (state4.u_.data_ - state3.u_.data_).norm(), 1e-10);

}


//////////////////////////////////////////////////////////////////////////////////////////////////
//                               Jacobian Tests
/////////////////////////////////////////////////////////////////////////////////////////////////


TYPED_TEST(BoxOTest, Jacobians) {


TypeParam state = TypeParam::Random(10);
typename TypeParam::Vec_SC tau = TypeParam::State::Log(state);
TypeParam state_test = TypeParam::State::Exp(tau);

ASSERT_LT( TypeParam::OMinus(state,state_test).norm(), 1e-9);

typename TypeParam::Mat_SC jacobian, Jr, Jr_inv, Jl, Jl_inv;
typename TypeParam::Vec_SC perturbation;
jacobian.setZero();
perturbation.setZero();
double dt = 1e-7;

Jr = TypeParam::Jr( tau );
Jr_inv = TypeParam::JrInv( tau);
Jl = TypeParam::Jl( -tau );             // Negate them so that Jl(- tau) = Jr(tau)
Jl_inv = TypeParam::JlInv( -tau );      // Negate them so that Jl(- tau) = Jr(tau)

// Right Jacobian
for (size_t ii = 0; ii < perturbation.rows(); ++ii) {
    perturbation.setZero();
    perturbation(ii,0) = dt;

    typename TypeParam::State state_perturbed = TypeParam::Exp(tau+perturbation);
    // std::cout << "sp: " << std::endl << state_perturbed.g_.data_ << std::endl;
    // std::cout << "s: " << std::endl << state.g_.data_ << std::endl;

    jacobian.block(0,ii,TypeParam::dim_,1) = TypeParam::State::OMinus( state_perturbed ,state) / dt;

}

// std::cout << "Jr: " << std::endl << Jr << std::endl;
// std::cout << "Jr est: " << std::endl << jacobian << std::endl;

ASSERT_LT( (Jr - jacobian).norm(), 1e-6);
ASSERT_LT( (Jl - jacobian).norm(), 1e-6);
ASSERT_LT( (Jr_inv - jacobian.inverse()).norm(), 1e-6);
ASSERT_LT( (Jl_inv - jacobian.inverse()).norm(), 1e-6);


}







} // namespace lie_groups