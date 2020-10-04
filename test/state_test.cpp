#include "state.h"
#include "gtest/gtest.h"

#include <Eigen/Dense>
#include <ctime>
#include <chrono>


namespace lie_groups {

using namespace std::chrono;

TEST(STATETEST, Constructors) {

// State<SO2,so2> x;

high_resolution_clock::time_point t1 = high_resolution_clock::now();

for (unsigned long int i(100000); i!=0; --i)
    SO2 x;

high_resolution_clock::time_point t2 = high_resolution_clock::now();

duration<double> time_span = duration_cast<duration<double>>(t2 - t1);

std::cout << "It took me " << time_span.count() << " seconds.";
std::cout << std::endl;

// std::cerr << x.Adjoint() << std::endl;

// State<SE2,se2> y;



}


}