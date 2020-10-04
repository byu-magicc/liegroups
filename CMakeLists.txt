cmake_minimum_required(VERSION 3.0)

project(lie_groups)

set(CMAKE_CXX_STANDARD 11)
# set(CMAKE_BUILD_TYPE DEBUG)
set(CMAKE_BUILD_TYPE RELEASE)
# set(CMAKE_CXX_FLAGS "-march=native")

find_package(Eigen3)


include_directories(include
${EIGEN3_INCLUDE_DIR})

add_subdirectory(test)