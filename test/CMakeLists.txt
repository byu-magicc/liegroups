cmake_minimum_required(VERSION 2.8.12)

########################################
# Install googletest in local directory
########################################

set(GTEST_DIR ${CMAKE_CURRENT_BINARY_DIR}/googletest-download)

# Download and install GoogleTest
configure_file(
	"${PROJECT_SOURCE_DIR}/cmake/googletest.cmake.in"
	"${GTEST_DIR}/CMakeLists.txt"
	IMMEDIATE @ONLY)

execute_process(COMMAND ${CMAKE_COMMAND} -G "${CMAKE_GENERATOR}" . WORKING_DIRECTORY ${GTEST_DIR} )
execute_process(COMMAND ${CMAKE_COMMAND} --build . WORKING_DIRECTORY ${GTEST_DIR})

# Prevent overriding the parent project's compiler/linker settings on Windows
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)

# Add googletest directly to the build. This defines the `gtest` and `gtest_main` targets. Don't add these
# targets to the ALL target so that `make install` will ignore googletest installation targets
add_subdirectory(${GTEST_DIR}/googletest-src ${GTEST_DIR}/googletest-build EXCLUDE_FROM_ALL)
 

include(CTest)
########################################
# Define Tests
########################################



# Lie Algebra Tests

add_executable(rn_test
lie_algebras/rn_test)
target_link_libraries(rn_test gtest_main)
add_test(NAME AllTestsInsrn_test COMMAND rn_test)


add_executable(so2_test
lie_algebras/so2_test
../src/lie_algebras/so2.cpp)
target_link_libraries(so2_test gtest_main)
add_test(NAME AllTestsInso2_test COMMAND so2_test)

add_executable(se2_test
lie_algebras/se2_test
../src/lie_algebras/se2.cpp)
target_link_libraries(se2_test gtest_main)
add_test(NAME AllTestsInse2_test COMMAND se2_test)

add_executable(so3_test
lie_algebras/so3_test
../src/lie_algebras/so3.cpp)
target_link_libraries(so3_test gtest_main)
add_test(NAME AllTestsInso3_test COMMAND so3_test)

add_executable(se3_test
lie_algebras/se3_test
../src/lie_algebras/se3.cpp
../src/lie_algebras/so3.cpp)
target_link_libraries(se3_test gtest_main)
add_test(NAME AllTestsInse3_test COMMAND se3_test)



# Lie group tests

add_executable(Rn_test
lie_groups/Rn_test)
target_link_libraries(Rn_test gtest_main)
add_test(NAME AllTestsInsRn_test COMMAND Rn_test)

add_executable(SO2_test
 lie_groups/SO2_test.cpp 
 ../src/lie_groups/SO2.cpp
 ../src/lie_algebras/so2.cpp)
target_link_libraries(SO2_test gtest_main)
add_test(NAME AllTestsInSO2_test COMMAND SO2_test)


add_executable(SO3_test
 lie_groups/SO3_test.cpp 
 ../src/lie_groups/SO3.cpp
 ../src/lie_algebras/so3.cpp)
target_link_libraries(SO3_test gtest_main)
add_test(NAME AllTestsInSO3_test COMMAND SO3_test)


add_executable(SE2_test
 lie_groups/SE2_test.cpp 
 ../src/lie_groups/SE2.cpp
 ../src/lie_algebras/se2.cpp)
target_link_libraries(SE2_test gtest_main)
add_test(NAME AllTestsInSE2_test COMMAND SE2_test)


add_executable(SE3_test
 lie_groups/SE3_test.cpp 
 ../src/lie_groups/SE3.cpp
 ../src/lie_algebras/se3.cpp
 ../src/lie_algebras/so3.cpp)
target_link_libraries(SE3_test gtest_main)
add_test(NAME AllTestsInSE3_test COMMAND SE3_test)

# State test

add_executable(STATE_test
state_test.cpp
../src/lie_algebras/so2.cpp
../src/lie_algebras/so3.cpp
../src/lie_algebras/se2.cpp
../src/lie_algebras/se3.cpp
../src/lie_groups/SO2.cpp
../src/lie_groups/SO3.cpp
../src/lie_groups/SE2.cpp
../src/lie_groups/SE3.cpp
)
target_link_libraries(STATE_test gtest_main)
add_test(NAME AllTestsInState_test COMMAND State_test)

