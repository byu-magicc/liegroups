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
  
########################################
# Define Tests
########################################



# Lie Algebra Tests
add_executable(so2_test
lie_algebras/so2_test
../src/lie_algebras/so2.cpp)
target_link_libraries(so2_test gtest_main)

add_executable(se2_test
lie_algebras/se2_test
../src/lie_algebras/se2.cpp)

target_link_libraries(se2_test gtest_main)

add_executable(so3_test
lie_algebras/so3_test
../src/lie_algebras/so3.cpp)

target_link_libraries(so3_test gtest_main)


# add_executable(SO2_test
#  lie_groups/SO2_test.cpp 
#  ../src/lie_groups/SO2.cpp)
# target_link_libraries(SO2_test gtest_main)

add_custom_target(gtests_SO2 SO2_test)