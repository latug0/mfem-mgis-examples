cmake_minimum_required(VERSION 3.13)
project(BimaterialExample)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

include(CTest)
# testing
set(CTEST_CONFIGURATION_TYPE "${JOB_BUILD_CONFIGURATION}")
# (must be placed *before* any add_subdirectory, cmake bug ?)
enable_testing()
if(CMAKE_CONFIGURATION_TYPES)
   if (NOT TARGET check)
      add_custom_target(check COMMAND 
         ${CMAKE_CTEST_COMMAND} -T test -C $<CONFIGURATION>)
   endif(NOT TARGET check)
else(CMAKE_CONFIGURATION_TYPES)
   if (NOT TARGET check)
      add_custom_target(check COMMAND 
      ${CMAKE_CTEST_COMMAND} -T test )
   endif(NOT TARGET check)
endif(CMAKE_CONFIGURATION_TYPES)

# Require MFEM-MGIS
if (NOT TARGET mfem-mgis::MFEMMGIS)
  find_package(MFEMMGIS REQUIRED)
endif(NOT TARGET mfem-mgis::MFEMMGIS)

# Define the executable to build for this example
add_executable(bimaterial bimaterial.cxx)
target_link_libraries(bimaterial mfem-mgis::MFEMMGIS)

# Copy needed file for execution
file(COPY
  ${PROJECT_SOURCE_DIR}/bimaterial.msh
  DESTINATION ${PROJECT_BINARY_DIR})

# Use mfront to generate libBehaviour.so
add_custom_command(
  OUTPUT "src/libBehaviour.so"
  COMMAND mfront --obuild --interface=generic
           ${PROJECT_SOURCE_DIR}/Elasticity.mfront
           ${PROJECT_SOURCE_DIR}/Plasticity.mfront
  DEPENDS ${PROJECT_SOURCE_DIR}/Elasticity.mfront
          ${PROJECT_SOURCE_DIR}/Plasticity.mfront
  COMMENT "Performing mfront code generation")

add_custom_target(generate-libBehaviour-ex9 ALL
  DEPENDS "src/libBehaviour.so")

add_dependencies(check bimaterial generate-libBehaviour-ex2)

# Define a test that will be called during the "make check"
set(test "bimaterial")
add_test(NAME ${test} COMMAND bimaterial)
