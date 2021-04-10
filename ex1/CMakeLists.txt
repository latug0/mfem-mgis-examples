cmake_minimum_required(VERSION 3.13)
project(UnilateralTensileTestExample)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Require MFEM
if (NOT TARGET mfem)
  find_package(MFEM REQUIRED)
endif(NOT TARGET mfem)

# Require MGIS
if (NOT TARGET mgis::MFrontGenericInterface)
  find_package(MFrontGenericInterface REQUIRED)
endif(NOT TARGET mgis::MFrontGenericInterface)

# Require MFEM-MGIS
if (NOT TARGET mfem-mgis::MFEMMGIS)
  find_package(MFEMMGIS REQUIRED)
endif(NOT TARGET mfem-mgis::MFEMMGIS)

# Specific setup to copy files from MFEM-MGIS main project
include(setup.cmake)

# Define the executable to build for this example

add_executable(UniaxialTensileTestEx UniaxialTensileTest.cxx)
target_link_libraries(UniaxialTensileTestEx mfem-mgis::MFEMMGIS)

# Use mfront to generate libBehaviour.so

add_custom_command(
  OUTPUT "src/libBehaviour.so"
  COMMAND mfront --obuild --interface=generic ${CMAKE_CURRENT_BINARY_DIR}/Plasticity.mfront
  DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/Plasticity.mfront
  COMMENT "Performing mfront code generation")

add_custom_target(generate-libBehaviour ALL
  DEPENDS "src/libBehaviour.so")

add_dependencies(check UniaxialTensileTestEx generate-libBehaviour)

# Define a test that will be called during the "make check"

function(add_uniaxial_tensile_test_ex)
  set(test "UniaxialTensileTestEx")
  add_test(NAME ${test} COMMAND UniaxialTensileTestEx)
endfunction(add_uniaxial_tensile_test_ex)

add_uniaxial_tensile_test_ex()




