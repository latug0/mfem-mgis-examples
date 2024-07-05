cmake_minimum_required(VERSION 3.13)
project(MixedOxideFuelsExample)

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

# Define a test that will be called during the "make check"

# Copy needed file for execution
file(COPY
  ${PROJECT_SOURCE_DIR}/matrixlaw.mfront
  ${PROJECT_SOURCE_DIR}/mesh
  ${PROJECT_SOURCE_DIR}/file_geo
  ${PROJECT_SOURCE_DIR}/script_merope
  DESTINATION ${PROJECT_BINARY_DIR})
  
  # Use mfront to generate libBehaviour.so
add_custom_command(
  OUTPUT "src/libBehaviour.so"
  COMMAND mfront --obuild --interface=generic ${PROJECT_BINARY_DIR}/matrixlaw.mfront
  DEPENDS ${PROJECT_BINARY_DIR}/matrixlaw.mfront
  COMMENT "Performing mfront code generation")
  
add_custom_target(generate-libBehaviour-ex7 ALL
  DEPENDS "src/libBehaviour.so")

add_executable(mox2 MixedOxideFuels.cxx)
target_link_libraries(mox2 mfem-mgis::MFEMMGIS)
target_include_directories(mox2
  SYSTEM
  PUBLIC "$<BUILD_INTERFACE:${MFEM_INCLUDE_DIRS}>"   
  PUBLIC "$<INSTALL_INTERFACE:${MFEM_INCLUDE_DIRS}>")
    
add_dependencies(check mox2 generate-libBehaviour-ex7)

set(test "MixedOxideFuelsTest")
add_test(NAME ${test} COMMAND mox2 --postprocessings 0 --order 1 --mesh mesh/OneSphere.msh)


