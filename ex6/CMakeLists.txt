cmake_minimum_required(VERSION 3.13)
project(RveNonLinearElasticExample)

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

# Define a test that will be called during the "make check"

# Copy needed file for execution
file(COPY
  ${PROJECT_SOURCE_DIR}/SaintVenantKirchhoffElasticity.mfront
  ${PROJECT_SOURCE_DIR}/cube_2mat_per.mesh
  ${PROJECT_SOURCE_DIR}/inclusions_49.geo
  DESTINATION ${PROJECT_BINARY_DIR})
  
  # Use mfront to generate libBehaviour.so
add_custom_command(
  OUTPUT "src/libBehaviour.so"
  COMMAND mfront --obuild --interface=generic ${PROJECT_BINARY_DIR}/SaintVenantKirchhoffElasticity.mfront
  DEPENDS ${PROJECT_BINARY_DIR}/SaintVenantKirchhoffElasticity.mfront
  COMMENT "Performing mfront code generation")
  
add_custom_target(generate-libBehaviour-ex6 ALL
  DEPENDS "src/libBehaviour.so")

add_executable(rve RveNonLinearElastic.cxx)
target_link_libraries(rve mfem-mgis::MFEMMGIS)
target_include_directories(rve
  SYSTEM
  PUBLIC "$<BUILD_INTERFACE:${MFEM_INCLUDE_DIRS}>"   
  PUBLIC "$<INSTALL_INTERFACE:${MFEM_INCLUDE_DIRS}>")
    
add_dependencies(check rve generate-libBehaviour-ex6)

set(test "RveNonLinearElasticTest")
add_test(NAME ${test} COMMAND rve --postprocessings 0)


