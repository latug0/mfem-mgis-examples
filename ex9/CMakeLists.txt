cmake_minimum_required(VERSION 3.13)
project(CermetElasticExample)

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
  ${PROJECT_SOURCE_DIR}/SaintVenantKirchhoffElasticity.mfront
  ${PROJECT_SOURCE_DIR}/geo
  ${PROJECT_SOURCE_DIR}/script
  DESTINATION ${PROJECT_BINARY_DIR})
  
  # Use mfront to generate libBehaviour.so
add_custom_command(
  OUTPUT "src/libBehaviour.so"
  COMMAND mfront --obuild --interface=generic ${PROJECT_BINARY_DIR}/SaintVenantKirchhoffElasticity.mfront
  DEPENDS ${PROJECT_BINARY_DIR}/SaintVenantKirchhoffElasticity.mfront
  COMMENT "Performing mfront code generation")
  
add_custom_target(generate-libBehaviour-ex9 ALL
  DEPENDS "src/libBehaviour.so")

add_executable(cermet_elasticity CermetElastic.cxx)
target_link_libraries(cermet_elasticity mfem-mgis::MFEMMGIS)
target_include_directories(cermet_elasticity
  SYSTEM
  PUBLIC "$<BUILD_INTERFACE:${MFEM_INCLUDE_DIRS}>"   
  PUBLIC "$<INSTALL_INTERFACE:${MFEM_INCLUDE_DIRS}>")
    
add_dependencies(check cermet_elasticity generate-libBehaviour-ex9)

set(test "CermetElasticTest")
add_test(NAME ${test} COMMAND cermet_elasticity --post-processing 0 --order 1 --mesh mesh/OneSphere.msh)

if(MFEM_USE_PETSC)
	add_test(NAME ${test}_PETSC COMMAND cermet_elasticity --post-processing 0 --order 1 --mesh mesh/OneSphere.msh --use-petsc true --petsc-configuration-file petscrc) 
endif()

