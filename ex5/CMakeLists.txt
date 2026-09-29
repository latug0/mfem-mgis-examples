cmake_minimum_required(VERSION 3.13)
project(SatohTestExample)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

include(CTest)
# testing
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

if(MFEM_USE_SUITESPARSE)
  # Copy the files needed at run time. Unlike file(COPY), configure_file updates
  # the copies when these files are modified.
  foreach(input IsotropicLinearThermoElasticity.mfront square.msh)
    configure_file(${PROJECT_SOURCE_DIR}/${input} ${PROJECT_BINARY_DIR}/${input}
                   COPYONLY)
  endforeach()
  
  # Use mfront to generate libBehaviour.so
  add_custom_command(
    OUTPUT "src/libBehaviour.so"
    COMMAND mfront --obuild --interface=generic ${PROJECT_BINARY_DIR}/IsotropicLinearThermoElasticity.mfront
    DEPENDS ${PROJECT_BINARY_DIR}/IsotropicLinearThermoElasticity.mfront
    COMMENT "Performing mfront code generation")
  
  add_custom_target(generate-libBehaviour-ex5 ALL
    DEPENDS "src/libBehaviour.so")

  add_executable(satoh SatohTest.cxx)
  target_link_libraries(satoh mfem-mgis::MFEMMGIS)
  target_include_directories(satoh
    SYSTEM
    PUBLIC "$<BUILD_INTERFACE:${MFEM_INCLUDE_DIRS}>"
    PUBLIC "$<INSTALL_INTERFACE:${MFEM_INCLUDE_DIRS}>")
    
  add_dependencies(check satoh generate-libBehaviour-ex5)

  # The test is the same in the full and restricted test modes, since it is a
  # single computation on a coarse mesh, whose reference reaction depends on
  # the mesh.
  set(test "satoh")
  add_test(NAME ${test} COMMAND ${test})
else()
  message(STATUS "Disabling satoh. Compile MFEM with SuiteSparse support to enable it")
endif()


