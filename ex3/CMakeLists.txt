cmake_minimum_required(VERSION 3.13)
project(TwoLayerCubeExample)

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

# The sources of this example are installed by mfem-mgis
set(MFEMMGIS_EX3_DIR "${MFEMMGIS_DIR}/../examples/ex3")

if(MFEM_USE_MPI)
  # Define the executable to build for this example
  add_executable(two_layer_cube ${MFEMMGIS_EX3_DIR}/TwoLayerCube.cxx)
  target_link_libraries(two_layer_cube mfem-mgis::MFEMMGIS)
  target_include_directories(two_layer_cube
      SYSTEM
      PUBLIC "$<BUILD_INTERFACE:${MFEM_INCLUDE_DIRS}>"
      PUBLIC "$<INSTALL_INTERFACE:${MFEM_INCLUDE_DIRS}>")

  # Copy the files needed at run time. Unlike file(COPY), configure_file updates
  # the copies when these files are modified.
  foreach(input IsotropicLinearElasticity.mfront cube_2mat_per.mesh Box.med Box.per)
    configure_file(${MFEMMGIS_EX3_DIR}/${input} ${PROJECT_BINARY_DIR}/${input}
                   COPYONLY)
  endforeach()

  # Use mfront to generate libBehaviour.so
  add_custom_command(
    OUTPUT "src/libBehaviour.so"
    COMMAND mfront --obuild --interface=generic ${PROJECT_BINARY_DIR}/IsotropicLinearElasticity.mfront
    DEPENDS ${PROJECT_BINARY_DIR}/IsotropicLinearElasticity.mfront
    COMMENT "Performing mfront code generation")

  add_custom_target(generate-libBehaviour-ex3 ALL
    DEPENDS "src/libBehaviour.so")

  add_dependencies(check two_layer_cube generate-libBehaviour-ex3)

  # Define the tests that will be called during the "make check", one per
  # loading case, each one being compared to its analytical solution. They are
  # the same in the full and restricted test modes, since each loading case is
  # a single computation on a coarse mesh.
  function(add_two_layer_cube_ex tcase)
    set(test "two_layer_cube-${tcase}")
    add_test(NAME ${test} COMMAND two_layer_cube --test-case ${tcase})
  endfunction(add_two_layer_cube_ex)

  foreach(tcase RANGE 0 5)
    add_two_layer_cube_ex(${tcase})
  endforeach()

endif()
