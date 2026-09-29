cmake_minimum_required(VERSION 3.13)
project(ssna303Example)

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

# Define the executable to build for this example
add_executable(ssna303 ssna303.cxx)
target_link_libraries(ssna303 mfem-mgis::MFEMMGIS)
target_include_directories(ssna303
    SYSTEM
    PUBLIC "$<BUILD_INTERFACE:${MFEM_INCLUDE_DIRS}>"
    PUBLIC "$<INSTALL_INTERFACE:${MFEM_INCLUDE_DIRS}>")

# Copy the files needed at run time. Unlike file(COPY), configure_file updates
# the copies when these files are modified.
foreach(input IsotropicLinearHardeningPlasticity.mfront ssna303.msh ssna303-force.ref
              ssna303-force-fbar.ref)
  configure_file(${PROJECT_SOURCE_DIR}/${input} ${PROJECT_BINARY_DIR}/${input}
                 COPYONLY)
endforeach()

# Use mfront to generate libBehaviour.so
add_custom_command(
  OUTPUT "src/libBehaviour.so"
  COMMAND mfront --obuild --interface=generic ${PROJECT_BINARY_DIR}/IsotropicLinearHardeningPlasticity.mfront
  DEPENDS ${PROJECT_BINARY_DIR}/IsotropicLinearHardeningPlasticity.mfront
  COMMENT "Performing mfront code generation")

add_custom_target(generate-libBehaviour-ex2 ALL
  DEPENDS "src/libBehaviour.so")

add_dependencies(check ssna303 generate-libBehaviour-ex2)

# Define the tests that will be called during the "make check". The optimised
# builds compute the complete loading. The debug and coverage builds, which are
# much slower, only compute its first two time steps, during which the plastic
# flow starts. In both cases, the results are compared to the first values of
# the reference files.
if(CMAKE_BUILD_TYPE MATCHES "^(Debug|Coverage)$")
  set(ssna303_loading --nbsteps 2 --end-time 0.05)
else()
  set(ssna303_loading --nbsteps 40 --end-time 1)
endif()

function(add_ssna303_test_ex)
  set(test "ssna303")
  add_test(NAME ${test} COMMAND ssna303 ${ssna303_loading}
                                   --reference-file ssna303-force.ref)
  if(MGIS_HAVE_TFEL)
    add_test(NAME ${test}-fbar COMMAND ssna303 --use-fbar ${ssna303_loading}
                                     --reference-file ssna303-force-fbar.ref
                                     --standard-reference-file ssna303-force.ref)
  endif(MGIS_HAVE_TFEL)
endfunction(add_ssna303_test_ex)

add_ssna303_test_ex()






