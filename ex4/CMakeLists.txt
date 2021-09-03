cmake_minimum_required(VERSION 3.13)
project(Ssna303_3dExample)


set(CMAKE_CXX_STANDARD 17)
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


# Copy needed file for execution
file(COPY
  ${PROJECT_SOURCE_DIR}/Plasticity.mfront
  ${PROJECT_SOURCE_DIR}/ssna303_3d.msh
  ${PROJECT_SOURCE_DIR}/rc_ex10p
  DESTINATION ${PROJECT_BINARY_DIR})

# Use mfront to generate libBehaviour.so
add_custom_command(
  OUTPUT "src/libBehaviour.so"
  COMMAND mfront --obuild --interface=generic ${PROJECT_BINARY_DIR}/Plasticity.mfront
  DEPENDS ${PROJECT_BINARY_DIR}/Plasticity.mfront
  COMMENT "Performing mfront code generation")

add_custom_target(generate-libBehaviour-ex4 ALL
  DEPENDS "src/libBehaviour.so")

# Define a test that will be called during the "make check"
function(add_ssna303_3d_test_ex solver)
 set(test "Ssna303_3d_${solver}")
 add_test(NAME ${test} COMMAND ${test})
endfunction(add_ssna303_3d_test_ex)

foreach(tg IN ITEMS mumps hypre petsc)
  # Define the executable to build for this example
  add_executable("Ssna303_3d_${tg}" ssna303_${tg}.cxx)
  target_link_libraries(Ssna303_3d_${tg} mfem-mgis::MFEMMGIS)
  target_include_directories(Ssna303_3d_${tg}
    SYSTEM
    PUBLIC "$<BUILD_INTERFACE:${MFEM_INCLUDE_DIRS}>"
    PUBLIC "$<INSTALL_INTERFACE:${MFEM_INCLUDE_DIRS}>")

  add_dependencies(check Ssna303_3d_${tg} generate-libBehaviour-ex4)
  add_ssna303_3d_test_ex(${tg})
endforeach()



