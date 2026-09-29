cmake_minimum_required(VERSION 3.13)
project(ssna303_3dExample)


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


# Copy the files needed at run time. Unlike file(COPY), configure_file updates
# the copies when these files are modified.
foreach(input Plasticity.mfront ssna303_3d.msh rc_ex10p ssna303_3d-force-1step.ref ssna303_3d-force-10steps.ref)
  configure_file(${PROJECT_SOURCE_DIR}/${input} ${PROJECT_BINARY_DIR}/${input}
                 COPYONLY)
endforeach()

# Use mfront to generate libBehaviour.so
add_custom_command(
  OUTPUT "src/libBehaviour.so"
  COMMAND mfront --obuild --interface=generic ${PROJECT_BINARY_DIR}/Plasticity.mfront ${PROJECT_BINARY_DIR}/FiniteStrainPlasticity.mfront
  DEPENDS ${PROJECT_BINARY_DIR}/Plasticity.mfront
  DEPENDS ${PROJECT_BINARY_DIR}/FiniteStrainPlasticity.mfront
  COMMENT "Performing mfront code generation")

add_custom_target(generate-libBehaviour-ex4 ALL
  DEPENDS "src/libBehaviour.so")

# Define a test that will be called during the "make check"
function(add_ssna303_3d_test_ex solver nsteps end_time reference_file use_fbar)
 set(test "ssna303_3d_${solver}")
 # the first time step of the default loading, which lasts several minutes
 add_test(NAME ${test} COMMAND ${test} --nbsteps ${nsteps} --end-time ${end_time}
                                  --reference-file ${reference_file})
 if(${use_fbar})
   if(MGIS_HAVE_TFEL)
     add_test(NAME ${test}-fbar COMMAND ${test} --use-fbar
                                                 --nbsteps ${nsteps} --end-time ${end_time}
                                                 --reference-file ${reference_file})
   endif(MGIS_HAVE_TFEL)
 endif()
endfunction(add_ssna303_3d_test_ex)

foreach(tg IN ITEMS mumps hypre petsc)
  # Define the executable to build for this example
  if(((NOT MFEM_USE_PETSC) AND (${tg} STREQUAL "petsc")) OR
      ((NOT MFEM_USE_MUMPS) AND (${tg} STREQUAL "mumps")))
  else()
    add_executable("ssna303_3d_${tg}" ssna303_${tg}.cxx)
    target_link_libraries(ssna303_3d_${tg} mfem-mgis::MFEMMGIS)
    target_include_directories(ssna303_3d_${tg}
      SYSTEM
      PUBLIC "$<BUILD_INTERFACE:${MFEM_INCLUDE_DIRS}>"
      PUBLIC "$<INSTALL_INTERFACE:${MFEM_INCLUDE_DIRS}>")
    add_dependencies(check ssna303_3d_${tg} generate-libBehaviour-ex4)
    if(${tg} STREQUAL "petsc")
      add_ssna303_3d_test_ex(${tg} 1 0.02 ssna303_3d-force-1step.ref OFF)
    else()
      if(${tg} STREQUAL "mumps")
        add_ssna303_3d_test_ex(${tg} 10 1 ssna303_3d-force-10steps.ref ON)
      else()
        add_ssna303_3d_test_ex(${tg} 10 1 ssna303_3d-force-10steps.ref OFF)
      endif()
    endif()
  endif()
endforeach()
