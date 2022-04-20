cmake_minimum_required(VERSION 3.13)
project(InclusionsExample)


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

# Require timers
if(MFEM_USE_MPI)
	add_compile_options(-D__MPI)
	set(TIMERS timer.cpp timer.hpp)
	add_library(timers SHARED ${TIMERS})
	target_include_directories(timers PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
	find_package(MPI REQUIRED)
	target_link_libraries(timers PUBLIC MPI::MPI_CXX)
endif()



if(MFEM_USE_MPI)
  # Define the executable to build for this example
  add_executable(perf perf.cxx)
  target_link_libraries(perf mfem-mgis::MFEMMGIS)
  target_link_libraries(perf timers)

  target_include_directories(perf
      SYSTEM
      PUBLIC "$<BUILD_INTERFACE:${MFEM_INCLUDE_DIRS}>"
      PUBLIC "$<INSTALL_INTERFACE:${MFEM_INCLUDE_DIRS}>"
      )

  add_executable(reader read_markdown.cxx)
  target_link_libraries(reader mfem-mgis::MFEMMGIS)
  target_link_libraries(reader timers)

  target_include_directories(reader
      SYSTEM
      PUBLIC "$<BUILD_INTERFACE:${MFEM_INCLUDE_DIRS}>"
      PUBLIC "$<INSTALL_INTERFACE:${MFEM_INCLUDE_DIRS}>"
      )
  
  add_executable(one_test one_test.cxx)
  target_link_libraries(one_test mfem-mgis::MFEMMGIS)
  target_link_libraries(one_test timers)
  
  target_include_directories(one_test
      SYSTEM
      PUBLIC "$<BUILD_INTERFACE:${MFEM_INCLUDE_DIRS}>"
      PUBLIC "$<INSTALL_INTERFACE:${MFEM_INCLUDE_DIRS}>"
      )

  # Copy needed file for execution
  file(COPY
    ${PROJECT_SOURCE_DIR}/Elasticity.mfront
    DESTINATION ${PROJECT_BINARY_DIR})

  ## mesh
  file(COPY
	${PROJECT_SOURCE_DIR}/cube_2mat_per.mesh
    	DESTINATION ${PROJECT_BINARY_DIR})

  # Use mfront to generate libBehaviour.so
  add_custom_command(
    OUTPUT "src/libBehaviour.so"
    COMMAND mfront --obuild --interface=generic ${PROJECT_BINARY_DIR}/Elasticity.mfront
    DEPENDS ${PROJECT_BINARY_DIR}/Elasticity.mfront
    COMMENT "Performing mfront code generation")
  
  add_custom_target(generate-libBehaviour-cas-cible-1 ALL
    DEPENDS "src/libBehaviour.so")
  
  add_dependencies(check perf generate-libBehaviour-cas-cible-1)
  
  # Define a test that will be called during the "make check"
  function(add_perf)
	  set(test "perf")
    add_test(NAME ${test}_cube_2mat_per COMMAND mpirun -n 3 ./perf --mesh cube_2mat_per.mesh) #see later why it fails with -n 4 and it works with n = 1,2,3,5,6,7,44
  endfunction(add_perf)

  function(add_solver_preconditioner)
	  #no foreach
    add_test(NAME periodic_hypreCG_HypreDiagScale_cube_2mat_per COMMAND mpirun -n 3 ./one_test -ls 3 --preconditioner 0 --mesh cube_2mat_per.mesh) 
    add_test(NAME periodic_hypreCG_HypreEuclid_cube_2mat_per COMMAND mpirun -n 3 ./one_test -ls 3 --preconditioner 2 --mesh cube_2mat_per.mesh) 
    add_test(NAME periodic_hypreCG_HypreBoomerAMG_cube_2mat_per COMMAND mpirun -n 3 ./one_test -ls 3 --preconditioner 3 --mesh cube_2mat_per.mesh) 
    add_test(NAME periodic_hypreCG_HypreILU_cube_2mat_per COMMAND mpirun -n 3 ./one_test -ls 3 --preconditioner 4 --mesh cube_2mat_per.mesh) 
    add_test(NAME periodic_hypreCG_any_cube_2mat_per COMMAND mpirun -n 3 ./one_test -ls 3 --preconditioner 5 --mesh cube_2mat_per.mesh) 
    add_test(NAME periodic_CGSolver_HypreDiagScale_cube_2mat_per COMMAND mpirun -n 3 ./one_test -ls 6 --preconditioner 0 --mesh cube_2mat_per.mesh) 
    add_test(NAME periodic_CGSolver_HypreParaSails_cube_2mat_per COMMAND mpirun -n 3 ./one_test -ls 6 --preconditioner 1 --mesh cube_2mat_per.mesh) 
    add_test(NAME periodic_CGSolver_any_cube_2mat_per COMMAND mpirun -n 3 ./one_test -ls 6 --preconditioner 5 --mesh cube_2mat_per.mesh) 
  endfunction(add_solver_preconditioner)

  add_perf()
  add_solver_preconditioner()

endif()




