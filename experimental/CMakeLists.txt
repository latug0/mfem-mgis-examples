cmake_minimum_required(VERSION 3.13)
project(InclusionsExample)


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

SET(HEARDERS include)



if(MFEM_USE_MPI)

	SET(HEARDERS "include")
	SET(utiles src/parameters/test_parameters.cpp src/solver/solver_name.cxx src/solver/precond_name.cxx src/common/data_gathering.cxx src/common/common.cxx src/parameters/catch_functions.cxx)
	add_library(labo-tools STATIC ${utiles})
	target_link_libraries(labo-tools mfem-mgis::MFEMMGIS)
	target_include_directories(labo-tools
		SYSTEM  
		PUBLIC "$<BUILD_INTERFACE:${MFEM_INCLUDE_DIRS}>"
		PUBLIC "$<INSTALL_INTERFACE:${MFEM_INCLUDE_DIRS}>"
		PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/include"
		)


	function(add_case exe name)
		target_link_libraries(${exe} ${name})
	endfunction()

	function(build_case name)
		add_library(${name} STATIC src/test-cases/${name}.cxx)
		target_link_libraries(${name} mfem-mgis::MFEMMGIS)
		target_include_directories(${name}
			SYSTEM  
			PUBLIC "$<BUILD_INTERFACE:${MFEM_INCLUDE_DIRS}>"
			PUBLIC "$<INSTALL_INTERFACE:${MFEM_INCLUDE_DIRS}>"
			PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/include"
			)
	endfunction()

	build_case(cas_cible_1)
	build_case(fissuration)
	build_case(check_manta_rev_elastic)


	# Define the executable to build for this example
	add_executable(perf perf.cxx)
	target_link_libraries(perf labo-tools)
	add_case(perf cas_cible_1)
	add_case(perf check_manta_rev_elastic)
	add_case(perf fissuration)

	add_executable(reader read_markdown.cxx)
	target_link_libraries(reader labo-tools)

	add_executable(one_test one_test.cxx)
	target_link_libraries(one_test labo-tools)

	add_case(one_test cas_cible_1)
	add_case(one_test check_manta_rev_elastic)
	add_case(one_test fissuration)

	add_executable(petsc_file_generation use_petsc.cxx)
	target_link_libraries(petsc_file_generation labo-tools)


	# Copy needed file for execution
	file(COPY
		${PROJECT_SOURCE_DIR}/mfront-law
		${PROJECT_SOURCE_DIR}/mesh
		DESTINATION ${PROJECT_BINARY_DIR})


	# Use mfront to generate libBehaviour.so
	add_custom_command(
		OUTPUT "src/libBehaviour-elasticity.so"
		COMMAND mfront --obuild --interface=generic ${PROJECT_BINARY_DIR}/mfront-law/Elasticity.mfront
		DEPENDS ${PROJECT_BINARY_DIR}/mfront-law/Elasticity.mfront
				COMMENT "Performing mfront code generation for cas cible 1")

	add_custom_command(
		OUTPUT "src/libBehaviour-SaintVenantKirchhoffElasticity.so"
		COMMAND mfront --obuild --interface=generic ${PROJECT_BINARY_DIR}/mfront-law/SaintVenantKirchhoffElasticity.mfront
		DEPENDS ${PROJECT_BINARY_DIR}/mfront-law/SaintVenantKirchhoffElasticity.mfront
				COMMENT "Performing mfront code generation for cas cible 2")
	
	add_custom_command(
		OUTPUT "src/libBehaviour-for-fissuration.so"
		COMMAND mfront --obuild --interface=generic ${PROJECT_BINARY_DIR}/mfront-law/AT1MicromorphicDamage.mfront ${PROJECT_BINARY_DIR}/mfront-law/MicromorphicDamageI.mfront
		DEPENDS ${PROJECT_BINARY_DIR}/mfront-law/AT1MicromorphicDamage.mfront ${PROJECT_BINARY_DIR}/mfront-law/MicromorphicDamageI.mfront
		COMMENT "Performing mfront code generation for fissuration")

		add_custom_target(generate-libBehaviour-cas-cible-1 ALL
				DEPENDS "src/libBehaviour-elasticity.so")

		add_custom_target(generate-libBehaviour-cas-cible-2 ALL
				DEPENDS "src/libBehaviour-SaintVenantKirchhoffElasticity.so")
	add_custom_target(generate-libBehaviour-fissuration ALL
		DEPENDS "src/libBehaviour-for-fissuration.so")
	
	add_dependencies(check perf generate-libBehaviour-fissuration generate-libBehaviour-cas-cible-1 generate-libBehaviour-cas-cible-2)

	# Define a test that will be called during the "make check"
	function(add_perf)
		set(test "perf")
		add_test(NAME ${test}_cube_2mat_per COMMAND mpirun -n 3 ./perf --mesh mesh/cube_2mat_per.mesh) #see later why it fails with -n 4 and it works with n = 1,2,3,5,6,7,44
	endfunction(add_perf)

	function(add_solver_preconditioner)
		#no foreach
		add_test(NAME periodic_hypreCG_HypreDiagScale_cube_2mat_per COMMAND mpirun -n 3 ./one_test -ls 3 --preconditioner 0 --mesh mesh/cube_2mat_per.mesh) 
		add_test(NAME periodic_hypreCG_HypreBoomerAMG_cube_2mat_per COMMAND mpirun -n 3 ./one_test -ls 3 --preconditioner 3 --mesh mesh/cube_2mat_per.mesh) 
		add_test(NAME periodic_hypreCG_HypreILU_cube_2mat_per COMMAND mpirun -n 3 ./one_test -ls 3 --preconditioner 4 --mesh mesh/cube_2mat_per.mesh) 
		add_test(NAME periodic_hypreCG_any_cube_2mat_per COMMAND mpirun -n 3 ./one_test -ls 3 --preconditioner 5 --mesh mesh/cube_2mat_per.mesh) 
		add_test(NAME periodic_CGSolver_HypreDiagScale_cube_2mat_per COMMAND mpirun -n 3 ./one_test -ls 6 --preconditioner 0 --mesh mesh/cube_2mat_per.mesh) 
		add_test(NAME periodic_CGSolver_HypreParaSails_cube_2mat_per COMMAND mpirun -n 3 ./one_test -ls 6 --preconditioner 1 --mesh mesh/cube_2mat_per.mesh) 
		add_test(NAME periodic_CGSolver_any_cube_2mat_per COMMAND mpirun -n 3 ./one_test -ls 6 --preconditioner 5 --mesh mesh/cube_2mat_per.mesh) 
	endfunction(add_solver_preconditioner)

	add_perf()
	add_solver_preconditioner()

endif()




