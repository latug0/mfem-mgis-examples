#include<string>
#include<iostream>
#include "MFEMMGIS/MFEMForward.hxx"
#include "MFEMMGIS/Material.hxx"
#include "MFEMMGIS/Profiler.hxx"
#include "MFEMMGIS/AnalyticalTests.hxx"
#include "MFEMMGIS/NonLinearEvolutionProblemImplementation.hxx"
#include "MFEMMGIS/PeriodicNonLinearEvolutionProblem.hxx"
#ifdef MFEM_USE_PETSC
#include "mfem/linalg/petsc.hpp"
#endif /* MFEM_USE_PETSC */

#include<MFEMMGIS/Profiler.hxx>
#include<solver/solver_name.hxx>
#include<solver/precond_name.hxx>

// TEST CASE
#include <test-cases/cas_cible_1.hxx>
#include <test-cases/fissuration.hxx>
#include <parameters/catch_functions.hxx>
#include <memory>
#include <common/memory_footprint.hxx>

// make a petsc file
//
using namespace configuration;
using namespace mfem_mgis;
int main(int argc, char* argv[]) 
{
        using std::function;
        typedef function<void(const TestParameters&, const bool, const solver_name, const precond_name, gather_information&)> regular_kernel;
        typedef function<int(const TestParameters&, const bool, const char*, gather_information&)> petsc_kernel;

	// mpi initialization here 
	mfem_mgis::initialize(argc, argv);

	// init timers
	Profiler::timers::init_timers();

	// get parameters
	auto p = parseCommandLineOptions_one_test(argc, argv);

	// add post processing
	const bool use_post_processing = (p.post_processing == 1);

	// get main kernel depending on the test case
	petsc_kernel kernel = get_petsc_kernel_with_file(p.tcase);
    
	// define petsc file name
	const char* petsc_file_name = "test_hypre.in"; 

	// class to gather data (is converged?, residu, runtime)
	gather_information my_info;

	// apply kernel if the solver matches with the preconditioner. 
	// some possibilities are removed depending on the refinement and the test case.
	int res = kernel(p, 
			use_post_processing,
		        petsc_file_name,	
			my_info
	);

	// print the memory footprint (mpi reduction)
	common::memory::print_memory_footprint();

	Profiler::timers::print_and_write_timers();
	return res;
}
