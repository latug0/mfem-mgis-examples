#include <memory>
#include <cstdlib>
#include <iostream>

#include <common/mfem_mgis_headers.hxx>
#include <parameters/test_parameters.hpp>
#include <common/timer.hpp>
#include <solver/solver_name.hxx>
#include <solver/precond_name.hxx>
#include <solver/config_solver.hxx>
#include <common/data_gathering.hxx>
#include <functional>
#include <common/memory_footprint.hxx>

// TEST CASE
#include <test-cases/cas_cible_1.hxx>
#include <test-cases/fissuration.hxx>
#include <parameters/catch_functions.hxx>


int main(int argc, char* argv[]) 
{
        using std::function;
        typedef function<void(const TestParameters&, const bool, const solver_name, const precond_name, gather_information&)> regular_kernel;

	// mpi initialization here 
	mfem_mgis::initialize(argc, argv);

	// init timers
	profiling::timers::init_timers();

	// get parameters
	auto p = parseCommandLineOptions_one_test(argc, argv);

	// add post processing
	const bool use_post_processing = (p.post_processing == 1);

	// define solver name
	const auto solver = solver_name(p.linearsolver);

	// define preconditioner name
	const auto preconditioner = precond_name(p.preconditioner);

	// get main kernel depending on the test case
	regular_kernel kernel = get_kernel<solver_name,precond_name>(p.tcase);

	// function to check if the solver matches with the preconditioner depending on the test case
	auto match = get_match(p.tcase);

	// class to gather data (is converged?, residu, runtime)
	gather_information my_info;

	int res=EXIT_SUCCESS;

	// apply kernel if the solver matches with the preconditioner. 
	// some possibilities are removed depending on the refinement and the test case.
	if(match(solver,preconditioner, p.refinement))
	{
		 kernel(p, 
				use_post_processing, 
				solver, 
				preconditioner,
				my_info
			    );
	}
	else res=EXIT_FAILURE;

	// print the memory footprint (mpi reduction)
	common::memory::print_memory_footprint();

	// print and write timetable
	profiling::timers::print_and_write_timers();
	return res;
}
