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

// TEST CASE
#include <test-cases/cas_cible_1.hxx>
#include <test-cases/fissuration.hxx>
#include <parameters/catch_functions.hxx>


int main(int argc, char* argv[]) 
{
        using std::function;
        typedef function<void(const TestParameters&, const bool, const solver_name, const precond_name, gather_information&)> regular_kernel;

	mfem_mgis::initialize(argc, argv);
	profiling::timers::init_timers();

	auto p = parseCommandLineOptions_one_test(argc, argv);
	const bool use_post_processing = (p.post_processing == 1);
	const auto solver = solver_name(p.linearsolver);
	const auto preconditioner = precond_name(p.preconditioner);
	regular_kernel kernel = get_kernel<solver_name,precond_name>(p.tcase);
	auto match = get_match(p.tcase);
	gather_information my_info;
	int res=EXIT_SUCCESS;

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

	profiling::timers::print_and_write_timers();
	return res;
}
