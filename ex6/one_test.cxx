/*!
 * \file   InclusionsEx.cxx
 * \brief
 * This example is modelling several inclusion within a periodic cube.
 *
 * Mechanical strain:
 *                 eps = E + grad_s v
 *
 *           with  E the given macrocoscopic strain
 *                 v the periodic displacement fluctuation
 * Displacement:
 *                   u = U + v
 *
 *           with  U the given displacement associated to E
 *                   E = grad_s U
 * The local microscopic strain is equal, on average, to the macroscopic strain:
 *           <eps> = <E>
 * \author Thomas Helfer, Guillaume Latu
 * \date   02/06/10/2021
 */

#include <memory>
#include <cstdlib>
#include <iostream>

#include<mfem_mgis_headers.hxx>

#include "timer.hpp"
#include "solver_name.hxx"
#include "precond_name.hxx"
#include "config_solver.hxx"
#include "data_gathering.hxx"
#include <test_parameters.hpp>

#include <cas_cible_1.hxx>
#include <fissuration.hxx>



int main(int argc, char* argv[]) 
{
        using std::function;
        typedef function<void(const TestParameters&, const bool, const solver_name, const precond_name, gather_information&)> regular_kernel;

	mfem_mgis::initialize(argc, argv);
	profiling::timers::init_timers();
	constexpr bool use_post_processing = false;

	auto p = parseCommandLineOptions_one_test(argc, argv);
	auto solver = solver_name(p.linearsolver);
	auto preconditioner = precond_name(p.preconditioner);
	regular_kernel kernel = get_kernel<solver_name,precond_name>(p.tcase);
	auto match = get_match(p.tcase);
	gather_information my_info;
	int res=EXIT_SUCCESS;

	if(match(solver,preconditioner))
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
