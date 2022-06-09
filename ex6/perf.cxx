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


#include "test_parameters.hpp"
#include "timer.hpp"
#include "solver_name.hxx"
#include "precond_name.hxx"
#include "config_solver.hxx"
#include "data_gathering.hxx"


// TEST CASE
#include <cas_cible_1.hxx>
#include <fissuration.hxx>


using namespace configuration;

auto get_test_case(int a_case)
{
	switch(a_case)
	{
		case 1: return cas_cible_1::kernel;
		case 2: return fissuration::kernel;
		default: return cas_cible_1::kernel;
	}
}

std::vector<solver_name> get_solvers(int a_case)
{
	switch(a_case)
	{
		case 1: return cas_cible_1::build_solvers_list();
		case 2: return fissuration::build_solvers_list();
		default: return cas_cible_1::build_solvers_list();
	}
}

std::vector<precond_name> get_pc(int a_case)
{
	switch(a_case)
	{
		case 1: return cas_cible_1::build_pc_list();
		case 2: return fissuration::build_pc_list();
		default: return cas_cible_1::build_pc_list();
	}
}

auto get_match(int a_case)
{
	switch(a_case)
	{
		case 1: return cas_cible_1::match;
		case 2: return fissuration::match;
		default: return cas_cible_1::match;
	}
}

#ifdef MFEM_USE_PETSC
auto get_solvers_with_petsc(int a_case)
{
	switch(a_case)
	{
		case 1: return cas_cible_1::build_solvers_list_with_petsc();
		case 2: return fissuration::build_solvers_list_with_petsc();
		default: return cas_cible_1::build_solvers_list_with_petsc();
	}
}
#endif

int try_several_solvers(TestParameters& p, const bool use_post_processing)
{
	if (mfem_mgis::getMPIrank() == 0)
		std::cout << "Number of processes: " << mfem_mgis::getMPIsize() << std::endl;

	gather_information data; // gather information (converged/residu/iterations) for all solver/precond used
	//constexpr bool all = false;
	constexpr bool all = false;
	//constexpr bool all = true;
	auto res = EXIT_SUCCESS;

	// build your solver / precond list
	// some of them are really slow without a precond such as X+GMRES
	// after a preliminary study, relevent solvers have been selected and stored in the "fast" traversal
	std::vector<solver_name> solver_traversal = get_solvers(p.tcase);
	std::vector<precond_name> pc_traversal = get_pc(p.tcase);
	auto fun = get_test_case(p.tcase);
	auto match = get_match(p.tcase);

	for(auto solver_iterator : solver_traversal)
	{
		START_TIMER(getName(solver_iterator));
		for(auto precond_iterator : pc_traversal)
		{
			if(match(solver_iterator, precond_iterator))
			{
				START_TIMER(getName(precond_iterator));
				fun(
						p,
						use_post_processing, 
						solver_iterator, 
						precond_iterator,
						data
				   );
			}
		}
	}

	data.write();
	data.writeMD();
	return res;
}


int main(int argc, char* argv[]) 
{
	mfem_mgis::initialize(argc, argv);
	profiling::timers::init_timers();
	constexpr bool use_post_processing = false;

	auto p = parseCommandLineOptions(argc, argv);
	auto res = try_several_solvers(p, use_post_processing);

	profiling::timers::print_and_write_timers();
	return res;
}
