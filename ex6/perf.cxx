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
#include <functional>

// TEST CASE
#include <cas_cible_1.hxx>
#include <fissuration.hxx>


int try_several_solvers(TestParameters& p, const bool use_post_processing)
{
	if (mfem_mgis::getMPIrank() == 0)
		std::cout << "Number of processes: " << mfem_mgis::getMPIsize() << std::endl;

	gather_information data; // gather information (converged/residu/iterations) for all solver/precond used
	//constexpr bool all = false;
	constexpr bool all = false;
	//constexpr bool all = true;
	auto res = EXIT_SUCCESS;

	using std::function;
	typedef function<void(const TestParameters&, const bool, const solver_name, const precond_name, gather_information&)> regular_kernel;

	// build your solver / precond list
	// some of them are really slow without a precond such as X+GMRES
	// after a preliminary study, relevent solvers have been selected and stored in the "fast" traversal
	std::vector<solver_name> solver_traversal = get_solvers(p.tcase);
	std::vector<precond_name> pc_traversal = get_pc(p.tcase);
	regular_kernel fun = get_kernel<solver_name, precond_name>(p.tcase);
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
