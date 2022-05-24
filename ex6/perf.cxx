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
#include "mfem/general/optparser.hpp"
#include "mfem/linalg/solvers.hpp"
#include "mfem/fem/datacollection.hpp"
#include "MFEMMGIS/MFEMForward.hxx"
#include "MFEMMGIS/Material.hxx"
#include "MFEMMGIS/Profiler.hxx"
#include "MFEMMGIS/AnalyticalTests.hxx"
#include "MFEMMGIS/NonLinearEvolutionProblemImplementation.hxx"
#include "MFEMMGIS/PeriodicNonLinearEvolutionProblem.hxx"

#ifdef MFEM_USE_PETSC
#include "mfem/linalg/petsc.hpp"
#endif /* MFEM_USE_PETSC */

#ifdef MFEM_USE_PETSC
#include "mfem/linalg/mumps.hpp"
#endif /* MFEM_USE_MUMPS */



constexpr double xmax = 1.;


struct TestParameters {
	const char* mesh_file = "cas-cible-1.mesh";
	const char* behaviour = "Elasticity";
	const char* library = "src/libBehaviour.so";
	const char* reference_file = "Elasticity.ref";
	int order = 1;
	int tcase = 1;
	int linearsolver = -1;
	double xmax = 1.;
	double ymax = 1.;
	double zmax = 1.;
	bool parallel = true;
};


#include "timer.hpp"
#include "solver_name.hxx"
#include "precond_name.hxx"
#include "config_solver.hxx"
#include "data_gathering.hxx"
#include "mfem_mgis_test.cxx"


// legacy
TestParameters parseCommandLineOptions(int& argc, char* argv[]) {
  START_TIMER("parse_command_line_options");
  TestParameters p;

  // options treatment
  mfem::OptionsParser args(argc, argv);
  args.AddOption(&p.mesh_file, "-m", "--mesh", "Mesh file to use.");
  args.AddOption(&p.library, "-l", "--library", "Material library.");
  args.AddOption(&p.order, "-o", "--order",
                 "Finite element order (polynomial degree).");
  args.AddOption(&p.tcase, "-t", "--test-case",
                 "identifier of the case : Exx->0, Eyy->1, Ezz->2, Exy->3, "
                 "Exz->4, Eyz->5");
  args.AddOption(
      &p.linearsolver, "-ls", "--linearsolver",
      "identifier of the linear solver: 0 -> GMRES, 1 -> CG, 2 -> UMFPack");
  args.Parse();
  if (!args.Good()) {
    if (mfem_mgis::getMPIrank() == 0)
      args.PrintUsage(std::cout);
    mfem_mgis::finalize();
    exit(0);
  }
  if (p.mesh_file == nullptr) {
    if (mfem_mgis::getMPIrank() == 0)
      std::cout << "ERROR: Mesh file missing" << std::endl;
    args.PrintUsage(std::cout);
    mfem_mgis::abort(EXIT_FAILURE);
  }
  if (mfem_mgis::getMPIrank() == 0)
    args.PrintOptions(std::cout);
  if ((p.tcase < 0) || (p.tcase > 5)) {
    std::cerr << "Invalid test case\n";
    mfem_mgis::abort(EXIT_FAILURE);
  }
  return p;
}

template<typename Func>
int try_several_solvers(TestParameters& p, const bool use_post_processing, Func & fun)
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
	std::vector<solver_name> solverTraversalAll = {
		solver_name::HyprePCG,
		solver_name::HypreFGMRES,
		solver_name::HypreGMRES,
		solver_name::MUMPSSolver,
//		solver_name::GMRESSolver,	// does not work with HypreBoomerAMG and does not converged without precond
		solver_name::CGSolver,
		solver_name::BiCGSTABSolver,
//		solver_name::UMFPackSolver,
		solver_name::MINRESSolver
	};

	std::vector<solver_name> solverTraversalFast = {
		solver_name::HypreFGMRES,
		solver_name::HypreGMRES,
		solver_name::HyprePCG,
//		solver_name::MUMPSSolver, // not really fast
		solver_name::CGSolver,
		solver_name::BiCGSTABSolver,
		solver_name::MINRESSolver
	};

	std::vector<precond_name> precondTraversal = 
	{
		precond_name::ANY,
		precond_name::HypreBoomerAMG,
		precond_name::HypreILU,
		precond_name::HypreEuclid,
	//	precond_name::HypreParaSails,
		precond_name::HypreDiagScale 
	};

	
	if constexpr (all)
	{
		for(auto solverIterator : solverTraversalAll)
		{
			START_TIMER(getName(solverIterator));
			for(auto precondIterator : precondTraversal)
			{
				if(match(solverIterator,precondIterator))
				{
					START_TIMER(getName(precondIterator));
					//executeMFEMMGISTest(
					fun(
						p,
						use_post_processing, 
						solverIterator, 
						precondIterator,
						data
					);
				}
			}
		}
	}
	else
	{
		for(auto solverIterator : solverTraversalFast)
		{
			START_TIMER(getName(solverIterator));
			for(auto precondIterator : precondTraversal)
			{
				if(match(solverIterator,precondIterator))
				{
					START_TIMER(getName(precondIterator));
					fun(
						p,
						use_post_processing, 
						solverIterator, 
						precondIterator,
						data
					);
				}
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
	//gather_information data;
	//auto res = executeMFEMMGISTest(p, use_post_processing, solver_name::HyprePCG, precond_name::ANY, data);
	auto res = try_several_solvers(p, use_post_processing, executeMFEMMGISTest);
	//auto res = try_several_solvers(p, use_post_processing, thomas_main_function);

	profiling::timers::print_and_write_timers();
	return res;
}
