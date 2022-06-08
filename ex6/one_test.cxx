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

#include "timer.hpp"
#include "solver_name.hxx"
#include "precond_name.hxx"
#include "config_solver.hxx"
#include "data_gathering.hxx"


constexpr double xmax = 1.;


struct TestParameters {
	const char* mesh_file = "cas-cible-1.mesh";
	const char* behaviour = "Elasticity";
	const char* library = "src/libBehaviour.so";
	const char* reference_file = "Elasticity.ref";
	int order = 1;
	int tcase = 1;
	int linearsolver = -1;
	int preconditioner = -1;
	double xmax = 1.;
	double ymax = 1.;
	double zmax = 1.;
	bool parallel = true;
};

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
      "identifier of the linear solver: 0 -> Default, 1 -> HypreFGMRES, 2 -> HypreGMRES, 3-> HyprePCG, 4 -> MUMPSSolver, 5 -> GMRESSolver, 6 -> CGSolver, 7 -> BiCGSTABSolver, 8 -> UMFPackSolver, 9 ->	MINRESSolver, 10 -> SLISolver");

  args.AddOption(
      &p.preconditioner, "-lp", "--preconditioner",
      " 0 -> HypreDiagScale, 1 -> HypreParaSails, 2 -> HypreEuclid, 3 -> HypreBoomerAMG, 4 -> HypreILU, 5 -> ANY");
  
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

int execute(const TestParameters& p, const bool use_post_processing, const solver_name a_solv, const precond_name a_precond) {
	constexpr const auto dim = mfem_mgis::size_type{3};
	
	//std::string string_solver = std::to_string(solv);
	std::string string_solver 	= getName(a_solv);
	std::string string_precond 	= getName(a_precond);
	std::string timer_name 		= "run_" + string_solver + "_" + string_precond;
	START_TIMER(timer_name);

	// creating the finite element workspace
	auto fed = std::make_shared<mfem_mgis::FiniteElementDiscretization>(
	mfem_mgis::Parameters{{"MeshFileName", p.mesh_file},
			    {"FiniteElementFamily", "H1"},
			    {"FiniteElementOrder", p.order},
			    {"UnknownsSize", dim},
			    {"NumberOfUniformRefinements", p.parallel ? 0 : 0},
			    {"Parallel", p.parallel}});

	mfem_mgis::PeriodicNonLinearEvolutionProblem problem(fed);

	{
		START_TIMER("set_mgis_stuff");
		problem.addBehaviourIntegrator("Mechanics", 1, p.library, "Elasticity");
		problem.addBehaviourIntegrator("Mechanics", 2, p.library, "Elasticity");
		// materials
		auto& m1 = problem.getMaterial(1);
		auto& m2 = problem.getMaterial(2);
		// setting the material properties
		auto set_properties = [](auto& m, const double l, const double mu) {
			mgis::behaviour::setMaterialProperty(m.s0, "FirstLameCoefficient", l);
			mgis::behaviour::setMaterialProperty(m.s0, "ShearModulus", mu);
			mgis::behaviour::setMaterialProperty(m.s1, "FirstLameCoefficient", l);
			mgis::behaviour::setMaterialProperty(m.s1, "ShearModulus", mu);
		};

		std::array<mfem_mgis::real,2> lambda({100, 200});
		std::array<mfem_mgis::real,2>     mu({75 , 150});
		set_properties(m1, lambda[0], mu[0]);
		set_properties(m2, lambda[1], mu[1]);
		//
		auto set_temperature = [](auto& m) {
			mgis::behaviour::setExternalStateVariable(m.s0, "Temperature", 293.15);
			mgis::behaviour::setExternalStateVariable(m.s1, "Temperature", 293.15);
		};
		set_temperature(m1);
		set_temperature(m2);

		// macroscopic strain
		std::vector<mfem_mgis::real> e(6, mfem_mgis::real{});
		if (p.tcase < 3) {
			e[p.tcase] = 1;
		} else {
			e[p.tcase] = 1.41421356237309504880 / 2;
		}
		problem.setMacroscopicGradientsEvolution([e](const double) { return e; });
	} // end timer set_mgis_stuff
	

	setLinearSolver(problem, a_solv, a_precond);
	setSolverParameters(problem);

	if(use_post_processing)
	{
		START_TIMER("add_postprocessing_and_outputs");
		problem.addPostProcessing(
			"ParaviewExportResults",
			{{"OutputFileName", "PeriodicTestOutput-" + std::to_string(p.tcase)}}	
		);
	} // end timer add_postprocessing_and_outputs

	mfem_mgis::NonLinearResolutionOutput solverStatistics;
	double measure = 0.;

	{
		START_TIMER("Solve");
		// solving the problem
		profiling::output::printMessage("Solving ...");
		measure = profiling::timers::chrono_section( [&](){
			solverStatistics = problem.solve(0, 1);
		});
		profiling::output::printMessage("... end solving");

		// check status
		if (!solverStatistics) {
			profiling::output::printMessage("INFO: ", string_solver,"+",string_precond," FAILED");
//			mfem_mgis::abort(EXIT_FAILURE);
		}
	} // end timer Solve

	if(use_post_processing)
	{
		START_TIMER("Postprocessing_step");
		problem.executePostProcessings(0, 1);
	} // en timer Postprocessing_step

	{
		START_TIMER("get_statistics");
		measure = profiling::output::reduce_max(measure);
		profiling::output::printMessage(" solver: "	, getName(a_solv), 
					" preconditioner: "	, getName(a_precond), 
					" ite: "		, solverStatistics.iterations, 
					" res: "		, solverStatistics.final_residual_norm, 
					" time: "		, measure);
	} // end timer get_statistics
	return(EXIT_SUCCESS);
}




int main(int argc, char* argv[]) 
{
	mfem_mgis::initialize(argc, argv);
	profiling::timers::init_timers();
	constexpr bool use_post_processing = false;

	auto p = parseCommandLineOptions(argc, argv);
	auto solver = solver_name(p.linearsolver);
	auto preconditioner = precond_name(p.preconditioner);
	int res;
	if(match(solver,preconditioner))
	{
		res = execute(p, 
				use_post_processing, 
				solver, 
				preconditioner);
	}
	else 
	{
		res = EXIT_FAILURE;
	}

	profiling::timers::print_and_write_timers();
	return res;
}
