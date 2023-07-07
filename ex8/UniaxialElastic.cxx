/*!
 * \file   Mono_UO2_CosH_Jaco.cxx
 * \brief
 * This example is modelling ...
 * \author Raphael Prat, Maxence Wangermez
 * \date   06/2023
 */

#include <memory>
#include <cstdlib>
#include <iostream>
#include <csignal>

#include "mfem/general/optparser.hpp"
#include "mfem/linalg/solvers.hpp"
#include "mfem/fem/datacollection.hpp"
#include "MFEMMGIS/MFEMForward.hxx"
#include "MFEMMGIS/Material.hxx"
#include "MFEMMGIS/Profiler.hxx"
#include "MFEMMGIS/AnalyticalTests.hxx"
#include "MFEMMGIS/NonLinearEvolutionProblemImplementation.hxx"
#include "MFEMMGIS/PeriodicNonLinearEvolutionProblem.hxx"
#include "MFEMMGIS/Config.hxx"
#include "MFEMMGIS/Parameters.hxx"
#include "MFEMMGIS/PartialQuadratureSpace.hxx"
#include "MFEMMGIS/ImposedDirichletBoundaryConditionAtClosestNode.hxx"
#include "MFEMMGIS/NonLinearEvolutionProblem.hxx"

#ifdef MFEM_USE_PETSC
#include "mfem/linalg/petsc.hpp"
#endif /* MFEM_USE_PETSC */

#ifdef MFEM_USE_PETSC
#include "mfem/linalg/mumps.hpp"
#endif /* MFEM_USE_MUMPS */

#include <MFEMMGIS/Profiler.hxx>
#include <functional>

/* 
		Problem :
		Parameters : 
*/


// We need this class for test case sources
struct TestParameters {
	const char* mesh_file = "periodic-cube.msh";
	const char* behaviour = "Elasticity";
	const char* library = "src/libBehaviour.so";
	int order = 1;
	bool parallel = true;
	int refinement = 0;
	int post_processing = 1; // default value : disabled
	int verbosity_level = 0; // default value : lower level
};

void common_parameters(mfem::OptionsParser& args, TestParameters& p)
{
	args.AddOption(&p.mesh_file, "-m", "--mesh", "Mesh file to use.");
	args.AddOption(&p.library, "-l", "--library", "Material library.");
	args.AddOption(&p.order, "-o", "--order", "Finite element order (polynomial degree).");
	args.AddOption(&p.refinement, "-r", "--refinement", "refinement level of the mesh, default = 0");
	args.AddOption(&p.post_processing, "-p", "--post-processing", "run post processing step");
	args.AddOption(&p.verbosity_level, "-v", "--verbosity-level", "choose the verbosity level");

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
}

	template<typename Problem>
void add_post_processings(Problem& p, std::string msg)
{
	p.addPostProcessing(
			"ParaviewExportResults",
			{{"OutputFileName", msg}}
			);
	p.addPostProcessing(
			"MeanThermodynamicForces",
			{{"OutputFileName", "avgStressElastic"}});
	p.addPostProcessing(
			"ParaviewExportIntegrationPointResultsAtNodes",
			{
				{"OutputFileName", "UniaxialElasticTestIntegrationPointOutput"},
				{"Results", {"Stress"}}
			}
			);
} // end timer add_postprocessing_and_outputs

	template<typename Problem>
void execute_post_processings(Problem& p, double start, double end)
{
	CatchTimeSection("common::post_processing_step");
	p.executePostProcessings(start, end);
}

void setup_properties(const TestParameters& p, mfem_mgis::PeriodicNonLinearEvolutionProblem& problem)
{
	using namespace mgis::behaviour;
	using real=mfem_mgis::real;

	CatchTimeSection("set_mgis_stuff");

	const int nMat = 1;

	for(int i = 0 ; i < nMat ; i++)
	{
		problem.addBehaviourIntegrator("Mechanics", i+1, p.library, p.behaviour);
	}
	// materials

//	auto& m1 = problem.getMaterial(1);
	auto set_properties = [](auto& m, const double yo, const double po) {
		setMaterialProperty(m.s0, "YoungModulus", yo);
		setMaterialProperty(m.s0, "PoissonRatio", po);
		setMaterialProperty(m.s1, "YoungModulus", yo);
		setMaterialProperty(m.s1, "PoissonRatio", po);
	};

	auto set_temperature = [](auto& m) {
		setExternalStateVariable(m.s0, "Temperature", 293.15);
		setExternalStateVariable(m.s1, "Temperature", 293.15);
	};

	for(int i = 0 ; i < nMat ; i++)
	{
		auto& mat = problem.getMaterial(i+1);
		set_properties(mat, 200.e9, 0.3);
		set_temperature(mat);
	}
	
	problem.setMacroscopicGradientsEvolution([](const double t) { 
		const auto eps = 0.02*t/200; // final eps = 0.02
		auto ret = std::vector<real>(6, real{0});
		ret[2] = eps;
		ret[1] = -0.3*eps;
		ret[0] = ret[1];
		return ret; 
	});

}


	template<typename Problem>		
static void setLinearSolver(Problem& p,
		const int verbosity = 0,
		const mfem_mgis::real Tol = 1e-12
		)
{
	CatchTimeSection("set_linear_solver");
	// pilote
	constexpr int defaultMaxNumOfIt	 	= 5000; 		// MaximumNumberOfIterations
	auto solverParameters = mfem_mgis::Parameters{};
	solverParameters.insert(mfem_mgis::Parameters{{"VerbosityLevel", verbosity}});
	solverParameters.insert(mfem_mgis::Parameters{{"MaximumNumberOfIterations", defaultMaxNumOfIt}});
	// solverParameters.insert(mfem_mgis::Parameters{{"AbsoluteTolerance", Tol}});
	// solverParameters.insert(mfem_mgis::Parameters{{"RelativeTolerance", Tol}});
	solverParameters.insert(mfem_mgis::Parameters{{"Tolerance", Tol}});


	// preconditionner hypreBoomerAMG
	auto options = mfem_mgis::Parameters{{"VerbosityLevel", verbosity}};
	auto preconditionner = mfem_mgis::Parameters{{"Name","HypreDiagScale"}, {"Options",options}};
	// auto preconditionner = mfem_mgis::Parameters{{"Name","HypreBoomerAMG"}, {"Options",options}};
	solverParameters.insert(mfem_mgis::Parameters{{"Preconditioner",preconditionner}});
	// solver HyprePCG
	p.setLinearSolver("HyprePCG", solverParameters);
	// p.setLinearSolver("CGSolver", solverParameters);
}

	template<typename Problem>
void run_solve(Problem& p, double start, double dt)
{
	CatchTimeSection("Solve");
	mfem_mgis::Profiler::Utils::Message("Solving problem from ",start,"to",start+dt);
	// solving the problem
	auto statistics = p.solve(start, dt);
	// check status
	if (!statistics.status) {
		mfem_mgis::Profiler::Utils::Message("INFO: SOLVE FAILED");
		std::abort();
	}
}

int main(int argc, char* argv[]) 
{
	// mpi initialization here 
	mfem_mgis::initialize(argc, argv);

	// init timers
	mfem_mgis::Profiler::timers::init_timers();

	// get parameters
	TestParameters p;
	mfem::OptionsParser args(argc, argv);
	common_parameters(args, p);

	// add post processing
	const bool use_post_processing = (p.post_processing == 1);

	// 3D
	constexpr const auto dim = mfem_mgis::size_type{3};

	// creating the finite element workspace
	auto fed = std::make_shared<mfem_mgis::FiniteElementDiscretization>(
			mfem_mgis::Parameters{{"MeshFileName", p.mesh_file},
			{"FiniteElementFamily", "H1"},
			{"FiniteElementOrder", p.order},
			{"UnknownsSize", dim},
			{"NumberOfUniformRefinements", p.parallel ? p.refinement : 0},
			{"Parallel", p.parallel}});
	
	mfem_mgis::PeriodicNonLinearEvolutionProblem problem(fed);

	// set problem
	setup_properties(p, problem);
	setLinearSolver(problem, p.verbosity_level);

	problem.setSolverParameters({{"VerbosityLevel", p.verbosity_level},
			{"RelativeTolerance", 1.e-4},
			{"AbsoluteTolerance", 1.e-4},
			{"MaximumNumberOfIterations", 15}});


	// add post processings
	if(use_post_processing) add_post_processings(problem, "OutputFile-Uniaxial-elastic");

	// main function here
	int nStep=1;
	double start=0;
	double end=200;
	const double dt = (end-start)/nStep;
	for(int i = 0 ; i < nStep ; i++)
	{
		run_solve(problem, i * dt, dt);
		if(use_post_processing)	execute_post_processings(problem, i * dt, dt);
		problem.update();
	}

	// print and write timetable
	mfem_mgis::Profiler::timers::print_and_write_timers();
	return(EXIT_SUCCESS);
}
