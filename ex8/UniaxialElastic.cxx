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

#include <common.hxx>

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


	// // preconditionner hypreBoomerAMG
	 auto options = mfem_mgis::Parameters{{"VerbosityLevel", verbosity}};
	 auto preconditionner = mfem_mgis::Parameters{{"Name","HypreDiagScale"}, {"Options",options}};
	// // auto preconditionner = mfem_mgis::Parameters{{"Name","HypreBoomerAMG"}, {"Options",options}};
	 solverParameters.insert(mfem_mgis::Parameters{{"Preconditioner",preconditionner}});
	// solver HyprePCG
	 p.setLinearSolver("HyprePCG", solverParameters);
	//p.setLinearSolver("MUMPSSolver", solverParameters);
}

void setup_properties(const TestParameters& p, mfem_mgis::PeriodicNonLinearEvolutionProblem& problem)
{
	using namespace mgis::behaviour;
	using real=mfem_mgis::real;

	CatchTimeSection("set_mgis_stuff");

	// const int nMat = 1;
	const int nMat = getMaterialsAttributes(*(problem.getFiniteElementDiscretizationPointer())).Max();
	mfem_mgis::Profiler::Utils::Message("Nombre de matériaux : ", nMat);

	for(int i = 0 ; i < nMat ; i++)
	{
		problem.addBehaviourIntegrator("Mechanics", i+1, p.library, p.behaviour);
	}
}

template <typename Problem>
void simulation(Problem &problem, double start, double end, double dt, bool pp)
{
	CatchTimeSection("simulation");

	const int nMat = getMaterialsAttributes(*(problem.getFiniteElementDiscretizationPointer())).Max();
	const double young = 200.e9;
	const double poisson = 0.3;

	// materials
	auto set_properties = [](auto &m, const double yo, const double po)
	{
		setMaterialProperty(m.s0, "YoungModulus", yo);
		setMaterialProperty(m.s0, "PoissonRatio", po);
		setMaterialProperty(m.s1, "YoungModulus", yo);
		setMaterialProperty(m.s1, "PoissonRatio", po);
	};

	auto set_temperature = [](auto &m)
	{
		setExternalStateVariable(m.s0, "Temperature", 293.15);
		setExternalStateVariable(m.s1, "Temperature", 293.15);
	};

	for (int i = 0; i < nMat; i++)
	{
		auto &mat = problem.getMaterial(i + 1);
		set_properties(mat, young, poisson);
		set_temperature(mat);
	}

	const int nStep = (end - start) / dt;
	for (int i = 0; i < nStep; i++)
	{
		// Transformation gradient components
		const double EZZ = 0.02 * ((i + 1.) * dt / end);
		double EXX = (-poisson) * EZZ;
		double EYY = (-poisson) * EZZ;

		// init macro cauchy stress components
		double SIGXX = 0.;
		double SIGYY = 0.;

		// init volume var
		double vol;

		// init MPI vars
		double SIG_local[2] = {SIGXX, SIGYY};
		double SIG_global[2] = {0., 0.};

		// Fixed-Point param
		double tolFP = 1e-4;
		double oldRes = 1.;
		double newRes = std::sqrt(SIGXX * SIGXX + SIGYY * SIGYY);
		int itFP = 0;
		int maxitFP = 20;
		while ((abs(oldRes - newRes) > tolFP) && (itFP <= maxitFP))
		{
			mfem_mgis::Profiler::Utils::Message("SIGXX = ", SIGXX);
			mfem_mgis::Profiler::Utils::Message("SIGYY = ", SIGYY);
			problem.setMacroscopicGradientsEvolution([&EXX, &EYY, &EZZ, &SIGXX, &SIGYY, young](const double t) { 
				auto ret = std::vector<mfem_mgis::real>(6, mfem_mgis::real{0});
				ret[2] = EZZ;
				ret[1] = EYY - (SIGYY)/young;
				ret[0] = EXX - (SIGXX)/young;
				EXX = ret[0];
				EYY = ret[1];
				mfem_mgis::Profiler::Utils::Message("EXX", EXX);
				mfem_mgis::Profiler::Utils::Message("EYY", EYY);
				mfem_mgis::Profiler::Utils::Message("EZZ", EZZ);
				return ret;
			});

			run_solve(problem, i * dt, dt);
			auto [tf_integrals, volumes] = mfem_mgis::computeMeanThermodynamicForcesValues<true>(problem.template getImplementation<true>());

			vol = 0.;
			for (const auto &mi : problem.getAssignedMaterialsIdentifiers())
			{
				vol += volumes[mi];
			}
			MPI_Allreduce(MPI_IN_PLACE, &vol, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

			SIGXX = 0.;
			SIGYY = 0.;
			for (const auto &mi : problem.getAssignedMaterialsIdentifiers())
			{
				// volumes[mi] ;
				SIGXX += tf_integrals[mi][0] / vol;
				SIGYY += tf_integrals[mi][1] / vol;
			}

			SIG_local[0] = SIGXX;
			SIG_local[1] = SIGYY;
			MPI_Allreduce(SIG_local, SIG_global, 2, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
			SIGXX = SIG_global[0];
			SIGYY = SIG_global[1];

			oldRes = newRes;
			newRes = std::sqrt(SIGXX * SIGXX + SIGYY * SIGYY);
			mfem_mgis::Profiler::Utils::Message("FIXED POINT iteration", itFP, ": |res| =", abs(oldRes - newRes));

			itFP += 1;
		}

		if (itFP >= maxitFP)
			mfem_mgis::Profiler::Utils::Message("warning: maximum number of iterations for the fixed-point algorithm attained, before the requested tolerance is reached");

		problem.update();		
		if (pp)
			execute_post_processings(problem, i * dt, dt);
		
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
	print_mesh_information(problem.getImplementation<true>());

	// set problem
	setup_properties(p, problem);
	setLinearSolver(problem, p.verbosity_level);

	problem.setSolverParameters({{"VerbosityLevel", p.verbosity_level},
			{"RelativeTolerance", 1.e-4},
			{"AbsoluteTolerance", 1.e-4},
			{"MaximumNumberOfIterations", 15}});


	// add post processings
	if(use_post_processing) add_post_processings(problem, "OutputFile-Uniaxial-elastic", "avgStressElastic", "UniaxialElasticOrthoTestIntegrationPointOutput");

	// main function here
	int nStep=2;
	double start=0;
	double end=200;
	const double dt = (end-start)/nStep;
	simulation(problem, start, end, dt, use_post_processing);

	// print and write timetable
	mfem_mgis::Profiler::timers::print_and_write_timers();
	return(EXIT_SUCCESS);
}
