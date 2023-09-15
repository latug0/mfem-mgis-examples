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
	const char* behaviour = "OrthotropicElastic";
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

auto norm(std::array<mfem_mgis::real, 3u> &u)
{
	mfem_mgis::real ret = 0.;
	for (auto val : u)
	{
		ret += val * val;
	}
	return (std::sqrt(ret));
}

std::vector<std::array<mfem_mgis::real, 3u>> readVectorsFromFile(const std::string &filename)
{
	// TODO: check if vectors size is equal to nMat
	std::vector<std::array<mfem_mgis::real, 3u>> vectors;
	std::ifstream inputFile(filename);

	auto checknorm = [](std::array<mfem_mgis::real, 3u> v)
	{
		mfem_mgis::real normv = norm(v);
		if (normv < 1.e-15)
		{
			throw std::invalid_argument("checknorm : vectors must not be null");
		}
		if (abs(normv - 1.) > 1.e-15)
		{
			std::cout << "checknorm : normalizing vector..." << std::endl;
			for (auto &comp : v)
			{
				comp *= 1. / normv;
			}
		}
		return v;
	};

	if (!inputFile.is_open())
	{
		std::cout << "Erreur lors de l'ouverture du fichier " << filename << std::endl;
		return vectors;
	}

	std::string line;
	while (std::getline(inputFile, line))
	{
		std::istringstream iss(line);
		std::array<mfem_mgis::real, 3u> vector;
		if (!(iss >> vector[0] >> vector[1] >> vector[2]))
		{
			std::cout << "Erreur de lecture du vecteur : " << line << std::endl;
			continue;
		}
		vectors.push_back(checknorm(vector));
	}

	inputFile.close();
	return vectors;
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
			{{"OutputFileName", "avgStressElasticOrtho"}});
	p.addPostProcessing(
			"ParaviewExportIntegrationPointResultsAtNodes",
			{
				{"OutputFileName", "UniaxialElasticOrthoTestIntegrationPointOutput"},
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

	// const int nMat = 1;
	const int nMat = getMaterialsAttributes(*(problem.getFiniteElementDiscretizationPointer())).Max();
	mfem_mgis::Profiler::Utils::Message("Nombre de matériaux : ", nMat);

	for(int i = 0 ; i < nMat ; i++)
	{
		problem.addBehaviourIntegrator("Mechanics", i+1, p.library, p.behaviour);
	}

	auto rot = [](std::array<mfem_mgis::real, 3u> u,
				  std::array<mfem_mgis::real, 3u> v)
	{
		if (std::inner_product(u.begin(), u.end(), v.begin(), 0.) > 1.e-15)
		{
			throw std::invalid_argument("rot : vectors are not orthogonals");
		}
		std::array<mfem_mgis::real, 9u> ret = {u[0], v[0], u[1] * v[2] - u[2] * v[1],
											   u[1], v[1], u[2] * v[0] - u[0] * v[2],
											   u[2], v[2], u[0] * v[1] - u[1] * v[0]};
		return ret;
	};

	for (int i = 0; i < nMat; i++)
	{
		auto &mat = problem.getMaterial(i + 1);
		if (mat.b.symmetry == mgis::behaviour::Behaviour::ORTHOTROPIC)
		{
			std::array<mfem_mgis::real, 9u> r = {0, -1, 0,  //
                                           		 1, 0, 0,  //
                                           		 0, 0, 1};
			mat.setRotationMatrix(mfem_mgis::RotationMatrix3D{r});
		}
	}

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
	// solverParameters.insert(mfem_mgis::Parameters{{"VerbosityLevel", verbosity}});
	// solverParameters.insert(mfem_mgis::Parameters{{"MaximumNumberOfIterations", defaultMaxNumOfIt}});
	// // solverParameters.insert(mfem_mgis::Parameters{{"AbsoluteTolerance", Tol}});
	// // solverParameters.insert(mfem_mgis::Parameters{{"RelativeTolerance", Tol}});
	// solverParameters.insert(mfem_mgis::Parameters{{"Tolerance", Tol}});


	// // preconditionner hypreBoomerAMG
	// auto options = mfem_mgis::Parameters{{"VerbosityLevel", verbosity}};
	// auto preconditionner = mfem_mgis::Parameters{{"Name","HypreDiagScale"}, {"Options",options}};
	// // auto preconditionner = mfem_mgis::Parameters{{"Name","HypreBoomerAMG"}, {"Options",options}};
	// solverParameters.insert(mfem_mgis::Parameters{{"Preconditioner",preconditionner}});
	// solver HyprePCG
	// p.setLinearSolver("HyprePCG", solverParameters);
	p.setLinearSolver("MUMPSSolver", solverParameters);
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

template <typename Problem>
void simulation(Problem &problem, double start, double end, double dt, bool pp)
{
	CatchTimeSection("simulation");

	const int nMat = getMaterialsAttributes(*(problem.getFiniteElementDiscretizationPointer())).Max();
	const double young1 = 10.e9;
	const double young2 = 10.e9;
	const double young3 = 230.e9;
	const double poisson12 = 0.3;
	const double poisson23 = 0.2;
	const double poisson13 = 0.2;
	const double shear12 = 3.e9;
	const double shear23 = 7.e9;
	const double shear13 = 7.e9;

	// materials
	auto set_properties = [](auto& m, 
			const double yo1 , const double yo2,	const double yo3,
			const double po12, const double po23, const double po13,
			const double sm12, const double sm23, const double sm13
			) {
		setMaterialProperty(m.s0, "YoungModulus1", yo1);
		setMaterialProperty(m.s0, "YoungModulus2", yo2);
		setMaterialProperty(m.s0, "YoungModulus3", yo3);
		setMaterialProperty(m.s0, "PoissonRatio12", po12);
		setMaterialProperty(m.s0, "PoissonRatio23", po23);
		setMaterialProperty(m.s0, "PoissonRatio13", po13);
		setMaterialProperty(m.s0, "ShearModulus12", sm12);
		setMaterialProperty(m.s0, "ShearModulus23", sm23);
		setMaterialProperty(m.s0, "ShearModulus13", sm13);

		setMaterialProperty(m.s1, "YoungModulus1", yo1);
		setMaterialProperty(m.s1, "YoungModulus2", yo2);
		setMaterialProperty(m.s1, "YoungModulus3", yo3);
		setMaterialProperty(m.s1, "PoissonRatio12", po12);
		setMaterialProperty(m.s1, "PoissonRatio23", po23);
		setMaterialProperty(m.s1, "PoissonRatio13", po13);
		setMaterialProperty(m.s1, "ShearModulus12", sm12);
		setMaterialProperty(m.s1, "ShearModulus23", sm23);
		setMaterialProperty(m.s1, "ShearModulus13", sm13);

	};

	auto set_temperature = [](auto &m)
	{
		setExternalStateVariable(m.s0, "Temperature", 293.15);
		setExternalStateVariable(m.s1, "Temperature", 293.15);
	};

	for (int i = 0; i < nMat; i++)
	{
		auto &mat = problem.getMaterial(i + 1);
		set_properties(mat,
					   young1, young2, young3,	  // young modulus
					   poisson12, poisson23, poisson13, // poisson ration
					   shear12, shear23, shear13		  // shear modulus
		);
		set_temperature(mat);
	}
	const int nStep = (end - start) / dt;
	for (int i = 0; i < nStep; i++)
	{
		// Transformation gradient components
		const double EZZ = 0.02 * ((i + 1.) * dt / end);
		double EXX = -((poisson12+1)*poisson13*young3)/(young1-poisson13*poisson13*young3)*EZZ;
		// double EYY = (poisson12*poisson12+poisson13)*young3/(young1-poisson12*poisson12*young3)*EZZ;

		// init macro cauchy stress components
		double SIGXX = 0.;

		// init volume var
		double vol;

		// init MPI vars
		double SIG_local[1] = {SIGXX};
		double SIG_global[1] = {0.};

		// Fixed-Point param
		double tolFP = 1e-6;
		double oldRes = 1.;
		double newRes = std::sqrt(SIGXX * SIGXX);
		int itFP = 0;
		int maxitFP = 100;
		while ((abs(oldRes - newRes) > tolFP) && (itFP <= maxitFP))
		{
			mfem_mgis::Profiler::Utils::Message("SIGXX = ", SIGXX);
			problem.setMacroscopicGradientsEvolution([&EXX, EZZ, &SIGXX, young1, young3, poisson13, poisson12](const double t) { 
				auto ret = std::vector<mfem_mgis::real>(6, mfem_mgis::real{0});
				ret[2] = EZZ;
				ret[0] = EXX + ((1+poisson12)*(2*poisson13*poisson13*young3+young1*(poisson12-1)))/(young1*(young1-poisson13*poisson13*young3))*SIGXX ;
				EXX = ret[0];
				mfem_mgis::Profiler::Utils::Message("EXX", EXX);
				mfem_mgis::Profiler::Utils::Message("EYY", ret[1]);
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
			for (const auto &mi : problem.getAssignedMaterialsIdentifiers())
			{
				// volumes[mi] ;
				mfem_mgis::Profiler::Utils::Message("tf_integrals[mi] = ", tf_integrals[mi][0]);
				SIGXX += tf_integrals[mi][0] / vol;
			}

			SIG_local[0] = SIGXX;
			MPI_Allreduce(SIG_local, SIG_global, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
			SIGXX = SIG_global[0];

			oldRes = newRes;
			newRes = std::sqrt(SIGXX * SIGXX);
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

	// set problem
	setup_properties(p, problem);
	setLinearSolver(problem, p.verbosity_level);

	problem.setSolverParameters({{"VerbosityLevel", p.verbosity_level},
			{"RelativeTolerance", 1.e-4},
			{"AbsoluteTolerance", 1.e-4},
			{"MaximumNumberOfIterations", 15}});


	// add post processings
	if(use_post_processing) add_post_processings(problem, "OutputFile-Uniaxial-elastic-ortho");

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
