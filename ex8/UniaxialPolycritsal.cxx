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
	const char* mesh_file = "cube.msh";
	const char* behaviour = "MonoCristal_UO2";
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
			{{"OutputFileName", "avgStress"}});
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

	const int nMat = 8;

	for(int i = 0 ; i < nMat ; i++)
	{
		problem.addBehaviourIntegrator("Mechanics", i+1, p.library, p.behaviour);
	}
	// materials

//	auto& m1 = problem.getMaterial(1);
	auto set_properties = [](auto& m, 
			const double yo1,
			const double yo2,
			const double yo3,
			const double po12,
			const double po23,
			const double po13,
			const double sm12,
			const double sm23,
			const double sm13,
			const double v1x,
			const double v1y,
			const double v1z,
			const double v2x,
			const double v2y,
			const double v2z
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

		// dunno why I have to comment these lines
/*
		setMaterialProperty(m.s0, "V1X", v1x);
		setMaterialProperty(m.s0, "V1Y", v1y);
		setMaterialProperty(m.s0, "V1Z", v1z);
		setMaterialProperty(m.s0, "V2X", v2x);
		setMaterialProperty(m.s0, "V2Y", v2y);
		setMaterialProperty(m.s0, "V2Z", v2z);
*/

		setMaterialProperty(m.s1, "YoungModulus1", yo1);
		setMaterialProperty(m.s1, "YoungModulus2", yo2);
		setMaterialProperty(m.s1, "YoungModulus3", yo3);
		setMaterialProperty(m.s1, "PoissonRatio12", po12);
		setMaterialProperty(m.s1, "PoissonRatio23", po23);
		setMaterialProperty(m.s1, "PoissonRatio13", po13);
		setMaterialProperty(m.s1, "ShearModulus12", sm12);
		setMaterialProperty(m.s1, "ShearModulus23", sm23);
		setMaterialProperty(m.s1, "ShearModulus13", sm13);

		/*
			 setMaterialProperty(m.s1, "V1X", v1x);
			 setMaterialProperty(m.s1, "V1Y", v1y);
			 setMaterialProperty(m.s1, "V1Z", v1z);
			 setMaterialProperty(m.s1, "V2X", v2x);
			 setMaterialProperty(m.s1, "V2Y", v2y);
			 setMaterialProperty(m.s1, "V2Z", v2z);
		 */
	};

	auto set_temperature = [](auto& m) {
		setExternalStateVariable(m.s0, "Temperature", 1600.);
		setExternalStateVariable(m.s1, "Temperature", 1600.);
	};

	for(int i = 0 ; i < nMat ; i++)
	{
		auto& mat = problem.getMaterial(i+1);
		set_properties(mat,     
				222.e9,              222.e9,              (222.e9, // youn modulus
				0.27,                0.27,                0.27, // poission ration
				54.e9,               54.e9,               54.e9, // shear modulus
				0.7071067811865475, -0.4086070447619255, -0.5770964243269279, // v1
				0.7071067811865475,  0.4086070447619256,  0.5770964243269281  // v2
				);
		set_temperature(mat);
	}


	// macroscopic strain
	std::vector<real> e(9, real{0});
	const int xx = 1;
	const int yy = 2;
	const int zz = 3;


	/* bar{E} = e33 *(-1/2 E1 x E1 + (-1/2) * E2 x E2 + E3 x E3)*/
	const double coef = 0.1;
	e[xx] = 1-0.5*coef;
	e[yy] = 1-0.5*coef;
	e[zz] = 1+1*coef;

	problem.setMacroscopicGradientsEvolution([e](const double dt) { 
			auto ret = e;
			for(auto& it : ret) it *= dt;
			return ret; 
			});

	// no idea of what is dis but I need to declare it
	//	if (m1.b.symmetry == mgis::behaviour::Behaviour::ISOTROPIC) {
	for(int i = 0 ; i < nMat ; i++)
	{
		auto& mat = problem.getMaterial(i+1);
		std::array<mfem_mgis::real, 9u> r = {1, 0, 0,  //
			0, 1, 0,  //
			0, 0, 1};
		mat.setRotationMatrix(mfem_mgis::RotationMatrix3D{r});
	}
	//	}
} 


	template<typename Problem>		
static void setLinearSolver(Problem& p,
		const int verbosity = 0,
		const mfem_mgis::real Tol = 1e-6
		)
{
	CatchTimeSection("set_linear_solver");
	// pilote
	constexpr int defaultMaxNumOfIt	 	= 500; 		// MaximumNumberOfIterations
	auto solverParameters = mfem_mgis::Parameters{};
	solverParameters.insert(mfem_mgis::Parameters{{"VerbosityLevel", verbosity}});
	solverParameters.insert(mfem_mgis::Parameters{{"MaximumNumberOfIterations", defaultMaxNumOfIt}});
	solverParameters.insert(mfem_mgis::Parameters{{"AbsoluteTolerance", Tol}});
	solverParameters.insert(mfem_mgis::Parameters{{"RelativeTolerance", Tol}});
	//solverParameters.insert(mfem_mgis::Parameters{{"Tolerance", Tol}});


	// preconditionner hypreBoomerAMG
	auto options = mfem_mgis::Parameters{{"VerbosityLevel", verbosity}};
	auto preconditionner = mfem_mgis::Parameters{{"Name","HypreBoomerAMG"}, {"Options",options}};
	solverParameters.insert(mfem_mgis::Parameters{{"Preconditioner",preconditionner}});
	// solver HyprePCG
	//p.setLinearSolver("HyprePCG", solverParameters);
	p.setLinearSolver("CGSolver", solverParameters);
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

	//	std::vector<mfem_mgis::real> corner1({0.,0.,0.});
	//	std::vector<mfem_mgis::real> corner2({1., 1., 1.});
	//	mfem_mgis::PeriodicNonLinearEvolutionProblem problem(fed, corner1, corner2);

	// set problem
	setup_properties(p, problem);
	setLinearSolver(problem, p.verbosity_level);

	problem.setSolverParameters({{"VerbosityLevel", p.verbosity_level},
			{"RelativeTolerance", 1e-6},
			{"AbsoluteTolerance", 0.},
			{"MaximumNumberOfIterations", 6}});


	// add post processings
	if(use_post_processing) add_post_processings(problem, "OutputFile-Uniaxial-polycristal");

	// main function here
	int nStep=1000;
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
