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
	const char* vect_file = "periodic-cube-vecteurs.txt";
	const char* behaviour = "MonoCristal_UO2";
	const char* library = "src/libBehaviour.so";
	int order = 2;
	bool parallel = true;
	int refinement = 0;
	int post_processing = 1; // default value : disabled
	int verbosity_level = 0; // default value : lower level
};

void common_parameters(mfem::OptionsParser& args, TestParameters& p)
{
	args.AddOption(&p.mesh_file, "-m", "--mesh", "Mesh file to use.");
	args.AddOption(&p.vect_file, "-f", "--vect", "Vector file to use.");
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

template<typename Implementation>
void print_mesh_information(Implementation& impl)
{
	CatchTimeSection("common::print_mesh_information");

	using mfem_mgis::Profiler::Utils::sum;
	using mfem_mgis::Profiler::Utils::Message;

	//getMesh
	auto mesh = impl.getFiniteElementSpace().GetMesh();

	//get the number of vertices
	int64_t numbers_of_vertices_local = mesh->GetNV();
	int64_t  numbers_of_vertices = sum(numbers_of_vertices_local);

	//get the number of elements
	int64_t numbers_of_elements_local = mesh->GetNE();
	int64_t numbers_of_elements = sum(numbers_of_elements_local);

	//get the element size
	double h = mesh->GetElementSize(0);

	// get n dofs
	auto& fespace = impl.getFiniteElementSpace();
	int64_t unknowns_local = fespace.GetTrueVSize(); 
	int64_t unknowns = sum(unknowns_local);

	Message("Info problem: number of vertices -> ", numbers_of_vertices);
	Message("Info problem: number of elements -> ", numbers_of_elements);
	Message("Info prolbem: element size -> ", h);
	Message("Info porblem: Number of finite element unknowns: " , unknowns);
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
			{{"OutputFileName", "avgStressPolycristal"}});
	// p.addPostProcessing(
	// 		"ParaviewExportIntegrationPointResultsAtNodes",
	// 		{{"OutputFileName", msg + "IntegrationPointOutputPKI"},
	// 		 {"Results", {"FirstPiolaKirchhoffStress"}}});
	// p.addPostProcessing(
	// 		"ParaviewExportIntegrationPointResultsAtNodes",
	// 		{{"OutputFileName", msg + "IntegrationPointOutputDG"},
	// 		 {"Results", {"DeformationGradient"}}});
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

	// const int nMat = 8;
	const int nMat = getMaterialsAttributes(*(problem.getFiniteElementDiscretizationPointer())).Max();
	mfem_mgis::Profiler::Utils::Message("Nombre de matériaux : ", nMat);

	for(int i = 0 ; i < nMat ; i++)
	{
		problem.addBehaviourIntegrator("Mechanics", i+1, p.library, p.behaviour);
	}
	
	std::vector<std::array<mfem_mgis::real, 3u>> vectors = readVectorsFromFile(p.vect_file);
	if (vectors.size()!=2*nMat) {
		throw std::invalid_argument("setup_properties : incorrect number of vectors in vector file");
	}

	std::array<mfem_mgis::MaterialAxis3D, 2u> r;
	for (int i = 0; i < nMat; i++)
	{
		auto &mat = problem.getMaterial(i + 1);
		if (mat.b.symmetry == mgis::behaviour::Behaviour::ORTHOTROPIC)
		{
			r[0] = vectors[2 * i];
			r[1] = vectors[2 * i + 1];
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
	solverParameters.insert(mfem_mgis::Parameters{{"VerbosityLevel", verbosity}});
	solverParameters.insert(mfem_mgis::Parameters{{"MaximumNumberOfIterations", defaultMaxNumOfIt}});
	// solverParameters.insert(mfem_mgis::Parameters{{"AbsoluteTolerance", Tol}});
	// solverParameters.insert(mfem_mgis::Parameters{{"RelativeTolerance", Tol}});
	solverParameters.insert(mfem_mgis::Parameters{{"Tolerance", Tol}});


	// preconditionner hypreBoomerAMG
	auto options = mfem_mgis::Parameters{{"VerbosityLevel", verbosity}};
	// auto preconditionner = mfem_mgis::Parameters{{"Name","HypreDiagScale"}, {"Options",options}};
	auto preconditionner = mfem_mgis::Parameters{{"Name","HypreBoomerAMG"}, {"Options",options}};
	solverParameters.insert(mfem_mgis::Parameters{{"Preconditioner",preconditionner}});
	// solver HyprePCG
	p.setLinearSolver("HyprePCG", solverParameters);
	// p.setLinearSolver("MUMPSSolver", solverParameters);
	// p.setLinearSolver("CGSolver", solverParameters);
}

	template<typename Problem>
void run_solve(Problem& p, double start, double dt)
{
	CatchTimeSection("Solve");
	// solving the problem
	auto statistics = p.solve(start, dt);
	// check status
	if (statistics.status == false) {
		mfem_mgis::Profiler::Utils::Message("INFO: SOLVE FAILED");
		std::abort();
	}
}

template <typename Problem>
void simulation(Problem &problem, double start, double end, double dt, bool pp)
{
	CatchTimeSection("simulation");

	const int nMat = getMaterialsAttributes(*(problem.getFiniteElementDiscretizationPointer())).Max();
	// cubic symmetry elasticity
	const double young1 = 222.e9;
	const double young2 = young1;
	const double young3 = young1;
	const double poisson12 = 0.27;
	const double poisson23 = poisson12;
	const double poisson13 = poisson12;
	const double shear12 = 54.e9;
	const double shear23 = shear12;
	const double shear13 = shear12;

	// ===========================
	// effective moduli estimations
	// ===========================
	const double C11 = young1 * (poisson12 - 1) / (2 * poisson12 * poisson12 + poisson12 - 1);
	const double C12 = -(poisson12 * young1) / (2 * poisson12 * poisson12 + poisson12 - 1);
	const double C44 = shear12;
	// ---- Bulk Modulus ----
	const double Keff = 1. / 3. * (C11 + 2 * C12);
	// ---- Shear Modulus ----
	// Reuss & Voigt bounds
	const double aniso = 2 * C44 / (C11 - C12);
	const double GReuss = 5 * C44 / (3 + 2 * aniso);
	const double GVoigt = (3 * aniso + 2) / (5 * aniso) * C44;
	// Hashin–Shtrikman bounds
	const double Gv = (C11 - C12) / 2;
	const double beta1 = -3 * (Keff + 2 * Gv) / (5 * Gv * (3 * Keff + 4 * Gv));
	const double beta2 = -3 * (Keff + 2 * C44) / (5 * C44 * (3 * Keff + 4 * C44));
	const double GHSl = Gv + 3 * ((C44 - Gv) - 4 * beta1) / 5;
	const double GHSu = C44 + 3 * ((Gv - C44) - 6 * beta2) / 5;
	// Kroener
	double m0, mus;
	double D1, sD1, D2;
	double GK = 0.;
	do
	{
		m0 = GK;
		D1 = Keff + 2 * GK;
		sD1 = 6 * D1;
		mus = GK * (9 * Keff + 8 * GK) / sD1;
		D2 = 5 * mus + 2 * C44 + 3 * Gv;
		GK = (3 * C44 * mus + 2 * Gv * mus + 5 * Gv * C44) / D2;
	} while (fabs(GK - m0) / Keff > 1e-10);

	// Shear modulus choice
	mfem_mgis::Profiler::Utils::Message("GReuss =", GReuss);
	mfem_mgis::Profiler::Utils::Message("GVoigt =", GVoigt);
	mfem_mgis::Profiler::Utils::Message("GHSl =", GHSl);
	mfem_mgis::Profiler::Utils::Message("GHSu =", GHSu);
	mfem_mgis::Profiler::Utils::Message("GK =", GK);

	const double Geff = GK; // GHSu
	// Young and Poisson effective moduli
	const double nueff = (3 * Keff - 2 * Geff) / (2 * (3 * Keff + Geff));
	const double Eeff = 9 * Keff * Geff / (3 * Keff + Geff);
	// ===========================

	// materials
	auto set_properties = [](auto &m,
							 const double yo1, const double yo2, const double yo3,
							 const double po12, const double po23, const double po13,
							 const double sm12, const double sm23, const double sm13)
	{
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
		setExternalStateVariable(m.s0, "Temperature", 1600.);
		setExternalStateVariable(m.s1, "Temperature", 1600.);
	};

	for (int i = 0; i < nMat; i++)
	{
		auto &mat = problem.getMaterial(i + 1);
		set_properties(mat,
					   young1, young2, young3,			// young modulus
					   poisson12, poisson23, poisson13, // poisson ration
					   shear12, shear23, shear13		// shear modulus
		);
		set_temperature(mat);
	}

	double FXX, FYY, FZZ;
	// init macro cauchy stress components
	double SIGXX = 0.;
	double SIGYY = 0.;
	problem.setMacroscopicGradientsEvolution([&FXX, &FYY, &FZZ, &SIGXX, &SIGYY, Eeff, nueff](const double t)
													 { 
				auto ret = std::vector<mfem_mgis::real>(9, mfem_mgis::real{});
				ret[0] = FXX;
				ret[1] = FYY;
				ret[2] = FZZ;
				return ret; });

	const int nStep = (end - start) / dt;
	for (int i = 0; i < nStep; i++)
	{
		mfem_mgis::Profiler::Utils::Message("\nSolving problem from ",i*dt,"to",(i+1)*dt);
		// Transformation gradient components
		FZZ = 1. + 0.1 * ((i + 1.) * dt / end);
		if (i == 0) {
			FXX = 1. / std::sqrt(FZZ);
			FYY = 1. / std::sqrt(FZZ);
		}

		// init volume var
		double vol;

		// init MPI vars
		double SIG_local[2] = {SIGXX, SIGYY};
		double SIG_global[2] = {0., 0.};

		// Fixed-Point param
		double tolFP = 1.e4;
		double oldRes = 1.e6;
		double newRes = std::sqrt(SIGXX * SIGXX + SIGYY * SIGYY);
		int itFP = 0;
		int maxitFP = 50;
		while ((abs(oldRes - newRes) > tolFP) && (itFP <= maxitFP))
		{
			//compute correction
			double Fcorr[2];
			Fcorr[0] = ((nueff*nueff-1)*SIGXX)/Eeff+(nueff*(nueff+1)*SIGYY)/Eeff;
			Fcorr[1] = ((nueff*nueff-1)*SIGYY)/Eeff+(nueff*(nueff+1)*SIGXX)/Eeff;

			// update load 
			FXX += Fcorr[0];
			FYY += Fcorr[1];

			// run solve
			run_solve(problem, i * dt, dt);
			problem.update();

			// post-processing
			// get MeanThermodynamicForcesValues and Volumes
			auto [tf_integrals, volumes] = mfem_mgis::computeMeanThermodynamicForcesValues<true>(problem.template getImplementation<true>());

			// compute volume on all Materials Identifiers
			vol = 0.;
			for (const auto &mi : problem.getAssignedMaterialsIdentifiers())
			{
				vol += volumes[mi];
			}
			MPI_Allreduce(MPI_IN_PLACE, &vol, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

			// compute mean value of Thermodynamic Forces on all Materials Identifiers
			SIGXX = 0.;
			SIGYY = 0.;
			for (const auto &mi : problem.getAssignedMaterialsIdentifiers())
			{
				SIGXX += tf_integrals[mi][0] / vol;
				SIGYY += tf_integrals[mi][1] / vol;
			}

			// PKI to VM
			SIGXX = SIGXX*FXX/(FXX*FYY*FZZ);
			SIGYY = SIGYY*FYY/(FXX*FYY*FZZ);

			SIG_local[0] = SIGXX;
			SIG_local[1] = SIGYY;
			MPI_Allreduce(SIG_local, SIG_global, 2, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
			SIGXX = SIG_global[0];
			SIGYY = SIG_global[1];

			oldRes = newRes;
			newRes = std::sqrt(SIGXX * SIGXX + SIGYY * SIGYY);

			mfem_mgis::Profiler::Utils::Message("Fixed Point iteration", itFP, ": |res| =", abs(oldRes - newRes));
			// mfem_mgis::Profiler::Utils::Message("SIGXX =", SIGXX);
			// mfem_mgis::Profiler::Utils::Message("SIGYY =", SIGYY);

			if (abs(oldRes - newRes) > tolFP) {
				problem.revert();
			}

			itFP += 1;
		}

		if (itFP >= maxitFP)
			mfem_mgis::Profiler::Utils::Message("warning: maximum number of iterations for the fixed-point algorithm attained, before the requested tolerance is reached");

		if (pp)
			execute_post_processings(problem, i * dt, dt);
		
		if (mfem_mgis::getMPIrank() == 0) {
			std::ofstream fichier("OutputFile-Uniaxial-polycristal-GradTrans.txt", std::ios::app);
			if (fichier.is_open()) {
				fichier << (i+1) * dt << " " << FXX << " " << FYY << " " << FZZ << "\n";
				fichier.close();
				std::cout << "Les valeurs ont été écrites dans le fichier." << std::endl;
			} else {
				std::cerr << "Erreur lors de l'ouverture du fichier." << std::endl;
			}
		}
		
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
			{"RelativeTolerance", 1e-4},
			{"AbsoluteTolerance", 1e-4},
			{"MaximumNumberOfIterations", 15}});


	// add post processings
	if(use_post_processing) add_post_processings(problem, "OutputFile-Uniaxial-polycristal");

	// main function here
	int nStep=600;
	double start=0;
	double end=200;
	const double dt = (end-start)/nStep;
	simulation(problem, start, end, dt, use_post_processing);

	// print and write timetable
	mfem_mgis::Profiler::timers::print_and_write_timers();
	return(EXIT_SUCCESS);
}
