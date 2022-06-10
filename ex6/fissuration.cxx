#include "test_parameters.hpp"
#include "solver_name.hxx"
#include "timer.hpp"
#include "precond_name.hxx"
#include "config_solver.hxx"
#include "data_gathering.hxx"
#include "MFEMMGIS/Config.hxx"
#include "MFEMMGIS/Parameters.hxx"
#include "MFEMMGIS/Material.hxx"
#include "MFEMMGIS/PartialQuadratureSpace.hxx"
#include "MFEMMGIS/ImposedDirichletBoundaryConditionAtClosestNode.hxx"
#include "MFEMMGIS/NonLinearEvolutionProblem.hxx"

namespace fissuration
{
	using namespace configuration;
	
	using solver=solver_name;
	using pc=precond_name;

	using petsc_solver=petsc_ksp_type;
	using petsc_pc=petsc_pc_type;

	using regular_solver_list=std::vector<solver>;
	using petsc_solver_list=std::vector<petsc_solver>;

	using regular_pc_list=std::vector<pc>;
	using petsc_pc_list=std::vector<petsc_pc>;

	/*
	 * brief return a list of solvers that work or do not run slowly for the test case 1
	 */
	regular_solver_list build_solvers_list()
	{
		regular_solver_list res = {
			solver_name::HypreFGMRES,
			solver_name::HypreGMRES,
			solver_name::HyprePCG,
			solver_name::MUMPSSolver,
			solver_name::GMRESSolver,
			solver_name::CGSolver,
			solver_name::BiCGSTABSolver,
			solver_name::UMFPackSolver,
			solver_name::MINRESSolver,
			solver_name::SLISolver,
		};
		return res;
	}

	/*
	 * brief return a list of preconditioners that work or do not run slowly for the test case 1
	 */
	regular_pc_list build_pc_list()
	{

		regular_pc_list res = {
			precond_name::HypreParaSails,
			precond_name::HypreBoomerAMG,
			precond_name::HypreILU,
			precond_name::HypreEuclid,
			precond_name::HypreDiagScale
		};
		return res;	
	}
	
	/*
	 * brief return a list of petsc solvers that work or do not run slowly for the test case 1
	 */
	petsc_solver_list build_solvers_list_with_petsc()
	{
		petsc_solver_list res = {
			petsc_solver::richardson,	petsc_solver::chebyshev,	petsc_solver::cg,
			petsc_solver::gmres,		petsc_solver::tcqmr,		petsc_solver::bcgs,
			petsc_solver::cgs,		petsc_solver::tfqmr,		petsc_solver::cr,
			petsc_solver::lsqr,		petsc_solver::bicg,		petsc_solver::preonly
		};
		return res;

	}

	petsc_pc_list build_pc_list_with_petsc()
	{
		petsc_pc_list res = {
			petsc_pc::jacobi,	petsc_pc::bjacobi,	petsc_pc::sor,
			petsc_pc::eisenstat,	petsc_pc::icc,		petsc_pc::ilu,
			petsc_pc::gasm,		petsc_pc::gamg,		petsc_pc::bddc,
			petsc_pc::ksp,		petsc_pc::lu,		petsc_pc::cholesky,
			petsc_pc::none
		};
		return res;
	}

	bool match(solver s, pc p)
	{
		return true;
	}

	bool match_with_petsc(petsc_solver s, petsc_pc p)
	{
		return true;
	}

	std::shared_ptr<mfem_mgis::NonLinearEvolutionProblem> buildMechanicalProblem(const mfem_mgis::Parameters& common_problem_parameters,
			const MeshParameters mesh_parameters) 
	{

		START_TIMER("fissuration::buildMechanicalProblem");
		auto lparameters = common_problem_parameters;
		lparameters.insert("UnknownsSize", 3);
		// building the non linear problem
		auto problem =
			std::make_shared<mfem_mgis::NonLinearEvolutionProblem>(lparameters);
		// materials
		problem->addBehaviourIntegrator("Mechanics", "Fuel", "src/libBehaviour.so",
				//                              "MicromorphicDamageI_DeviatoricSplit");
				//                              "MicromorphicDamageI_SpectralSplit");
			"MicromorphicDamageI");
		auto& m1 = problem->getMaterial("Fuel");
		// material properties at the beginning and the end of the time step
		for (auto& s : {&m1.s0, &m1.s1}) {
			mgis::behaviour::setMaterialProperty(*s, "YoungModulus", 150e9);
			mgis::behaviour::setMaterialProperty(*s, "PoissonRatio", 0.3);
		}
		// temperature
		mgis::behaviour::setExternalStateVariable(m1.s0, "Temperature", 293.15);
		mgis::behaviour::setExternalStateVariable(m1.s0, "Damage", 0);
		// boundary conditions
		std::array<mfem_mgis::real, 3u> p1{0, 0, mesh_parameters.dh};
		std::array<mfem_mgis::real, 3u> p2{0, 0, mesh_parameters.hp - mesh_parameters.dh};
		std::array<mfem_mgis::real, 3u> p3{mesh_parameters.rp, 0, 0};
		problem->addBoundaryCondition(
				std::make_unique<
				mfem_mgis::ImposedDirichletBoundaryConditionAtClosestNode>(
					problem->getFiniteElementDiscretizationPointer(), p1, 0));
		problem->addBoundaryCondition(
				std::make_unique<
				mfem_mgis::ImposedDirichletBoundaryConditionAtClosestNode>(
					problem->getFiniteElementDiscretizationPointer(), p1, 1));
		problem->addBoundaryCondition(
				std::make_unique<
				mfem_mgis::ImposedDirichletBoundaryConditionAtClosestNode>(
					problem->getFiniteElementDiscretizationPointer(), p1, 2));
		problem->addBoundaryCondition(
				std::make_unique<
				mfem_mgis::ImposedDirichletBoundaryConditionAtClosestNode>(
					problem->getFiniteElementDiscretizationPointer(), p2, 0));
		problem->addBoundaryCondition(
				std::make_unique<
				mfem_mgis::ImposedDirichletBoundaryConditionAtClosestNode>(
					problem->getFiniteElementDiscretizationPointer(), p2, 1));
		problem->addBoundaryCondition(
				std::make_unique<
				mfem_mgis::ImposedDirichletBoundaryConditionAtClosestNode>(
					problem->getFiniteElementDiscretizationPointer(), p3, 1));
		return problem;
	}


	std::shared_ptr<mfem_mgis::NonLinearEvolutionProblem> buildMicromorphicProblem(
			const mfem_mgis::Parameters& common_problem_parameters) {

		START_TIMER("fissuration::buildMicromorphicProblem");
		//  1/(2*E)*smax**2*lc = Gc
		constexpr auto Gc = mfem_mgis::real{5};
		constexpr auto DGc = mfem_mgis::real{0.5};
		//  constexpr auto l0 = mfem_mgis::real{60e-5};
		constexpr auto l0 = mfem_mgis::real{30e-5};
		constexpr auto beta = mfem_mgis::real{300};
		auto lparameters = common_problem_parameters;
		lparameters.insert({{"UnknownsSize", 1}});
		auto problem =
			std::make_shared<mfem_mgis::NonLinearEvolutionProblem>(lparameters);
		problem->addBehaviourIntegrator("MicromorphicDamage", "Fuel",
				"src/libBehaviour.so",
				"AT1MicromorphicDamage");
		auto& m = problem->getMaterial("Fuel");
		// material properties
		for (const auto& mp :
				std::map<std::string, double>{{"CharacteristicLength", l0},
				{"PenalisationFactor", beta}}) {
								       mgis::behaviour::setMaterialProperty(m.s0, mp.first, mp.second);
								       mgis::behaviour::setMaterialProperty(m.s1, mp.first, mp.second);
							       }
		std::default_random_engine Gc_generator;
		std::normal_distribution<double> Gc_distribution(Gc,DGc);
		auto Gc_mp = mfem_mgis::PartialQuadratureFunction::evaluate(
				problem->getMaterial("Fuel").getPartialQuadratureSpacePointer(),
				[&Gc_generator, &Gc_distribution](const mfem_mgis::real,
					const mfem_mgis::real,
					const mfem_mgis::real) {

				return std::max(Gc_distribution(Gc_generator), Gc - 3 * DGc);
				});
		mgis::behaviour::setMaterialProperty(m.s0, "FractureEnergy", Gc_mp->getValues());
		mgis::behaviour::setMaterialProperty(m.s1, "FractureEnergy", Gc_mp->getValues());
		mgis::behaviour::setExternalStateVariable(m.s0, "Temperature", 293.15);
		mgis::behaviour::setExternalStateVariable(m.s1, "Temperature", 293.15);
		mgis::behaviour::setExternalStateVariable(m.s0, "EnergyReleaseRate", 0);
		return problem;
	}
}
