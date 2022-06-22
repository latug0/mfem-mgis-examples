#pragma once
// thomas example
#include <common/timer.hpp>
#include <solver/solver_name.hxx>
#include <solver/precond_name.hxx>
#include <solver/config_solver.hxx>
#include <parameters/test_parameters.hpp>

namespace fissuration
{
	using namespace configuration;
	using solver=solver_name;
	using pc=precond_name;
	using regular_solver_list=std::vector<solver>;
	using regular_pc_list=std::vector<pc>;
	regular_solver_list 	build_solvers_list();
	regular_pc_list 	build_pc_list();
	bool match(solver, pc, int);

	using petsc_solver=petsc_ksp_type;
	using petsc_pc=petsc_pc_type;
	using petsc_solver_list=std::vector<petsc_solver>;
	using petsc_pc_list=std::vector<petsc_pc>;

	petsc_solver_list 	build_solvers_list_with_petsc();
	petsc_pc_list 		build_pc_list_with_petsc();
	bool match_with_petsc(petsc_solver, petsc_pc, int r);
	
	// specific functions
	struct MeshParameters {
		// fuel pellet radius
		const mfem_mgis::real rp = mfem_mgis::real{4.085e-3};
		// fuel pellet height
		const mfem_mgis::real hp = mfem_mgis::real{13.5e-3};
		// dishing height
		const mfem_mgis::real dh = mfem_mgis::real{3.2e-4};
	};

	static std::shared_ptr<mfem_mgis::NonLinearEvolutionProblem> buildMechanicalProblem(const mfem_mgis::Parameters& common_problem_parameters,
			const MeshParameters mesh_parameters) 
	{

		START_TIMER("fissuration::buildMechanicalProblem");
		auto lparameters = common_problem_parameters;
		lparameters.insert("UnknownsSize", 3);
		// building the non linear problem
		auto problem =
			std::make_shared<mfem_mgis::NonLinearEvolutionProblem>(lparameters);
	
		
		// materials
		problem->addBehaviourIntegrator("Mechanics", 
				"Fuel", 
				"src/libBehaviour.so",
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
		common::print_mesh_information(problem->getImplementation<true>());
		return problem;
	}

	static std::shared_ptr<mfem_mgis::NonLinearEvolutionProblem> buildMicromorphicProblem(
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
		common::print_mesh_information(problem->getImplementation<true>());
		return problem;
	}


	template<typename Solver, typename Pc>
		int kernel(const TestParameters& p, const bool use_post_processing, const solver_name a_solv, const precond_name a_precond, gather_information& a_info)
		{
			constexpr auto iter_max = mfem_mgis::size_type{5000};
			static constexpr const auto parallel = true;
			auto success = true;
			// building the non linear problems
			const auto common_problem_parameters = mfem_mgis::Parameters{
				{"MeshFileName", "Rod3D_Maillage_2.50_.msh"},
					{"Materials", mfem_mgis::Parameters{{"Fuel", 4}}},
					{"FiniteElementFamily", "H1"},
					{"FiniteElementOrder", 1},
					{"NumberOfUniformRefinements", parallel ? p.refinement : 0},
					{"Hypothesis", "Tridimensional"},
					{"Parallel", parallel}};
			const auto mesh_parameters = MeshParameters{};
			// mechanical problem
			auto mechanical_problem = buildMechanicalProblem(common_problem_parameters,
					mesh_parameters);
			// micromorphic problem
			auto micromorphic_problem =
				buildMicromorphicProblem(common_problem_parameters);
			// quadrature functions used to transfer information from one problem to the
			// other
			mfem_mgis::PartialQuadratureFunction Y(
					micromorphic_problem->getMaterial("Fuel")
					.getPartialQuadratureSpacePointer(),
					1u);
			mfem_mgis::PartialQuadratureFunction d(
					micromorphic_problem->getMaterial("Fuel")
					.getPartialQuadratureSpacePointer(),
					1u);
			// using external storage allows to directly modify the values of the
			// quadrature functions Y and d
			mgis::behaviour::setExternalStateVariable(
					mechanical_problem->getMaterial("Fuel").s1, "Damage", d.getValues(),
					mgis::behaviour::MaterialStateManager::EXTERNAL_STORAGE);
			mgis::behaviour::setExternalStateVariable(
					micromorphic_problem->getMaterial("Fuel").s1, "EnergyReleaseRate",
					Y.getValues(), mgis::behaviour::MaterialStateManager::EXTERNAL_STORAGE);
			// solving the problem on 1 time step
			constexpr auto n = mfem_mgis::size_type{50};
			constexpr auto tb = mfem_mgis::real{0};
			constexpr auto te = mfem_mgis::real{1};
			auto times = [] {
				const auto dt = (te - tb) / static_cast<mfem_mgis::real>(n);
				auto r = std::vector<mfem_mgis::real>{};
				auto t = tb;
				r.push_back(t);
				for (mfem_mgis::size_type i = 0; i != n; ++i) {
					t += dt;
					r.push_back(t);
				}
				return r;
			}();
			for (mfem_mgis::size_type i = 0; i != n; ++i) 
			{
				START_TIMER("fissuration::main_loop");
				const auto t0 = times[i];
				const auto t1 = times[i + 1];
				const auto dt = t1 - t0;
				if(mfem_mgis::getMPIrank() == 0) {
					std::cout << "Time step " << i << " from " << t0 << " to " << t1 << '\n';
				}

				// update temperature
				auto Tg = mfem_mgis::PartialQuadratureFunction::evaluate(
						mechanical_problem->getMaterial("Fuel").getPartialQuadratureSpacePointer(),
						[&t1, &mesh_parameters](const mfem_mgis::real x,
							const mfem_mgis::real y,
							const mfem_mgis::real) {
						constexpr auto Tref = mfem_mgis::real{600.15};
						constexpr auto dT = (mfem_mgis::real{1500} - Tref);
						const auto rp = mesh_parameters.rp;
						const auto dTe = dT * (t1 - tb) / (te - tb);
						const auto r = std::sqrt(x * x + y * y) / rp;
						const auto T = dTe * (1 - r * r) + Tref;
						return T;
						});
				mgis::behaviour::setExternalStateVariable(mechanical_problem->getMaterial("Fuel").s1,
						"Temperature", Tg->getValues());
				auto converged = false;
				auto iter = mfem_mgis::size_type{};
				auto mechanical_problem_initial_residual = mfem_mgis::real{};
				auto micromorphic_problem_initial_residual = mfem_mgis::real{};

				// set solver and preconditionner
				setLinearSolver(*mechanical_problem, a_solv, a_precond, p.verbosity_level, mfem_mgis::real{1e-8}, mfem_mgis::real{1e-6});
				setLinearSolver(*micromorphic_problem, a_solv, a_precond, p.verbosity_level, mfem_mgis::real{0}, mfem_mgis::real{1e-6});

				// alternate miminisation algorithm
				while (!converged) {
					constexpr auto reps = mfem_mgis::real{1e-4};
					if(mfem_mgis::getMPIrank() == 0) {
						std::cout << "time step " << i  //
							<< ", alternate minimisation iteration, " << iter << '\n';
					}
					if (iter == 0) {
						mechanical_problem->setSolverParameters({{"AbsoluteTolerance", 1e-10},
								{"RelativeTolerance", reps}});
						micromorphic_problem->setSolverParameters({{"AbsoluteTolerance", 1e-10},
								{"RelativeTolerance", reps}});
					} else {
						mechanical_problem->setSolverParameters(
								{{"AbsoluteTolerance",
								mechanical_problem_initial_residual * reps},
								{"RelativeTolerance", reps}});
						micromorphic_problem->setSolverParameters(
								{{"AbsoluteTolerance",
								micromorphic_problem_initial_residual * reps},
								{"RelativeTolerance", reps}});
					}

					// solving the mechanical problem
					auto mechanical_output = mechanical_problem->solve(t0, dt);
					if (!mechanical_output.status) {
						mfem_mgis::raise("non convergence of the mechanical problem");
					}
					// passing the energy release rate to the micromorphic problem
					Y = mfem_mgis::getInternalStateVariable(
							mechanical_problem->getMaterial("Fuel"), "EnergyReleaseRate");
					// solving the micromorphic problem
					auto micromorphic_output = micromorphic_problem->solve(t0, dt);
					if (!micromorphic_output.status) {
						mfem_mgis::raise("non convergence of the micromorphic problem");
					}
					// passing the damage to the mechanical problem
					d = mfem_mgis::getInternalStateVariable(
							micromorphic_problem->getMaterial("Fuel"), "Damage");
					if (iter == 0) {
						mechanical_problem_initial_residual =
							mechanical_output.initial_residual_norm;
						micromorphic_problem_initial_residual =
							micromorphic_output.initial_residual_norm;
					} else {
						converged = (mechanical_output.iterations == 0) &&
							(micromorphic_output.iterations == 0);
					}
					++iter;
					// check convergence
					if ((iter == iter_max) && (!converged)) {
						mfem_mgis::raise("non convergence of the fixed-point problem");
					}
				}
				mechanical_problem->executePostProcessings(t0, dt);
				micromorphic_problem->executePostProcessings(t0, dt);
				mechanical_problem->update();
				micromorphic_problem->update();
			}
			//
			return success ? EXIT_SUCCESS : EXIT_FAILURE;
		}
};
