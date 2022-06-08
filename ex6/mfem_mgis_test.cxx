#pragma once




// first example (from ex3 / Inclusion)
int executeMFEMMGISTest(const TestParameters& p, const bool use_post_processing, const solver_name a_solv, const precond_name a_precond, gather_information& a_info) {
	constexpr const auto dim = mfem_mgis::size_type{3};
	
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
			    {"NumberOfUniformRefinements", p.parallel ? 2 : 0},
			    {"Parallel", p.parallel}});

	mfem_mgis::PeriodicNonLinearEvolutionProblem problem(fed);

	//getMesh
	auto mesh = problem.getImplementation<true>().getFiniteElementSpace().GetMesh();
	//get the number of vertices
	int numbers_of_vertices = mesh->GetNV();
	//get the number of elements
	int numbers_of_elements = mesh->GetNE();
	//get the element size
	double h = mesh->GetElementSize(0);

	profiling::output::printMessage("INFO: number of vertices -> ", numbers_of_vertices);
	profiling::output::printMessage("INFO: number of elements -> ", numbers_of_elements);
	profiling::output::printMessage("INFO: element size -> ", h);

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
		measure = profiling::timers::chrono_section( [&](){
			solverStatistics = problem.solve(0, 1);
		});

		// check status
		if (!solverStatistics) {
			profiling::output::printMessage("INFO: ", string_solver,"+",string_precond," FAILED");
		}
	} // end timer Solve

	if(use_post_processing)
	{
		START_TIMER("Postprocessing_step");
		problem.executePostProcessings(0, 1);
	} // end timer Postprocessing_step

	{
		START_TIMER("get_statistics");
		measure = profiling::output::reduce_max(measure);
		// fill info data.
		a_info.add(
				info{
					a_solv, 
					a_precond,
					solverStatistics.status, 
					solverStatistics.iterations, 
					solverStatistics.final_residual_norm,
					measure
					}
			);
		// print information
                profiling::output::printMessage(
				"INFO: ", 
				string_solver,
				"+",
				string_precond,
				" works correctly, elapsed time: ", 
				measure
			);

	} // end timer get_statistics
	
	profiling::output::printMessage("");
	
	return(EXIT_SUCCESS);
}




// thomas example
#include "MFEMMGIS/Config.hxx"
#include "MFEMMGIS/Parameters.hxx"
#include "MFEMMGIS/Material.hxx"
#include "MFEMMGIS/PartialQuadratureSpace.hxx"
#include "MFEMMGIS/ImposedDirichletBoundaryConditionAtClosestNode.hxx"
#include "MFEMMGIS/NonLinearEvolutionProblem.hxx"

struct MeshParameters {
	// fuel pellet radius
	const mfem_mgis::real rp = mfem_mgis::real{4.085e-3};
	// fuel pellet height
	const mfem_mgis::real hp = mfem_mgis::real{13.5e-3};
	// dishing height
	const mfem_mgis::real dh = mfem_mgis::real{3.2e-4};
};

static std::shared_ptr<mfem_mgis::NonLinearEvolutionProblem>
buildMechanicalProblem(const mfem_mgis::Parameters& common_problem_parameters,
                       const MeshParameters mesh_parameters) {

	START_TIMER("buildMechanicalProblem");
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



static std::shared_ptr<mfem_mgis::NonLinearEvolutionProblem>
buildMicromorphicProblem(
    const mfem_mgis::Parameters& common_problem_parameters) {
 
  START_TIMER("buildMicromorphicProblem");
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
                                  "MicromorphicDamageII");
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


int thomas_main_function(const TestParameters& p, const bool use_post_processing, const solver_name a_solv, const precond_name a_precond, gather_information& a_info)
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
      {"NumberOfUniformRefinements", parallel ? 2 : 0},
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

  for (mfem_mgis::size_type i = 0; i != n; ++i) {
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
    setLinearSolver(*micromorphic_problem, a_solv, a_precond);
    setLinearSolver(*mechanical_problem, a_solv, a_precond);

    // alternate miminisation algorithm
    while (!converged) {
      constexpr auto reps = mfem_mgis::real{1e-4};
      if(mfem_mgis::getMPIrank() == 0) {
        std::cout << "time step " << i  //
                  << ", alternate minimisation iteration, " << iter << '\n';
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
