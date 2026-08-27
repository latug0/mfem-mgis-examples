#include <cstdlib>
#include <functional>
#include <iostream>
#include <memory>
#include <sys/resource.h>
#include <sys/time.h>

#include <chrono>
#include <string>
#include <map>
#include <algorithm>
#include <vector>

#include "mfem/fem/datacollection.hpp"
#include "mfem/general/optparser.hpp"
#include "mfem/linalg/solvers.hpp"

#include "MGIS/Model/Model.hxx"
#include "MGIS/Behaviour/Integrate.hxx"

#include "MFEMMGIS/Config.hxx"
#include "MFEMMGIS/MFEMForward.hxx"
#include "MFEMMGIS/Material.hxx"
#include "MFEMMGIS/NonLinearEvolutionProblem.hxx"
#include "MFEMMGIS/NonLinearEvolutionProblemImplementation.hxx"
#include "MFEMMGIS/Parameters.hxx"
#include "MFEMMGIS/Profiler.hxx"
#include "MFEMMGIS/UniformDirichletBoundaryCondition.hxx"
#include "MFEMMGIS/UniformHeatSourceBoundaryCondition.hxx"
#include "MFEMMGIS/ParaviewExportIntegrationPointResultsAtNodes.hxx"
#include "MFEMMGIS/PhysicalSystem.hxx"
#include "MFEMMGIS/PointWiseModel.hxx"
#include "MFEMMGIS/NonLinearModel.hxx"
#include "MFEMMGIS/IterativeCouplingScheme.hxx"
#include "MFEMMGIS/FirstIterationConvergenceCriterion.hxx"
#include "MFEMMGIS/LoopCouplingScheme.hxx"
#include "MFEMMGIS/Simulation.hxx"
#include "MFEMMGIS/TimeStep.hxx"
#include "MFEMMGIS/BehaviourIntegratorBase.hxx"

#ifdef MFEM_USE_PETSC
#include "mfem/linalg/petsc.hpp"
#endif /* MFEM_USE_PETSC */

#ifdef MFEM_USE_PETSC
#include "mfem/linalg/mumps.hpp"
#endif /* MFEM_USE_MUMPS */

/* Project specific includes */
#include "headers/BoundaryConditions.hxx"
#include "headers/Setup.hxx"
#include "headers/Utils.hxx"
#include "headers/Debug.hxx"

void common_parameters(mfem::OptionsParser &args, TestParameters &p) {
  args.AddOption(&p.mesh_file, "-m", "--mesh", "Mesh file to use.");
  args.AddOption(&p.libraryU3SI2, "-lU", "--libraryU3SI2",
                 "Material library for said material.");
  args.AddOption(&p.libraryALFENI, "-lA", "--libraryALFENI",
                 "Material library for said material.");
  args.AddOption(&p.solver_thermo, "-svTh", "--solverTh",
                 "Solver for heat_transfer.");
  args.AddOption(&p.precond_thermo, "-pcTh", "--preconditionnerTh",
                 "Preconditionner for heat transfer.");
  args.AddOption(&p.solver_meca, "-svMc", "--solverMc",
                 "Solver for mechanics.");
  args.AddOption(&p.precond_meca, "-pcMc", "--preconditionnerMc",
                 "Preconditionner for mechanics.");
  args.AddOption(&p.order, "-o", "--order",
                 "Finite element order (polynomial degree).");
  args.AddOption(&p.refinement, "-r", "--refinement",
                 "refinement level of the mesh, default = 0");
  args.AddOption(&p.post_processing, "-p", "--post-processing",
                 "run post processing step");
  args.AddOption(&p.verbosity_level, "-v", "--verbosity-level",
                 "choose the verbosity level");
  args.AddOption(&p.debug, "-d", "--debug", "-nd", "--nodebug",
                 "Enable physics statistics debug output.");

  args.Parse();

  if (!args.Good()) {
    if (mfem_mgis::getMPIrank() == 0) args.PrintUsage(std::cout);
    mfem_mgis::finalize();
    exit(0);
  }
  if (mfem_mgis::getMPIrank() == 0) args.PrintOptions(std::cout);
  mfem_mgis::declareDefaultOptions(args);
}

int main(int argc, char *argv[]) {
  using namespace mfem_mgis;
  using namespace mfem;
  initialize(argc, argv);

  auto ctx = mgis::Context{};
  ctx.enableProfiling(true);
  auto or_die = ctx.getFatalFailureHandler();

  TestParameters p;

  OptionsParser args(argc, argv);
  common_parameters(args, p);

  auto power_history = [p](const double t) {
    const double t_ramp = 1e5;
    return (t <= t_ramp) ? p.source * (t / t_ramp) : p.source;
  };

  auto mesh =
      construct<MeshDiscretization>(
          ctx, ctx,
          Parameters{
              {"MeshFileName", p.mesh_file},
              {"Materials",
               Parameters{{"comb", 1}, {"gaine", 2}, {"stiffeners", 3}}},
              {"NumberOfUniformRefinements", p.parallel ? p.refinement : 0},
              {"Parallel", p.parallel}}) |
      or_die;

  auto heat_transfer_model =
      make_shared<NonLinearModel>(ctx, mesh,
                                  Parameters{{"FiniteElementFamily", "H1"},
                                             {"FiniteElementOrder", p.order},
                                             {"Hypothesis", "Tridimensional"},
                                             {"UnknownsSize", 1},
                                             {"Name", "Thermal"}}) |
      or_die;

  auto mechanics_model =
      make_shared<NonLinearModel>(ctx, mesh,
                                  Parameters{{"FiniteElementFamily", "H1"},
                                             {"FiniteElementOrder", p.order},
                                             {"Hypothesis", "Tridimensional"},
                                             {"UnknownsSize", 3},
                                             {"Name", "Mechanics"}}) |
      or_die;

  auto &heat_transfer = heat_transfer_model->getProblem();
  auto &mechanics = mechanics_model->getProblem();

  heat_transfer.setSolverParameters({{"VerbosityLevel", 2},
                                     {"RelativeTolerance", 1e-6},
                                     {"AbsoluteTolerance", 1e-6},
                                     {"MaximumNumberOfIterations", 10}});

  mechanics.setSolverParameters({{"VerbosityLevel", 2},
                                 {"RelativeTolerance", 1e-6},
                                 {"AbsoluteTolerance", 1e-6},
                                 {"MaximumNumberOfIterations", 10}});

  print_mesh_information(heat_transfer.getImplementation<true>());
  print_mesh_information(mechanics.getImplementation<true>());
  print_memory_footprint("After_problem_creation:");

  /*Test de récupération du déplacement pour l'envoyer à Robin*/
  auto mechanics_fed = mechanics.getFiniteElementDiscretizationPointer();
#ifdef MFEM_USE_MPI
  auto &mech_fes = mechanics_fed->getFiniteElementSpace<true>();
#else
  auto &mech_fes = mechanics_fed->getFiniteElementSpace<false>();
#endif
  double *u_data = mechanics.getUnknowns(mfem_mgis::ets).GetData();
  mfem::GridFunction u_mech(&mech_fes, u_data);
  /*Fin de test*/

  const auto setup =
      setup_properties(ctx, p, heat_transfer, mechanics, power_history);

  apply_boundary_conditions(heat_transfer, mechanics, p, power_history,
                            &u_mech);

  setLinearSolver(ctx, heat_transfer, "heat_transfer", p, p.verbosity_level);
  setLinearSolver(ctx, mechanics, "mechanics", p, p.verbosity_level);

  if (p.post_processing == 1) {
    add_post_processings(mechanics, "Results/Mechanics");
    add_post_processings(heat_transfer, "Results/Thermal");

    // // Exportation du swelling + déformation plastique
    // mfem_mgis::Parameters params_plast = {
    //     {"Results", "EquivalentPlasticStrain"},
    //     {"OutputFileName", "Resultats/MFront_Plasticite"}
    // };
    // mechanics.addPostProcessing("ParaviewExportIntegrationPointResultsAtNodes",
    // params_plast);

    mfem_mgis::Parameters params_swell;
    params_swell.insert("Results", "SwellingExport");
    params_swell.insert("OutputFileName", "Results/Swelling");

    std::vector<mfem_mgis::Parameter> mat_filter = {"comb"};
    params_swell.insert("Materials", mat_filter);

    mechanics.addPostProcessing("ParaviewExportIntegrationPointResultsAtNodes",
                                params_swell);
  }

  auto ps = construct<PhysicalSystem>(ctx, mesh) | or_die;

  auto c = std::make_shared<IterativeCouplingScheme>(ctx, mesh) | or_die;
  auto criterion = std::make_shared<FirstIterationConvergenceCriterion>();

  c->setMaximumNumberOfIterations(ctx, 10);
  c->addConvergenceCriterion(ctx, criterion);

  auto updater_model = std::make_shared<FieldUpdaterModel>(
      ctx, mesh, setup.fields[0].Pow_s0_sw, setup.fields[0].Pow_s1_sw,
      power_history);

  c->addModel(ctx, heat_transfer_model) | or_die;
  c->addModel(ctx, updater_model) | or_die;
  c->addModel(ctx, setup.swelling_model) | or_die;
  c->addModel(ctx, mechanics_model) | or_die;
  ps.setCouplingScheme(ctx, c) | or_die;

  // declaring the simulation
  const auto times =
      construct<Simulation::TimesDescription>(ctx, 0, p.duree, p.nbsteps) |
      or_die;
  auto s = construct<Simulation>(ctx, ctx, ps, times) |
           or_die;  // mgis::construct could be implemented such that it
                    // uses/transfers the first context we give
  // running the simulation
  const auto [status, output] = s.run(ctx);
  if (status != ExitStatus::success) {
    std::cerr << "simulation failed: " << ctx.getErrorMessage() << '\n';
    print_memory_footprint("After Solving:");
    return EXIT_FAILURE;
  }
  print_memory_footprint("After Solving:");

  if (p.debug) {
    debug_print_physics_stats(ctx, heat_transfer, mechanics, p.parallel, setup);
  }

  Profiler::OutputManager::printTimeTable(ctx);
  return EXIT_SUCCESS;
}