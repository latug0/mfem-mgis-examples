#include <cmath>
#include <string>
#include <vector>
#include <cstdlib>
#include <string_view>

#include "mfem/general/optparser.hpp"

#include "MFEMMGIS/Config.hxx"
#include "MFEMMGIS/Parameters.hxx"
#include "MFEMMGIS/Profiler.hxx"
#include "MFEMMGIS/MeshDiscretization.hxx"
#include "MFEMMGIS/NonLinearEvolutionProblemImplementation.hxx"
#include "MFEMMGIS/NonLinearModel.hxx"
#include "MFEMMGIS/PhysicalSystem.hxx"
#include "MFEMMGIS/IterativeCouplingScheme.hxx"
#include "MFEMMGIS/FirstIterationConvergenceCriterion.hxx"
#include "MFEMMGIS/Simulation.hxx"

/* Project specific includes */
#include "headers/BoundaryConditions.hxx"
#include "headers/Setup.hxx"
#include "headers/Utils.hxx"
#include "headers/Debug.hxx"

// @see Setup.hxx
void common_parameters(mfem::OptionsParser& args, TestParameters& p) {
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
  args.AddOption(&p.post_processing, "-p", "--post-processing", "-no-p",
                 "--no-post-processing", "Export the results to Paraview.");
  args.AddOption(&p.verbosity_level, "-v", "--verbosity-level",
                 "Verbosity level of the linear solvers.");
  args.AddOption(&p.debug, "-d", "--debug", "-nd", "--nodebug",
                 "Print the statistics of the fields.");
  args.AddOption(&p.reference_file, "-rf", "--reference-file",
                 "Reference statistics of the fields, no comparison if "
                 "empty.");
  args.AddOption(&p.duree, "-dur", "--duree",
                 "Total simulation duration, default = 1e5");
  args.AddOption(&p.nbsteps, "-ns", "--nbsteps",
                 "Number of time steps, default = 1");
  args.AddOption(&p.t_ramp, "-tr", "--t-ramp",
                 "Duration of the power ramp (0 disables it), default = 1e5");
  args.AddOption(&p.h_conv, "-hc", "--h-conv",
                 "Thermal convection coefficient, default = 5e4");
  args.AddOption(&p.water_pressure, "-wp", "--water-pressure",
                 "Coolant pressure, default = 1e6");

  mfem_mgis::declareDefaultOptions(args);
  args.Parse();

  if (args.Help()) {
    args.PrintUsage(mfem_mgis::getOutputStream());
    mfem_mgis::finalize();
    std::exit(EXIT_SUCCESS);
  }
  if (!args.Good()) {
    args.PrintUsage(mfem_mgis::getOutputStream());
    mfem_mgis::abort(EXIT_FAILURE);
  }
  args.PrintOptions(mfem_mgis::getOutputStream());
}

int main(int argc, char* argv[]) {
  using namespace mfem_mgis;
  using namespace mfem;
  initialize(argc, argv);

  auto ctx = mgis::Context{};
  ctx.enableProfiling(true);
  auto or_die = ctx.getFatalFailureHandler();

  TestParameters p;

  OptionsParser args(argc, argv);
  common_parameters(args, p);

  if (p.t_ramp < 0) {
    ctx.log() << "the duration of the power ramp must not be negative\n";
    finalize();
    return EXIT_FAILURE;
  }
  const auto ramp_steps = p.t_ramp * p.nbsteps / p.duree;
  if ((p.t_ramp < p.duree) &&
      (std::abs(ramp_steps - std::round(ramp_steps)) > 1e-9)) {
    ctx.log() << "the end of the power ramp (t = " << p.t_ramp
              << " s) must be a time step boundary\n";
    finalize();
    return EXIT_FAILURE;
  }

  auto power_history = [p](const double t) {
    return (t < p.t_ramp) ? p.source * (t / p.t_ramp) : p.source;
  };

  auto mesh =
      construct<MeshDiscretization>(
          ctx, ctx,
          Parameters{{"MeshFileName", p.mesh_file},
                     {"Materials",
                      Parameters{{"comb", 1}, {"gaine", 2}, {"stiffeners", 3}}},
                     {"NumberOfUniformRefinements", p.refinement},
                     {"Parallel", true}}) |
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

  auto& heat_transfer = heat_transfer_model->getProblem();
  auto& mechanics = mechanics_model->getProblem();

  heat_transfer.setSolverParameters(ctx, {{"VerbosityLevel", 2},
                                          {"RelativeTolerance", 1e-6},
                                          {"AbsoluteTolerance", 1e-6},
                                          {"MaximumNumberOfIterations", 10}}) |
      or_die;

  mechanics.setSolverParameters(ctx, {{"VerbosityLevel", 2},
                                      {"RelativeTolerance", 1e-6},
                                      {"AbsoluteTolerance", 1e-6},
                                      {"MaximumNumberOfIterations", 10}}) |
      or_die;

  print_mesh_information(heat_transfer.getImplementation<true>());
  print_mesh_information(mechanics.getImplementation<true>());
  print_memory_footprint("After_problem_creation:");

  auto mechanics_fed = mechanics.getFiniteElementDiscretizationPointer();
  mfem::ParGridFunction u_mech(&mechanics_fed->getFiniteElementSpace<true>());
  u_mech = 0.0;

  const auto setup = setup_properties(mfem_mgis::may_abort, ctx, p,
                                      heat_transfer, mechanics, power_history);

  apply_boundary_conditions(mfem_mgis::may_abort, ctx, heat_transfer, mechanics,
                            p, power_history, &u_mech);

  setLinearSolver(mfem_mgis::may_abort, ctx, heat_transfer, "heat_transfer", p,
                  p.verbosity_level);
  setLinearSolver(mfem_mgis::may_abort, ctx, mechanics, "mechanics", p,
                  p.verbosity_level);

  if (p.post_processing) {
    add_post_processings(mfem_mgis::may_abort, ctx, mechanics,
                         "Results/Mechanics", "Displacement");
    add_post_processings(mfem_mgis::may_abort, ctx, heat_transfer,
                         "Results/Thermal", "Temperature");
    // swelling in the fuel
    mechanics.addPostProcessing(
        ctx, "ParaviewExportIntegrationPointResultsAtNodes",
        {{"Results", "SwellingExport"},
         {"OutputFileName", "Results/Swelling"},
         {"Materials", std::vector<mfem_mgis::Parameter>{"comb"}}}) |
        or_die;
  }

  auto ps = mfem_mgis::construct<PhysicalSystem>(ctx, mesh) | or_die;

  auto c = mfem_mgis::make_shared<IterativeCouplingScheme>(ctx, mesh) | or_die;
  auto criterion =
      mfem_mgis::make_shared<FirstIterationConvergenceCriterion>(ctx);

  c->setMaximumNumberOfIterations(ctx, 10) | or_die;
  c->addConvergenceCriterion(ctx, criterion) | or_die;

  auto updater_model = std::make_shared<FieldUpdaterModel>(
      ctx, mesh, setup.fields[0].Pow_s0_sw, setup.fields[0].Pow_s1_sw,
      power_history, &u_mech, &mechanics.getUnknowns(mfem_mgis::ets));

  c->addModel(ctx, updater_model) | or_die;
  c->addModel(ctx, heat_transfer_model) | or_die;
  c->addModel(ctx, setup.swelling_model) | or_die;
  c->addModel(ctx, mechanics_model) | or_die;
  ps.setCouplingScheme(ctx, c) | or_die;

  // declaring the simulation
  const auto times =
      construct<Simulation::TimesDescription>(ctx, 0, p.duree, p.nbsteps) |
      or_die;
  auto s = construct<Simulation>(ctx, ctx, ps, times) | or_die;
  // running the simulation
  const auto [status, output] = s.run(ctx);
  if (status != ExitStatus::success) {
    getErrorStream() << "simulation failed: " << ctx.getErrorMessage() << '\n';
    print_memory_footprint("After Solving:");
    return EXIT_FAILURE;
  }
  print_memory_footprint("After Solving:");

  const auto stats = computePhysicsStatistics(heat_transfer, mechanics, setup);
  if (p.debug) {
    printPhysicsStatistics(getOutputStream(), stats);
  }
  Profiler::OutputManager::printTimeTable(ctx);
  // the swelling is always compared to its exact value
  auto success = checkSwelling(stats.at("Swelling"), p);
  if (!std::string_view{p.reference_file}.empty()) {
    success = checkPhysicsStatistics(stats, p.reference_file) && success;
  }
  return success ? EXIT_SUCCESS : EXIT_FAILURE;
}
