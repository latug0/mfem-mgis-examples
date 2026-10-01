#include <memory>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include "mfem/general/optparser.hpp"
#include "MFEMMGIS/Simulation.hxx"
#include "MFEMMGIS/LinearSolverFactory.hxx"
#include "MFEMMGIS/Material.hxx"
#include "MFEMMGIS/MechanicalPostProcessings.hxx"
#include "MFEMMGIS/NonLinearEvolutionProblem.hxx"
#include "MFEMMGIS/Config.hxx"
#include "MFEMMGIS/Profiler.hxx"
#include "ThirdMediumUtils.hxx"

int main(int argc, char **argv) {
  thirdmedium_utils::TestParameters p;
  p.mesh_file = "point_load_approx_full.msh";
  mfem::OptionsParser args = parse_options(p, argc, argv);

  constexpr const auto dim = mfem_mgis::size_type{2};
  constexpr auto nsteps = mfem_mgis::size_type{100};
  constexpr auto H_disp = mfem_mgis::real{0.95};  // Total downward displacement
  constexpr auto te = mfem_mgis::real{1};
  auto factor = mfem_mgis::real{p.gamma};

  // Material parameters
  constexpr auto Ks = mfem_mgis::real{200e9};
  constexpr auto Gs = mfem_mgis::real{200e9};
  constexpr auto E = mfem_mgis::real{210e9};
  constexpr auto nu = mfem_mgis::real{0.3};

  mfem_mgis::initialize(argc, argv);

  auto ctx = mfem_mgis::Context{};
  auto or_die = ctx.getFatalFailureHandler();
  ctx.enableProfiling(true);

  const char *behaviour1 = "Elasticity";
  const char *behaviour2 = "ThirdMediumContactHyperElasticBehaviour";
  const char *library = "src/libBehaviour.so";

#if defined(MFEM_USE_MUMPS) && defined(MFEM_USE_MPI)
  constexpr bool parallel = true;
#else
  constexpr bool parallel = false;
#endif

  if (mfem_mgis::getMPIrank() == 0) {
    args.PrintOptions(std::cout);
    std::cout << "Parallel : " << p.parallel << '\n';
  }

  mfem_mgis::NonLinearEvolutionProblem mechanics(
      ctx, {{"MeshFileName", p.mesh_file},
            {"FiniteElementFamily", "H1"},
            {"FiniteElementOrder", p.order},
            {"UnknownsSize", dim},
            {"NumberOfUniformRefinements", p.parallel ? p.refinement : 0},
            {"Hypothesis", "PlaneStrain"},
            {"Parallel", p.parallel}});

  thirdmedium_utils::print_mesh_information(
      mechanics.getImplementation<parallel>());

  auto faltus_parameters = mfem_mgis::Parameters{
      {"Regularization",
       mfem_mgis::Parameters{
           {"Faltus2026",
            mfem_mgis::Parameters{{"PenalizationCoefficient", p.alpha}}}}}};

  // Map behaviours strictly to the new .geo surfaces[cite: 1]
  mechanics.addBehaviourIntegrator(ctx, "Mechanics", "VOID", library,
                                   behaviour2, faltus_parameters) |
      or_die;
  mechanics.addBehaviourIntegrator(ctx, "Mechanics", "SOLID", library,
                                   behaviour1) |
      or_die;

  // Material properties for VOID[cite: 1]
  auto &mV = mechanics.getMaterial(ctx, "VOID", 0) | or_die;
  setMaterialProperty(ctx, mV.s0, "BulkModulus", factor * Ks) | or_die;
  setMaterialProperty(ctx, mV.s1, "BulkModulus", factor * Ks) | or_die;
  setMaterialProperty(ctx, mV.s0, "ShearModulus", factor * Gs) | or_die;
  setMaterialProperty(ctx, mV.s1, "ShearModulus", factor * Gs) | or_die;

  // Material properties for SOLID[cite: 1]
  auto &mS = mechanics.getMaterial(ctx, "SOLID", 0) | or_die;
  setMaterialProperty(ctx, mS.s0, "YoungModulus", E) | or_die;
  setMaterialProperty(ctx, mS.s1, "YoungModulus", E) | or_die;
  setMaterialProperty(ctx, mS.s0, "PoissonRatio", nu) | or_die;
  setMaterialProperty(ctx, mS.s1, "PoissonRatio", nu) | or_die;

  // Boundary Conditions setup
  mechanics.addUniformDirichletBoundaryCondition(
      ctx, {{"Boundary", "LEFT_WALL"}, {"Component", 0}}) |
      or_die;
  mechanics.addUniformDirichletBoundaryCondition(
      ctx, {{"Boundary", "LEFT_WALL"}, {"Component", 1}}) |
      or_die;

  // Displacement driven loading acting on Gamma_D
  mechanics.addUniformDirichletBoundaryCondition(
      ctx,
      {{"Boundary", "GAMMA_D"},
       {"Component", 1},
       {"LoadingEvolution", [](const auto t) { return -H_disp * (t / te); }}}) |
      or_die;

  mechanics.setPredictionPolicy(
      {.strategy =
           mfem_mgis::PredictionStrategy::BEGINNING_OF_TIME_STEP_PREDICTION});

  mechanics.setSolverParameters(
      ctx, {{"VerbosityLevel", p.verbosity_level},
            {"RelativeTolerance", 1e-4},
            {"AbsoluteTolerance", 0.},
            {"MaximumNumberOfIterations", p.newton_iterations}}) |
      or_die;

  auto solverParameters = thirdmedium_utils::set_solver_parameters(p);
  mfem_mgis::Parameters prec =
      thirdmedium_utils::set_preconditioner_parameters(p);

  if constexpr (parallel) {
    if (p.preconditioner != -1)
      solverParameters.insert(mfem_mgis::may_throw,
                              mfem_mgis::Parameters{{"Preconditioner", prec}});
    mechanics.setLinearSolver(ctx, p.solver_name, solverParameters) | or_die;
  } else {
    mechanics.setLinearSolver(ctx, "UMFPackSolver", {}) | or_die;
  }

  mechanics.addPostProcessing(
      ctx, "ParaviewExportResults",
      {{"OutputFileName", std::string("c_shape_results")}}) |
      or_die;

  // Reduced iteration check array based on the new single void body
  std::vector<std::string> materials = {"VOID"};
  thirdmedium_utils::setup_additional_convergence_criteria(
      mechanics, ctx, p, solverParameters, factor * Ks, factor * Gs,
      std::move(materials));

  auto [exit_status, sim_output] =
      thirdmedium_utils::run_solve(ctx, mechanics, 0, te, nsteps);

  if (mfem_mgis::getMPIrank() == 0) {
    std::cout << "Exit status : "
              << (exit_status == mfem_mgis::ExitStatus::success ? "Success\n"
                                                                : "Failed\n");
    std::ofstream status("run.status");
    status << (exit_status == mfem_mgis::ExitStatus::success ? "SUCCESS\n"
                                                             : "FAILED\n");

    if (sim_output.has_value()) {
      thirdmedium_utils::check_and_save_convergence_info(sim_output.value(),
                                                         exit_status);
    }
  }

  thirdmedium_utils::print_memory_footprint("After Solving:");
  mfem_mgis::Profiler::OutputManager::printTimeTable(ctx);

  return EXIT_SUCCESS;
}
