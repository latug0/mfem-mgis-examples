#include <cstdlib>
#include <fstream>
#include <iostream>
#include <cmath>
#include "mfem/general/optparser.hpp"
#include "MFEMMGIS/LinearSolverFactory.hxx"
#include "MFEMMGIS/Material.hxx"
#include "MFEMMGIS/NonLinearEvolutionProblem.hxx"
#include "MFEMMGIS/Config.hxx"
#include "MFEMMGIS/Profiler.hxx"
#include "ThirdMediumUtils.hxx"

int main(int argc, char **argv) {
  thirdmedium_utils::TestParameters p;
  p.mesh_file = "hertz2D.msh";
  mfem::OptionsParser args = parse_options(p, argc, argv);

  constexpr const auto dim = mfem_mgis::size_type{2};
  constexpr auto nsteps = mfem_mgis::size_type{50};
  constexpr auto H = mfem_mgis::real{0.2}; // Applied downward displacement
  constexpr auto te = mfem_mgis::real{1.0};
  
  auto factor = mfem_mgis::real{p.gamma};
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
  std::string library = p.library;

  #if defined(MFEM_USE_MUMPS) && defined(MFEM_USE_MPI)
    constexpr bool parallel = true;
  #else
    constexpr bool parallel = false;
  #endif

  if (mfem_mgis::getMPIrank() == 0) {
    args.PrintOptions(std::cout);
  }

  mfem_mgis::NonLinearEvolutionProblem mechanics(ctx, {
      {"MeshFileName", p.mesh_file},
      {"FiniteElementFamily", "H1"},
      {"FiniteElementOrder", p.order},
      {"UnknownsSize", dim},
      {"NumberOfUniformRefinements", p.parallel ? p.refinement : 0},
      {"Hypothesis", "PlaneStrain"},
      {"Parallel", p.parallel}
  });

    // Mesh and memory info
  thirdmedium_utils::print_mesh_information(mechanics.getImplementation<parallel>());
  thirdmedium_utils::print_memory_footprint("After_problem:");
  thirdmedium_utils::save_unknowns_info(mechanics.getImplementation<parallel>(),p.output_dir);


  auto faltus_parameters = mfem_mgis::Parameters{
      {"Regularization", mfem_mgis::Parameters{{"Faltus2026", mfem_mgis::Parameters{{"PenalizationCoefficient", p.alpha}}}}}};

  // Assign Behaviours
  mechanics.addBehaviourIntegrator(ctx, "Mechanics", "VOID", library, behaviour2, faltus_parameters) | or_die;
  mechanics.addBehaviourIntegrator(ctx, "Mechanics", "FOUNDATION", library, behaviour1) | or_die;
  mechanics.addBehaviourIntegrator(ctx, "Mechanics", "INDENTER", library, behaviour1) | or_die;

  // Material Properties
  auto &m_void = mechanics.getMaterial(ctx, "VOID", 0) | or_die;
  setMaterialProperty(ctx, m_void.s0, "BulkModulus", factor * Ks) | or_die;
  setMaterialProperty(ctx, m_void.s1, "BulkModulus", factor * Ks) | or_die;
  setMaterialProperty(ctx, m_void.s0, "ShearModulus", factor * Gs) | or_die;
  setMaterialProperty(ctx, m_void.s1, "ShearModulus", factor * Gs) | or_die;

  for (const auto &n : {"FOUNDATION", "INDENTER"}) {
    auto &m = mechanics.getMaterial(ctx, n, 0) | or_die;
    setMaterialProperty(ctx, m.s0, "YoungModulus", E) | or_die;
    setMaterialProperty(ctx, m.s1, "YoungModulus", E) | or_die;
    setMaterialProperty(ctx, m.s0, "PoissonRatio", nu) | or_die;
    setMaterialProperty(ctx, m.s1, "PoissonRatio", nu) | or_die;
  }

  // Boundary Conditions
  mechanics.addUniformDirichletBoundaryCondition(ctx,{{"Boundary", "FOUND_BOT"}, {"Component", 0}})|or_die;
  mechanics.addUniformDirichletBoundaryCondition(ctx,{{"Boundary", "FOUND_BOT"}, {"Component", 1}})|or_die;
  mechanics.addUniformDirichletBoundaryCondition(ctx,{{"Boundary", "SYM_FOUND"}, {"Component", 0}})|or_die;
  mechanics.addUniformDirichletBoundaryCondition(ctx,{{"Boundary", "SYM_IND"}, {"Component", 0}})|or_die;
  
  mechanics.addUniformDirichletBoundaryCondition(ctx,
    {
      {"Boundary", "IND_TOP"},
      {"Component", 1},
      {"LoadingEvolution", [](const auto t) { return -H * (t / te); }}
  })|or_die;

  mechanics.setPredictionPolicy({.strategy = mfem_mgis::PredictionStrategy::BEGINNING_OF_TIME_STEP_PREDICTION});
  mechanics.setSolverParameters(ctx, {{"VerbosityLevel", p.verbosity_level}, {"RelativeTolerance", 1e-4}, {"AbsoluteTolerance", 0.}, {"MaximumNumberOfIterations", p.newton_iterations}}) | or_die;

  auto solverParameters = thirdmedium_utils::set_solver_parameters(p);
  mfem_mgis::Parameters prec = thirdmedium_utils::set_preconditioner_parameters(p); 

  if constexpr (parallel) {
    if (p.preconditioner != -1 && p.solver != 0)
      solverParameters.insert(mfem_mgis::may_throw, mfem_mgis::Parameters{{"Preconditioner", prec}});
    mechanics.setLinearSolver(ctx, p.solver_name, solverParameters) | or_die; 
  } else {
    mechanics.setLinearSolver(ctx, "UMFPackSolver", {}) | or_die;
  }

  if (p.post_processing) {
    mechanics.addPostProcessing(ctx, "ParaviewExportResults", {{"OutputFileName", p.output_dir + std::string("hertz2D_disp")}}) | or_die;
    auto results = std::vector<mfem_mgis::Parameter>{"StressExport","JacobiExport"};
    auto third_medium_materials = std::vector<mfem_mgis::Parameter>{"VOID"};
    //mechanics.addPostProcessing(ctx, "ParaviewExportIntegrationPointResultsAtNodes", {{"OutputFileName", p.output_dir + std::string("ResultsAtNodes")}, {"Materials", third_medium_materials}, {"Results", results}}) | or_die;
  }

  std::vector<std::string> materials = {"VOID"};
  thirdmedium_utils::setup_additional_convergence_criteria(mechanics, ctx, p, solverParameters, factor * Ks, factor * Gs, std::move(materials));

  auto [exit_status, sim_output] = thirdmedium_utils::run_solve(ctx, mechanics, 0, te, nsteps);

  if (mfem_mgis::getMPIrank() == 0) {
    std::cout << "Exit status : " << (exit_status == mfem_mgis::ExitStatus::success ? "Success\n" : "Failed\n");
    if (sim_output.has_value()) {
        thirdmedium_utils::check_and_save_convergence_info(sim_output.value(), exit_status);
    }
  }

  mfem_mgis::Profiler::OutputManager::printTimeTable(ctx);
  return EXIT_SUCCESS;
}