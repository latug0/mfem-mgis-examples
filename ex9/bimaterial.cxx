/*!
 * \file   bimaterial.cxx
 * \brief
 * \author th202608
 * \date   27/08/2026
 */

#include <memory>
#include <cstdlib>
#include <iostream>
#include "MFEMMGIS/Profiler.hxx"
#include "MFEMMGIS/Parameters.hxx"
#include "MFEMMGIS/Material.hxx"
#include "MFEMMGIS/UniformDirichletBoundaryCondition.hxx"
#include "MFEMMGIS/NonLinearEvolutionProblem.hxx"
#include "MFEMMGIS/NonLinearEvolutionProblemImplementation.hxx"
#include "MFEMMGIS/MechanicalPostProcessings.hxx"
#include "MFEMMGIS/ParaviewExportIntegrationPointResultsAtNodes.hxx"
#include "MFEMMGIS/Simulation.hxx"

int main(int argc, char** argv) {
  constexpr auto parallel = true;
  auto ctx = mgis::Context{};
  auto or_die = ctx.getFatalFailureHandler();
  constexpr const auto dim = mfem_mgis::size_type{2};
  // options treatment
  mfem_mgis::initialize(argc, argv);
  auto success = true;
  // building the non linear problem
  mfem_mgis::NonLinearEvolutionProblem problem(
      ctx, {{"MeshFileName", "bimaterial.msh"},
            {"FiniteElementFamily", "H1"},
            {"FiniteElementOrder", 3},
            {"UnknownsSize", dim},
            {"NumberOfUniformRefinements", 0},  // faster for testing
            //{"NumberOfUniformRefinements", parameters.parallel ? 1 : 0},
            {"Hypothesis", "PlaneStrain"},
            {"Parallel", parallel}});
  // solving the problem without petsc
  if (!mfem_mgis::usePETSc()) {
    problem.setPredictionPolicy(
        {.strategy =
             mfem_mgis::PredictionStrategy::BEGINNING_OF_TIME_STEP_PREDICTION});
    problem.setSolverParameters(ctx, {{"VerbosityLevel", 2},
                                      {"RelativeTolerance", 1e-6},
                                      {"AbsoluteTolerance", 0.},
                                      {"MaximumNumberOfIterations", 10}}) |
        or_die;
    if (parallel) {
      problem.setLinearSolver(ctx, "MUMPSSolver", {}) | or_die;
    } else {
      problem.setLinearSolver(ctx, "UMFPackSolver", {}) | or_die;
    }
  }
  // materials
  const auto& b1 =
      problem.addBehaviourIntegrator(ctx, "Mechanics", "FirstMaterial",
                                     "src/libBehaviour.so", "Elasticity") |
      or_die;
  const auto& b2 =
      problem.addBehaviourIntegrator(ctx, "Mechanics", "SecondMaterial",
                                     "src/libBehaviour.so", "Plasticity") |
      or_die;
  for (const auto& material : {"FirstMaterial", "SecondMaterial"}) {
    auto& m = problem.getMaterial(ctx, material, 0) | or_die;
    mgis::behaviour::setExternalStateVariable(ctx, m.s0, "Temperature",
                                              293.15) |
        or_die;
    mgis::behaviour::setExternalStateVariable(ctx, m.s1, "Temperature",
                                              293.15) |
        or_die;
  }
  problem.addBoundaryCondition(
      ctx, std::make_unique<mfem_mgis::UniformDirichletBoundaryCondition>(
               problem.getFiniteElementDiscretizationPointer(), "left", 0)) |
      or_die;
  problem.addBoundaryCondition(
      ctx, std::make_unique<mfem_mgis::UniformDirichletBoundaryCondition>(
               problem.getFiniteElementDiscretizationPointer(), "bottom", 1)) |
      or_die;
  problem.addBoundaryCondition(
      ctx, std::make_unique<mfem_mgis::UniformDirichletBoundaryCondition>(
               problem.getFiniteElementDiscretizationPointer(), "right", 0,
               [](const auto t) noexcept { return 20e-4 * t; })) |
      or_die;
  // post-processings
  problem.addPostProcessing(
      ctx, "ComputeResultantForceOnBoundary",
      {{"Boundary", "right"}, {"OutputFileName", "force.txt"}}) |
      or_die;
  problem.addPostProcessing(ctx, "ParaviewExportResults",
                            {{"OutputFileName", "Displacements"}}) |
      or_die;
  problem.addPostProcessing(ctx, "ParaviewExportIntegrationPointResultsAtNodes",
                            {{{"Results", "Stress"},
                              {"OutputFileName", "Stress"}}}) |
      or_die;
  problem.addPostProcessing(ctx, "ParaviewExportIntegrationPointResultsAtNodes",
                            {{{"Materials", "SecondMaterial"},
                              {"Results",
                               std::vector<mfem_mgis::Parameter>{
                                   "Stress", "EquivalentPlasticStrain"}},  //
                              {"OutputFileName", "Test"}}}) |
      or_die;
  problem.getImplementation<parallel>().addPostProcessing(
      ctx,
      std::make_unique<
          mfem_mgis::
              ParaviewExportIntegrationPointPostProcessingsResultsAtNodes<
                  true>>(
          ctx, problem.getImplementation<true>(), "EquivalentPlasticStrain",
          std::vector<int>{
              problem.getMaterialIdentifier(ctx, "SecondMaterial") | or_die},
          1,
          [&problem, &or_die](mfem_mgis::Context& ctx,
                              mfem_mgis::PartialQuadratureFunction& f) {
            const auto mid = f.getPartialQuadratureSpace().getId();
            const auto& bi =
                problem.getBehaviourIntegrator(ctx, mid, 0) | or_die;
            const auto& m = bi.getMaterial(ctx) | or_die;
            assign_values(ctx, f,
                          mfem_mgis::getInternalStateVariable(
                              ctx, m, "EquivalentPlasticStrain",
                              mfem_mgis::Material::END_OF_TIME_STEP) |
                              or_die) |
                or_die;
            return true;
          },
          "EquivalentPlasticStrainOutput")) |
      or_die;
  problem.getImplementation<parallel>().addPostProcessing(
      ctx, std::make_unique<
               mfem_mgis::
                   ParaviewExportIntegrationPointPostProcessingsResultsAtNodes<
                       true>>(
               ctx, problem.getImplementation<true>(), "vonMisesStress",
               problem.getAssignedMaterialsIdentifiers(), 1,
               [&problem, &or_die](mfem_mgis::Context& ctx,
                                   mfem_mgis::PartialQuadratureFunction& f) {
                 const auto mid = f.getPartialQuadratureSpace().getId();
                 const auto& bi =
                     problem.getBehaviourIntegrator(ctx, mid, 0) | or_die;
                 const auto& m = bi.getMaterial(ctx) | or_die;
                 return mfem_mgis::computeVonMisesEquivalentStress(
                     ctx, f, m, mfem_mgis::Material::END_OF_TIME_STEP);
               },
               "vonMisesStressOutput")) |
      or_die;
  problem.getImplementation<parallel>().addPostProcessing(
      ctx, std::make_unique<
               mfem_mgis::
                   ParaviewExportIntegrationPointPostProcessingsResultsAtNodes<
                       true>>(
               ctx, problem.getImplementation<true>(), "vonMisesStress",
               std::vector<int>{
                   problem.getMaterialIdentifier(ctx, "SecondMaterial") |
                   or_die},
               1,
               [&problem, &or_die](mfem_mgis::Context& ctx,
                                   mfem_mgis::PartialQuadratureFunction& f) {
                 const auto mid = f.getPartialQuadratureSpace().getId();
                 const auto& bi =
                     problem.getBehaviourIntegrator(ctx, mid, 0) | or_die;
                 const auto& m = bi.getMaterial(ctx) | or_die;
                 return mfem_mgis::computeVonMisesEquivalentStress(
                     ctx, f, m, mfem_mgis::Material::END_OF_TIME_STEP);
               },
               "vonMisesStressOutput-2")) |
      or_die;
  // running the simulation
  const auto times = mfem_mgis::Simulation::TimesDescription{0, 1, 30};
  auto s = mfem_mgis::Simulation{ctx, problem, times};
  auto [status, ooutput] = s.run(ctx);
  success = status.shallContinue();
  if (!success) {
    ctx.log() << ctx.getErrorMessage() << '\n';
  }
  //
  mfem_mgis::Profiler::OutputManager::printTimeTable(ctx);
  return success ? EXIT_SUCCESS : EXIT_FAILURE;
}
