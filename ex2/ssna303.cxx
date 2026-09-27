/*!
 * \file   ssna3030.cxx
 * \brief
 * \author Thomas Helfer, Guillaume Latu
 * \date   06/04/2021
 */

#include <memory>
#include <cstdlib>
#include <iostream>
#include "mfem/general/optparser.hpp"
#include "MFEMMGIS/Material.hxx"
#include "MFEMMGIS/NonLinearEvolutionProblem.hxx"

int main(int argc, char** argv) {
  auto ctx = mgis::Context{};
  auto or_die = ctx.getFatalFailureHandler();
  constexpr const auto dim = mfem_mgis::size_type{2};
  // Initialize mfem_mgis (it includes a call to MPI_Init)
  mfem_mgis::initialize(argc, argv);

  const char* mesh_file = "ssna303.msh";
  const char* behaviour = "Plasticity";
  const char* library = "src/libBehaviour.so";
#if defined(MFEM_USE_MUMPS) && defined(MFEM_USE_MPI)
  bool parallel = true;
#else
  bool parallel = false;
#endif
  auto order = 1;
  // options treatment
  mfem::OptionsParser args(argc, argv);
  mfem_mgis::declareDefaultOptions(args);
  args.AddOption(&order, "-o", "--order",
                 "Finite element order (polynomial degree).");
  args.AddOption(&parallel, "-p", "--parallel", "-no-p", "--no-parallel",
                 "Perform parallel computations.");
  args.Parse();
  if (!args.Good()) {
    args.PrintUsage(std::cout);
    mfem_mgis::abort(EXIT_FAILURE);
  }
  args.PrintOptions(std::cout);
  // the non linear problem
  auto problem = mfem_mgis::construct<mfem_mgis::NonLinearEvolutionProblem>(ctx,
									    mfem_mgis::Parameters{{"MeshFileName", mesh_file},
                                                {"FiniteElementFamily", "H1"},
                                                {"FiniteElementOrder", order},
                                                {"UnknownsSize", dim},
                                                {"Hypothesis", "PlaneStrain"},
                                                {"Parallel", parallel}})|or_die;
  problem.setMaterialsNames(ctx, {{1, "NotchedBeam"}})|or_die;
  problem.setBoundariesNames(ctx,
      {{3, "LowerBoundary"}, {4, "SymmetryAxis"}, {2, "UpperBoundary"}})|or_die;
  problem.addBehaviourIntegrator(ctx, "Mechanics", "NotchedBeam", library,
                                 behaviour)|or_die;
  // materials
  auto& m1 = problem.getMaterial(ctx, "NotchedBeam", 0)|or_die;
  mgis::behaviour::setExternalStateVariable(ctx, m1.s0, "Temperature", 293.15)|or_die;
  mgis::behaviour::setExternalStateVariable(ctx, m1.s1, "Temperature", 293.15)|or_die;
  // boundary conditions
  problem.addUniformDirichletBoundaryCondition(ctx,
      {{"Boundary", "LowerBoundary"}, {"Component", 1}})|or_die;
  problem.addUniformDirichletBoundaryCondition(ctx,
      {{"Boundary", "SymmetryAxis"}, {"Component", 0}})|or_die;
  problem.addUniformDirichletBoundaryCondition(ctx,
      {{"Boundary", "UpperBoundary"},
       {"Component", 1},
       {"LoadingEvolution", [](const auto t) {
          const auto u = 6e-3 * t;
          return u;
        }}})|or_die;
  // solving the problem
  if (!mfem_mgis::usePETSc()) {
    problem.setSolverParameters(ctx,{{"VerbosityLevel", 0},
                                 {"RelativeTolerance", 1e-6},
                                 {"AbsoluteTolerance", 0.},
                                 {"MaximumNumberOfIterations", 10}})|or_die;
    // selection of the linear solver
    if (parallel) {
      problem.setLinearSolver(ctx, "MUMPSSolver", {})|or_die;
    } else {
      problem.setLinearSolver(ctx, "UMFPackSolver", {})|or_die;
    }
  }
  // post-processings
  problem.addPostProcessing(ctx, "ComputeResultantForceOnBoundary",
                            {{"Boundary", 2}, {"OutputFileName", "force.txt"}})|or_die;

  problem.addPostProcessing(ctx, "ParaviewExportResults",
                            {{"OutputFileName", "ssna303-displacements"}})|or_die;
  problem.addPostProcessing(ctx,
      "ParaviewExportIntegrationPointResultsAtNodes",
      {{{"Results", "Stress"}, {"OutputFileName", "ssna303-stress"}}})|or_die;
  problem.addPostProcessing(ctx,
      "ParaviewExportIntegrationPointResultsAtNodes",
      {{{"Results", "EquivalentPlasticStrain"},
        {"OutputFileName", "ssna303-equivalent-plastic-strain"}}})|or_die;
  // loop over time step
  const auto nsteps = mfem_mgis::size_type{50};
  const auto dt = mfem_mgis::real{1} / nsteps;
  auto t = mfem_mgis::real{0};
  auto iteration = mfem_mgis::size_type{};
  for (mfem_mgis::size_type i = 0; i != nsteps; ++i) {
    std::cout << "iteration " << iteration << " from " << t << " to " << t + dt
              << '\n';
    // resolution
    auto ct = t;
    auto dt2 = dt;
    auto nsteps = mfem_mgis::size_type{1};
    auto nsubsteps = mfem_mgis::size_type{0};
    while (nsteps != 0) {
      auto converged = problem.solve(ctx, ct, dt2);
      if (converged) {
        --nsteps;
        ct += dt2;
        problem.update(ctx)|or_die;
      } else {
        std::cout << "\nsubstep: " << nsubsteps << '\n';
        nsteps *= 2;
        dt2 /= 2;
        ++nsubsteps;
        problem.revert(ctx) | or_die;
        if (nsubsteps == 10) {
          mfem_mgis::abort("maximum number of substeps");
        }
      }
    }
    problem.executePostProcessings(ctx, t, dt) | or_die;
    t += dt;
    ++iteration;
    std::cout << '\n';
  }
  return EXIT_SUCCESS;
}
