/*!
 * \file   ssna303_mumps.cxx
 * \brief
 * \author Thomas Helfer
 * \date   14/12/2020
 */

#include <memory>
#include <string_view>
#include <cstdlib>
#include <iostream>
#include "mfem/general/optparser.hpp"
#include "mfem/linalg/solvers.hpp"
#include "mfem/linalg/hypre.hpp"

#ifdef MFEM_USE_PETSC
#include "mfem/linalg/petsc.hpp"
#endif /* MFEM_USE_PETSC */

#ifdef MFEM_USE_MUMPS
#include "mfem/linalg/mumps.hpp"
#endif /* MFEM_USE_MUMPS */

#include "mfem/fem/datacollection.hpp"
#include "MGIS/Raise.hxx"
#include "MFEMMGIS/Material.hxx"
#include "MFEMMGIS/Profiler.hxx"
#include "MFEMMGIS/UniformDirichletBoundaryCondition.hxx"
#include "MFEMMGIS/NonLinearEvolutionProblem.hxx"
#include "MFEMMGIS/NonLinearEvolutionProblemImplementation.hxx"
#include "MFEMMGIS/LinearSolverFactory.hxx"
#include "CheckResultantForce.hxx"

int main(int argc, char** argv) {
  auto ctx = mgis::Context{};
  auto or_die = ctx.getFatalFailureHandler();
  // ctx.enableProfiling(true);
  mfem_mgis::initialize(argc, argv);
  constexpr const auto dim = mfem_mgis::size_type{3};
  const char* mesh_file = "ssna303_3d.msh";
  const char* behaviour = "Plasticity";
  const char* library = "src/libBehaviour.so";
  // not null, since mfem::OptionsParser::PrintUsage stops at the first null
  // string
  const char* reference_file = "";
  auto parallel = int{1};
  auto order = 1;
  auto nbsteps = 50;
  auto end_time = mfem_mgis::real{1};

  // options treatment
  mfem::OptionsParser args(argc, argv);
  mfem_mgis::declareDefaultOptions(args);
  args.AddOption(&parallel, "-p", "--parallel",
                 "Perform parallel computations.");
  args.AddOption(&order, "-o", "--order",
                 "Finite element order (polynomial degree).");
  args.AddOption(&nbsteps, "-ns", "--nbsteps", "Number of time steps.");
  args.AddOption(
      &end_time, "-et", "--end-time",
      "End time. The displacement of the upper boundary is 6e-3 * t.");
  args.AddOption(&reference_file, "-r", "--reference-file",
                 "Reference values of the resultant force on the upper "
                 "boundary, no comparison if empty.");
  args.Parse();
  if (args.Help()) {
    args.PrintUsage(std::cout);
    return EXIT_SUCCESS;
  }
  if (!args.Good()) {
    args.PrintUsage(std::cout);
    return EXIT_FAILURE;
  }
  args.PrintOptions(std::cout);

  // loading the mesh
  auto problem =
      mfem_mgis::construct<mfem_mgis::NonLinearEvolutionProblem>(
          ctx, mfem_mgis::Parameters{{"MeshFileName", mesh_file},
                                     {"FiniteElementFamily", "H1"},
                                     {"FiniteElementOrder", order},
                                     {"UnknownsSize", dim},
                                     {"Hypothesis", "Tridimensional"},
                                     {"Parallel", bool(parallel)}}) |
      or_die;

  // 2 1 "Volume"
  problem.addBehaviourIntegrator(ctx, "Mechanics", 1, library, behaviour) |
      or_die;
  // materials
  auto& m1 = problem.getMaterial(ctx, 1, 0) | or_die;
  mgis::behaviour::setExternalStateVariable(ctx, m1.s0, "Temperature", 293.15) |
      or_die;
  mgis::behaviour::setExternalStateVariable(ctx, m1.s1, "Temperature", 293.15) |
      or_die;
  // boundary conditions

  // 3 LowerBoundary
  problem.addBoundaryCondition(
      ctx, make_unique<mfem_mgis::UniformDirichletBoundaryCondition>(
               ctx, problem.getFiniteElementDiscretizationPointer(), 3, 1) |
               or_die) |
      or_die;
  // 4 SymmetryPlane1
  problem.addBoundaryCondition(
      ctx, make_unique<mfem_mgis::UniformDirichletBoundaryCondition>(
               ctx, problem.getFiniteElementDiscretizationPointer(), 4, 0) |
               or_die) |
      or_die;
  // 5 SymmetryPlane2
  problem.addBoundaryCondition(
      ctx, make_unique<mfem_mgis::UniformDirichletBoundaryCondition>(
               ctx, problem.getFiniteElementDiscretizationPointer(), 5, 2) |
               or_die) |
      or_die;
  // 2 UpperBoundary
  problem.addBoundaryCondition(
      ctx, make_unique<mfem_mgis::UniformDirichletBoundaryCondition>(
               ctx, problem.getFiniteElementDiscretizationPointer(), 2, 1,
               [](const auto t) {
                 const auto u = 6e-3 * t;
                 return u;
               }) |
               or_die) |
      or_die;

  // the default prediction concentrates the increment of the imposed
  // displacement in the elements next to the upper boundary
  problem.setPredictionPolicy(
      {.strategy =
           mfem_mgis::PredictionStrategy::BEGINNING_OF_TIME_STEP_PREDICTION});

  // solving the problem without petsc
  if (!mfem_mgis::usePETSc()) {
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

  // vtk export
  problem.addPostProcessing(
      ctx, "ParaviewExportResults",
      {{"OutputFileName", std::string("ssna303-displacements")}}) |
      or_die;
  problem.addPostProcessing(
      ctx, "ComputeResultantForceOnBoundary",
      {{"Boundary", 2}, {"OutputFileName", "force.txt"}}) |
      or_die;

  // loop over time step
  const auto nsteps = mfem_mgis::size_type(nbsteps);
  const auto dt = end_time / nsteps;
  auto t = mfem_mgis::real{0};
  auto iteration = mfem_mgis::size_type{};
  for (mfem_mgis::size_type i = 0; i != nsteps; ++i) {
    std::cout << "iteration " << iteration << " from " << t << " to " << t + dt
              << '\n';
    // resolution
    auto ct = t;
    auto dt2 = dt;
    auto nsteps = mfem_mgis::size_type{1};
    auto niter = mfem_mgis::size_type{0};
    while (nsteps != 0) {
      bool converged = problem.solve(ctx, ct, dt2);
      if (converged) {
        --nsteps;
        ct += dt2;
        if (nsteps == 0) {
          problem.executePostProcessings(ctx, t, dt) | or_die;
        }
        problem.update(ctx) | or_die;
      } else {
        std::cout << "\nsubstep: " << niter << '\n';
        nsteps *= 2;
        dt2 /= 2;
        ++niter;
        problem.revert(ctx) | or_die;
        if (niter == 10) {
          mgis::abort("maximum number of substeps");
        }
      }
    }
    t += dt;
    ++iteration;
    std::cout << '\n';
  }
  // comparison to the reference values, only on the process writing the
  // resultant force
  if ((!std::string_view{reference_file}.empty()) &&
      (mfem_mgis::isMainProcess(problem.getFiniteElementDiscretization()))) {
    if (!checkVerticalForce("force.txt", reference_file)) {
      return EXIT_FAILURE;
    }
  }
  // mfem_mgis::Profiler::OutputManager::printTimeTable(ctx);
  return EXIT_SUCCESS;
}
