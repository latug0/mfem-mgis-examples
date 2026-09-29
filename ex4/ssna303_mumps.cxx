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
  using namespace mfem_mgis;
  auto ctx = Context{};
  auto or_die = ctx.getFatalFailureHandler();
  // ctx.enableProfiling(true);
  initialize(argc, argv);
  constexpr const auto dim = size_type{3};
  const char* mesh_file = "ssna303_3d.msh";
  const char* behaviour = "IsotropicLinearHardeningPlasticity";
  const char* library = "src/libBehaviour.so";
  bool use_fbar = false;
  // not null, since mfem::OptionsParser::PrintUsage stops at the first null
  // string
  const char* reference_file = "";
  const char* standard_reference_file = "";
  auto parallel = true;
  auto order = 1;
  auto nbsteps = 50;
  auto end_time = mfem_mgis::real{1};

  // options treatment
  mfem::OptionsParser args(argc, argv);
  declareDefaultOptions(args);
  args.AddOption(&parallel, "-p", "--parallel", "-no-p", "--no-parallel",
                 "Perform parallel computations.");
  args.AddOption(&order, "-o", "--order",
                 "Finite element order (polynomial degree).");
  args.AddOption(&nbsteps, "-ns", "--nbsteps", "Number of time steps.");
  args.AddOption(
      &end_time, "-et", "--end-time",
      "End time. The displacement of the upper boundary is 6e-3 * t.");
  args.AddOption(&reference_file, "-rf", "--reference-file",
                 "Reference values of the resultant force on the upper "
                 "boundary, no comparison if empty.");
#ifdef MGIS_HAVE_TFEL
  args.AddOption(&use_fbar, "-fb", "--use-fbar", "-no-fb", "--no-use-fbar",
                 "Use the FBar formulation.");
  args.AddOption(&standard_reference_file, "-srf", "--standard-reference-file",
                 "Reference values of the resultant force on the upper "
                 "boundary computed without FBar, compared with a larger "
                 "tolerance, no comparison if empty.");
#endif /* MGIS_HAVE_TFEL */
  args.Parse();
  if (args.Help()) {
    args.PrintUsage(mfem_mgis::getOutputStream());
    mfem_mgis::finalize();
    return EXIT_SUCCESS;
  }
  if (!args.Good()) {
    args.PrintUsage(mfem_mgis::getOutputStream());
    mfem_mgis::abort(EXIT_FAILURE);
  }
  args.PrintOptions(mfem_mgis::getOutputStream());
  const auto* const output_file = use_fbar ? "force-fbar.txt" : "force.txt";
  // the non linear problem
  auto problem = construct<NonLinearEvolutionProblem>(
                     ctx, dict{{"MeshFileName", mesh_file},
                               {"Materials", dict{{"NotchedBeam", 1}}},
                               {"FiniteElementFamily", "H1"},
                               {"FiniteElementOrder", order},
                               {"UnknownsSize", dim},
                               {"Hypothesis", "Tridimensional"},
                               {"Parallel", parallel}}) |
                 or_die;

  // 2 1 "Volume"
#ifdef MGIS_HAVE_TFEL
  if (use_fbar) {
    problem.addBehaviourIntegrator(
        ctx, "Mechanics", "NotchedBeam", library, behaviour,
        {{"Regularization", dict{{"FBar", dict{}}}}}) |
        or_die;
  } else {
    problem.addBehaviourIntegrator(ctx, "Mechanics", "NotchedBeam", library,
                                   behaviour) |
        or_die;
  }
#else  /* MGIS_HAVE_TFEL */
  problem.addBehaviourIntegrator(ctx, "Mechanics", "NotchedBeam", library,
                                 behaviour) |
      or_die;
#endif /* MGIS_HAVE_TFEL */
  // materials
  auto& m1 = problem.getMaterial(ctx, "NotchedBeam", 0) | or_die;
  mgis::behaviour::setExternalStateVariable(ctx, m1.s0, "Temperature", 293.15) |
      or_die;
  mgis::behaviour::setExternalStateVariable(ctx, m1.s1, "Temperature", 293.15) |
      or_die;

  // boundary conditions

  // 3 LowerBoundary
  problem.addBoundaryCondition(
      ctx, make_unique<UniformDirichletBoundaryCondition>(
               ctx, problem.getFiniteElementDiscretizationPointer(), 3, 1) |
               or_die) |
      or_die;
  // 4 SymmetryPlane1
  problem.addBoundaryCondition(
      ctx, make_unique<UniformDirichletBoundaryCondition>(
               ctx, problem.getFiniteElementDiscretizationPointer(), 4, 0) |
               or_die) |
      or_die;
  // 5 SymmetryPlane2
  problem.addBoundaryCondition(
      ctx, make_unique<UniformDirichletBoundaryCondition>(
               ctx, problem.getFiniteElementDiscretizationPointer(), 5, 2) |
               or_die) |
      or_die;
  // 2 UpperBoundary
  problem.addBoundaryCondition(
      ctx, make_unique<UniformDirichletBoundaryCondition>(
               ctx, problem.getFiniteElementDiscretizationPointer(), 2, 1,
               [](const auto t) {
                 const auto u = 6e-3 * t;
                 return u;
               }) |
               or_die) |
      or_die;

  // solving the problem without petsc
  if (!usePETSc()) {
    // the default prediction concentrates the increment of the imposed
    // displacement in the elements next to the upper boundary
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

  // vtk export
  problem.addPostProcessing(
      ctx, "ParaviewExportResults",
      {{"OutputFileName", std::string("ssna303-displacements")}}) |
      or_die;
  problem.addPostProcessing(
      ctx, "ComputeResultantForceOnBoundary",
      {{"Boundary", 2}, {"OutputFileName", output_file}}) |
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
    auto nsteps = size_type{1};
    auto niter = size_type{0};
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
  // resultant force. The tolerance is above the rounding of the forces, which
  // are written with 6 significant digits. It is larger for the reference
  // values computed without FBar, which are only close to the results with
  // FBar.
  if (mfem_mgis::isMainProcess(problem.getFiniteElementDiscretization())) {
    if ((!std::string_view{reference_file}.empty()) &&
        (!checkVerticalForce(output_file, reference_file, 1e-4))) {
      return EXIT_FAILURE;
    }
    if ((!std::string_view{standard_reference_file}.empty()) &&
        (!checkVerticalForce(output_file, standard_reference_file, 1e-3))) {
      return EXIT_FAILURE;
    }
  }
  // mfem_mgis::Profiler::OutputManager::printTimeTable(ctx);
  return EXIT_SUCCESS;
}
