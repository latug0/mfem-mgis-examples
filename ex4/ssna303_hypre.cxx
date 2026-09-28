/*!
 * \file   ssna303.cxx
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
#include "mfem/linalg/petsc.hpp"
#include "mfem/fem/datacollection.hpp"
#include "MGIS/Raise.hxx"
#include "MFEMMGIS/Material.hxx"
#include <MFEMMGIS/Profiler.hxx>
#include "MFEMMGIS/UniformDirichletBoundaryCondition.hxx"
#include "MFEMMGIS/NonLinearEvolutionProblem.hxx"
#include "MFEMMGIS/NonLinearEvolutionProblemImplementation.hxx"
#include "MFEMMGIS/LinearSolverFactory.hxx"
#include "CheckResultantForce.hxx"

#define PRINT_DEBUG (std::cout << __FILE__ << ":" << __LINE__ << std::endl)

int main(int argc, char** argv) {
  auto ctx = mgis::Context{};
  auto or_die = ctx.getFatalFailureHandler();
  // ctx.enableProfiling(true);
  mfem_mgis::initialize(argc, argv);
  bool parallel = true;
  constexpr const auto dim = mfem_mgis::size_type{3};
  const char* mesh_file = "ssna303_3d.msh";
  const char* behaviour = "Plasticity";
  const char* library = "src/libBehaviour.so";
  // not null, since mfem::OptionsParser::PrintUsage stops at the first null
  // string
  const char* reference_file = "";
  auto solver = "HypreFGMRES";
  auto preconditioner = "HypreBoomerAMG";  //"";//
  auto ref_para = 0;
  auto ref_seq = 0;
  auto order = 1;
  auto nbsteps = 50;
  auto end_time = mfem_mgis::real{1};

  // file creation
  std::string const myFile("test.txt");
  std::ofstream out(myFile.c_str());

  // options treatment
  mfem::OptionsParser args(argc, argv);
  args.AddOption(&order, "-o", "--order",
                 "Finite element order (polynomial degree).");
  args.AddOption(&nbsteps, "-ns", "--nbsteps", "Number of time steps.");
  args.AddOption(
      &end_time, "-et", "--end-time",
      "End time. The displacement of the upper boundary is 6e-3 * t.");
  args.AddOption(&reference_file, "-rf", "--reference-file",
                 "Reference values of the resultant force on the upper "
                 "boundary, no comparison if empty.");
  args.AddOption(&solver, "-s", "--solver", "Solver of the Problem.");
  args.AddOption(&preconditioner, "-p", "--preconditioner",
                 "Preconditioner for the Problem.");
  args.AddOption(&ref_para, "-rp", "--refinement_parallel",
                 "Number of Refinement for parallel call.");
  args.AddOption(&ref_seq, "-rs", "--refinement_sequential",
                 "Number of Refinement for sequential call.");
  args.Parse();
  if (!args.Good()) {
    args.PrintUsage(std::cout);
    return EXIT_FAILURE;
  }
  args.PrintOptions(std::cout);

  // loading the mesh
  {
    auto problem =
        mfem_mgis::construct<mfem_mgis::NonLinearEvolutionProblem>(
            ctx, mfem_mgis::Parameters{{"MeshFileName", mesh_file},
                                       {"FiniteElementFamily", "H1"},
                                       {"FiniteElementOrder", order},
                                       {"UnknownsSize", dim},
                                       {"NumberOfUniformRefinements",
                                        parallel ? ref_para : ref_seq},
                                       {"Hypothesis", "Tridimensional"},
                                       {"Parallel", true}}) |
        or_die;

    auto mesh =
        problem.getImplementation<true>().getFiniteElementSpace().GetMesh();
    // get the number of vertices
    int numbers_of_vertices = mesh->GetNV();
    // get the number of elements
    int numbers_of_elements = mesh->GetNE();
    // get the element size
    double h = mesh->GetElementSize(0);

    // 2 1 "Volume"
    problem.addBehaviourIntegrator(ctx, "Mechanics", 1, library, behaviour) |
        or_die;
    // materials
    auto& m1 = problem.getMaterial(ctx, 1, 0) | or_die;
    mgis::behaviour::setExternalStateVariable(ctx, m1.s0, "Temperature",
                                              293.15) |
        or_die;
    mgis::behaviour::setExternalStateVariable(ctx, m1.s1, "Temperature",
                                              293.15) |
        or_die;
    // boundary conditions

    // 3 LowerBoundary
    problem.addBoundaryCondition(
        ctx,
        mfem_mgis::make_unique<mfem_mgis::UniformDirichletBoundaryCondition>(
            ctx, problem.getFiniteElementDiscretizationPointer(), 3, 1) |
            or_die) |
        or_die;
    // 4 SymmetryPlane1
    problem.addBoundaryCondition(
        ctx,
        mfem_mgis::make_unique<mfem_mgis::UniformDirichletBoundaryCondition>(
            ctx, problem.getFiniteElementDiscretizationPointer(), 4, 0) |
            or_die) |
        or_die;
    // 5 SymmetryPlane2
    problem.addBoundaryCondition(
        ctx,
        mfem_mgis::make_unique<mfem_mgis::UniformDirichletBoundaryCondition>(
            ctx, problem.getFiniteElementDiscretizationPointer(), 5, 2) |
            or_die) |
        or_die;
    // 2 UpperBoundary
    problem.addBoundaryCondition(
        ctx,
        mfem_mgis::make_unique<mfem_mgis::UniformDirichletBoundaryCondition>(
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

    // solving the problem
    problem.setSolverParameters(ctx, {{"VerbosityLevel", 0},
                                      {"RelativeTolerance", 1e-6},
                                      {"AbsoluteTolerance", 0.},
                                      {"MaximumNumberOfIterations", 20}}) |
        or_die;

    // selection of the linear solver without preconditioner
    if (std::string_view{solver}.empty()) {
      return EXIT_FAILURE;
    }
    if (std::string_view{preconditioner}.empty()) {
      problem.setLinearSolver(ctx, solver,
                              {{"VerbosityLevel", 0},
                               //{"AbsoluteTolerance", 1e-12},
                               //{"KDim", 3},
                               {"Tolerance", 1e-12},
                               {"MaximumNumberOfIterations", 1000}}) |
          or_die;
    } else {
      // with the HypreBoomerAMG preconditioner
      //  auto prec_none = mfem_mgis::Parameters{{"Name", "None"}};
      auto prec_boomer = mfem_mgis::Parameters{
          {"Name", preconditioner},  //
          {"Options",
           mfem_mgis::Parameters{
               //                           {"Strategy", "Elasticity"},
               {"VerbosityLevel", 0}}}};

      // the number of iterations increases with the plastic strain: more
      // than 300 iterations are needed at the end of the loading
      problem.setLinearSolver(ctx, solver,
                              {{"VerbosityLevel", 0},
                               //{"AbsoluteTolerance", 1e-12},
                               //{"RelativeTolerance", 1e-12},
                               //{"Tolerance", 1e-12},
                               {"MaximumNumberOfIterations", 1000},
                               {"Preconditioner", prec_boomer}}) |
          or_die;
    }
    // print on file
    out << " SetLinearSolver" << std::endl;
    out << " VerbosityLevel = " << 0 << std::endl;
    out << " RelativeTolerance = " << 1e-12 << std::endl;
    out << " MaximumNumberOfIterations = " << 1000 << std::endl;
    out << " Preconditioner = " << preconditioner << std::endl;
    out << " taille_maille = " << h << std::endl;
    out << " 1/h = " << 1 / h << std::endl;
    out << " nbr_ref_parallel = " << ref_para << std::endl;
    out << " nbr_ref_sequential = " << ref_seq << std::endl;
    out << " numbers_of_vertices = " << numbers_of_vertices << std::endl;
    out << " numbers_of_elements = " << numbers_of_elements << std::endl;

    // vtk export
    problem.addPostProcessing(
        ctx, "ParaviewExportResults",
        {{"OutputFileName",
          std::string("ssna303-displacements-HFGMRES_WS_1")}}) |
        or_die;
    problem.addPostProcessing(
        ctx, "ComputeResultantForceOnBoundary",
        {{"Boundary", 2}, {"OutputFileName", "force_HFGMRES_WS_1.txt"}}) |
        or_die;

    // loop over time step
    const auto nsteps = mfem_mgis::size_type(nbsteps);
    const auto dt = end_time / nsteps;
    auto t = mfem_mgis::real{0};
    auto iteration = mfem_mgis::size_type{};
    for (mfem_mgis::size_type i = 0; i != nsteps; ++i) {
      using namespace mfem_mgis;
      CatchTimeSection(ctx, "time_loop");
      std::cout << "iteration " << iteration << " from " << t << " to "
                << t + dt << '\n';
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
      if (!checkVerticalForce("force_HFGMRES_WS_1.txt", reference_file)) {
        return EXIT_FAILURE;
      }
    }
  }
  // mfem_mgis::Profiler::OutputManager::printTimeTable(ctx);
  return EXIT_SUCCESS;
}
