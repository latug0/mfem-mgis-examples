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
  mfem_mgis::NonLinearEvolutionProblem problem({{"MeshFileName", mesh_file},
                                                {"FiniteElementFamily", "H1"},
                                                {"FiniteElementOrder", order},
                                                {"UnknownsSize", dim},
                                                {"Hypothesis", "PlaneStrain"},
                                                {"Parallel", parallel}});
  problem.setMaterialsNames({{1, "NotchedBeam"}});
  problem.setBoundariesNames(
      {{3, "LowerBoundary"}, {4, "SymmetryAxis"}, {2, "UpperBoundary"}});
  problem.addBehaviourIntegrator("Mechanics", "NotchedBeam", library,
                                 behaviour);
  // materials
  auto& m1 = problem.getMaterial("NotchedBeam");
  mgis::behaviour::setExternalStateVariable(m1.s0, "Temperature", 293.15);
  mgis::behaviour::setExternalStateVariable(m1.s1, "Temperature", 293.15);
  // boundary conditions
  problem.addUniformDirichletBoundaryCondition(
      {{"Boundary", "LowerBoundary"}, {"Component", 1}});
  problem.addUniformDirichletBoundaryCondition(
      {{"Boundary", "SymmetryAxis"}, {"Component", 0}});
  problem.addUniformDirichletBoundaryCondition(
      {{"Boundary", "UpperBoundary"},
       {"Component", 1},
       {"LoadingEvolution", [](const auto t) {
            const auto u = 6e-3 * t;
            return u;
        }}});
  // solving the problem
  if (!mfem_mgis::usePETSc()) {
    problem.setSolverParameters({{"VerbosityLevel", 0},
                                 {"RelativeTolerance", 1e-6},
                                 {"AbsoluteTolerance", 0.},
                                 {"MaximumNumberOfIterations", 10}});
  // selection of the linear solver
    if (parallel) {
      problem.setLinearSolver("MUMPSSolver", {});
    } else {
      problem.setLinearSolver("UMFPackSolver", {});
    }
  }
  // post-processings
  problem.addPostProcessing("ComputeResultantForceOnBoundary",
                            {{"Boundary", 2}, {"OutputFileName", "force.txt"}});

  problem.addPostProcessing("ParaviewExportResults",
                            {{"OutputFileName", "ssna303-displacements"}});
  problem.addPostProcessing("ParaviewExportIntegrationPointResultsAtNodes",
                            {{{"Results", "Stress"},
                              {"OutputFileName", "ssna303-stress"}}});
  problem.addPostProcessing(
      "ParaviewExportIntegrationPointResultsAtNodes",
      {{{"Results", "EquivalentPlasticStrain"},
        {"OutputFileName", "ssna303-equivalent-plastic-strain"}}});
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
      auto converged = problem.solve(ct, dt2);
      if (converged) {
        --nsteps;
        ct += dt2;
        problem.update();
      } else {
        std::cout << "\nsubstep: " << nsubsteps << '\n';
        nsteps *= 2;
        dt2 /= 2;
        ++nsubsteps;
        problem.revert();
        if (nsubsteps == 10) {
          mfem_mgis::raise("maximum number of substeps");
        }
      }
    }
    problem.executePostProcessings(t, dt);
    t += dt;
    ++iteration;
    std::cout << '\n';
  }
  return EXIT_SUCCESS;
}
