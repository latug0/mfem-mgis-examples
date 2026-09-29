/*!
 * \file   ssna303.cxx
 * \brief
 * \author Thomas Helfer, Guillaume Latu
 * \date   06/04/2021
 */

#include <cmath>
#include <memory>
#include <string>
#include <vector>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <iostream>
#include <string_view>
#include "mfem/general/optparser.hpp"
#include "MFEMMGIS/Material.hxx"
#include "MFEMMGIS/NonLinearEvolutionProblem.hxx"

/*!
 * \return the vertical component of the resultant force written by the
 * `ComputeResultantForceOnBoundary` post-processing
 * \param[in] f: file name
 */
static std::vector<mfem_mgis::real> readVerticalForce(const std::string& f) {
  auto fy = std::vector<mfem_mgis::real>{};
  auto in = std::ifstream(f);
  auto line = std::string{};
  while (std::getline(in, line)) {
    if ((line.empty()) || (line[0] == '#')) {
      continue;
    }
    auto t = mfem_mgis::real{};
    auto fx = mfem_mgis::real{};
    auto v = mfem_mgis::real{};
    std::istringstream(line) >> t >> fx >> v;
    fy.push_back(v);
  }
  return fy;
}  // end of readVerticalForce

/*!
 * \return if the vertical component of the resultant force matches the
 * reference values
 * \param[in] f: file written by the `ComputeResultantForceOnBoundary`
 * post-processing
 * \param[in] r: reference file
 */
static bool checkVerticalForce(const std::string& f, const std::string& r) {
  // relative tolerance, above the rounding of the forces which are written
  // with 6 significant digits
  constexpr auto eps = mfem_mgis::real{1e-3};
  const auto values = readVerticalForce(f);
  const auto references = readVerticalForce(r);
  if ((references.empty()) || (values.size() != references.size())) {
    std::cerr << "'" << f << "' and '" << r
              << "' do not have the same number of values\n";
    return false;
  }
  for (std::size_t i = 0; i != values.size(); ++i) {
    if (std::abs(values[i] - references[i]) > eps * std::abs(references[i])) {
      std::cerr << "invalid vertical force at time step " << i + 1 << " ("
                << values[i] << " vs " << references[i] << ")\n";
      return false;
    }
  }
  return true;
}  // end of checkVerticalForce

int main(int argc, char** argv) {
  using namespace mfem_mgis;
  constexpr const auto dim = size_type{2};
  auto ctx = Context{};
  auto or_die = ctx.getFatalFailureHandler();
  // Initialize mfem_mgis (it includes a call to MPI_Init)
  initialize(argc, argv);

  const char* mesh_file = "ssna303.msh";
  const char* behaviour = "Plasticity";
  const char* library = "src/libBehaviour.so";
  bool use_fbar = false;
  // not null, since mfem::OptionsParser::PrintUsage stops at the first null
  // string
  const char* reference_file = "";
#if defined(MFEM_USE_MUMPS) && defined(MFEM_USE_MPI)
  bool parallel = true;
#else
  bool parallel = false;
#endif
  auto order = 1;
  auto nbsteps = 50;
  auto end_time = mfem_mgis::real{1};
  // options treatment
  mfem::OptionsParser args(argc, argv);
  declareDefaultOptions(args);
  args.AddOption(&order, "-o", "--order",
                 "Finite element order (polynomial degree).");
  args.AddOption(&nbsteps, "-ns", "--nbsteps", "Number of time steps.");
  args.AddOption(
      &end_time, "-et", "--end-time",
      "End time. The displacement of the upper boundary is 6e-3 * t.");
  args.AddOption(&nbsteps, "-ns", "--nbsteps", "Number of time steps.");
  args.AddOption(&reference_file, "-rf", "--reference-file",
                 "Reference values of the resultant force on the upper "
                 "boundary, no comparison if empty.");
  args.AddOption(&parallel, "-p", "--parallel", "-no-p", "--no-parallel",
                 "Perform parallel computations.");
#ifdef MGIS_HAVE_TFEL
  args.AddOption(&use_fbar, "", "--use-fbar", "", "--no-use-fbar",
                 "Use Fbar formulation.");
#endif /* MGIS_HAVE_TFEL */
  args.Parse();
  if (args.Help()) {
    args.PrintUsage(mfem_mgis::getOutputStream());
    mfem_mgis::finalize();
    return EXIT_SUCCESS;
  }
  if (!args.Good()) {
    args.PrintUsage(std::cout);
    abort(EXIT_FAILURE);
  }
  args.PrintOptions(mfem_mgis::getOutputStream());
//
  const auto* const output_file = use_fbar ? "force-fbar.txt" : "force.txt";
  // the non linear problem
  auto problem = construct<NonLinearEvolutionProblem>(
      ctx, dict{{"MeshFileName", mesh_file},
                {"FiniteElementFamily", "H1"},
                {"FiniteElementOrder", order},
                {"UnknownsSize", dim},
                {"Materials", dict{{"NotchedBeam", 1}}},
                {"Boundaries", dict{{"LowerBoundary", 3},
				    {"SymmetryAxis", 4},
				    {"UpperBoundary", 2}}},
                {"Hypothesis", "PlaneStrain"},
                {"Parallel", parallel}}) |or_die; 
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
#else /* MGIS_HAVE_TFEL */
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
  problem.addUniformDirichletBoundaryCondition(
      ctx, {{"Boundary", "LowerBoundary"}, {"Component", 1}}) |
      or_die;
  problem.addUniformDirichletBoundaryCondition(
      ctx, {{"Boundary", "SymmetryAxis"}, {"Component", 0}}) |
      or_die;
  problem.addUniformDirichletBoundaryCondition(ctx,
                                               {{"Boundary", "UpperBoundary"},
                                                {"Component", 1},
                                                {"LoadingEvolution",
                                                 [](const auto t) {
                                                   const auto u = 6e-3 * t;
                                                   return u;
                                                 }}}) |
      or_die;
  // solving the problem
  if (!usePETSc()) {
    problem.setPredictionPolicy(
        {.strategy =
             mfem_mgis::PredictionStrategy::BEGINNING_OF_TIME_STEP_PREDICTION});
    problem.setSolverParameters(ctx, {{"VerbosityLevel", 0},
                                      {"RelativeTolerance", 1e-6},
                                      {"AbsoluteTolerance", 0.},
                                      {"MaximumNumberOfIterations", 10}}) |
        or_die;
    // selection of the linear solver
    if (parallel) {
      problem.setLinearSolver(ctx, "MUMPSSolver", {}) | or_die;
    } else {
      problem.setLinearSolver(ctx, "UMFPackSolver", {}) | or_die;
    }
  }
  // post-processings
  problem.addPostProcessing(
      ctx, "ComputeResultantForceOnBoundary",
      {{"Boundary", 2}, {"OutputFileName", output_file}}) |
      or_die;
  problem.addPostProcessing(ctx, "ParaviewExportResults",
                            {{"OutputFileName", "ssna303-displacements"}}) |
      or_die;
  problem.addPostProcessing(ctx, "ParaviewExportIntegrationPointResultsAtNodes",
                            {{{"Results", "FirstPiolaKirchhoffStress"},
                              {"OutputFileName", "ssna303-stress"}}}) |
      or_die;
  problem.addPostProcessing(
      ctx, "ParaviewExportIntegrationPointResultsAtNodes",
      {{{"Results", "EquivalentPlasticStrain"},
        {"OutputFileName", "ssna303-equivalent-plastic-strain"}}}) |
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
    auto nsubsteps = size_type{0};
    while (nsteps != 0) {
      auto converged = problem.solve(ctx, ct, dt2);
      if (converged) {
        --nsteps;
        ct += dt2;
        problem.update(ctx) | or_die;
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
  // comparison to the reference values, only on the process writing the
  // resultant force
  if ((!std::string_view{reference_file}.empty()) &&
      (mfem_mgis::isMainProcess(problem.getFiniteElementDiscretization()))) {
    if (!checkVerticalForce(output_file, reference_file)) {
      return EXIT_FAILURE;
    }
  }
  return EXIT_SUCCESS;
}
