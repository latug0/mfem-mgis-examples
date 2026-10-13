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

//! \brief vertical component of the resultant force at a given time
struct VerticalForce {
  //! \brief time
  mfem_mgis::real t;
  //! \brief vertical component of the resultant force
  mfem_mgis::real fy;
};

/*!
 * \return the vertical component of the resultant force written by the
 * `ComputeResultantForceOnBoundary` post-processing at each time, or an empty
 * vector if the file can't be read
 * \param[in] f: file name
 */
static std::vector<VerticalForce> readVerticalForce(const std::string& f) {
  auto forces = std::vector<VerticalForce>{};
  auto in = std::ifstream(f);
  auto line = std::string{};
  while (std::getline(in, line)) {
    if ((line.empty()) || (line[0] == '#')) {
      continue;
    }
    auto fx = mfem_mgis::real{};
    auto v = VerticalForce{};
    if (!(std::istringstream(line) >> v.t >> fx >> v.fy)) {
      std::cerr << "invalid line '" << line << "' in '" << f << "'\n";
      return {};
    }
    forces.push_back(v);
  }
  return forces;
}  // end of readVerticalForce

/*!
 * \return if the vertical component of the resultant force matches the
 * reference values. The computed times must be the first times of the
 * reference file, so that the beginning of a loading can be compared to the
 * reference values of the complete loading.
 * \param[in] f: file written by the `ComputeResultantForceOnBoundary`
 * post-processing
 * \param[in] r: reference file
 * \param[in] eps: relative tolerance
 */
static bool checkVerticalForce(const std::string& f,
                               const std::string& r,
                               const mfem_mgis::real eps) {
  const auto values = readVerticalForce(f);
  const auto references = readVerticalForce(r);
  if (values.empty()) {
    std::cerr << "no value read in '" << f << "'\n";
    return false;
  }
  if (references.empty()) {
    std::cerr << "no value read in '" << r << "'\n";
    return false;
  }
  if (values.size() > references.size()) {
    std::cerr << "'" << f << "' has more values than '" << r << "'\n";
    return false;
  }
  for (std::size_t i = 0; i != values.size(); ++i) {
    // the times are written with 6 significant digits
    if (std::abs(values[i].t - references[i].t) >
        1e-6 * std::abs(references[i].t)) {
      std::cerr << "invalid time at time step " << i + 1 << " (" << values[i].t
                << " vs " << references[i].t << ")\n";
      return false;
    }
    if (std::abs(values[i].fy - references[i].fy) >
        eps * std::abs(references[i].fy)) {
      std::cerr << "invalid vertical force at time step " << i + 1 << " ("
                << values[i].fy << " vs " << references[i].fy << ")\n";
      return false;
    }
  }
  return true;
}  // end of checkVerticalForce

int main(int argc, char** argv) {
  constexpr const auto dim = mfem_mgis::size_type{2};
  auto ctx = mfem_mgis::Context{};
  auto or_die = ctx.getFatalFailureHandler();
  // Initialize mfem_mgis (it includes a call to MPI_Init)
  mfem_mgis::initialize(argc, argv);

  const char* mesh_file = "ssna303.msh";
  const char* behaviour = "IsotropicLinearHardeningPlasticity";
  const char* library = "src/libBehaviour.so";
  bool use_fbar = false;
  // not null, since mfem::OptionsParser::PrintUsage stops at the first null
  // string
  const char* reference_file = "";
  const char* standard_reference_file = "";
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
  mfem_mgis::declareDefaultOptions(args);
  args.AddOption(&order, "-o", "--order",
                 "Finite element order (polynomial degree).");
  args.AddOption(&nbsteps, "-ns", "--nbsteps", "Number of time steps.");
  args.AddOption(
      &end_time, "-et", "--end-time",
      "End time. The displacement of the upper boundary is 6e-3 * t.");
  args.AddOption(&reference_file, "-rf", "--reference-file",
                 "Reference values of the resultant force on the upper "
                 "boundary, no comparison if empty.");
  args.AddOption(&parallel, "-p", "--parallel", "-no-p", "--no-parallel",
                 "Perform parallel computations.");
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
  auto problem =
      mfem_mgis::construct<mfem_mgis::NonLinearEvolutionProblem>(
          ctx,
          mfem_mgis::Parameters{
              {"MeshFileName", mesh_file},
              {"FiniteElementFamily", "H1"},
              {"FiniteElementOrder", order},
              {"UnknownsSize", dim},
              {"Materials", mfem_mgis::Parameters{{"NotchedBeam", 1}}},
              {"Boundaries", mfem_mgis::Parameters{{"LowerBoundary", 3},
                                                   {"SymmetryAxis", 4},
                                                   {"UpperBoundary", 2}}},
              {"Hypothesis", "PlaneStrain"},
              {"Parallel", parallel}}) |
      or_die;
#ifdef MGIS_HAVE_TFEL
  if (use_fbar) {
    problem.addBehaviourIntegrator(
        ctx, "Mechanics", "NotchedBeam", library, behaviour,
        {{"Regularization",
          mfem_mgis::Parameters{{"FBar", mfem_mgis::Parameters{}}}}}) |
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
  if (!mfem_mgis::usePETSc()) {
    // the default prediction concentrates the increment of the imposed
    // displacement in the elements next to the upper boundary
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
    auto nsteps = mfem_mgis::size_type{1};
    auto nsubsteps = mfem_mgis::size_type{0};
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
  return EXIT_SUCCESS;
}
