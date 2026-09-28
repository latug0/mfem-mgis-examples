/*!
 * \file   RveNonLinearElastic.cxx
 * \brief
 * This example models a periodic Representative Volume Element (RVE) made of
 * two materials described by the Saint Venant-Kirchhoff hyperelastic
 * behaviour, under an imposed macroscopic deformation gradient.
 *
 * The default mesh, cube_2mat_per.mesh, is a unit cube made of two layers
 * split at x = 0.5. In this case, the solution is compared to the analytical
 * solution (see checkSolution). A RVE with spherical inclusions can be meshed
 * from inclusions_49.geo and computed without this comparison:
 *
 *   gmsh -3 inclusions_49.geo
 *   ./rve --mesh inclusions_49.msh --no-check
 */

#include <array>
#include <cmath>
#include <cstdlib>
#include <algorithm>
#include "mfem/general/optparser.hpp"
#include "MFEMMGIS/Profiler.hxx"
#include "MFEMMGIS/Material.hxx"
#include "MFEMMGIS/PartialQuadratureSpace.hxx"
#include "MFEMMGIS/PeriodicNonLinearEvolutionProblem.hxx"

// Young moduli of the two materials
constexpr auto young1 = mfem_mgis::real{2e11};
constexpr auto young2 = mfem_mgis::real{8e11};
// Poisson ratio of the two materials
constexpr auto poisson = mfem_mgis::real{0.3};
// imposed macroscopic value of Fxx, the other components of the deformation
// gradient being those of the identity
constexpr auto Fxx = mfem_mgis::real{1.1};

// command line options
struct TestParameters {
  const char* mesh_file = "cube_2mat_per.mesh";
  const char* behaviour = "SaintVenantKirchhoffElasticity";
  const char* library = "src/libBehaviour.so";
  int order = 1;
  int refinement = 0;
  bool post_processing = true;
  bool check = true;
  int verbosity_level = 1;
};

void common_parameters(mfem::OptionsParser& args, TestParameters& p) {
  args.AddOption(&p.mesh_file, "-m", "--mesh", "Mesh file to use.");
  args.AddOption(&p.library, "-l", "--library", "Material library.");
  args.AddOption(&p.order, "-o", "--order",
                 "Finite element order (polynomial degree).");
  args.AddOption(&p.refinement, "-r", "--refinement",
                 "Number of uniform refinements of the mesh.");
  args.AddOption(&p.post_processing, "-pp", "--post-processing", "-no-pp",
                 "--no-post-processing", "Export the results to Paraview.");
  args.AddOption(&p.check, "-c", "--check", "-no-c", "--no-check",
                 "Compare the solution to the analytical solution of the "
                 "two-layer cube, only valid for cube_2mat_per.mesh.");
  args.AddOption(&p.verbosity_level, "-v", "--verbosity-level",
                 "Verbosity level of the solvers.");

  args.Parse();

  if (args.Help()) {
    args.PrintUsage(mfem_mgis::getOutputStream());
    mfem_mgis::finalize();
    std::exit(EXIT_SUCCESS);
  }
  if (!args.Good()) {
    args.PrintUsage(mfem_mgis::getOutputStream());
    mfem_mgis::abort(EXIT_FAILURE);
  }
  args.PrintOptions(mfem_mgis::getOutputStream());
}

template <typename Problem>
void add_post_processings(mfem_mgis::attributes::MayAbort,
                          mfem_mgis::Context& ctx,
                          Problem& p,
                          std::string msg) {
  auto or_die = ctx.getFatalFailureHandler();
  p.addPostProcessing(ctx, "ParaviewExportResults", {{"OutputFileName", msg}}) |
      or_die;
}  // end of add_post_processings

template <typename Problem>
void execute_post_processings(mfem_mgis::attributes::MayAbort,
                              mfem_mgis::Context& ctx,
                              Problem& p,
                              double start,
                              double end) {
  CatchTimeSection(ctx, "common::post_processing_step");
  auto or_die = ctx.getFatalFailureHandler();
  p.executePostProcessings(ctx, start, end) | or_die;
}

void setup_properties(mfem_mgis::attributes::MayAbort,
                      mfem_mgis::Context& ctx,
                      const TestParameters& p,
                      mfem_mgis::PeriodicNonLinearEvolutionProblem& problem) {
  using namespace mgis::behaviour;
  using real = mfem_mgis::real;
  auto or_die = ctx.getFatalFailureHandler();
  CatchTimeSection(ctx, "set_mgis_stuff");
  problem.addBehaviourIntegrator(ctx, "Mechanics", 1, p.library, p.behaviour) |
      or_die;
  problem.addBehaviourIntegrator(ctx, "Mechanics", 2, p.library, p.behaviour) |
      or_die;
  // materials
  auto& m1 = problem.getMaterial(ctx, 1, 0) | or_die;
  auto& m2 = problem.getMaterial(ctx, 2, 0) | or_die;
  auto set_properties = [&ctx, &or_die](auto& m, const double yo,
                                        const double po) {
    setMaterialProperty(ctx, m.s0, "YoungModulus", yo) | or_die;
    setMaterialProperty(ctx, m.s0, "PoissonRatio", po) | or_die;
    setMaterialProperty(ctx, m.s1, "YoungModulus", yo) | or_die;
    setMaterialProperty(ctx, m.s1, "PoissonRatio", po) | or_die;
  };

  set_properties(m1, young1, poisson);
  set_properties(m2, young2, poisson);

  //
  auto set_temperature = [&ctx, &or_die](auto& m) {
    setExternalStateVariable(ctx, m.s0, "Temperature", 293.15) | or_die;
    setExternalStateVariable(ctx, m.s1, "Temperature", 293.15) | or_die;
  };
  set_temperature(m1);
  set_temperature(m2);

  // macroscopic deformation gradient
  std::vector<real> e(9, real{0});
  e[0] = Fxx;
  e[1] = 1.0;
  e[2] = 1.0;
  problem.setMacroscopicGradientsEvolution([e](const double) { return e; });
}

template <typename Problem>
static void setLinearSolver(mfem_mgis::attributes::MayAbort,
                            mfem_mgis::Context& ctx,
                            Problem& p,
                            const int verbosity = 0,
                            const mfem_mgis::real Tol = 1e-12) {
  CatchTimeSection(ctx, "set_linear_solver");
  auto or_die = ctx.getFatalFailureHandler();
  // preconditioner hypreBoomerAMG
  auto options = mfem_mgis::Parameters{{"VerbosityLevel", verbosity}};
  auto preconditioner =
      mfem_mgis::Parameters{{"Name", "HypreBoomerAMG"}, {"Options", options}};
  // solver HyprePCG
  p.setLinearSolver(ctx, "HyprePCG",
                    {{"VerbosityLevel", verbosity},
                     {"MaximumNumberOfIterations", 5000},
                     {"Tolerance", Tol},
                     {"Preconditioner", preconditioner}}) |
      or_die;
}

template <typename Problem>
void run_solve(mfem_mgis::attributes::MayAbort,
               mfem_mgis::Context& ctx,
               Problem& p,
               double start,
               double end) {
  CatchTimeSection(ctx, "Solve");
  auto or_die = ctx.getFatalFailureHandler();
  p.solve(ctx, start, end) | or_die;
}

/*!
 * \brief compare the solution to the analytical solution of the two-layer
 * cube.
 *
 * The two layers are normal to the x-axis and have the same thickness. The
 * deformation gradient is uniform in each layer, F = diag(l_i, 1, 1), and the
 * mean value of the l_i is the imposed value Fxx. The first Piola-Kirchhoff
 * stress Pxx = M_i l_i (l_i^2 - 1) / 2, with
 * M_i = E_i (1 - nu) / ((1 + nu) (1 - 2 nu)), is the same in both layers,
 * which gives l_1 by Newton's method.
 */
[[nodiscard]] static bool checkSolution(
    mfem_mgis::Context& ctx, mfem_mgis::PeriodicNonLinearEvolutionProblem& p) {
  using real = mfem_mgis::real;
  auto or_die = ctx.getFatalFailureHandler();
  const auto M1 = young1 * (1 - poisson) / ((1 + poisson) * (1 - 2 * poisson));
  const auto M2 = young2 * (1 - poisson) / ((1 + poisson) * (1 - 2 * poisson));
  auto l1 = Fxx;
  for (int i = 0; i != 20; ++i) {
    const auto l2 = 2 * Fxx - l1;
    const auto r = M1 * l1 * (l1 * l1 - 1) - M2 * l2 * (l2 * l2 - 1);
    const auto dr = M1 * (3 * l1 * l1 - 1) + M2 * (3 * l2 * l2 - 1);
    l1 -= r / dr;
  }
  const auto l = std::array<real, 2>{l1, 2 * Fxx - l1};
  const auto Pxx = M1 * l1 * (l1 * l1 - 1) / 2;
  // maximum error on the deformation gradient and relative error on Pxx
  auto error = real{};
  for (const auto m : {1, 2}) {
    const auto& material = p.getMaterial(ctx, m, 0) | or_die;
    const auto F =
        mfem_mgis::getGradient(ctx, material, "DeformationGradient") | or_die;
    const auto P = mfem_mgis::getThermodynamicForce(
                       ctx, material, "FirstPiolaKirchhoffStress") |
                   or_die;
    // components Fxx, Fyy, Fzz, Fxy, Fyx, Fxz, Fzx, Fyz, Fzy
    const auto Fe = std::array<real, 9>{l[m - 1], 1, 1, 0, 0, 0, 0, 0, 0};
    const auto n = F.getPartialQuadratureSpace().getNumberOfIntegrationPoints();
    for (mfem_mgis::size_type i = 0; i != n; ++i) {
      const auto Fi = F.getIntegrationPointValues(i);
      for (std::size_t c = 0; c != Fe.size(); ++c) {
        error = std::max(error, std::abs(Fi[c] - Fe[c]));
      }
      const auto Pi = P.getIntegrationPointValues(i);
      error = std::max(error, std::abs(Pi[0] - Pxx) / Pxx);
    }
  }
  MPI_Allreduce(MPI_IN_PLACE, &error, 1, MPI_DOUBLE, MPI_MAX,
                mfem_mgis::getMPICommunicator(p));
  if (error > 1e-10) {
    mfem_mgis::getErrorStream()
        << "the solution does not match the analytical solution (error: "
        << error << ")\n";
    return false;
  }
  return true;
}  // end of checkSolution

int main(int argc, char* argv[]) {
  auto ctx = mfem_mgis::Context{};
  auto or_die = ctx.getFatalFailureHandler();
  ctx.enableProfiling(true);
  // mpi initialization here
  mfem_mgis::initialize(argc, argv);

  // get parameters
  TestParameters p;
  mfem::OptionsParser args(argc, argv);
  common_parameters(args, p);

  // 3D
  constexpr const auto dim = mfem_mgis::size_type{3};

  // creating the finite element workspace
  auto fed =
      mfem_mgis::make_shared<mfem_mgis::FiniteElementDiscretization>(
          ctx,
          mfem_mgis::Parameters{{"MeshFileName", p.mesh_file},
                                {"FiniteElementFamily", "H1"},
                                {"FiniteElementOrder", p.order},
                                {"UnknownsSize", dim},
                                {"NumberOfUniformRefinements", p.refinement},
                                {"Parallel", true}}) |
      or_die;
  auto problem =
      mfem_mgis::construct<mfem_mgis::PeriodicNonLinearEvolutionProblem>(ctx,
                                                                         fed) |
      or_die;

  // set problem
  setup_properties(mfem_mgis::may_abort, ctx, p, problem);
  setLinearSolver(mfem_mgis::may_abort, ctx, problem, p.verbosity_level);
  problem.setSolverParameters(ctx, {{"VerbosityLevel", p.verbosity_level},
                                    {"RelativeTolerance", 1e-12},
                                    {"AbsoluteTolerance", 0.},
                                    {"MaximumNumberOfIterations", 10}}) |
      or_die;

  // add post processings
  if (p.post_processing) {
    add_post_processings(mfem_mgis::may_abort, ctx, problem,
                         "OutputFile-rve-non-linear-elastic");
  }
  // main function here
  run_solve(mfem_mgis::may_abort, ctx, problem, 0, 1);

  if (p.post_processing) {
    execute_post_processings(mfem_mgis::may_abort, ctx, problem, 0, 1);
  }
  const auto success = p.check ? checkSolution(ctx, problem) : true;

  // print and write timetable
  mfem_mgis::Profiler::OutputManager::printTimeTable(ctx);
  return success ? EXIT_SUCCESS : EXIT_FAILURE;
}
