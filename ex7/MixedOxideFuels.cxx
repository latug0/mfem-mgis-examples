#include <cmath>
#include <string>
#include <vector>
#include <cassert>
#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <string_view>
#include <sys/resource.h>

#include "mfem/general/optparser.hpp"
#include "MFEMMGIS/Profiler.hxx"
#include "MFEMMGIS/Material.hxx"
#include "MFEMMGIS/NonLinearEvolutionProblemImplementation.hxx"
#include "MFEMMGIS/PeriodicNonLinearEvolutionProblem.hxx"

/*
Problem : Rve mox2 phases with a viscoplastic behavior law

Parameters :

start time = 0
end time = 5s
number of time step = 40

Imposed strain : eps = a * t, with a = 0.012 s^-1
[ - a / 2 ,         0 ,   0 ]
[ 0       ,   - a / 2 ,   0 ]
[ 0       ,         0 ,   a ]

Solver : HypreGMRES
Preconditioner : HypreBoomerAMG

Behavior law parameters : NortonViscoplasticityWithThreshold
[ parameters       , matrix   , inclusions ]
[ Young Modulus    , 8.182e9  , 2*8.182e9  ];
[ Poisson Ratio    , 0.364    , 0.364      ];
[ Stress Threshold , 100.0e6  , 100.0e12   ];
[ Norton Exponent  , 3.333333 , 3.333333   ];
[ Temperature      , 293.15   , 293.15     ];

Element :

Family H1
Order 2
*/

// command line options
struct TestParameters {
  const char* mesh_file = "mesh/OneSphere.msh";
  const char* behaviour = "NortonViscoplasticityWithThreshold";
  const char* library = "src/libBehaviour.so";
  const char* reference_file = "";
  int order = 2;
  int refinement = 0;
  int nbsteps = 40;
  bool post_processing = true;
  int verbosity_level = 0;
};

void common_parameters(mfem::OptionsParser& args, TestParameters& p) {
  args.AddOption(&p.mesh_file, "-m", "--mesh", "Mesh file to use.");
  args.AddOption(&p.library, "-l", "--library", "Material library.");
  args.AddOption(&p.order, "-o", "--order",
                 "Finite element order (polynomial degree).");
  args.AddOption(&p.refinement, "-r", "--refinement",
                 "Number of uniform refinements of the mesh.");
  args.AddOption(&p.nbsteps, "-ns", "--nbsteps",
                 "Number of time steps, the end time being 5 s.");
  args.AddOption(&p.post_processing, "-pp", "--post-processing", "-no-pp",
                 "--no-post-processing", "Export the results to Paraview.");
  args.AddOption(&p.reference_file, "-rf", "--reference-file",
                 "Reference values of the mean stresses in each material, "
                 "no comparison if empty.");
  args.AddOption(&p.verbosity_level, "-v", "--verbosity-level",
                 "Verbosity level of the linear solvers.");
  // PETSc options, handled by mfem_mgis::initialize
  mfem_mgis::declareDefaultOptions(args);

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

template <typename Implementation>
void print_mesh_information(mfem_mgis::attributes::MayAbort,
                            mfem_mgis::Context& ctx,
                            Implementation& impl) {
  using mfem_mgis::Profiler::Utils::sum;
  mfem::out << "INFO: print_mesh_information\n";

  // getMesh
  auto mesh = impl.getFiniteElementSpace().GetMesh();

  // get the number of elements
  int64_t numbers_of_elements_local = mesh->GetNE();
  int64_t numbers_of_elements = sum(numbers_of_elements_local);

  // get n dofs
  auto& fespace = impl.getFiniteElementSpace();
  int64_t unknowns_local = fespace.GetTrueVSize();
  int64_t unknowns = sum(unknowns_local);

  mfem::out << "INFO: number of elements -> " << numbers_of_elements << '\n'
            << "INFO: Number of finite element unknowns: " << unknowns << '\n';
}

long get_memory_checkpoint() {
  rusage obj;
  [[maybe_unused]] const auto r = getrusage(RUSAGE_SELF, &obj);
  assert((r == 0) && "error: getrusage has failed");
  long res = 0;
  MPI_Reduce(&(obj.ru_maxrss), &(res), 1, MPI_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
  return res;
}

void print_memory_footprint(mfem_mgis::attributes::MayAbort,
                            mfem_mgis::Context& ctx,
                            std::string msg) {
  long mem = get_memory_checkpoint();
  double m = double(mem) * 1e-6;  // conversion kb to Gb
  mfem::out << msg << " memory footprint: " << m << " GB\n";
}

template <typename Problem>
void add_post_processings(mfem_mgis::attributes::MayAbort,
                          mfem_mgis::Context& ctx,
                          Problem& p,
                          const bool paraview) {
  auto or_die = ctx.getFatalFailureHandler();
  if (paraview) {
    p.addPostProcessing(ctx, "ParaviewExportResults",
                        {{"OutputFileName", "OutputFile-mixed-oxide-fuels"}}) |
        or_die;
  }
  p.addPostProcessing(ctx, "MeanThermodynamicForces",
                      {{"OutputFileName", "avgStress"}}) |
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
  // materials: 1 is the matrix, 2 the inclusions
  auto& m1 = problem.getMaterial(ctx, 1, 0) | or_die;
  auto& m2 = problem.getMaterial(ctx, 2, 0) | or_die;
  auto set_properties = [&ctx, &or_die](auto& m, const double yo,
                                        const double po, const double st,
                                        const double no) {
    setMaterialProperty(ctx, m.s0, "YoungModulus", yo) | or_die;
    setMaterialProperty(ctx, m.s0, "PoissonRatio", po) | or_die;
    setMaterialProperty(ctx, m.s0, "StressThreshold", st) | or_die;
    setMaterialProperty(ctx, m.s0, "NortonExponent", no) | or_die;

    setMaterialProperty(ctx, m.s1, "YoungModulus", yo) | or_die;
    setMaterialProperty(ctx, m.s1, "PoissonRatio", po) | or_die;
    setMaterialProperty(ctx, m.s1, "StressThreshold", st) | or_die;
    setMaterialProperty(ctx, m.s1, "NortonExponent", no) | or_die;
  };

  set_properties(m1, 8.182e9, 0.364, 100.0e6, 3.333333);
  set_properties(m2, 2 * 8.182e9, 0.364, 100.0e12, 3.333333);

  //
  auto set_temperature = [&ctx, &or_die](auto& m) {
    setExternalStateVariable(ctx, m.s0, "Temperature", 293.15) | or_die;
    setExternalStateVariable(ctx, m.s1, "Temperature", 293.15) | or_die;
  };
  set_temperature(m1);
  set_temperature(m2);

  // macroscopic strain
  std::vector<real> e(6, real{0});
  const int xx = 0;
  const int yy = 1;
  const int zz = 2;

  /* bar{E} = a * t * (-1/2 E1 x E1 + (-1/2) * E2 x E2 + E3 x E3)*/
  const double a = 0.012;
  e[xx] = -0.5 * a;
  e[yy] = -0.5 * a;
  e[zz] = a;
  problem.setMacroscopicGradientsEvolution([e](const double t) {
    auto ret = e;
    for (auto& it : ret) it *= t;
    return ret;
  });
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
  // solver HypreGMRES
  p.setLinearSolver(ctx, "HypreGMRES",
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
               double dt) {
  CatchTimeSection(ctx, "Solve");
  auto or_die = ctx.getFatalFailureHandler();
  p.solve(ctx, start, dt) | or_die;
}

/*!
 * \return the values written by the `MeanThermodynamicForces`
 * post-processing, one row per time step
 * \param[in] f: file name
 */
static std::vector<std::vector<mfem_mgis::real>> readMeanStresses(
    const std::string& f) {
  auto rows = std::vector<std::vector<mfem_mgis::real>>{};
  auto in = std::ifstream(f);
  auto line = std::string{};
  while (std::getline(in, line)) {
    if ((line.empty()) || (line[0] == '#')) {
      continue;
    }
    auto& row = rows.emplace_back();
    auto is = std::istringstream(line);
    auto v = mfem_mgis::real{};
    while (is >> v) {
      row.push_back(v);
    }
  }
  return rows;
}  // end of readMeanStresses

/*!
 * \return if the mean stresses match the reference values
 * \param[in] f: file written by the `MeanThermodynamicForces` post-processing
 * \param[in] r: reference file
 */
static bool checkMeanStresses(const std::string& f, const std::string& r) {
  // relative tolerance, above the rounding of the values which are written
  // with 6 significant digits
  constexpr auto eps = mfem_mgis::real{1e-4};
  const auto values = readMeanStresses(f);
  const auto references = readMeanStresses(r);
  if (references.empty()) {
    mfem_mgis::getErrorStream() << "no value read in '" << r << "'\n";
    return false;
  }
  if (values.size() != references.size()) {
    mfem_mgis::getErrorStream() << "'" << f << "' and '" << r
                                << "' do not have the same number of rows\n";
    return false;
  }
  for (std::size_t i = 0; i != values.size(); ++i) {
    const auto& v = values[i];
    const auto& ref = references[i];
    // the first column is the time, the stresses are compared to the largest
    // one
    auto s = mfem_mgis::real{};
    for (std::size_t j = 1; j < ref.size(); ++j) {
      s = std::max(s, std::abs(ref[j]));
    }
    auto ok = (!ref.empty()) && (v.size() == ref.size()) &&
              (std::abs(v[0] - ref[0]) <= eps * std::abs(ref[0]));
    for (std::size_t j = 1; ok && (j != ref.size()); ++j) {
      ok = std::abs(v[j] - ref[j]) <= eps * s;
    }
    if (!ok) {
      mfem_mgis::getErrorStream()
          << "invalid mean stresses at row " << i + 1 << " of '" << f << "'\n";
      return false;
    }
  }
  return true;
}  // end of checkMeanStresses

int main(int argc, char* argv[]) {
  auto ctx = mgis::Context{};
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
                                {"MeshReadMode", "FromScratch"},
                                {"UnknownsSize", dim},
                                {"NumberOfUniformRefinements", p.refinement},
                                {"Parallel", true}}) |
      or_die;
  auto problem =
      mfem_mgis::construct<mfem_mgis::PeriodicNonLinearEvolutionProblem>(ctx,
                                                                         fed) |
      or_die;
  print_mesh_information(mfem_mgis::may_abort, ctx,
                         problem.getImplementation<true>());
  print_memory_footprint(mfem_mgis::may_abort, ctx, "After_problem:");

  // set problem
  setup_properties(mfem_mgis::may_abort, ctx, p, problem);

  if (!mfem_mgis::usePETSc()) {
    setLinearSolver(mfem_mgis::may_abort, ctx, problem, p.verbosity_level);
  }

  problem.setSolverParameters(ctx, {{"VerbosityLevel", 1},
                                    {"RelativeTolerance", 1e-6},
                                    {"AbsoluteTolerance", 0.},
                                    {"MaximumNumberOfIterations", 6}}) |
      or_die;

  // add post processings: the mean stresses in each material are always
  // computed
  add_post_processings(mfem_mgis::may_abort, ctx, problem, p.post_processing);

  // main function here
  const int nStep = p.nbsteps;
  double start = 0;
  double end = 5;
  const double dt = (end - start) / nStep;
  for (int i = 0; i < nStep; i++) {
    mfem::out << "Solving: from " << i * dt << " to " << (i + 1) * dt << '\n';
    run_solve(mfem_mgis::may_abort, ctx, problem, i * dt, dt);
    execute_post_processings(mfem_mgis::may_abort, ctx, problem, i * dt, dt);
    problem.update(ctx) | or_die;
  }

  // print and write timetable
  print_memory_footprint(mfem_mgis::may_abort, ctx, "After Solving:");
  mfem_mgis::Profiler::OutputManager::printTimeTable(ctx);
  // comparison to the reference values, only on the process writing the mean
  // stresses
  if ((!std::string_view{p.reference_file}.empty()) &&
      (mfem_mgis::isMainProcess(problem.getFiniteElementDiscretization()))) {
    if (!checkMeanStresses("avgStress", p.reference_file)) {
      return EXIT_FAILURE;
    }
  }
  return EXIT_SUCCESS;
}
