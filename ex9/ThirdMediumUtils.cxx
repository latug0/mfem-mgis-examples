#include "ThirdMediumUtils.hxx"
#include <cstdlib>
#include <iostream>
#include "mfem/general/optparser.hpp"
#include "MFEMMGIS/Simulation.hxx"
#include "MFEMMGIS/LinearSolverFactory.hxx"

#include "MFEMMGIS/Config.hxx"

#include "MFEMMGIS/Profiler.hxx"

#include <sys/time.h>

#ifdef MFEM_MGIS_HAS_ADDITIONAL_CONVERGENCE_CRITERION
#include "MFEMMGIS/AbstractAdditionalConvergenceCriterion.hxx"

namespace {

  // Concrete implementation of AbstractAdditionalConvergenceCriterion to
  // manipulate linear solver tolerance
  struct AdditionalConvergenceCriterionTolerance final
      : mfem_mgis::nonlinear_solver::AbstractAdditionalConvergenceCriterion {
    mfem_mgis::NonLinearEvolutionProblem& problem;
    mfem_mgis::Context& ctx;
    mfem_mgis::real lin_sol_tolerance_0;
    mfem_mgis::real ls_tolerance;
    mfem_mgis::Parameters params;
    thirdmedium_utils::TestParameters& p;

    AdditionalConvergenceCriterionTolerance(
        mfem_mgis::NonLinearEvolutionProblem& prob,
        mfem_mgis::Context& c,
        mfem_mgis::real l_0,
        mfem_mgis::real lst,
        mfem_mgis::Parameters par,
        thirdmedium_utils::TestParameters& tp)
        : problem(prob),
          ctx(c),
          lin_sol_tolerance_0(l_0),
          ls_tolerance(lst),
          params(std::move(par)),
          p(tp) {}

    ~AdditionalConvergenceCriterionTolerance() override = default;

    void updateLinearSolverTolerance() {
      auto or_die = ctx.getFatalFailureHandler();
      if (p.solver > 0 && p.solver <= 8) {
        if (p.solver <= 5) {
          params.replaceOrInsert("AbsoluteTolerance", ls_tolerance);
          params.replaceOrInsert("RelativeTolerance", ls_tolerance);
        } else {
          params.replaceOrInsert("Tolerance", ls_tolerance);
        }
        std::cout << "Updating Tolerance value of " << p.solver_name << " with "
                  << ls_tolerance << "\n";
        problem.setLinearSolver(ctx, p.solver_name, params) | or_die;
      }
    }

    void helper() noexcept override {
      ls_tolerance = lin_sol_tolerance_0;
      updateLinearSolverTolerance();
      reset();
    }

    void reset() noexcept override {
      constexpr double factor = 100.0;
      ls_tolerance = factor * lin_sol_tolerance_0;
      updateLinearSolverTolerance();
    }

    [[nodiscard]] std::optional<bool> check(
        mfem_mgis::Context&, const CheckArguments& carg) noexcept override {
      if (!carg.converged) {
        return false;
      }
      if (ls_tolerance <= lin_sol_tolerance_0) {
        return true;
      }
      constexpr double divisor = 10.0;
      ls_tolerance = (ls_tolerance / divisor <= lin_sol_tolerance_0)
                         ? lin_sol_tolerance_0
                         : ls_tolerance / divisor;
      updateLinearSolverTolerance();
      return false;
    }
  };

  // Concrete implementation of AbstractAdditionalConvergenceCriterion to
  // manipulate material parameters
  struct AdditionalConvergenceCriterionGamma final
      : mfem_mgis::nonlinear_solver::AbstractAdditionalConvergenceCriterion {
    mfem_mgis::NonLinearEvolutionProblem& problem;
    mfem_mgis::Context& ctx;
    mfem_mgis::real bulk_modulus_0;
    mfem_mgis::real shear_modulus_0;
    mfem_mgis::real bulkm;
    mfem_mgis::real shearm;
    std::vector<std::string> materials;

    AdditionalConvergenceCriterionGamma(mfem_mgis::NonLinearEvolutionProblem& p,
                                        mfem_mgis::Context& c,
                                        mfem_mgis::real b_0,
                                        mfem_mgis::real s_0,
                                        mfem_mgis::real bm,
                                        mfem_mgis::real sm,
                                        std::vector<std::string> mats)
        : problem(p),
          ctx(c),
          bulk_modulus_0(b_0),
          shear_modulus_0(s_0),
          bulkm(bm),
          shearm(sm),
          materials(std::move(mats)) {}

    ~AdditionalConvergenceCriterionGamma() override = default;

    void setMaterialProperties() {
      auto or_die = ctx.getFatalFailureHandler();
      std::cout << "Value setMatProp: " << bulkm << " and " << shearm << "\n";
      for (const auto& n : materials) {
        auto& m = problem.getMaterial(ctx, n, 0) | or_die;
        mgis::behaviour::setMaterialProperty(ctx, m.s0, "BulkModulus", bulkm) |
            or_die;
        mgis::behaviour::setMaterialProperty(ctx, m.s1, "BulkModulus", bulkm) |
            or_die;
        mgis::behaviour::setMaterialProperty(ctx, m.s0, "ShearModulus",
                                             shearm) |
            or_die;
        mgis::behaviour::setMaterialProperty(ctx, m.s1, "ShearModulus",
                                             shearm) |
            or_die;
      }
    }

    void helper() noexcept override {
      bulkm = bulk_modulus_0;
      shearm = shear_modulus_0;
      std::cout << "Value helper : " << bulkm << " and " << shearm << "\n";
      setMaterialProperties();
    }

    void reset() noexcept override {
      constexpr int factor = 3000;
      bulkm = factor * bulk_modulus_0;
      shearm = factor * shear_modulus_0;
      std::cout << "Value reset : " << bulkm << " and " << shearm << "\n";
      setMaterialProperties();
    }

    [[nodiscard]] std::optional<bool> check(
        mfem_mgis::Context&, const CheckArguments& carg) noexcept override {
      if (!carg.converged) {
        return false;
      }
      if (bulkm <= bulk_modulus_0) {
        return true;
      }
      constexpr double divisor = 50.0;
      bulkm = (bulkm / divisor <= bulk_modulus_0) ? bulk_modulus_0
                                                  : bulkm / divisor;
      shearm = (shearm / divisor <= shear_modulus_0) ? shear_modulus_0
                                                     : shearm / divisor;
      setMaterialProperties();
      return false;
    }
  };

}  // anonymous namespace
#endif

namespace thirdmedium_utils {

  void setup_additional_convergence_criteria(
      mfem_mgis::NonLinearEvolutionProblem& mechanics,
      mfem_mgis::Context& ctx,
      TestParameters& p,
      mfem_mgis::Parameters& solver_parameters,
      double bulk_modulus,
      double shear_modulus,
      std::vector<std::string> materials) {
#ifdef MFEM_MGIS_HAS_ADDITIONAL_CONVERGENCE_CRITERION
    auto& solv = mechanics.getSolver();

    if (p.variable_gamma != 0) {
      auto cv_crit_gamma =
          std::make_unique<AdditionalConvergenceCriterionGamma>(
              mechanics, ctx, bulk_modulus, shear_modulus, bulk_modulus,
              shear_modulus, std::move(materials));
      solv.addAdditionalConvergenceCheck(std::move(cv_crit_gamma));
    }

    if (p.variable_tolerance != 0) {
      auto cv_crit_tolerance =
          std::make_unique<AdditionalConvergenceCriterionTolerance>(
              mechanics, ctx, p.tolerance, p.tolerance, solver_parameters, p);
      solv.addAdditionalConvergenceCheck(std::move(cv_crit_tolerance));
    }
#else
    if (p.variable_gamma != 0 || p.variable_tolerance != 0) {
      std::cerr << "WARNING: Variable gamma and tolerance features are not "
                   "implemented "
                   "in this version of the library. Proceeding without "
                   "additional convergence criteria.\n";
    }
#endif
  }

  long* get_memory_checkpoint() {
    rusage obj;
    int who = 0;
    [[maybe_unused]] auto test = getrusage(who, &obj);
    assert((test = -1) && "error: getrusage has failed");
    long* res = new long[3];
    MPI_Reduce(&(obj.ru_maxrss), &res[0], 1, MPI_LONG, MPI_SUM, 0,
               MPI_COMM_WORLD);
    MPI_Reduce(&(obj.ru_maxrss), &res[1], 1, MPI_LONG, MPI_MAX, 0,
               MPI_COMM_WORLD);
    MPI_Reduce(&(obj.ru_maxrss), &res[2], 1, MPI_LONG, MPI_MIN, 0,
               MPI_COMM_WORLD);

    return res;
  };

  void print_memory_footprint(std::string msg) {
    long* mem = get_memory_checkpoint();
    int nprocs;
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    double sum = double(mem[0]) * 1e-6;  // conversion kb to Gb
    double max = double(mem[1]) * 1e-6;  // conversion kb to Gb
    double min = double(mem[2]) * 1e-6;  // conversion kb to Gb
    double mean = sum / nprocs;

    mfem_mgis::Profiler::Utils::Message(
        msg, " memory footprint: sum = ", sum, " GB | ", "max = ", max,
        " GB | ", "min = ", min, " GB | ", "mean = ", mean, " GB | ");
  }

  mfem::OptionsParser parse_options(TestParameters& p, int argc, char** argv) {
    mfem::OptionsParser args(argc, argv);
    // mfem_mgis::declareDefaultOptions(args);
    args.AddOption(&p.mesh_file, "-m", "--mesh", "Mesh file to use.");
    args.AddOption(&p.order, "-o", "--order",
                   "Finite element order (polynomial degree).");
    args.AddOption(&p.refinement, "-r", "--refinement",
                   "refinement level of the mesh, default = 0");
    args.AddOption(
        &p.solver, "-s", "--solver",
        "Number of the solver to be used. Default : 0 = MUMPSSolver. Options "
        ":\n- Direct : 0 = MUMPSSolver\n- Iterative : 1 = CGSolver, 2 = "
        "GMRESSolver,3 =  BiCGSTABSolver,4 = MINRESSolver, 5 = SLISolver, 6 = "
        "HyprePCG, 7 = HypreGMRES and 8 = HypreFGMRES (require Hypre "
        "preconditioner)");
    args.AddOption(&p.preconditioner, "-p", "--preconditioner",
                   "Number of the preconditioner to be used. Default : -1 = "
                   "none .Options :\n- 0 = HypreBoomerAMG, 1 = HypreEuclid, 2 "
                   "= HypreILU, 3 = HypreParaSails, 4 = HypreDiagScale");
    args.AddOption(&p.k_dim, "-k", "--krylov-dimension",
                   "Maximum dimension for the iterations of the GMRES type "
                   "linear solvers");
    args.AddOption(
        &p.ilu_level_of_fill, "-ilu", "--ilu-level-of-fill",
        "Level of fill for the ILU preconditioner. Default : 1. Values should "
        "be greater than 0. Higher values require more memory and computation");
    args.AddOption(&p.newton_iterations, "-n", "--newton-iterations",
                   "Maximum number of iterations for the non-linear solver "
                   "(Newton algorithm)");
    args.AddOption(&p.linear_iterations, "-i", "--linear-iterations",
                   "Maximum number of iterations for the linear solver "
                   "(applicable for iterative solvers)");
    args.AddOption(
        &p.gamma, "-g", "--gamma",
        "Coefficient used for the third medium behavior. Default 5e-7.");
    args.AddOption(&p.alpha, "-a", "--alpha",
                   "Coefficient used for the Faltus 2026 regularization "
                   "penalization. Default 1e11.");
    args.AddOption(
        &p.verbosity_level, "-v", "--verbosity",
        "Value that defines the verbosity level of the solvers. Default value "
        "0 (low verbosity). Values : 0 and 1 (higher verbosity)");
    args.AddOption(&p.tolerance, "-t", "--tolerance",
                   "Value used for the tolerance threshold in the linear "
                   "solvers. Default : 1e-8.");
    args.AddOption(
        &p.variable_gamma, "-vg", "--variable-gamma",
        "Flag used to enable or disable the variation of the gamma parameter. "
        "Default : 0 (disabled). Use a non-zero value to enable");
    args.AddOption(&p.variable_tolerance, "-vt", "--variable-tolerance",
                   "Flag used to enable or disable the variation of the linear "
                   "solver tolerance parameter. Default : 0 (disabled). Use a "
                   "non-zero value to enable");
    args.AddOption(&p.output_dir, "-d", "--output-dir",
                   "Output directory for the post-processing paraview files "
                   "and the performance files");
    args.AddOption(&p.library, "-l", "--library",
                   "Path to the behaviour library.");
    args.AddOption(
        &p.post_processing, "-pp", "--post-processing",
        "Enable or disable post-processing (Paraview output). Default : "
        "disabled/0. Values : O to disable, non-zero to enable");
    args.Parse();

    //  add verification that the output dir has a final slash
    if (!p.output_dir.empty() && p.output_dir.back() != '/') {
      p.output_dir += '/';
    }
    // std::cout << std::setprecision(17) << p.tolerance << "\n";
    if (!args.Good()) {
      args.PrintUsage(std::cout);
      mfem_mgis::abort(EXIT_FAILURE);
    }
    return args;
  }

  // Standard stream insertion operator overload for the struct
  std::ostream& operator<<(std::ostream& os, const TestParameters& p) {
    const int label_width = 35;

    os << "===================================================================="
          "==\n"
       << "                           Test Parameters                          "
          "  \n"
       << "===================================================================="
          "==\n"
       << std::left << std::setw(label_width)
       << "Mesh file (-m):" << p.mesh_file << "\n"
       << std::setw(label_width) << "Library:" << p.library << "\n"
       << std::setw(label_width) << "Finite element order (-o):" << p.order
       << "\n"
       << std::setw(label_width) << "Mesh refinement (-r):" << p.refinement
       << "\n"
       << std::setw(label_width) << "Solver (-s):" << p.solver << " ("
       << (p.solver_name ? p.solver_name : "none") << ")\n"
       << std::setw(label_width) << "Preconditioner (-p):" << p.preconditioner
       << " (" << (p.preconditioner_name != "" ? p.preconditioner_name : "none")
       << ")\n"
       << std::setw(label_width) << "KDim (-k):" << p.k_dim << "\n"
       << std::setw(label_width)
       << "ILULevelOfFill (-ilu):" << p.ilu_level_of_fill << "\n"
       << std::setw(label_width)
       << "Linear iterations (-i):" << p.linear_iterations << "\n"
       << std::setw(label_width)
       << "Newton iterations (-n):" << p.newton_iterations << "\n"
       << std::setw(label_width)
       << "Tolerance of the linear solver (-t):" << p.tolerance << "\n"
       << std::setw(label_width) << "Gamma (-g):" << p.gamma << "\n"
       << std::setw(label_width) << "Alpha (-a):" << p.alpha << "\n"
       << std::setw(label_width)
       << "Parallel execution:" << (p.parallel ? "true" : "false") << "\n"
       << std::setw(label_width) << "Post-processing (-p):" << p.post_processing
       << (p.post_processing != 0 ? " (enabled)" : " (disabled)") << "\n"
       << std::setw(label_width) << "Verbosity level (-v):" << p.verbosity_level
       << "\n"
       << std::setw(label_width) << "Variable gamma parameter (-vg): "
       << (p.variable_gamma ? "enabled" : "disabled") << "\n"
       << std::setw(label_width) << "Variable tolerance parameter (-vt): "
       << (p.variable_gamma ? "enabled" : "disabled") << "\n"
       << std::setw(label_width) << "Output directory (-d): " << p.output_dir
       << "\n"
       << "===================================================================="
          "==\n";

    return os;
  }

  mfem_mgis::Parameters set_solver_parameters(TestParameters& p) {
    mfem_mgis::Parameters solverParameters = mfem_mgis::Parameters{};
    // In all cases (except MUMPS), set the verbosity
    if (p.solver != 0)
      solverParameters.insert(
          mfem_mgis::may_throw,
          mfem_mgis::Parameters{{"VerbosityLevel", p.verbosity_level}});
    // Change the default iteration number according to the solver chosen, if
    // option not specified. Values taken from
    // https://thelfer.github.io/mfem-mgis/user_guide/glossary.html#setup-your-preconditioner
    if (p.parallel) {
      switch (p.solver) {
        case 0:
#if defined(MFEM_USE_MUMPS)
          p.solver_name = "MUMPSSolver";
          p.linear_iterations = 0;  // Not needed
#else
          std::cerr
              << "Error : MUMPSSolver requested but MFEM is not using MUMPS\n";
          mfem_mgis::abort(EXIT_FAILURE);
#endif

          break;

        case 1:
          p.solver_name = "CGSolver";
          p.linear_iterations = 5000;
          solverParameters.insert(mfem_mgis::may_throw, "AbsoluteTolerance",
                                  p.tolerance);
          solverParameters.insert(mfem_mgis::may_throw, "RelativeTolerance",
                                  p.tolerance);
          solverParameters.insert(mfem_mgis::may_throw,
                                  "MaximumNumberOfIterations",
                                  p.linear_iterations);

          break;

        case 2:
          p.solver_name = "GMRESSolver";
          p.linear_iterations = 10000;
          solverParameters.insert(mfem_mgis::may_throw, "AbsoluteTolerance",
                                  p.tolerance);
          solverParameters.insert(mfem_mgis::may_throw, "RelativeTolerance",
                                  p.tolerance);
          solverParameters.insert(mfem_mgis::may_throw,
                                  "MaximumNumberOfIterations",
                                  p.linear_iterations);
          solverParameters.insert(mfem_mgis::may_throw, "KDim", p.k_dim);
          // Could also add "KDim" >=1 parameter.
          break;

        case 3:
          p.solver_name = "BiCGSTABSolver";
          p.linear_iterations = 1000;
          solverParameters.insert(mfem_mgis::may_throw, "AbsoluteTolerance",
                                  p.tolerance);
          solverParameters.insert(mfem_mgis::may_throw, "RelativeTolerance",
                                  p.tolerance);
          solverParameters.insert(mfem_mgis::may_throw,
                                  "MaximumNumberOfIterations",
                                  p.linear_iterations);

          break;

        case 4:
          p.solver_name = "MINRESSolver";
          p.linear_iterations = 1000;
          solverParameters.insert(mfem_mgis::may_throw, "AbsoluteTolerance",
                                  p.tolerance);
          solverParameters.insert(mfem_mgis::may_throw, "RelativeTolerance",
                                  p.tolerance);
          solverParameters.insert(mfem_mgis::may_throw,
                                  "MaximumNumberOfIterations",
                                  p.linear_iterations);

          break;
        case 5:
          p.solver_name = "SLISolver";
          p.linear_iterations = 1000;
          solverParameters.insert(mfem_mgis::may_throw, "AbsoluteTolerance",
                                  p.tolerance);
          solverParameters.insert(mfem_mgis::may_throw, "RelativeTolerance",
                                  p.tolerance);
          solverParameters.insert(mfem_mgis::may_throw,
                                  "MaximumNumberOfIterations",
                                  p.linear_iterations);

          break;
        case 6:
          p.solver_name = "HyprePCG";
          p.linear_iterations = 5000;
          solverParameters.insert(mfem_mgis::may_throw, "Tolerance",
                                  p.tolerance);
          solverParameters.insert(mfem_mgis::may_throw,
                                  "MaximumNumberOfIterations",
                                  p.linear_iterations);

          break;
        case 7:
          p.solver_name = "HypreGMRES";
          p.linear_iterations = 10000;
          solverParameters.insert(mfem_mgis::may_throw, "Tolerance",
                                  p.tolerance);
          solverParameters.insert(mfem_mgis::may_throw,
                                  "MaximumNumberOfIterations",
                                  p.linear_iterations);
          solverParameters.insert(mfem_mgis::may_throw, "KDim", p.k_dim);
          // Could also add "KDim" >=1 parameter.
          break;
        case 8:
          p.solver_name = "HypreFGMRES";
          p.linear_iterations = 10000;
          solverParameters.insert(mfem_mgis::may_throw, "Tolerance",
                                  p.tolerance);
          solverParameters.insert(mfem_mgis::may_throw,
                                  "MaximumNumberOfIterations",
                                  p.linear_iterations);
          solverParameters.insert(mfem_mgis::may_throw, "KDim", p.k_dim);
          // Could also add "KDim" >= 1 parameter
          break;
        default:
          std::cerr << "Provided linear solver option not recognised\n";
          mfem_mgis::abort(EXIT_FAILURE);
      }

    } else {  // Default option if not parallel
      p.solver_name = "UMFPackSolver";
      p.linear_iterations = 0;  // Not needed.
    }

    return solverParameters;
  }

  mfem_mgis::Parameters set_preconditioner_parameters(TestParameters& p) {
    mfem_mgis::Parameters options;
    mfem_mgis::Parameters prec = mfem_mgis::Parameters{};
    if (p.parallel) {
      // Change the preconditioner according to the option given. If none given,
      // no preconditioner is used.
      switch (p.preconditioner) {
        case -1:
          p.preconditioner_name = "";  // No preconditioner.
          prec = {};
          break;
        case 0:
          p.preconditioner_name = "HypreBoomerAMG";
          options =
              mfem_mgis::Parameters{{"VerbosityLevel", p.verbosity_level}};

          prec.insert(mfem_mgis::may_throw, "Name", p.preconditioner_name);
          prec.insert(mfem_mgis::may_throw, "Options", options);
          break;

        case 1:
          options =
              mfem_mgis::Parameters{{"VerbosityLevel", p.verbosity_level}};

          p.preconditioner_name = "HypreEuclid";
          prec.insert(mfem_mgis::may_throw, "Name", p.preconditioner_name);
          prec.insert(mfem_mgis::may_throw, "Options", options);
          break;

        case 2:
          options = mfem_mgis::Parameters{
              {"VerbosityLevel", p.verbosity_level},
              {"HypreILULevelOfFill", p.ilu_level_of_fill}};

          p.preconditioner_name = "HypreILU";
          prec.insert(mfem_mgis::may_throw, "Name", p.preconditioner_name);
          prec.insert(mfem_mgis::may_throw, "Options", options);
          break;

        case 3:
          options =
              mfem_mgis::Parameters{{"VerbosityLevel", p.verbosity_level}};

          p.preconditioner_name = "HypreParaSails";
          prec.insert(mfem_mgis::may_throw, "Name", p.preconditioner_name);
          prec.insert(mfem_mgis::may_throw, "Options", options);
          break;

        case 4:
          options =
              mfem_mgis::Parameters{{"VerbosityLevel", p.verbosity_level}};

          p.preconditioner_name = "HypreDiagScale";

          prec.insert(mfem_mgis::may_throw, "Name", p.preconditioner_name);
          prec.insert(mfem_mgis::may_throw, "Options", options);
          break;

        default:
          std::cerr << "Provided preconditioner option not recognised\n";
          mfem_mgis::abort(EXIT_FAILURE);
      }
    }  // end of if (parallel) for preconditionner
    return prec;
  }

  void check_and_save_convergence_info(mfem_mgis::SimulationOutput& out_params,
                                       mfem_mgis::ExitStatus exit_status) {
    // Retrieve the array of all successfully computed time steps
    auto step_outputs = mfem_mgis::get<std::vector<mfem_mgis::Parameter>>(
        mfem_mgis::may_throw, out_params, "TimeStepOutputs");
    if (exit_status != mfem_mgis::ExitStatus::success) {
      // Simulation failed: Find the exact timestep that caused the divergence
      if (!step_outputs.empty()) {
        auto last_step = mfem_mgis::get<mfem_mgis::Parameters>(
            mfem_mgis::may_throw, step_outputs.back());
        double last_time = mfem_mgis::get<mfem_mgis::real>(
            mfem_mgis::may_throw, last_step, "EndOfTimeStep");

        std::cout << "Simulation diverged after successfully computing time: "
                  << last_time << '\n';
      } else {
        std::cout << "Simulation diverged at the very first time step.\n";
      }

    } else {
      // Simulation succeeded: Quantify linear solver convergence
      for (const auto& step_param : step_outputs) {
        auto step = mfem_mgis::get<mfem_mgis::Parameters>(mfem_mgis::may_throw,
                                                          step_param);

        if (step.contains("ComputeNextStateOutput")) {
          // Access the inner Parameters object containing the solver metrics
          auto next_state = mfem_mgis::get<mfem_mgis::Parameters>(
              mfem_mgis::may_throw, step, "ComputeNextStateOutput");

          // Extract metrics using the exact string keys defined in
          // NonLinearResolutionOutput.cxx
          auto iterations = mfem_mgis::get<mfem_mgis::size_type>(
              mfem_mgis::may_throw, next_state, "Iterations");
          auto initial_res = mfem_mgis::get<mfem_mgis::real>(
              mfem_mgis::may_throw, next_state, "InitialResidualNorm");
          auto final_res = mfem_mgis::get<mfem_mgis::real>(
              mfem_mgis::may_throw, next_state, "FinalResidualNorm");

          // Print or store the convergence data
          std::cout << "Iterations: " << iterations
                    << " | Initial Residual: " << initial_res
                    << " | Final Residual: " << final_res << '\n';
        }
      }
    }
  }

}  // namespace thirdmedium_utils
