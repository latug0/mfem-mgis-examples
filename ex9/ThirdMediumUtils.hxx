#pragma once
#include <cstdlib>
#include <fstream>
#include <iostream>
#include "mfem/general/optparser.hpp"
#include "MFEMMGIS/Simulation.hxx"
#include "MFEMMGIS/LinearSolverFactory.hxx"
#include "MFEMMGIS/ParaviewExportIntegrationPointResultsAtNodes.hxx"

#include "MFEMMGIS/Config.hxx"

#include "MFEMMGIS/Profiler.hxx"

#include <sys/time.h>
#include <sys/resource.h>

namespace thirdmedium_utils {

  struct TestParameters {
    std::string mesh_file = "third_medium.msh";
    std::string library = "src/libBehaviour.so";
    int solver = 0;  // Default MUMPS
    char const* solver_name = "MUMPSSolver";
    int preconditioner = -1;  // Default no preconditioner
    std::string preconditioner_name = "";
    int k_dim = 200;
    int ilu_level_of_fill = 1;
    int newton_iterations = 100;  // Non-linear solver iterations
    int linear_iterations = 500;  // Linear solver iterations (default will be
                                  // set according to the solver)
    double tolerance = 1e-8;
    double gamma = 5e-7;
    double alpha = 1e11;
    int order = 1;
    bool parallel = true;
    int refinement = 0;
    int post_processing = 0;     // default value : disabled
    int verbosity_level = 0;     // default value : lower level
    int variable_gamma = 0;      // default value : disabled
    int variable_tolerance = 0;  // default value : disabled
    std::string output_dir = "output/";
  };

  void setup_additional_convergence_criteria(
      mfem_mgis::NonLinearEvolutionProblem& mechanics,
      mfem_mgis::Context& ctx,
      TestParameters& p,
      mfem_mgis::Parameters& solver_parameters,
      double bulk_modulus,
      double shear_modulus,
      std::vector<std::string> materials);

  long* memory_checkpoint();

  void print_memory_footprint(std::string msg);

  mfem::OptionsParser parse_options(TestParameters& p, int argc, char** argv);
  mfem_mgis::Parameters set_solver_parameters(TestParameters& p);
  mfem_mgis::Parameters set_preconditioner_parameters(TestParameters& p);
  std::ostream& operator<<(std::ostream& os, const TestParameters& p);

  void check_and_save_convergence_info(mfem_mgis::SimulationOutput& out_params,
                                       mfem_mgis::ExitStatus exit_status);

  template <typename Implementation>
  void print_mesh_information(Implementation& impl) {
    using mfem_mgis::Profiler::Utils::Message;
    using mfem_mgis::Profiler::Utils::sum;
    Message("INFO: print_mesh_information");
    // getMesh
    auto mesh = impl.getFiniteElementSpace().GetMesh();
    std::cout << std::flush;
    // get the number of vertices
    //  was int64_t
    int numbers_of_vertices_local = mesh->GetNV();
    int64_t numbers_of_vertices = sum(numbers_of_vertices_local);
    // get the number of elements
    int64_t numbers_of_elements_local = mesh->GetNE();
    int64_t numbers_of_elements = sum(numbers_of_elements_local);
    // get the element size
    double h = (mesh->GetNE() > 0) ? mesh->GetElementSize(0) : 0.0;
    // double h = mesh->GetElementSize(0);
    //  get n dofs
    auto& fespace = impl.getFiniteElementSpace();
    int64_t unknowns_local = fespace.GetTrueVSize();
    int64_t unknowns = sum(unknowns_local);

    Message("INFO: number of vertices -> ", numbers_of_vertices);
    Message("INFO: number of elements -> ", numbers_of_elements);
    Message("INFO: element size -> ", h);
    Message("INFO: Number of finite element unknowns: ", unknowns);
  }

  template <typename Implementation>
  void save_unknowns_info(Implementation& impl, const std::string& output_dir) {
    using mfem_mgis::Profiler::Utils::sum;
    // get n dofs
    auto& fespace = impl.getFiniteElementSpace();
    int64_t unknowns_local = fespace.GetTrueVSize();
    int64_t unknowns = sum(unknowns_local);

    if (mfem_mgis::getMPIrank() == 0) {
      std::ofstream file(output_dir + "unknowns.status");
      file << unknowns;
      file.close();
    }
  }

  /*
      template<typename Problem>
      mfem_mgis::ExitStatus run_solve(mfem_mgis::Context& ctx, Problem&
     mechanics, double ts, double te, int nsteps)
      {
        CatchTimeSection("Solve");
        // loop over time steps
        const auto times = mfem_mgis::Simulation::TimesDescription{ts, te,
     nsteps}; mfem_mgis::Parameters params; params.insert(mfem_mgis::may_throw,
     "KeepOutputs", true); params.insert(mfem_mgis::may_throw, "Times", times);
        auto s = mfem_mgis::Simulation{mechanics, times};
        std::pair run_out = s.run(ctx);
        std::cout << ctx.getErrorMessage() << '\n';

        return run_out.get(0); // Gets the exit status of the simulation
      }
  */
  template <typename Problem>
  std::pair<mfem_mgis::ExitStatus, std::optional<mfem_mgis::SimulationOutput>>
  run_solve(mfem_mgis::Context& ctx,
            Problem& mechanics,
            double ts,
            double te,
            int nsteps) {
    CatchTimeSection(ctx, "Solve");
    auto or_die = ctx.getFatalFailureHandler();

    // Setup Parameters to keep simulation outputs
    mfem_mgis::Parameters params;
    params.insert(mfem_mgis::may_throw, "KeepOutputs", true);

    // Build the time sequence manually for the Parameters object
    std::vector<mfem_mgis::Parameter> times_vec;
    double dt = (te - ts) / nsteps;
    for (int i = 0; i <= nsteps; ++i) {
      times_vec.push_back(ts + i * dt);
    }
    params.insert(mfem_mgis::may_throw, "Times", times_vec);

    // Initialize simulation with parameters
    auto s =
        mfem_mgis::construct<mfem_mgis::Simulation>(ctx, mechanics, params) |
        or_die;
    const auto run_out = s.run(ctx);
    if (run_out.first.shallStop()) {
      mfem::out << ctx.getErrorMessage() << '\n';
    }
    // Return the full pair to access both status and output data
    return run_out;
  }

}  // namespace thirdmedium_utils
