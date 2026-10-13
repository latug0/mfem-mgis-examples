#pragma once

#include <string>
#include <iostream>
#include <cstdlib>
#include <cstdint>
#include <vector>
#include <cassert>
#include <algorithm>

#include <sys/resource.h>
#include <mpi.h>

#include "mfem.hpp"

#include "MFEMMGIS/Parameter.hxx"

#include "Setup.hxx"

inline bool contains(const std::vector<std::string>& list,
                     const std::string& value) {
  return std::find(list.begin(), list.end(), value) != list.end();
}

// UMFPackSolver is not listed since it is only available for sequential
// computations
inline const std::vector<std::string> direct_solvers = {"MUMPSSolver"};

inline const std::vector<std::string> iterative_solvers = {
    "CGSolver",  "GMRESSolver", "BiCGSTABSolver", "MINRESSolver",
    "SLISolver", "HyprePCG",    "HypreGMRES",     "HypreFGMRES"};

inline const std::vector<std::string> preconditioners = {
    "HypreBoomerAMG", "HypreDiagScale", "HypreEuclid", "HypreILU",
    "HypreParaSails"};

inline long get_memory_checkpoint() {
  rusage obj;
  [[maybe_unused]] const auto r = getrusage(RUSAGE_SELF, &obj);
  assert((r == 0) && "error: getrusage has failed");
  long res = 0;
  MPI_Reduce(&(obj.ru_maxrss), &(res), 1, MPI_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
  return res;
}

inline void print_memory_footprint(std::string msg) {
  long mem = get_memory_checkpoint();
  double m = double(mem) * 1e-6;
  mfem_mgis::Profiler::Utils::Message(msg, " memory footprint: ", m, " GB");
}

template <typename Implementation>
inline void print_mesh_information(Implementation& impl) {
  using mfem_mgis::Profiler::Utils::Message;
  using mfem_mgis::Profiler::Utils::sum;
  Message("INFO: print_mesh_information");
  auto mesh = impl.getFiniteElementSpace().GetMesh();
  int64_t numbers_of_elements_local = mesh->GetNE();
  int64_t numbers_of_elements = sum(numbers_of_elements_local);
  auto& fespace = impl.getFiniteElementSpace();
  int64_t unknowns_local = fespace.GetTrueVSize();
  int64_t unknowns = sum(unknowns_local);
  Message("INFO: number of elements -> ", numbers_of_elements);
  Message("INFO: Number of finite element unknowns: ", unknowns);
}

template <typename Problem>
inline void add_post_processings(mfem_mgis::attributes::MayAbort,
                                 mfem_mgis::Context& ctx,
                                 Problem& p,
                                 std::string msg,
                                 std::string field_name) {
  auto or_die = ctx.getFatalFailureHandler();
  p.addPostProcessing(
      ctx, "ParaviewExportResults",
      {{"OutputFileName", msg}, {"OutputFieldName", field_name}}) |
      or_die;
}

template <typename Problem>
inline static void setLinearSolver(mfem_mgis::attributes::MayAbort,
                                   mfem_mgis::Context& ctx,
                                   Problem& p,
                                   const std::string& physics_type,
                                   const TestParameters& param,
                                   const int verbosity = 0,
                                   const mfem_mgis::real Tol = 1e-9) {
  CatchTimeSection(ctx, "set_linear_solver");
  auto or_die = ctx.getFatalFailureHandler();
  std::string solver;
  std::string precond;

  if (physics_type == "heat_transfer") {
    solver = param.solver_thermo;
    precond = param.precond_thermo;
  } else if (physics_type == "mechanics") {
    solver = param.solver_meca;
    precond = param.precond_meca;
  } else {
    mfem_mgis::abort("unknown physics (" + physics_type + ")");
  }

  if (contains(iterative_solvers, solver)) {
    auto solverParameters = mfem_mgis::Parameters{
        {"VerbosityLevel", verbosity}, {"MaximumNumberOfIterations", 10000}};
    if (solver.starts_with("Hypre")) {
      // relative tolerance of the Hypre solvers
      solverParameters.insert(ctx, "Tolerance", Tol) | or_die;
    } else {
      // the tolerances of the MFEM solvers are zero by default
      solverParameters.insert(ctx, "RelativeTolerance", Tol) | or_die;
      solverParameters.insert(ctx, "AbsoluteTolerance", 0.) | or_die;
    }

    if (!precond.empty()) {
      if (!contains(preconditioners, precond)) {
        mfem_mgis::abort("Invalid preconditioner: " + precond);
      }

      auto options = mfem_mgis::Parameters{{"VerbosityLevel", verbosity}};

      auto preconditioner =
          mfem_mgis::Parameters{{"Name", precond}, {"Options", options}};

      solverParameters.insert(ctx, "Preconditioner", preconditioner) | or_die;
    }

    p.setLinearSolver(ctx, solver, solverParameters) | or_die;
  }

  else if (contains(direct_solvers, solver)) {
    p.setLinearSolver(ctx, solver, mfem_mgis::Parameters{}) | or_die;
  }

  else {
    mfem_mgis::abort("Unknown solver type: " + solver);
  }
}
