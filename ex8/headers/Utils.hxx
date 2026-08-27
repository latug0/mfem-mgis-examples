#pragma once

#include <string>
#include <iostream>
#include <cstdlib>
#include <cstdint>
#include <vector>
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

inline const std::vector<std::string> direct_solvers = {"MUMPSSolver",
                                                        "UMFPackSolver"};

inline const std::vector<std::string> iterative_solvers = {
    "CGSolver",  "GMRESSolver", "BiCGSTABSolver", "MINRESSolver",
    "SLISolver", "HyprePCG",    "HypreGMRES",     "HypreFGMRES"};

inline const std::vector<std::string> preconditionners = {
    "HypreBoomerAMG", "HypreDiagScale", "HypreEuclid", "HypreILU",
    "HypreParaSails"};

inline long get_memory_checkpoint() {
  rusage obj;
  int who = 0;
  [[maybe_unused]] auto test = getrusage(who, &obj);
  long res;
  MPI_Reduce(&(obj.ru_maxrss), &(res), 1, MPI_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
  return res;
};

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
  int64_t numbers_of_vertices_local = mesh->GetNV();
  int64_t numbers_of_vertices = sum(numbers_of_vertices_local);
  int64_t numbers_of_elements_local = mesh->GetNE();
  int64_t numbers_of_elements = sum(numbers_of_elements_local);
  double h = mesh->GetElementSize(0);
  auto& fespace = impl.getFiniteElementSpace();
  int64_t unknowns_local = fespace.GetTrueVSize();
  int64_t unknowns = sum(unknowns_local);
  Message("INFO: number of vertices -> ", numbers_of_vertices);
  Message("INFO: number of elements -> ", numbers_of_elements);
  Message("INFO: element size -> ", h);
  Message("INFO: Number of finite element unknowns: ", unknowns);
}

template <typename Problem>
inline void add_post_processings(Problem& p, std::string msg) {
  p.addPostProcessing(
      "ParaviewExportResults",
      {{"OutputDirectory", "Resultats"}, {"OutputFileName", msg}});
}

template <typename Problem>
inline void execute_post_processings(mgis::Context& ctx,
                                     Problem& p,
                                     double start,
                                     double end) {
  CatchTimeSection(ctx, "common::post_processing_step");
  p.executePostProcessings(start, end);
}

template <typename Problem>
inline static void setLinearSolver(mgis::Context& ctx,
                                   Problem& p,
                                   const std::string& physics_type,
                                   const TestParameters& param,
                                   const int verbosity = 0,
                                   const mfem_mgis::real Tol = 1e-9) {
  CatchTimeSection(ctx, "set_linear_solver");

  std::string solver;
  std::string precond;

  if (physics_type == "heat_transfer") {
    solver = param.solver_thermo;
    precond = param.precond_thermo;
  } else if (physics_type == "mechanics") {
    solver = param.solver_meca;
    precond = param.precond_meca;
  } else {
    std::cerr << "Erreur : Physique inconnue (" << physics_type << ")"
              << std::endl;
    std::abort();
  }

  if (contains(iterative_solvers, solver)) {
    constexpr int defaultMaxNumOfIt = 10e3;

    auto solverParameters = mfem_mgis::Parameters{};
    solverParameters.insert(
        mfem_mgis::Parameters{{"VerbosityLevel", verbosity}});
    solverParameters.insert(mfem_mgis::Parameters{
        {"MaximumNumberOfIterations", defaultMaxNumOfIt}});

    if (solver == "MINRESSolver" || solver == "BiCGSTABSolver" ||
        solver == "CGSolver" || solver == "GMRESSolver") {
      solverParameters.insert(
          mfem_mgis::Parameters{{"AbsoluteTolerance", Tol}});
    } else {
      solverParameters.insert(mfem_mgis::Parameters{{"Tolerance", Tol}});
    }

    if (!precond.empty()) {
      if (!contains(preconditionners, precond)) {
        std::cerr << "Invalid preconditioner: " << precond << std::endl;
        std::abort();
      }

      auto options = mfem_mgis::Parameters{{"VerbosityLevel", verbosity}};

      auto preconditioner =
          mfem_mgis::Parameters{{"Name", precond}, {"Options", options}};

      solverParameters.insert(
          mfem_mgis::Parameters{{"Preconditioner", preconditioner}});
    }

    p.setLinearSolver(solver, solverParameters);
  }

  else {
    std::cerr << "Unknown solver type: " << solver << std::endl;
    std::abort();
  }
}