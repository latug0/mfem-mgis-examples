/*!
 * \file   ssna303.cxx
 * \brief
 * \author Thomas Helfer
 * \date   14/12/2020
 */

#include <memory>
#include <cstdlib>
#include <iostream>
#include "mfem/general/optparser.hpp"
#include "mfem/linalg/solvers.hpp"
#include "mfem/linalg/hypre.hpp"

#ifdef MFEM_USE_PETSC
#include "mfem/linalg/petsc.hpp"
#endif /* MFEM_USE_PETSC */

#ifdef MFEM_USE_MUMPS
#include "mfem/linalg/mumps.hpp"
#endif /* MFEM_USE_MUMPS */

#include "mfem/fem/datacollection.hpp"
#include "MGIS/Raise.hxx"
#include "MFEMMGIS/Material.hxx"
#include "MFEMMGIS/Profiler.hxx"
#include "MFEMMGIS/UniformDirichletBoundaryCondition.hxx"
#include "MFEMMGIS/NonLinearEvolutionProblem.hxx"
#include "MFEMMGIS/NonLinearEvolutionProblemImplementation.hxx"
#include "MFEMMGIS/LinearSolverFactory.hxx"

#define PRINT_DEBUG (std::cout <<  __FILE__ << ":" <<  __LINE__ << std::endl)

#define USE_PROFILER 1
#ifndef USE_PROFILER
#define USE_PROFILER 0
#endif
#if USE_PROFILER == 1
#define LIB_PROFILER_IMPLEMENTATION
#define LIB_PROFILER_PRINTF MpiPrintf
#include "libProfiler.h"
#else
#define LogProfiler()
#define PROFILER_ENABLE
#define PROFILER_DISABLE
#define PROFILER_START(x)
#define PROFILER_END()
#endif


int main(int argc, char** argv) {
  PROFILER_ENABLE;
  mfem_mgis::initialize(argc, argv);
  PROFILER_START(0_total);
  PROFILER_START(1_initialize);
  constexpr const auto dim = mfem_mgis::size_type{3};
  const char* mesh_file = "ssna303_3d.msh";
  const char* behaviour = "Plasticity";
  const char* library = "src/libBehaviour.so";
  const char* petscrc_file = "rc_ex10p";
  int parallel = 1;
  int ref_para = 2;
  int ref_seq = 0;
  int order = 1;
  
  //file creation 
  std::string const myFile("./data.txt");
  std::ofstream out(myFile.c_str());

  // options treatment
  mfem::OptionsParser args(argc, argv);
  mfem_mgis::declareDefaultOptions(args);// PETSc Initialize 
  if (!mfem_mgis::usePETSc())  
    mfem_mgis::setPETSc(petscrc_file);
  args.AddOption(&order, "-o", "--order",
                 "Finite element order (polynomial degree).");
  args.AddOption(&parallel, "-p", "--parallel",
                 "Perform parallel computations.");
  args.AddOption(&order, "-o", "--order",
                 "Finite element order (polynomial degree).");
  args.AddOption(&ref_para,"-rp", "--refinement_parallel",
                 "Number of Refinement for parallel call.");
  args.AddOption(&ref_seq,"-rs", "--refinement_sequential",
                 "Number of Refinement for sequential call.");
  args.Parse();
  if (!args.Good()) {
    args.PrintUsage(mfem::out);
    return EXIT_FAILURE;
  }
  args.PrintOptions(mfem::out);

  {
    // loading the mesh and timer
    //  const auto main_timer = mfem_mgis::getTimer("main_timer");
    mfem_mgis::NonLinearEvolutionProblem problem(
                                                 {{"MeshFileName", mesh_file},
                                                     {"FiniteElementFamily", "H1"},
                                                       {"FiniteElementOrder", order},
                                                         {"UnknownsSize", dim},
                                                           {"Hypothesis", "Tridimensional"},
                                                             {"Parallel", parallel !=0},
                                                               {"NumberOfUniformRefinements", parallel ? ref_para : ref_seq}});

    auto mesh = problem.getImplementation<true>().getFiniteElementSpace().GetMesh();
    //get the number of vertices
    int numbers_of_vertices = mesh->GetNV();
    //get the number of elements
    int numbers_of_elements = mesh->GetNE();
    //get the element size
    double h = mesh->GetElementSize(0);
    PROFILER_END(); PROFILER_START(2_mfem_mgis_settings);

    // 2 1 "Volume"
    problem.addBehaviourIntegrator("Mechanics", 1, library, behaviour);
    // materials
    auto& m1 = problem.getMaterial(1);
    mgis::behaviour::setExternalStateVariable(m1.s0, "Temperature", 293.15);
    mgis::behaviour::setExternalStateVariable(m1.s1, "Temperature", 293.15);
    // boundary conditions

    // 3 LowerBoundary
    problem.addBoundaryCondition(
                                 std::make_unique<mfem_mgis::UniformDirichletBoundaryCondition>(
                                                                                                problem.getFiniteElementDiscretizationPointer(), 3, 1));
    // 4 SymmetryPlane1
    problem.addBoundaryCondition(
                                 std::make_unique<mfem_mgis::UniformDirichletBoundaryCondition>(
                                                                                                problem.getFiniteElementDiscretizationPointer(), 4, 0));
    // 5 SymmetryPlane2
    problem.addBoundaryCondition(
                                 std::make_unique<mfem_mgis::UniformDirichletBoundaryCondition>(
                                                                                                problem.getFiniteElementDiscretizationPointer(), 5, 2));
    // 2 UpperBoundary 
    problem.addBoundaryCondition(
                                 std::make_unique<mfem_mgis::UniformDirichletBoundaryCondition>(
                                                                                                problem.getFiniteElementDiscretizationPointer(), 2, 1,
                                                                                                [](const auto t) {
                                                                                                  const auto u = 6e-3 * t;
                                                                                                  return u;
                                                                                                }));

    // solving the problem without petsc
    if (!mfem_mgis::usePETSc()) {
      problem.setSolverParameters({{"VerbosityLevel", 0},
            {"RelativeTolerance", 1e-8},
              {"AbsoluteTolerance", 0.},
                {"MaximumNumberOfIterations", 20}});
      if (parallel) {
        mfem::out << "MUMPS" << std::endl;
        problem.setLinearSolver("MUMPSSolver", {});
      } else {
        mfem::out << "UMFSolver" << std::endl;
        problem.setLinearSolver("UMFPackSolver", {});
      }
    }

    // print on file
    out << " USE_PETSc = " << mfem_mgis::usePETSc() << std::endl;
    out << " taille_maille = " << h << std::endl;
    out << " 1/h = " << 1/h << std::endl;
    out << " numbers_of_vertices = " << numbers_of_vertices << std::endl;
    out << " numbers_of_elements = " << numbers_of_elements << std::endl;

    // vtk export
    //  problem.addPostProcessing("ParaviewExportResults",
    //                            {{"OutputFileName", std::string("ssna303-displacements")}});
    problem.addPostProcessing("ComputeResultantForceOnBoundary",
                              {{"Boundary", 2}, {"OutputFileName", "force.txt"}});
    PROFILER_END(); PROFILER_START(3_time_loop);
  
    // loop over time step
    //  const auto nsteps = mfem_mgis::size_type{50};
    //  const auto dt = mfem_mgis::real{1} / nsteps;
    const auto nsteps = 8;
    const auto dt = mfem_mgis::real{4} / 50;
    auto t = mfem_mgis::real{0};
    auto iteration = mfem_mgis::size_type{};
    for (mfem_mgis::size_type i = 0; i != nsteps; ++i) {
      mfem::out << "iteration " << iteration << " from " << t << " to " << t + dt
                << '\n';
      // resolution
      PROFILER_START(3.1_solve);
      auto ct = t;
      auto dt2 = dt;
      auto nsteps = mfem_mgis::size_type{1};
      auto niter  = mfem_mgis::size_type{0};
      while (nsteps != 0) {
        bool converged = problem.solve(ct, dt2);
        if (converged) {
          --nsteps;
          ct += dt2;
          problem.update();
        } else {
          mfem::out << "\nsubstep: " << niter << '\n';
          nsteps *= 2;
          dt2 /= 2;
          ++niter;
          problem.revert();
          if (niter == 10) {
            mgis::raise("maximum number of substeps");
          }
        }
      }
      PROFILER_END(); PROFILER_START(3.2_postproc);
      problem.executePostProcessings(t, dt);
      PROFILER_END(); 
      t += dt;
      ++iteration;
      mfem::out << '\n';
    }
  }
  PROFILER_END(); 
  PROFILER_END(); 
  if (mfem_mgis::getMPIrank() == 0) 
    LogProfiler();
  PROFILER_DISABLE;
  //  mfem_mgis::Profiler::getProfiler().print(out);
  return EXIT_SUCCESS;
}
