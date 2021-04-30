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
#include "mfem/fem/datacollection.hpp"
#include "MGIS/Raise.hxx"
#include "MFEMMGIS/Material.hxx"
#include "MFEMMGIS/Profiler.hxx"
#include "MFEMMGIS/UniformDirichletBoundaryCondition.hxx"
#include "MFEMMGIS/NonLinearEvolutionProblem.hxx"
#include "MFEMMGIS/NonLinearEvolutionProblemImplementation.hxx"
#include "MFEMMGIS/LinearSolverFactory.hxx"

#define PRINT_DEBUG (std::cout <<  __FILE__ << ":" <<  __LINE__ << std::endl)

int main(int argc, char** argv) {
	
  mfem_mgis::initialize(argc, argv);
  bool parallel = true;
  constexpr const auto dim = mfem_mgis::size_type{3};
  const char* mesh_file = "ssna303_3d.msh";
  const char* behaviour = "Plasticity";
  const char* library = "src/libBehaviour.so";
  auto order = 1;

  //file creation 
https://meet.jit.si/cea-cad
  std::string const myFile("/home/hc265945/spack_codes/mfem-mgis/ssna303/ssna303-3D/Test_Ssna303.txt");
  std::ofstream out(myFile.c_str());

  // options treatment
  mfem::OptionsParser args(argc, argv);
  args.AddOption(&order, "-o", "--order",
                 "Finite element order (polynomial degree).");
  args.Parse();
  if (!args.Good()) {
    args.PrintUsage(std::cout);
    return EXIT_FAILURE;
  }
  args.PrintOptions(std::cout);
  // loading the mesh
  int nbr_ref_parallel = 1;
  int nbr_ref_sequential = 1;
  {
  const auto main_timer = mfem_mgis::getTimer("main_timer");
mfem_mgis::NonLinearEvolutionProblem problem(
      {{"MeshFileName", mesh_file},
       {"FiniteElementFamily", "H1"},
       {"FiniteElementOrder", order},
       {"UnknownsSize", dim},
       {"NumberOfUniformRefinements", parallel ? nbr_ref_parallel : nbr_ref_sequential},
       {"Hypothesis", "Tridimensional"},
       {"Parallel", true}});

  auto mesh = problem.getImplementation<true>().getFiniteElementSpace().GetMesh();

//get the number of vertices
  int numbers_of_vertices = mesh->GetNV();
//get the number of elements
  int numbers_of_elements = mesh->GetNE();
//get the element size
  double h;
  h = mesh->GetElementSize(0);
  double t_m = 1/h;

/*  mfem_mgis::NonLinearEvolutionProblem problem(
      {{"MeshFileName", mesh_file},
       {"FiniteElementFamily", "H1"},
       {"FiniteElementOrder", order},
       {"UnknownsSize", dim},
//       {"NumberOfUniformRefinements", parallel ? 2 : 0},
       {"Hypothesis", "Tridimensional"},
       {"Parallel", true}});*/

/*  // building the non linear problem
  mfem_mgis::NonLinearEvolutionProblem problem(
      std::make_shared<mfem_mgis::FiniteElementDiscretization>(
          mesh, std::make_shared<mfem::H1_FECollection>(order, dim), dim),
      mgis::behaviour::Hypothesis::TRIDIMENSIONAL); */
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

  // solving the problem
  
  problem.setSolverParameters({{"VerbosityLevel", 0},
                               {"RelativeTolerance", 1e-6},
                               {"AbsoluteTolerance", 0.},
                               {"MaximumNumberOfIterations", 10}});

  // selection of the linear solver


    problem.setLinearSolver("CGSolver", {{"VerbosityLevel", 0},
                                             {"RelativeTolerance", 1e-12},
                                             {"MaximumNumberOfIterations",
                                             300}});
	// print on file
  out << "SetLinearSolver" << std::endl;
  out << "VerbosityLevel = " << 0 << std::endl;
  out << "RelativeTolerance = " << 1e-12 << std::endl;
  out << "MaximumNumberOfIterations = " << 300 << std::endl;
  out << "Preconditioner = " << false << std::endl;
  out << "taille_maille = " << h << std::endl;
  out << "1/h = " << 1/h << std::endl;
  out << " nbr_ref_parallel = " << nbr_ref_parallel << std::endl;
  out << " nbr_ref_sequential = " << nbr_ref_sequential << std::endl;
  out << " numbers_of_vertices = " << numbers_of_vertices << std::endl;
  out << " numbers_of_elements = " << numbers_of_elements << std::endl;
  exit(0);
    //auto prec_none = mfem_mgis::Parameters{{"Name", "None"}};

/*  auto prec_boomer =
      mfem_mgis::Parameters{{"Name", "HypreBoomerAMG"},  //
                           {"Options", mfem_mgis::Parameters{
//                                            {"Strategy", "Elasticity"},  //
                                            {"VerbosityLevel", 0}        //
                                        }}};*/

/*problem.setLinearSolver("CGSolver", {{"VerbosityLevel", 1},
                                          {"AbsoluteTolerance", 1e-12},
                                          {"RelativeTolerance", 1e-12},
                                          {"MaximumNumberOfIterations", 300},
                                          {"Preconditioner", prec_boomer}});*/

/*problem.setLinearSolver("BiCGSTABSolver", {{"VerbosityLevel", 0},
                                          {"AbsoluteTolerance", 1e-12},
                                          {"RelativeTolerance", 1e-12},
                                          {"MaximumNumberOfIterations", 300},
                                          {"Preconditioner", prec_boomer}});*/

/*problem.setLinearSolver("MINRESSolver", {{"VerbosityLevel", 1},
                                          {"AbsoluteTolerance", 1e-12},
                                          {"RelativeTolerance", 1e-12},
                                          {"MaximumNumberOfIterations", 300},
                                          {"Preconditioner", prec_boomer}});*/

/*  problem.setLinearSolver("GMRESSolver", {{"VerbosityLevel", 1},
                                          {"AbsoluteTolerance", 1e-12},
                                          {"RelativeTolerance", 1e-12},
                                          {"MaximumNumberOfIterations", 300},
                                          {"Preconditioner", prec_boomer}});*/

	// print on file
/*  out << "SetLinearSolver" << std::endl;
  out << "VerbosityLevel = " << 0 << std::endl;
  out << "RelativeTolerance = " << 1e-12 << std::endl;
  out << "MaximumNumberOfIterations = " << 300 << std::endl;
  out << "Preconditioner = " << false << std::endl;
  out << "taille_maille = " << h << std::endl;
  out << "1/h = " << t_m << std::endl;
  out << " nbr_ref_parallel = " << nbr_ref_parallel << std::endl;
  out << " nbr_ref_sequential = " << nbr_ref_sequential << std::endl;
  out << " numbers_of_vertices = " << numbers_of_vertices << std::endl;
  out << " numbers_of_elements = " << numbers_of_elements << std::endl;*/

/*     constexpr auto parallel = true;
     auto& lsf =
     mfem_mgis::LinearSolverFactory<parallel>::getFactory();
     auto lsh = lsf.generate("GMRESSolver",
     problem.getImplementation<parallel>(),
                             {{"VerbosityLevel", 1},
                              {"AbsoluteTolerance", 1e-12},
                              {"RelativeTolerance", 1e-12},
                              {"MaximumNumberOfIterations",
                              300}});
     if constexpr (parallel) {
       auto amg = std::make_unique<mfem::HypreTriSolve>();//remplacer par les autres préconditionneur
       //       amg->SetElasticityOptions(
       //
       &(problem.getImplementation<parallel>().getFiniteElementSpace());
       //    amg->SetSystemsOptions(dim);
       amg->SetPrintLevel(0);
       lsh.preconditioner = std::move(amg);//
     }
     problem.getImplementation<parallel>().updateLinearSolver(std::move(lsh));*/

//version Guillaume pour les hypreboomer

/*     constexpr auto parallel = true;
     auto& lsf =
     mfem_mgis::LinearSolverFactory<parallel>::getFactory();
     auto lsh = lsf.generate("GMRESSolver",
			     problem.getImplementation<parallel>(),
                             {{"VerbosityLevel", 0},
                              {"AbsoluteTolerance", 1e-12},
                              {"RelativeTolerance", 1e-12},
                              {"MaximumNumberOfIterations",
                              300}});
     if constexpr (parallel) {
	 //       auto amg = std::make_unique<mfem::HypreBoomerAMG>();
       //       amg->SetElasticityOptions(
       //
	 auto amg = std::make_unique<mfem::HypreSmoother>();
	 amg->SetType(mfem::HypreSmoother::l1GStr,3);
       &(problem.getImplementation<parallel>().getFiniteElementSpace());
       //    amg->SetSystemsOptions(dim);
       //       amg->SetPrintLevel(0);
       lsh.preconditioner = std::move(amg);
     }
     problem.getImplementation<parallel>().updateLinearSolver(std::move(lsh));*/

//print file preconditionner
/*  out << "SetLinearSolver" << std::endl;
  out << " " << std::endl;
  out << "VerbosityLevel = " << 1 << std::endl;
  out << "RelativeTolerance = " << 1e-12 << std::endl;
  out << "RelativeTolerance = " << 1e-12 << std::endl;
  out << "MaximumNumberOfIterations = " << 300 << std::endl;
  out << "Preconditioner = " << true << std::endl;
  out << "Preconditioner's name = " << "HypreBoomerAMG" << std::endl;
  //out << "Number of relaxation = " << 0 <<std::endl; */

	// print on file
/*  out << "SetLinearSolver" << std::endl;
  out << "VerbosityLevel = " << 0 << std::endl;
  out << "RelativeTolerance = " << 1e-12 << std::endl;
  out << "MaximumNumberOfIterations = " << 300 << std::endl;
  out << "Preconditioner = " << "HypreSmoother" << std::endl;
  out << "taille_maille = " << h << std::endl;
  out << "1/h = " << t_m << std::endl;
  out << " nbr_ref_parallel = " << nbr_ref_parallel << std::endl;
  out << " nbr_ref_sequential = " << nbr_ref_sequential << std::endl;
  out << " numbers_of_vertices = " << numbers_of_vertices << std::endl;
  out << " numbers_of_elements = " << numbers_of_elements << std::endl;*/

// selection of the linear solver
//#ifdef MFEM_USE_SUITESPARSE
//  problem.setLinearSolver("UMFPackSolver", {});
//#else
//  problem.setLinearSolver("CGSolver", {{"VerbosityLevel", 1},
//                                       {"RelativeTolerance", 1e-12},
 //                                      {"MaximumNumberOfIterations", 300}});
//#endif

  // vtk export
  problem.addPostProcessing("ParaviewExportResults",
                            {{"OutputFileName", std::string("ssna303-displacements")}});
  problem.addPostProcessing("ComputeResultantForceOnBoundary",
                            {{"Boundary", 2}, {"OutputFileName", "force.txt"}});

  // loop over time step
  const auto nsteps = mfem_mgis::size_type{50};
  const auto dt = mfem_mgis::real{1} / nsteps;
  auto t = mfem_mgis::real{0};
  auto iteration = mfem_mgis::size_type{};
  for (mfem_mgis::size_type i = 0; i != nsteps; ++i) {
    std::cout << "iteration " << iteration << " from " << t << " to " << t + dt
              << '\n';
//PRINT_DEBUG;
    // resolution
    auto ct = t;
    auto dt2 = dt;
    auto nsteps = mfem_mgis::size_type{1};
    auto niter  = mfem_mgis::size_type{0};
//PRINT_DEBUG;
    while (nsteps != 0) {
      bool converged = true;
      try {
//PRINT_DEBUG;
        problem.solve(ct, dt2);
//PRINT_DEBUG;
      } catch (std::runtime_error&) {
        converged = false;
      }
      if (converged) {
//PRINT_DEBUG;
        --nsteps;
        ct += dt2;
      } else {
//PRINT_DEBUG;
        std::cout << "\nsubstep: " << niter << '\n';
        nsteps *= 2;
        dt2 /= 2;
        ++niter;
        problem.revert();
        if (niter == 10) {
          mgis::raise("maximum number of substeps");
        }
      }
    }
    problem.executePostProcessings(t, dt);
    problem.update();
    t += dt;
    ++iteration;
    std::cout << '\n';
  }
  }
  mfem_mgis::Profiler::getProfiler().print(out);
  return EXIT_SUCCESS;
}
