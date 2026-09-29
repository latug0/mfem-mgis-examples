/*!
 * \file   BendingTest.cxx
 * \brief
 * \author Thomas Helfer, Guillaume Latu
 * \date   06/04/2021
 */

// #include <memory>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include "mfem/general/optparser.hpp"
// #include "mfem/fem/datacollection.hpp"
// #include "MFEMMGIS/Simulation.hxx"
// #include "MFEMMGIS/L2Projection.hxx"
#include "MFEMMGIS/LinearSolverFactory.hxx"
#include "MFEMMGIS/Material.hxx"
// #include "MFEMMGIS/MechanicalPostProcessings.hxx"
#include "MFEMMGIS/NonLinearEvolutionProblem.hxx"
// #include "MFEMMGIS/ParaviewExportIntegrationPointResultsAtNodes.hxx"

#include "MFEMMGIS/Config.hxx"


#include "MFEMMGIS/Profiler.hxx"

#include <sys/time.h>
#include <sys/resource.h>


#include "ThirdMediumUtils.hxx"



int main(int argc, char **argv) {
///feenableexcept(FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW); // DEBUG 
   // options treatment (contains default values for the parameters. Default values are in the header ThirdMediumUtils.hxx)
  thirdmedium_utils::TestParameters p;
  p.mesh_file="bending.msh";
  mfem::OptionsParser args = parse_options(p,argc,argv);

    
  constexpr const auto dim = mfem_mgis::size_type{2};
  constexpr auto nsteps = mfem_mgis::size_type{50};
  constexpr auto H = mfem_mgis::real{0.6};
  constexpr auto te = mfem_mgis::real{1};
  auto factor = mfem_mgis::real{p.gamma};
  constexpr auto Ks = mfem_mgis::real{200e9};
  constexpr auto Gs = mfem_mgis::real{200e9};
  constexpr auto E = mfem_mgis::real{210e9};
  constexpr auto nu = mfem_mgis::real{0.3};
  // Initialize mfem_mgis (it includes a call to MPI_Init)
  mfem_mgis::initialize(argc, argv);
  //
  // init timers (deprecated)
  //mfem_mgis::Profiler::timers::init_timers();

  auto ctx = mfem_mgis::Context{};
  auto or_die = ctx.getFatalFailureHandler();

  // init timers
  ctx.enableProfiling(true);

  const char *behaviour1 = "Elasticity";
  const char *behaviour2 = "ThirdMediumContactHyperElasticBehaviour";
  std::string library = p.library;
 #if defined(MFEM_USE_MUMPS) && defined(MFEM_USE_MPI)
     constexpr bool parallel = true;
   #else
  constexpr bool parallel = false;
  #endif
  
   if (mfem_mgis::getMPIrank() == 0){
  args.PrintOptions(std::cout);
   std::cout << "Parallel : " << p.parallel << '\n';
  }



  // the non linear problem
  mfem_mgis::NonLinearEvolutionProblem mechanics(ctx,{{"MeshFileName", p.mesh_file},
                                                  {"FiniteElementFamily", "H1"},
                                                  {"FiniteElementOrder", p.order},
                                                  {"UnknownsSize", dim},
                                                  {"NumberOfUniformRefinements", p.parallel ? p.refinement : 0},
                                                  {"Hypothesis", "PlaneStrain"},
                                                  {"Parallel", p.parallel}});

    // Mesh and memory info
  thirdmedium_utils::print_mesh_information(mechanics.getImplementation<parallel>());
  thirdmedium_utils::print_memory_footprint("After_problem:");
  thirdmedium_utils::save_unknowns_info(mechanics.getImplementation<parallel>(),p.output_dir);

  //
  auto faltus_parameters = mfem_mgis::Parameters{
      {"Regularization",
       mfem_mgis::Parameters{
           {"Faltus2026",   
            mfem_mgis::Parameters{{"PenalizationCoefficient", p.alpha}}}}}};


//  mechanics.addBehaviourIntegrator(ctx, "Mechanics", "VOID1", library, behaviour2) | or_die;
//  mechanics.addBehaviourIntegrator(ctx, "Mechanics", "VOID2", library, behaviour2) | or_die;

  mechanics.addBehaviourIntegrator(ctx, "Mechanics", "VOID1", library,
                                   behaviour2, faltus_parameters) |
      or_die;
  mechanics.addBehaviourIntegrator(ctx, "Mechanics", "VOID2", library,
                                   behaviour2, faltus_parameters) |
      or_die;

  mechanics.addBehaviourIntegrator(ctx, "Mechanics", "POUTRE", library,
                                   behaviour1) |
      or_die;

  // material proprerties
  for (const auto &n : {"VOID1", "VOID2"}) {
    auto &m = mechanics.getMaterial(ctx, n, 0) | or_die;
    setMaterialProperty(ctx, m.s0, "BulkModulus", factor * Ks) | or_die;
    setMaterialProperty(ctx, m.s1, "BulkModulus", factor * Ks) | or_die;
    setMaterialProperty(ctx, m.s0, "ShearModulus", factor * Gs) | or_die;
    setMaterialProperty(ctx, m.s1, "ShearModulus", factor * Gs) | or_die;
  }
  auto &m = mechanics.getMaterial(ctx, "POUTRE", 0) | or_die;
  setMaterialProperty(ctx, m.s0, "YoungModulus", E) | or_die;
  setMaterialProperty(ctx, m.s1, "YoungModulus", E) | or_die;
  setMaterialProperty(ctx, m.s0, "PoissonRatio", nu) | or_die;
  setMaterialProperty(ctx, m.s1, "PoissonRatio", nu) | or_die;

  // boundary conditions
  for (const auto &n : {"LPOUG", "LVOI1G", "LVOI1D", "LBOD1D"}) {
    mechanics.addUniformDirichletBoundaryCondition(ctx,
        {{"Boundary", n}, {"Component", 0}})|or_die;
  }
  mechanics.addUniformDirichletBoundaryCondition(ctx,
      {{"Boundary", "LBOD1D"},
       {"Component", 1},
       {"LoadingEvolution", [](const auto t) {
          const auto u = -H * (t / te);
          return u;
        }}})|or_die;
  mechanics.addUniformDirichletBoundaryCondition(ctx,
      {{"Boundary", "LBOD2H"}, {"Component", 0}})|or_die;
  mechanics.addUniformDirichletBoundaryCondition(ctx,
      {{"Boundary", "LBOD2H"}, {"Component", 1}})|or_die;

  // solving the problem
  mechanics.setPredictionPolicy(
      {.strategy =
       mfem_mgis::PredictionStrategy::BEGINNING_OF_TIME_STEP_PREDICTION});

       // mfem_mgis::PredictionStrategy::
       //     CONSTANT_GRADIENTS_INTEGRATION_PREDICTION});

  mechanics.setSolverParameters(ctx, {{"VerbosityLevel", p.verbosity_level},
                                      {"RelativeTolerance", 1e-4},
                                      {"AbsoluteTolerance", 0.},
                                      {"MaximumNumberOfIterations", p.newton_iterations}}) |
      or_die;

  // returns the solver parameters and sets the solver name string and number of iterations in the struct p.  
  auto solverParameters = thirdmedium_utils::set_solver_parameters(p); 

  // returns the preconditioner parameters and sets the preconditioner string in the struct p. 
  mfem_mgis::Parameters prec = thirdmedium_utils::set_preconditioner_parameters(p); 

  // Print out TestParameters once the solver has been selectioned
  if (mfem_mgis::getMPIrank()==0) std::cout << p << '\n';
 
      // selection of the linear solver
  if constexpr (parallel) {
    if (p.preconditioner != -1 && p.solver !=0) // Solver 0 (MUMPS) doesn't need a preconditioner (error is thrown otherwise)
      solverParameters.insert(mfem_mgis::may_throw,mfem_mgis::Parameters{{"Preconditioner",prec}});
        
    mechanics.setLinearSolver(ctx, p.solver_name, solverParameters) | or_die; 

      //mechanics.setLinearSolver(ctx, p.solver_name, solverParameters) | or_die; 
  } else {
    mechanics.setLinearSolver(ctx, "UMFPackSolver", {}) | or_die;
  }
  //
  // Enable or disable post-processing for benchmarks
  if (p.post_processing){
    mechanics.addPostProcessing(
      ctx, "ParaviewExportResults",
      {{"OutputFileName", p.output_dir + std::string("displacements")}}) |
      or_die;
  }
  

std::vector<std::string> materials = {"VOID1", "VOID2"};
thirdmedium_utils::setup_additional_convergence_criteria(
    mechanics, ctx, p, solverParameters, factor * Ks, factor * Gs, std::move(materials));
  // loop over time steps
 // const auto times = mfem_mgis::Simulation::TimesDescription{0, te, nsteps};
 // auto s = mfem_mgis::Simulation{mechanics, times};
 // std::ignore = s.run(ctx);
 // std::cout << ctx.getErrorMessage() << '\n';

  auto [exit_status, sim_output] = thirdmedium_utils::run_solve(ctx, mechanics, 0, te, nsteps);
    if (mfem_mgis::getMPIrank()==0){
        std::cout << "Exit status : " << (exit_status == mfem_mgis::ExitStatus::success ? "Success\n" : "Failed\n");
        
        std::ofstream status(p.output_dir + "run.status");

        if (exit_status == mfem_mgis::ExitStatus::success)
        {
                status << "SUCCESS\n";
        }
        else
        {
                status << "FAILED\n";
        }
                
        if (sim_output.has_value()) {
            thirdmedium_utils::check_and_save_convergence_info(sim_output.value(), exit_status);
        }else{
            std::cout<< "No value !!!!!!\n" ;
        }
    }


  // print and write timetable
  thirdmedium_utils::print_memory_footprint("After Solving:");
  mfem_mgis::Profiler::OutputManager::printTimeTable(ctx);
  mfem_mgis::Profiler::OutputManager::writeFile(ctx,p.output_dir + mfem_mgis::Profiler::OutputManager::build_name());
  std::flush(std::cout);
  std::flush(std::cerr);
  return EXIT_SUCCESS;
}   
