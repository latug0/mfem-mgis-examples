/*!
 * \file   InclusionsEx.cxx
 * \brief
 * This example is modelling several inclusion within a periodic cube.
 *
 * Mechanical strain:
 *                 eps = E + grad_s v
 *
 *           with  E the given macrocoscopic strain
 *                 v the periodic displacement fluctuation
 * Displacement:
 *                   u = U + v
 *
 *           with  U the given displacement associated to E
 *                   E = grad_s U
 * The local microscopic strain is equal, on average, to the macroscopic strain:
 *           <eps> = <E>
 * \author Thomas Helfer, Guillaume Latu
 * \date   02/06/10/2021
 */

#include <memory>
#include <cstdlib>
#include <iostream>
#include "mfem/general/optparser.hpp"
#include "mfem/linalg/solvers.hpp"
#include "mfem/fem/datacollection.hpp"
#include <MFEMMGIS/Profiler.hxx>
#include "MFEMMGIS/MFEMForward.hxx"
#include "MFEMMGIS/Material.hxx"
#include "MFEMMGIS/AnalyticalTests.hxx"
#include "MFEMMGIS/NonLinearEvolutionProblemImplementation.hxx"
#include "MFEMMGIS/PeriodicNonLinearEvolutionProblem.hxx"

constexpr double xmax = 1.;

void (*getSolution(const std::size_t i))(mfem::Vector&, const mfem::Vector&) {
  constexpr const auto xthr = xmax / 2.;
  std::array<void (*)(mfem::Vector&, const mfem::Vector&), 6u> solutions = {
      +[](mfem::Vector& u, const mfem::Vector& x) {
        constexpr const auto gradx = mfem_mgis::real(1) / 3;
        u = mfem_mgis::real{};
        if (x(0) < xthr) {
          u(0) = gradx * x(0);
        } else {
          u(0) = gradx * xthr - gradx * (x(0) - xthr);
        }
      },
      +[](mfem::Vector& u, const mfem::Vector& x) {
        constexpr const auto gradx = mfem_mgis::real(4) / 30;
        u = mfem_mgis::real{};
        if (x(0) < xthr) {
          u(0) = gradx * x(0);
        } else {
          u(0) = gradx * xthr - gradx * (x(0) - xthr);
        }
      },
      +[](mfem::Vector& u, const mfem::Vector& x) {
        constexpr const auto gradx = mfem_mgis::real(4) / 30;
        u = mfem_mgis::real{};
        if (x(0) < xthr) {
          u(0) = gradx * x(0);
        } else {
          u(0) = gradx * xthr - gradx * (x(0) - xthr);
        }
      },
      +[](mfem::Vector& u, const mfem::Vector& x) {
        constexpr const auto gradx = mfem_mgis::real(1) / 3;
        u = mfem_mgis::real{};
        if (x(0) < xthr) {
          u(1) = gradx * x(0);
        } else {
          u(1) = gradx * xthr - gradx * (x(0) - xthr);
        }
      },
      +[](mfem::Vector& u, const mfem::Vector& x) {
        constexpr const auto gradx = mfem_mgis::real(1) / 3;
        u = mfem_mgis::real{};
        if (x(0) < xthr) {
          u(2) = gradx * x(0);
        } else {
          u(2) = gradx * xthr - gradx * (x(0) - xthr);
        }
      },
      +[](mfem::Vector& u, const mfem::Vector&) { u = mfem_mgis::real{}; }};
  return solutions[i];
}

static void setLinearSolver(mfem_mgis::AbstractNonLinearEvolutionProblem& p,
                            const std::size_t i) {
  if (i == 0) {
    p.setLinearSolver("GMRESSolver", {{"VerbosityLevel", 1},
                                      {"AbsoluteTolerance", 1e-12},
                                      {"RelativeTolerance", 1e-12},
                                      {"MaximumNumberOfIterations", 5000}});
  } else if (i == 1) {
    p.setLinearSolver("CGSolver", {{"VerbosityLevel", 1},
                                   {"AbsoluteTolerance", 1e-12},
                                   {"RelativeTolerance", 1e-12},
                                   {"MaximumNumberOfIterations", 5000}});
#ifdef MFEM_USE_SUITESPARSE
  } else if (i == 2) {
    p.setLinearSolver("UMFPackSolver", {});
#endif
#ifdef MFEM_USE_MUMPS
  } else if (i == 3) {
    p.setLinearSolver("MUMPSSolver",
                      {{"Symmetric", true}, {"PositiveDefinite", true}});
#endif
  } else {
    std::cerr << "unsupported linear solver\n";
    mfem_mgis::abort(EXIT_FAILURE);
  }
}

static void setSolverParameters(
    mfem_mgis::AbstractNonLinearEvolutionProblem& problem) {
  problem.setSolverParameters({{"VerbosityLevel", 0},
                               {"RelativeTolerance", 1e-12},
                               {"AbsoluteTolerance", 1e-12},
                               {"MaximumNumberOfIterations", 10}});
}  // end of setSolverParmeters

bool checkSolution(mfem_mgis::NonLinearEvolutionProblem& problem,
                   const std::size_t i) {
  const auto b = mfem_mgis::compareToAnalyticalSolution(
      problem, getSolution(i), {{"CriterionThreshold", 1e-7}});
  if (!b) {
    if (mfem_mgis::getMPIrank() == 0)
      std::cerr << "Error is greater than threshold\n";
    return false;
  }
  if (mfem_mgis::getMPIrank() == 0)
    std::cerr << "Error is lower than threshold\n";
  return true;
}

struct TestParameters {
  const char* mesh_file = "cube_2mat_per.mesh";
  const char* behaviour = "Elasticity";
  const char* library = "src/libBehaviour.so";
  const char* reference_file = "Elasticity.ref";
  int order = 1;
  int tcase = 1;
  int linearsolver = 1;
  double xmax = 1.;
  double ymax = 1.;
  double zmax = 1.;
  bool parallel = true;
};

TestParameters parseCommandLineOptions(int& argc, char* argv[]) {
  TestParameters p;

  // options treatment
  mfem::OptionsParser args(argc, argv);
  args.AddOption(&p.mesh_file, "-m", "--mesh", "Mesh file to use.");
  args.AddOption(&p.library, "-l", "--library", "Material library.");
  args.AddOption(&p.order, "-o", "--order",
                 "Finite element order (polynomial degree).");
  args.AddOption(&p.xmax, "-xm", "--xmax",
                 "Corner, coordinate x direction.");
  args.AddOption(&p.ymax, "-ym", "--ymax",
                 "Corner coordinate y direction.");
  args.AddOption(&p.zmax, "-zm", "--zmax",
                 "Corner coordinate z direction.");
  args.AddOption(&p.tcase, "-t", "--test-case",
                 "identifier of the case : Exx->0, Eyy->1, Ezz->2, Exy->3, "
                 "Exz->4, Eyz->5");
  args.AddOption(
      &p.linearsolver, "-ls", "--linearsolver",
      "identifier of the linear solver: 0 -> GMRES, 1 -> CG, 2 -> UMFPack");
  args.Parse();
  if (!args.Good()) {
    if (mfem_mgis::getMPIrank() == 0)
      args.PrintUsage(std::cout);
    mfem_mgis::finalize();
    exit(0);
  }
  if (p.mesh_file == nullptr) {
    if (mfem_mgis::getMPIrank() == 0)
      std::cout << "ERROR: Mesh file missing" << std::endl;
    args.PrintUsage(std::cout);
    mfem_mgis::abort(EXIT_FAILURE);
  }
  if (mfem_mgis::getMPIrank() == 0)
    args.PrintOptions(std::cout);
  if ((p.tcase < 0) || (p.tcase > 5)) {
    std::cerr << "Invalid test case\n";
    mfem_mgis::abort(EXIT_FAILURE);
  }
  return p;
}

int executeMFEMMGISTest(const TestParameters& p) {
  constexpr const auto dim = mfem_mgis::size_type{3};
  // creating the finite element workspace

  auto fed = std::make_shared<mfem_mgis::FiniteElementDiscretization>(
      mfem_mgis::Parameters{{"MeshFileName", p.mesh_file},
                            {"FiniteElementFamily", "H1"},
                            {"FiniteElementOrder", p.order},
                            {"UnknownsSize", dim},
                            {"NumberOfUniformRefinements", p.parallel ? 0 : 0},
                            {"Parallel", p.parallel}});

  {
    if (mfem_mgis::getMPIrank() == 0)
      std::cout << "Number of processes: " << mfem_mgis::getMPIsize() << std::endl;
    // building the non linear problem

    std::vector<mfem_mgis::real> corner1({0.,0.,0.});
    std::vector<mfem_mgis::real> corner2({p.xmax, p.ymax, p.zmax});
    mfem_mgis::PeriodicNonLinearEvolutionProblem problem(fed, corner1, corner2);


    //    const mfem::Mesh &m = fed->getMesh<true>();

    
    problem.addBehaviourIntegrator("Mechanics", 1, p.library, "Elasticity");
    problem.addBehaviourIntegrator("Mechanics", 2, p.library, "Elasticity");
    // materials
    auto& m1 = problem.getMaterial(1);
    auto& m2 = problem.getMaterial(2);
    // setting the material properties
    auto set_properties = [](auto& m, const double l, const double mu) {
      mgis::behaviour::setMaterialProperty(m.s0, "FirstLameCoefficient", l);
      mgis::behaviour::setMaterialProperty(m.s0, "ShearModulus", mu);
      mgis::behaviour::setMaterialProperty(m.s1, "FirstLameCoefficient", l);
      mgis::behaviour::setMaterialProperty(m.s1, "ShearModulus", mu);
    };

    std::array<mfem_mgis::real,2> lambda({100, 200});
    std::array<mfem_mgis::real,2>     mu({75 , 150});
    set_properties(m1, lambda[0], mu[0]);
    set_properties(m2, lambda[1], mu[1]);
    //
    auto set_temperature = [](auto& m) {
      mgis::behaviour::setExternalStateVariable(m.s0, "Temperature", 293.15);
      mgis::behaviour::setExternalStateVariable(m.s1, "Temperature", 293.15);
    };
    set_temperature(m1);
    set_temperature(m2);

    // macroscopic strain
    std::vector<mfem_mgis::real> e(6, mfem_mgis::real{});
    if (p.tcase < 3) {
      e[p.tcase] = 1;
    } else {
      e[p.tcase] = 1.41421356237309504880 / 2;
    }
    problem.setMacroscopicGradientsEvolution([e](const double) { return e; });
    //
    setLinearSolver(problem, p.linearsolver);
    setSolverParameters(problem);

    // Add postprocessing and outputs
    problem.addPostProcessing(
        "ParaviewExportResults",
        {{"OutputFileName", "PeriodicTestOutput-" + std::to_string(p.tcase)}});
    // solving the problem
    if (!problem.solve(0, 1)) {
      mfem_mgis::abort(EXIT_FAILURE);
    }
    problem.executePostProcessings(0, 1);
    //
    if (!checkSolution(problem, p.tcase)) {
      return(EXIT_FAILURE);
    }
    mfem_mgis::Profiler::timers::print_and_write_timers();
    return(EXIT_SUCCESS);
  }
}

int main(int argc, char* argv[]) {
  mfem_mgis::initialize(argc, argv);
  const auto p = parseCommandLineOptions(argc, argv);
  return(executeMFEMMGISTest(p));
}
