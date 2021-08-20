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
#include <limits>
#include <cstdlib>
#include <iostream>
#include "mfem/general/optparser.hpp"
#include "mfem/linalg/solvers.hpp"
#include "mfem/fem/datacollection.hpp"
#include "MFEMMGIS/MFEMForward.hxx"
#include "MFEMMGIS/Material.hxx"
#include "MFEMMGIS/Parameter.hxx"
#include "MFEMMGIS/AnalyticalTests.hxx"
#include "MFEMMGIS/NonLinearEvolutionProblemImplementation.hxx"
#include "MFEMMGIS/PeriodicNonLinearEvolutionProblem.hxx"


constexpr const double xmax = 1.;
constexpr const double xthr = xmax / 2.;

void (*getSolution(const std::size_t i))(mfem::Vector&, const mfem::Vector&) {
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
                   const std::size_t i, const mfem_mgis::Parameters &params) {
  const auto b = mfem_mgis::compareToAnalyticalSolution(
      problem, getSolution(i),
      params);
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
  const char* mesh_file = "n8-id1.msh";
  const char* behaviour = "Elasticity";
  const char* library = "src/libBehaviour.so";
  const char* reference_file = "Elasticity.ref";
  int order = 1;
  int tcase = 1;
  int linearsolver = 1;
  static constexpr bool parallel = true;
};

TestParameters parseCommandLineOptions(int& argc, char* argv[]) {
  TestParameters p;

  // options treatment
  mfem::OptionsParser args(argc, argv);
  args.AddOption(&p.mesh_file, "-m", "--mesh", "Mesh file to use.");
  args.AddOption(&p.library, "-l", "--library", "Material library.");
  args.AddOption(&p.order, "-o", "--order",
                 "Finite element order (polynomial degree).");
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


template <bool parallel>
std::vector<double> findMinPoint(mfem_mgis::FiniteElementDiscretization &fed) {
  const mfem_mgis::FiniteElementSpace<parallel> & fes = fed.getFiniteElementSpace<parallel>();
  bool bynodes = fes.GetOrdering() == mfem::Ordering::byNODES;
  const auto* const mesh = fes.GetMesh();
  mfem::GridFunction nodes(&(fed.getFiniteElementSpace<parallel>()));
  mesh->GetNodes(nodes);
  const size_t nbnodes = static_cast<size_t>(nodes.Size());
  const size_t dim = static_cast<size_t>(mesh->Dimension()); //fes.GetVDim()
  std::vector<double> minCoord(dim);
  const size_t tuples = nbnodes /dim;
  MFEM_VERIFY(3 == dim, "Number of dim is not 3");
  for (int j = 0; j < dim; ++j) {
    minCoord[j] = std::numeric_limits<double>::max();
  }
  for (int i = 0; i < tuples; ++i) {
    double coord[dim];  // coordinates of a node
    for (int j = 0; j < dim; ++j) {
      if (bynodes)
	coord[j] = (nodes)[j * tuples + i];
      else
	coord[j] = (nodes)[i * dim + j];
    }
    bool isInferior = (coord[0] < minCoord[0]);
    if (coord[0] == minCoord[0]) isInferior = (coord[1] < minCoord[1]);
    if (coord[0] == minCoord[0] && coord[1] == minCoord[1]) isInferior = (coord[2] < minCoord[2]);
    if (isInferior) {
      for (int j = 0; j < dim; ++j) {
	minCoord[j] = coord[j];
      }
    }
  }
  return (minCoord);
}

template <bool parallel>
void translateMesh(mfem_mgis::FiniteElementDiscretization &fed, std::vector<double> vectorTranslate) {
  const mfem_mgis::FiniteElementSpace<parallel> & fes = fed.getFiniteElementSpace<parallel>();
  bool bynodes = fes.GetOrdering() == mfem::Ordering::byNODES;
  auto* mesh = fes.GetMesh();
  mfem::GridFunction nodes(&(fed.getFiniteElementSpace<parallel>()));
  mesh->GetNodes(nodes);
  const size_t nbnodes = static_cast<size_t>(nodes.Size());
  const size_t dim = static_cast<size_t>(mesh->Dimension()); //fes.GetVDim()
  std::vector<double> minCoord(dim);
  const size_t tuples = nbnodes /dim;
  MFEM_VERIFY(3 == dim, "Number of dim is not 3");
  MFEM_VERIFY(bynodes, "This ordering is not supported");
  mfem::Vector disp(nbnodes);
  for (int i = 0; i < tuples; ++i) {
    for (int j = 0; j < dim; ++j) {
      (disp)[j * tuples + i] = vectorTranslate[j];
    }
  }
  mesh->MoveNodes(disp);
}


// std::ifstream is("numbers.txt");
//  std::istream_iterator<double> start(is), end;
//  std::vector<double> numbers(start, end);
//  std::cout << "Read " << numbers.size() << " numbers" << std::endl;
//
//  // print the numbers to stdout
//  std::cout << "numbers read in:\n";
//  std::copy(numbers.begin(), numbers.end(), 
//            std::ostream_iterator<double>(std::cout, " "));
//  std::cout << std::endl;


int executeMFEMMGISTest(const TestParameters& p) {
  constexpr const auto dim = mfem_mgis::size_type{3};
  const auto params = mfem_mgis::Parameters {{"MeshFileName", p.mesh_file},
                            {"FiniteElementFamily", "H1"},
                            {"FiniteElementOrder", p.order},
                            {"UnknownsSize", dim},
                            {"NumberOfUniformRefinements", p.parallel ? 0 : 0},
                            {"Parallel", p.parallel},
			    {"GeneralVerbosityLevel", 1}};
  // creating the finite element workspace
  auto fed = std::make_shared<mfem_mgis::FiniteElementDiscretization>(params);
  {
    if (mfem_mgis::getMPIrank() == 0)
      std::cout << "Number of processes: " << mfem_mgis::getMPIsize() << std::endl;

    const auto verbosity = mfem_mgis::get_if<int>(params, "GeneralVerbosityLevel", 0);
    // building the non linear problem

    std::vector<mfem_mgis::real> minCoord ={ 0., 0., 0.};
    minCoord = findMinPoint<p.parallel>(*fed.get());
    for (int d=0; d < dim; d++) minCoord[d] *= -1.;
    translateMesh<p.parallel>(*fed.get(), minCoord);
    minCoord = findMinPoint<p.parallel>(*fed.get());
    std::cout << "minCoord " << minCoord[0] << " " <<
      minCoord[1] << " " << minCoord[2] << "\n";
    
    std::vector<mfem_mgis::real> corner1({0.,0.,0.});
    std::vector<mfem_mgis::real> corner2({1., 1., 1.});
    mfem_mgis::PeriodicNonLinearEvolutionProblem problem(fed, corner1, corner2);

    const mfem::ParMesh& mesh = fed->getMesh<p.parallel>();

    int nrelem = mesh.GetNE();
    const int nbgrains = mesh.attributes.Size();
    if (verbosity > 1) std::cout << "nbgrains " << nbgrains << "\n";
    std::vector<std::array<double,dim>> barycenter(nbgrains, {0., 0., 0.});
    std::vector<int> barycenter_nb(nbgrains, 0);
    
    for (int iel = 0; iel < nrelem; ++iel) {
      const mfem::Element *el = mesh.GetElement(iel);
      double sum_coords[dim] ={};
      int attr = el->GetAttribute()-1;
      MFEM_VERIFY(attr < nbgrains, "Some elements of the mesh has a number exceeding the number of grains");
      mfem::Array<int> vertices;
      el->GetVertices(vertices);
      int nrvert = vertices.Size();
      for (int iv = 0; iv < nrvert; ++iv)  {
         int vert_idx = vertices[iv];
         const double *coords = mesh.GetVertex(vert_idx);
	 for (int  d=0; d<dim; d++)
	   sum_coords[d] += coords[d];
	 barycenter_nb[attr] += 1;
      }
      for (int  d=0; d<dim; d++) {
	barycenter[attr][d] += sum_coords[d];
	
      }
    }
    // TODO: add mpi_allreduce for parallel execution
    for (int ig = 0; ig < nbgrains; ++ig)  {
      double coord[dim];
      for (int  d=0; d<dim; d++) {
	coord[d] =  std::fmod(barycenter[ig][d] / barycenter_nb[ig], xmax);
	barycenter[ig][d] = coord[d];
      }
//      std::cout << mfem_mgis::getMPIrank() << " " << mesh.GetGlobalElementNum(ig) <<
//	" attr " << ig << " bary " << barycenter[ig][0]<<" " << barycenter[ig][1]<< " " << barycenter[ig][2]<<" bool "<<
//	(barycenter[ig][0] < xthr) <<std::endl;
    }
    
    std::array<mfem_mgis::real,2> lambda({100, 200});
    std::array<mfem_mgis::real,2>     mu({75 , 150});
    for (int i=1; i<=nbgrains; i++) {
      problem.addBehaviourIntegrator("Mechanics", i, p.library, "Elasticity");
      // materials
      auto& m1 = problem.getMaterial(i);
      // setting the material properties
      auto set_properties = [](auto& m, const double l, const double mu) {
	mgis::behaviour::setMaterialProperty(m.s0, "FirstLameCoefficient", l);
	mgis::behaviour::setMaterialProperty(m.s0, "ShearModulus", mu);
	mgis::behaviour::setMaterialProperty(m.s1, "FirstLameCoefficient", l);
	mgis::behaviour::setMaterialProperty(m.s1, "ShearModulus", mu);
      };
      if (barycenter[i-1][0] < xthr)
	set_properties(m1, lambda[0], mu[0]);
      else
	set_properties(m1, lambda[1], mu[1]);
      auto set_temperature = [](auto& m) {
	mgis::behaviour::setExternalStateVariable(m.s0, "Temperature", 293.15);
	mgis::behaviour::setExternalStateVariable(m.s1, "Temperature", 293.15);
      };
      set_temperature(m1);
    }

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
#define POSTCHECK
#ifdef POSTCHECK
    if (!checkSolution(problem, p.tcase,
		       mfem_mgis::Parameters(params).insert("CriterionThreshold", 1e-4))) {
      return(EXIT_FAILURE);
    }
#endif
    return(EXIT_SUCCESS);
  }
}

int main(int argc, char* argv[]) {
  mfem_mgis::initialize(argc, argv);
  const auto p = parseCommandLineOptions(argc, argv);
  return(executeMFEMMGISTest(p));
}
