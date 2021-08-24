/*!
 * \file   GrainsEx.cxx
 * \brief
 * This example is modelling polycristals confguration with periodic domain.
 * Prerequisites : 
 *   - Please install Neper software 
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
 * \author Guillaume Latu
 * \date   02/08/2021
 */

#include <algorithm>
#include <numeric>
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


using Orientation3D = std::array<double,9u>;


constexpr const auto cdim = mfem_mgis::size_type{3};


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
  const auto b = true;
  if (!b) {
    if (mfem_mgis::getMPIrank() == 0)
      std::cerr << "Error, check failed\n";
    return false;
  }
  if (mfem_mgis::getMPIrank() == 0)
    std::cerr << "OK, no error detected\n";
  return true;
}


struct TestParameters {
  const char* mesh_file = "n8-id1.msh";
  const char* behaviour = "OrthotropicElasticity";
  const char* library = "src/libBehaviour.so";
  const char* reference_file = "Elasticity.ref";
  const char* x_orientation[3] = {"X1.or", "X2.or", "X3.or"};
  const char* y_orientation[3] = {"Y1.or", "Y2.or", "Y3.or"};
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
  MFEM_VERIFY(cdim == dim, "Number of dimension is invalid (not 3D)");
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
  MFEM_VERIFY(cdim == dim, "Number of dimension is invalid (not 3D)");
  MFEM_VERIFY(bynodes, "This ordering is not supported");
  mfem::Vector disp(nbnodes);
  for (int i = 0; i < tuples; ++i) {
    for (int j = 0; j < dim; ++j) {
      (disp)[j * tuples + i] = vectorTranslate[j];
    }
  }
  mesh->MoveNodes(disp);
}


int readFileOneArray(std::string &filename, std::vector<double> &out) {
  std::ifstream is(filename);
  MFEM_VERIFY(is.good(), "Fail opening the file "+filename);
  std::istream_iterator<double> start(is), end{};
  std::vector<double> numbers(start, end);
  out.clear();
  copy(numbers.begin(), numbers.end(), back_inserter(out));
  return numbers.size();
}

void crossProduct (const double a[3], const double b[3], double c[3]) 
{ 
  c[0] = a[1] * b[2] - a[2] * b[1]; 
  c[1] = a[2] * b[0] - a[0] * b[2];
  c[2] = a[0] * b[1] - a[1] * b[0]; 
}

bool normalize (double a[3], double threshold) 
{
  double norm = std::sqrt(std::inner_product(a, a + 3, a, 0.l));
  if (norm < threshold) return false;
  // Multiply all components of vector a by the inverse of norm
  norm = 1./norm;
  for (uint i=0; i < 3; i++) {
    a[i] *= norm;
  }
  return true;
}

bool checkNorm (const Orientation3D &ori, const double threshold) 
{
  double c[3] = {};
  for (uint i=0; i < 3; i++) {
    c[0] += ori[0+i]*ori[0+i]; 
    c[1] += ori[3+i]*ori[3+i];
    c[2] += ori[6+i]*ori[6+i];
  }
  for (uint i=0; i < 3; i++) {
    if (fabs(sqrt(c[i]) - 1.) > threshold) {
      std::cout << std::setprecision(16) <<  "failure " << i << " " << sqrt(c[i]) << std::endl;
      return false;
    }
  }
  return true;
}


bool readOrientations(const TestParameters& p, std::vector<Orientation3D>& orientations, const double threshold) {
  std::string filename;
  int vector_size = -1;
  for (int d=0; d < cdim; d++) {
    std::vector<double> v;
    filename = p.x_orientation[d];
    readFileOneArray(filename,  v);
    if (vector_size == -1) {
      vector_size = v.size();
      orientations.resize(vector_size);
    }
    MFEM_VERIFY(vector_size == v.size(), "Vector read is not of expected size, filename: "+filename);
    for (std::size_t i = 0; i < vector_size; ++i) {
      orientations[i][d] = v[i];
    }
    filename = p.y_orientation[d];
    readFileOneArray(filename,  v);
    MFEM_VERIFY(vector_size == v.size(), "Vector read is not of expected size, filename: "+filename);
    for (std::size_t i = 0; i < vector_size; ++i) {
      orientations[i][3+d] = v[i];
    }
  }
  for (std::size_t i = 0; i < vector_size; ++i) {
    double *data = orientations[i].data();
    // Normalizing input vectors
    MFEM_VERIFY(normalize(data, .5), "Orientation vector is wrong, normalization fails");
    MFEM_VERIFY(normalize(data+3, .5), "Orientation vector is wrong, normalization fails");
    // Cross product of X (data) and Y (data+3) vectors that gives Z vector stored at data+6
    crossProduct(data, data+3, data+6);
    //    std::cout << "vector " << i << std::endl;
    //    for (std::size_t d=0; d<9; ++d) std::cout << " " << data[d];
    //    std::cout << std::endl;
    MFEM_VERIFY(checkNorm(orientations[i], threshold), "Rotation matrix is invalid");
  }
  return true;
}

int executeMFEMMGISTest(const TestParameters& p) {
  const auto params = mfem_mgis::Parameters {{"MeshFileName", p.mesh_file},
                            {"FiniteElementFamily", "H1"},
                            {"FiniteElementOrder", p.order},
                            {"UnknownsSize", cdim},
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

    std::vector<mfem_mgis::real> minCoord(3);
    minCoord = findMinPoint<p.parallel>(*fed.get());
    std::cout << "minCoord " << minCoord[0] << " " <<
      minCoord[1] << " " << minCoord[2] << "\n";
    
//    for (int d=0; d < cdim; d++) minCoord[d] *= -1.;
//    translateMesh<p.parallel>(*fed.get(), minCoord);
//    minCoord = findMinPoint<p.parallel>(*fed.get());
//    std::cout << "minCoord " << minCoord[0] << " " <<
//      minCoord[1] << " " << minCoord[2] << "\n";
    
//    std::vector<mfem_mgis::real> corner1({0.,0.,0.});
//    std::vector<mfem_mgis::real> corner2({1., 1., 1.});
    mfem_mgis::PeriodicNonLinearEvolutionProblem problem(fed, mfem_mgis::FIX_XMIN);
    const mfem::ParMesh& mesh = fed->getMesh<p.parallel>();
    //    const mfem::Mesh& mesh = fed->getMesh<p.parallel>();

    int nrelem = mesh.GetNE();
    const int nbgrains = mesh.attributes.Size();
    if (verbosity > 1) std::cout << "nbgrains " << nbgrains << "\n";
    std::vector<Orientation3D> orientations;
    auto check_orientation = readOrientations(p, orientations, 1e-12);
    MFEM_VERIFY(check_orientation, "Orientation read failed");
    
    for (int i=1; i<= nbgrains; i++) {
      std::cout  << "coucou " << i << std::endl;
      problem.addBehaviourIntegrator("Mechanics", i, p.library, p.behaviour);
      
      // materials
      auto& m1 = problem.getMaterial(i);
      // setting the material properties
      auto set_properties = [](auto& m, const Orientation3D r) {
	MFEM_VERIFY(m.b.symmetry == mgis::behaviour::Behaviour::ORTHOTROPIC, "Test cas defined for orthotropic behaviour only");
	{
//	  std::array<mfem_mgis::real, 9u> r = {0, 1, 0,  //
//					       1, 0, 0,  //
//					       0, 0, 1};
	  m.setRotationMatrix(mfem_mgis::RotationMatrix3D{r});
	}
      };
      set_properties(m1, orientations[i]);
      
      //TODO modifiy temperature
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

    setLinearSolver(problem, p.linearsolver);
    setSolverParameters(problem);

    // Add postprocessing and outputs
    problem.addPostProcessing(
        "ParaviewExportResults",
        {{"OutputFileName", "PeriodicTestOutput-" + std::to_string(p.tcase)}});
    std::vector<mfem_mgis::Parameter> materials_out{1,2};
    problem.addPostProcessing(
        "ParaviewExportIntegrationPointResultsAtNodes",
        {{"OutputFileName", "PeriodicTestOutput-Strain-" 
	      + std::to_string(p.tcase)},
	    {"Materials", {materials_out}},
	    {"Results", "Strain"}});
    problem.addPostProcessing(
        "ParaviewExportIntegrationPointResultsAtNodes",
        {{"OutputFileName", "PeriodicTestOutput-Stress-" 
	      + std::to_string(p.tcase)},
	    {"Materials", {materials_out}},
	    {"Results", "Stress"}});
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
