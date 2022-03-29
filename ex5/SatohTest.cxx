/*!
 * \file   SatohTest.cxx
 * \brief
 * \author Thomas Helfer
 * \date   28/03/2022
 */

#include <memory>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include "mfem/linalg/vector.hpp"
#include "mfem/fem/fespace.hpp"
#include "MFEMMGIS/Parameters.hxx"
#include "MFEMMGIS/Material.hxx"
#include "MFEMMGIS/PartialQuadratureSpace.hxx"
#include "MFEMMGIS/UniformDirichletBoundaryCondition.hxx"
#include "MFEMMGIS/NonLinearEvolutionProblem.hxx"

template <bool parallel>
static void dumpPartialQuadratureFunction(
    std::ostream& os,
    const mfem_mgis::ImmutablePartialQuadratureFunctionView& f) {
  const auto& s = f.getPartialQuadratureSpace();
  const auto& fed = s.getFiniteElementDiscretization();
  const auto& fespace = fed.getFiniteElementSpace<parallel>();
  const auto m = s.getId();
  for (mfem_mgis::size_type i = 0; i != fespace.GetNE(); ++i) {
    if (fespace.GetAttribute(i) != m) {
      continue;
    }
    const auto& fe = *(fespace.GetFE(i));
    auto& tr = *(fespace.GetElementTransformation(i));
    const auto& ir = s.getIntegrationRule(fe, tr);
    for (mfem_mgis::size_type g = 0; g != ir.GetNPoints(); ++g) {
      // get the gradients of the shape functions
      mfem::Vector p;
      const auto& ip = ir.IntPoint(g);
      tr.SetIntPoint(&ip);
      tr.Transform(tr.GetIntPoint(), p);
      os << p[0] << " " << p[1] << " " << f.getIntegrationPointValue(i, g)
         << '\n';
    }
  }
}  // end of dumpPartialQuadratureFunction

/*
 * This test models a 2D plate of lenght 1 in plane strain clamped on the left
 * and right boundaries and submitted to a parabolic thermal gradient along the
 * x-axis:
 *
 * - the temperature profile is minimal on the left and right boundaries
 * - the temperature profile is maximal for x = 0.5
 *
 * This example shows how to define an external state variable using an
 * analytical profile.
 */
int main(int argc, char** argv) {
  //
  static constexpr const auto parallel = false;
  // options treatment
  mfem_mgis::initialize(argc, argv);
  auto success = true;
  // building the non linear problem
  mfem_mgis::NonLinearEvolutionProblem problem(
      {{"MeshFileName", "./cube.msh"},
       {"Materials", mfem_mgis::Parameters{{"plate", 1}}},
       {"Boundaries", mfem_mgis::Parameters{{"left", 2}, {"right", 4}}},
       {"FiniteElementFamily", "H1"},
       {"FiniteElementOrder", 2},
       {"UnknownsSize", 2},
       {"NumberOfUniformRefinements", parallel ? 2 : 0},
       {"Hypothesis", "PlaneStrain"},
       {"Parallel", parallel}});
  // materials
  problem.addBehaviourIntegrator("Mechanics", "plate", "./src/libBehaviour.so",
                                 "Elasticity");
  auto& m1 = problem.getMaterial("plate");
  // material properties at the beginning and the end of the time step
  for (auto& s : {&m1.s0, &m1.s1}) {
    mgis::behaviour::setMaterialProperty(*s, "YoungModulus", 200e9);
    mgis::behaviour::setMaterialProperty(*s, "PoissonRatio", 0.499);
  }
  // temperature
  mgis::behaviour::setExternalStateVariable(m1.s0, "Temperature", 293.15);
  auto Tg = mfem_mgis::PartialQuadratureFunction::evaluate(
      m1.getPartialQuadratureSpacePointer(),
      [](const mfem_mgis::real x, const mfem_mgis::real) {
        constexpr auto Tref = mfem_mgis::real{293.15};
        constexpr auto dT = mfem_mgis::real{2000} - Tref;
        return 4 * dT * x * (1 - x) + Tref;
      });
  mgis::behaviour::setExternalStateVariable(m1.s1, "Temperature",
                                            Tg->getValues());
  // boundary conditions
  for (const auto boundary : {"left", "right"}) {
    for (const auto dof : {0, 1}) {
      problem.addBoundaryCondition(
          std::make_unique<mfem_mgis::UniformDirichletBoundaryCondition>(
              problem.getFiniteElementDiscretizationPointer(), boundary, dof));
    }
  }
  // set the solver parameters
  problem.setLinearSolver("UMFPackSolver", {});
  problem.setSolverParameters({{"VerbosityLevel", 0},
                               {"RelativeTolerance", 1e-12},
                               {"AbsoluteTolerance", 0.},
                               {"MaximumNumberOfIterations", 10}});
  // vtk export
  problem.addPostProcessing("ParaviewExportResults",
                            {{"OutputFileName", "SatohTestOutput"}});
  auto results = std::vector<mfem_mgis::Parameter>{
      "Stress", "ImposedTemperature", "HydrostaticPressure"};
  problem.addPostProcessing(
      "ParaviewExportIntegrationPointResultsAtNodes",
      {{"OutputFileName", "SatohTestIntegrationPointOutput"},
       {"Materials", {"plate"}},
       {"Results", results}});
  // solving the problem on 1 time step
  auto r = problem.solve(0, 1);
  problem.executePostProcessings(0, 1);
  //
  std::ofstream output("HydrostaticPressure.txt");
  const auto pr = getInternalStateVariable(
      static_cast<const mfem_mgis::Material&>(m1), "HydrostaticPressure");
  dumpPartialQuadratureFunction<parallel>(output, pr);
  //
  return success ? EXIT_SUCCESS : EXIT_FAILURE;
}
