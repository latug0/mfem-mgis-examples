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
  auto ctx = mgis::Context{};
  auto or_die = ctx.getFatalFailureHandler();
  //
  static constexpr const auto parallel = false;
  // options treatment
  mfem_mgis::initialize(argc, argv);
  auto success = true;
  // building the non linear problem
  auto problem =
      mfem_mgis::construct<mfem_mgis::NonLinearEvolutionProblem>(
          ctx,
          mfem_mgis::Parameters{
              {"MeshFileName", "./cube.msh"},
              {"Materials", mfem_mgis::Parameters{{"plate", 1}}},
              {"Boundaries", mfem_mgis::Parameters{{"left", 2}, {"right", 4}}},
              {"FiniteElementFamily", "H1"},
              {"FiniteElementOrder", 2},
              {"UnknownsSize", 2},
              {"NumberOfUniformRefinements", parallel ? 2 : 0},
              {"Hypothesis", "PlaneStrain"},
              {"Parallel", parallel}}) |
      or_die;
  // materials
  problem.addBehaviourIntegrator(ctx, "Mechanics", "plate",
                                 "./src/libBehaviour.so", "Elasticity") |
      or_die;
  auto& m1 = problem.getMaterial(ctx, "plate", 0) | or_die;
  // material properties at the beginning and the end of the time step
  for (auto& s : {&m1.s0, &m1.s1}) {
    mgis::behaviour::setMaterialProperty(ctx, *s, "YoungModulus", 150e9) |
        or_die;
    mgis::behaviour::setMaterialProperty(ctx, *s, "PoissonRatio", 0.3) | or_die;
  }
  // temperature
  mgis::behaviour::setExternalStateVariable(ctx, m1.s0, "Temperature", 293.15) |
      or_die;
  auto Tg = mfem_mgis::PartialQuadratureFunction::evaluate(
      m1.getPartialQuadratureSpacePointer(),
      [](const mfem_mgis::real x, const mfem_mgis::real) {
        constexpr auto Tref = mfem_mgis::real{293.15};
        constexpr auto dT = mfem_mgis::real{2000} - Tref;
        return 4 * dT * x * (1 - x) + Tref;
      });
  mgis::behaviour::setExternalStateVariable(ctx, m1.s1, "Temperature",
                                            Tg->getValues()) |
      or_die;
  // boundary conditions
  for (const auto boundary : {"left", "right"}) {
    for (const auto dof : {0, 1}) {
      problem.addBoundaryCondition(
          ctx,
          mfem_mgis::make_unique<mfem_mgis::UniformDirichletBoundaryCondition>(
              ctx, problem.getFiniteElementDiscretizationPointer(), boundary,
              dof) |
              or_die) |
          or_die;
    }
  }
  // set the solver parameters
  problem.setLinearSolver(ctx, "UMFPackSolver", {}) | or_die;
  problem.setSolverParameters(ctx, {{"VerbosityLevel", 0},
                                    {"RelativeTolerance", 1e-12},
                                    {"AbsoluteTolerance", 0.},
                                    {"MaximumNumberOfIterations", 10}}) |
      or_die;
  // vtk export
  problem.addPostProcessing(ctx, "ParaviewExportResults",
                            {{"OutputFileName", "SatohTestOutput"}}) |
      or_die;
  auto results = std::vector<mfem_mgis::Parameter>{
      "Stress", "ImposedTemperature", "HydrostaticPressure"};
  problem.addPostProcessing(
      ctx, "ParaviewExportIntegrationPointResultsAtNodes",
      {{"OutputFileName", "SatohTestIntegrationPointOutput"},
       {"Materials", {"plate"}},
       {"Results", results}}) |
      or_die;
  // solving the problem on 1 time step
  auto r = problem.solve(ctx, 0, 1);
  if (r) {
    problem.executePostProcessings(ctx, 0, 1) | or_die;
    //
    std::ofstream output("HydrostaticPressure.txt");
    const auto pr = getInternalStateVariable(
                        ctx, static_cast<const mfem_mgis::Material&>(m1),
                        "HydrostaticPressure") |
                    or_die;
    dumpPartialQuadratureFunction<parallel>(output, pr);
  }
  //
  return success ? EXIT_SUCCESS : EXIT_FAILURE;
}
