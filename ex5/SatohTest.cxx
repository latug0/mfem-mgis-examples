/*!
 * \file   SatohTest.cxx
 * \brief
 * \author Thomas Helfer
 * \date   28/03/2022
 */

#include <cmath>
#include <memory>
#include <cstdlib>
#include <algorithm>
#include <fstream>
#include <iostream>
#include "mfem/linalg/vector.hpp"
#include "mfem/fem/fespace.hpp"
#include "MFEMMGIS/Parameters.hxx"
#include "MFEMMGIS/Material.hxx"
#include "MFEMMGIS/PartialQuadratureSpace.hxx"
#include "MFEMMGIS/UniformDirichletBoundaryCondition.hxx"
#include "MFEMMGIS/NonLinearEvolutionProblem.hxx"

/*!
 * \brief call a function at each integration point of a partial quadrature
 * function, with the coordinates of the point and the value of the function
 * \param[in] f: partial quadrature function
 * \param[in] c: function called
 */
template <bool parallel, typename Callable>
static void forEachIntegrationPoint(
    const mfem_mgis::ImmutablePartialQuadratureFunctionView& f, Callable&& c) {
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
      mfem::Vector p;
      const auto& ip = ir.IntPoint(g);
      tr.SetIntPoint(&ip);
      tr.Transform(tr.GetIntPoint(), p);
      c(p[0], p[1], f.getIntegrationPointValue(i, g));
    }
  }
}  // end of forEachIntegrationPoint

//! \return the temperature imposed at the end of the time step
static mfem_mgis::real getTemperature(const mfem_mgis::real x) noexcept {
  constexpr auto Tref = mfem_mgis::real{293.15};
  constexpr auto dT = mfem_mgis::real{2000} - Tref;
  return 4 * dT * x * (1 - x) + Tref;
}  // end of getTemperature

/*!
 * \return if the temperature seen by the behaviour at each integration point
 * is the imposed one
 * \param[in] T: temperature at the integration points
 */
template <bool parallel>
static bool checkImposedTemperature(
    const mfem_mgis::ImmutablePartialQuadratureFunctionView& T) {
  auto e = mfem_mgis::real{0};
  forEachIntegrationPoint<parallel>(
      T, [&e](const mfem_mgis::real x, const mfem_mgis::real,
              const mfem_mgis::real v) {
        e = std::max(e, std::abs(v - getTemperature(x)));
      });
  // the temperatures are expressed in Kelvin
  if (e > 1e-9) {
    mfem_mgis::getErrorStream()
        << "invalid imposed temperature (maximal error: " << e << ")\n";
    return false;
  }
  return true;
}  // end of checkImposedTemperature

/*!
 * \return if the horizontal reactions of the left and right boundaries are
 * opposite and match the reference value
 * \param[in, out] ctx: execution context
 * \param[in] p: non linear evolution problem
 */
static bool checkReactions(mfem_mgis::Context& ctx,
                           mfem_mgis::NonLinearEvolutionProblem& p) {
  auto or_die = ctx.getFatalFailureHandler();
  // horizontal reaction of a boundary
  auto reaction = [&ctx, &or_die, &p](const char* const b) {
    const auto id = p.getBoundaryIdentifier(ctx, b) | or_die;
    const auto dofs = mfem_mgis::getElementsDegreesOfFreedomOnBoundary(p, id);
    auto F = mfem::Vector{};
    mfem_mgis::computeResultantForceOnBoundary(ctx, F, p, dofs) | or_die;
    return F[0];
  };
  const auto Fl = reaction("left");
  const auto Fr = reaction("right");
  if (std::abs(Fl + Fr) > 1e-10 * std::abs(Fl)) {
    mfem_mgis::getErrorStream()
        << "the reactions of the left and right "
        << "boundaries are not opposite (" << Fl << " vs " << Fr << ")\n";
    return false;
  }
  // reference value, in N/m. The converged value, computed on finer meshes,
  // is 3e-4 lower. It lies between the values given by 1D models of the
  // plate, with free (2.44e9 N/m) and blocked (4.27e9 N/m) upper and lower
  // boundaries.
  constexpr auto Fref = mfem_mgis::real{2.743179149139e9};
  if (std::abs(Fl - Fref) > 1e-8 * Fref) {
    mfem_mgis::getErrorStream() << "invalid reaction of the left boundary ("
                                << Fl << " vs " << Fref << ")\n";
    return false;
  }
  return true;
}  // end of checkReactions

/*
 * This test models a 2D plate of length 1 in plane strain clamped on the left
 * and right boundaries and submitted to a parabolic temperature profile along
 * the x-axis:
 *
 * - the temperature profile is minimal on the left and right boundaries
 * - the temperature profile is maximal for x = 0.5
 *
 * This example shows how to define an external state variable using an
 * analytical profile.
 *
 * The test checks that:
 *
 * - the temperature seen by the behaviour at each integration point is the
 *   imposed profile,
 * - the horizontal reactions of the left and right boundaries are opposite,
 *   as imposed by the equilibrium,
 * - the horizontal reaction of the left boundary matches its reference value.
 */
int main(int argc, char** argv) {
  auto ctx = mgis::Context{};
  auto or_die = ctx.getFatalFailureHandler();
  //
  static constexpr const auto parallel = false;
  // options treatment
  mfem_mgis::initialize(argc, argv);
  // building the non linear problem
  auto problem =
      mfem_mgis::construct<mfem_mgis::NonLinearEvolutionProblem>(
          ctx,
          mfem_mgis::Parameters{
              {"MeshFileName", "./square.msh"},
              {"Materials", mfem_mgis::Parameters{{"plate", 1}}},
              {"Boundaries", mfem_mgis::Parameters{{"left", 2}, {"right", 4}}},
              {"FiniteElementFamily", "H1"},
              {"FiniteElementOrder", 2},
              {"UnknownsSize", 2},
              {"Hypothesis", "PlaneStrain"},
              {"Parallel", parallel}}) |
      or_die;
  // materials
  problem.addBehaviourIntegrator(ctx, "Mechanics", "plate",
                                 "./src/libBehaviour.so",
                                 "IsotropicLinearThermoElasticity") |
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
        return getTemperature(x);
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
  problem.solve(ctx, 0, 1) | or_die;
  problem.executePostProcessings(ctx, 0, 1) | or_die;
  // hydrostatic pressure at the integration points
  const auto& cm1 = static_cast<const mfem_mgis::Material&>(m1);
  const auto pr =
      getInternalStateVariable(ctx, cm1, "HydrostaticPressure") | or_die;
  auto output = std::ofstream("HydrostaticPressure.txt");
  forEachIntegrationPoint<parallel>(
      pr, [&output](const mfem_mgis::real x, const mfem_mgis::real y,
                    const mfem_mgis::real v) {
        output << x << " " << y << " " << v << '\n';
      });
  // checks
  const auto T =
      getInternalStateVariable(ctx, cm1, "ImposedTemperature") | or_die;
  const auto temperature = checkImposedTemperature<parallel>(T);
  const auto reactions = checkReactions(ctx, problem);
  return (temperature && reactions) ? EXIT_SUCCESS : EXIT_FAILURE;
}
