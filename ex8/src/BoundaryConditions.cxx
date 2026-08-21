#include "../headers/BoundaryConditions.hxx"
#include <memory>
#include <vector>

void apply_boundary_conditions(
    mfem_mgis::NonLinearEvolutionProblem &heat_transfer,
    mfem_mgis::NonLinearEvolutionProblem &mechanics,
    const TestParameters& p,
    const std::function<double(double)>& power_history,
    mfem::GridFunction* u_mech) {

  // Constrain X and Y displacements on the upper surface of the stiffeners.
  // Only the Z displacement remains free.
  for (int j = 0; j <= 1; j++) {
    mechanics.addBoundaryCondition(std::make_unique<mfem_mgis::UniformDirichletBoundaryCondition>(
        mechanics.getFiniteElementDiscretizationPointer(), 9, j, [](const auto) noexcept { return 0.0; }));

    mechanics.addBoundaryCondition(std::make_unique<mfem_mgis::UniformDirichletBoundaryCondition>(
        mechanics.getFiniteElementDiscretizationPointer(), 10, j, [](const auto) noexcept { return 0.0; }));
  }

  mechanics.addBoundaryCondition(std::make_unique<mfem_mgis::UniformDirichletBoundaryCondition>(
      mechanics.getFiniteElementDiscretizationPointer(), 9, 2, [](const auto) noexcept { return 0.0; }));

  for (const int surface_id : {5, 7, 8}) {
    // Apply the coolant pressure.
    mechanics.addBoundaryCondition(std::make_unique<mfem_mgis::UniformImposedPressureBoundaryCondition>(
        mechanics.getFiniteElementDiscretizationPointer(), surface_id,
        [p](const auto) noexcept { return p.water_pressure; }));

    // Convective heat transfer at the coolant interface.
    heat_transfer.addBoundaryCondition(std::make_unique<mfem_mgis::RobinBC>(
        heat_transfer.getFiniteElementDiscretizationPointer(),
        surface_id, p.h_conv, p.Te, u_mech));
  }

  // Volumetric heat generation in the fuel.
  heat_transfer.addBoundaryCondition(std::make_unique<mfem_mgis::UniformHeatSourceBoundaryCondition>(
      heat_transfer.getFiniteElementDiscretizationPointer(), 1,
      [power_history](const auto t) { return power_history(t); }));
}