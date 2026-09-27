#include "../headers/BoundaryConditions.hxx"
#include <memory>
#include <vector>

void apply_boundary_conditions(
			       mfem_mgis::attributes::MayAbort,
			       mfem_mgis::Context& ctx,
    mfem_mgis::NonLinearEvolutionProblem& heat_transfer,
    mfem_mgis::NonLinearEvolutionProblem& mechanics,
    const TestParameters& p,
    const std::function<double(double)>& power_history,
    mfem::GridFunction* u_mech) {
  auto or_die = ctx.getFatalFailureHandler();
  // Constrain X and Y displacements on the upper surface of the stiffeners.
  // Only the Z displacement remains free.
  for (int j = 0; j <= 1; j++) {
    mechanics.addBoundaryCondition(ctx,
				   mfem_mgis::make_unique<mfem_mgis::UniformDirichletBoundaryCondition>(ctx,
            mechanics.getFiniteElementDiscretizationPointer(), 9, j,
            [](const auto) noexcept { return 0.0; })|or_die)|or_die;

    mechanics.addBoundaryCondition(ctx,
				   mfem_mgis::make_unique<mfem_mgis::UniformDirichletBoundaryCondition>(ctx,
            mechanics.getFiniteElementDiscretizationPointer(), 10, j,
            [](const auto) noexcept { return 0.0; })|or_die)|or_die;
  }

  mechanics.addBoundaryCondition(ctx,
				 mfem_mgis::make_unique<mfem_mgis::UniformDirichletBoundaryCondition>(ctx,
          mechanics.getFiniteElementDiscretizationPointer(), 9, 2,
          [](const auto) noexcept { return 0.0; })|or_die)|or_die;

  for (const int surface_id : {5, 7, 8}) {
    // Apply the coolant pressure.
    mechanics.addBoundaryCondition(ctx,
				   mfem_mgis::make_unique<mfem_mgis::UniformImposedPressureBoundaryCondition>(ctx,
            mechanics.getFiniteElementDiscretizationPointer(), surface_id,
            [p](const auto) noexcept { return p.water_pressure; })|or_die)|or_die;

    // Convective heat transfer at the coolant interface.
    heat_transfer.addBoundaryCondition(ctx,
				       mfem_mgis::make_unique<mfem_mgis::RobinBC>(ctx,
										  heat_transfer.getFiniteElementDiscretizationPointer(), surface_id,
        p.h_conv, p.Te, u_mech)|or_die)|or_die;
  }

  // Volumetric heat generation in the fuel.
  heat_transfer.addBoundaryCondition(ctx, 
				     mfem_mgis::make_unique<mfem_mgis::UniformHeatSourceBoundaryCondition>(ctx,
          heat_transfer.getFiniteElementDiscretizationPointer(), 1,
          [power_history](const auto t) { return power_history(t); })|or_die)|or_die;
}
