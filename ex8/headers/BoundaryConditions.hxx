#pragma once
#include <functional>
#include <vector>
#include <memory>
#include "MFEMMGIS/ModelBase.hxx"
#include "MFEMMGIS/TimeStep.hxx"
#include "MFEMMGIS/NonLinearEvolutionProblem.hxx"
#include "MFEMMGIS/UniformDirichletBoundaryCondition.hxx"
#include "MFEMMGIS/UniformHeatSourceBoundaryCondition.hxx"
#include "MFEMMGIS/UniformImposedPressureBoundaryCondition.hxx"

#include "Setup.hxx"
#include "RobinBC.hxx"

void apply_boundary_conditions(
    mfem_mgis::NonLinearEvolutionProblem& heat_transfer,
    mfem_mgis::NonLinearEvolutionProblem& mechanics,
    const TestParameters& p,
    const std::function<double(double)>& power_history,
    mfem::GridFunction* u_mech
);

// Model used to manually update the source term before passing it to the
// constitutive model responsible for swelling.
class FieldUpdaterModel : public mfem_mgis::ModelBase {
private:
    std::shared_ptr<std::vector<mgis::real>> Pow_s1;
    std::function<double(double)> power_history;
    std::string name;

public:
    FieldUpdaterModel(const mfem_mgis::MeshDiscretization& mesh,
                      std::shared_ptr<std::vector<mgis::real>> pow_array,
                      std::function<double(double)> history_func)
        : mfem_mgis::ModelBase(mesh), 
          Pow_s1(pow_array), 
          power_history(history_func), 
          name("FieldUpdater")
    {}

    [[nodiscard]] std::string getName() const noexcept override {
        return name;
    }

    std::pair<mfem_mgis::ExitStatus, std::optional<mfem_mgis::ComputeNextStateOutput>>
    computeNextState(mfem_mgis::Context& ctx, const mfem_mgis::TimeStep& ts) noexcept override {
        
        double current_time = ts.end;

        double current_power = power_history(current_time);
        
        for (auto& val : *Pow_s1) {
            val = current_power;
        }

        return {mfem_mgis::ExitStatus::success, mfem_mgis::ComputeNextStateOutput{}};
    }
};