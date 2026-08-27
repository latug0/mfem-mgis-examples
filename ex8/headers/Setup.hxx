#pragma once

#include <vector>
#include <memory>
#include <functional>

//#include "MFEMMGIS/Context.hxx"
#include "MFEMMGIS/NonLinearEvolutionProblem.hxx"
#include "MFEMMGIS/PointWiseModel.hxx"
#include "MFEMMGIS/Profiler.hxx"
#include "MGIS/Behaviour/MaterialStateManager.hxx"

struct TestParameters {
  const char *mesh_file = "../mesh/assemblage_hexa.msh";
  const char *behaviour = "ConductionThermique";
  const char *libraryALFENI = "src/libALFENI-generic.so";
  const char *libraryU3SI2 = "src/libU3SI2-generic.so";
  const char *solver_thermo = "HypreGMRES";
  const char *precond_thermo = "HypreBoomerAMG";
  const char *solver_meca = "HyprePCG";
  const char *precond_meca = "HypreBoomerAMG";
  int order = 1;
  bool parallel = true;
  bool debug = true;
  int refinement = 0;
  int post_processing = 1;  // default value : disabled
  int verbosity_level = 0;  // default value : lower level

  // Physical properties
  double Ti = 293.15;
  double Te = 315.0;
  double source = 1e10;
  double water_pressure = 0.0;
  double duree = 1e5;
  int nbsteps = 1;
  double h_conv = 5e4;
};

struct GaussFieldStorage {
  std::shared_ptr<std::vector<double>> T_s0;
  std::shared_ptr<std::vector<double>> T_s1;

  std::shared_ptr<std::vector<double>> Pow_s0_sw;
  std::shared_ptr<std::vector<double>> Pow_s1_sw;

  std::shared_ptr<std::vector<double>> Pow_s0_mmc;
  std::shared_ptr<std::vector<double>> Pow_s1_mmc;

  std::shared_ptr<std::vector<double>> Pow_s0_th;
  std::shared_ptr<std::vector<double>> Pow_s1_th;
};

struct SetupPropertiesResult {
  std::vector<GaussFieldStorage> fields;
  std::shared_ptr<mfem_mgis::PointWiseModel> swelling_model;
};

/*!
 * \brief Configures materials, models, and field storages
 */
inline SetupPropertiesResult setup_properties(
    mgis::Context &ctx,
    const TestParameters &p,
    mfem_mgis::NonLinearEvolutionProblem &heat_transfer,
    mfem_mgis::NonLinearEvolutionProblem &mechanics,
    const std::function<double(double)> &power_history) {
  using namespace mfem_mgis;
  using namespace mgis::behaviour;
  // using namespace mgis::model;
  using real = mfem_mgis::real;

  CatchTimeSection(ctx, "set_mgis_stuff");
  auto or_die = ctx.getFatalFailureHandler();

  SetupPropertiesResult result;

  for (auto ts : {bts, ets}) {
    mechanics.getUnknowns(ts) = real{0};
    heat_transfer.getUnknowns(ts) = p.Ti;
  }

  mechanics.addBehaviourIntegrator("Mechanics", 1, p.libraryU3SI2,
                                   "U3SI2_NortonPRQ");
  mechanics.addBehaviourIntegrator("Mechanics", 2, p.libraryALFENI,
                                   "ALFENI_IsotropicLinearHardeningPlasticity");
  mechanics.addBehaviourIntegrator("Mechanics", 3, p.libraryALFENI,
                                   "ALFENI_IsotropicLinearHardeningPlasticity");

  heat_transfer.addBehaviourIntegrator("StationaryNonLinearHeatTransfer", 1,
                                       p.libraryU3SI2,
                                       "U3SI2_ThermiqueCouplee");
  heat_transfer.addBehaviourIntegrator("StationaryNonLinearHeatTransfer", 2,
                                       p.libraryALFENI,
                                       "ALFENI_ThermiqueCouplee");
  heat_transfer.addBehaviourIntegrator("StationaryNonLinearHeatTransfer", 3,
                                       p.libraryALFENI,
                                       "ALFENI_ThermiqueCouplee");

  for (const int mat_id : {1, 2, 3}) {
    auto &m_th = heat_transfer.getMaterial(mat_id);
    auto &m_mc = mechanics.getMaterial(mat_id);

    GaussFieldStorage storage;
    storage.T_s0 = std::make_shared<std::vector<mgis::real>>(m_mc.n, p.Ti);
    storage.T_s1 = std::make_shared<std::vector<mgis::real>>(m_mc.n, p.Ti);

    setExternalStateVariable(m_th.s0, "Temperature", *storage.T_s0,
                             MaterialStateManager::EXTERNAL_STORAGE,
                             MaterialStateManager::UPDATE);
    setExternalStateVariable(m_th.s1, "Temperature", *storage.T_s1,
                             MaterialStateManager::EXTERNAL_STORAGE,
                             MaterialStateManager::NOUPDATE);
    setExternalStateVariable(m_mc.s0, "Temperature", *storage.T_s0,
                             MaterialStateManager::EXTERNAL_STORAGE,
                             MaterialStateManager::NOUPDATE);
    setExternalStateVariable(m_mc.s1, "Temperature", *storage.T_s1,
                             MaterialStateManager::EXTERNAL_STORAGE,
                             MaterialStateManager::NOUPDATE);

    mgis::behaviour::setExternalStateVariable(
        m_th.s0, "DeformationGradient", m_mc.s0.gradients,
        MaterialStateManager::EXTERNAL_STORAGE, MaterialStateManager::NOUPDATE);
    mgis::behaviour::setExternalStateVariable(
        m_th.s1, "DeformationGradient", m_mc.s1.gradients,
        MaterialStateManager::EXTERNAL_STORAGE, MaterialStateManager::NOUPDATE);

    if (mat_id == 1) {
      auto sw_model = make_shared<PointWiseModel>(
                          ctx, m_mc.getPartialQuadratureSpacePointer(),
                          Parameters{{"Library", p.libraryU3SI2},
                                     {"Model", "U3SI2_SolidSwelling"},
                                     {"Hypothesis", "Tridimensional"}}) |
                      or_die;

      auto &m_sw = sw_model->getMaterial();
      const double initial_power = power_history(0.0);
      storage.Pow_s0_sw =
          std::make_shared<std::vector<mgis::real>>(m_sw.n, initial_power);
      storage.Pow_s1_sw =
          std::make_shared<std::vector<mgis::real>>(m_sw.n, initial_power);

      // Point all variables to the same array
      storage.Pow_s0_mmc = storage.Pow_s0_sw;
      storage.Pow_s1_mmc = storage.Pow_s1_sw;

      setExternalStateVariable(m_sw.s0, "Temperature", *storage.T_s0,
                               MaterialStateManager::EXTERNAL_STORAGE,
                               MaterialStateManager::NOUPDATE);
      setExternalStateVariable(m_sw.s1, "Temperature", *storage.T_s1,
                               MaterialStateManager::EXTERNAL_STORAGE,
                               MaterialStateManager::NOUPDATE);

      setExternalStateVariable(m_sw.s0, "PowerDensity", *storage.Pow_s0_sw,
                               MaterialStateManager::EXTERNAL_STORAGE,
                               MaterialStateManager::UPDATE);
      setExternalStateVariable(m_sw.s1, "PowerDensity", *storage.Pow_s1_sw,
                               MaterialStateManager::EXTERNAL_STORAGE,
                               MaterialStateManager::NOUPDATE);

      setExternalStateVariable(m_mc.s0, "Swelling",
                               m_sw.s0.internal_state_variables,
                               MaterialStateManager::EXTERNAL_STORAGE,
                               MaterialStateManager::NOUPDATE);
      setExternalStateVariable(m_mc.s1, "Swelling",
                               m_sw.s1.internal_state_variables,
                               MaterialStateManager::EXTERNAL_STORAGE,
                               MaterialStateManager::NOUPDATE);

      setExternalStateVariable(m_mc.s0, "PowerDensity", *storage.Pow_s0_mmc,
                               MaterialStateManager::EXTERNAL_STORAGE,
                               MaterialStateManager::NOUPDATE);
      setExternalStateVariable(m_mc.s1, "PowerDensity", *storage.Pow_s1_mmc,
                               MaterialStateManager::EXTERNAL_STORAGE,
                               MaterialStateManager::NOUPDATE);

      result.swelling_model = sw_model;
    }
    result.fields.push_back(std::move(storage));
  }
  return result;
}