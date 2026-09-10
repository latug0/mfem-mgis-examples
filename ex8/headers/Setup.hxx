#pragma once

#include <vector>
#include <memory>
#include <functional>

// #include "MFEMMGIS/Context.hxx"
#include "MFEMMGIS/NonLinearEvolutionProblem.hxx"
#include "MFEMMGIS/PointWiseModel.hxx"
#include "MFEMMGIS/Profiler.hxx"
#include "MGIS/Behaviour/MaterialStateManager.hxx"

struct TestParameters {
  const char* mesh_file =
      "../mesh/assemblage_hexa.msh";  // path to the mesh file
  const char* behaviour =
      "ConductionThermique";  // default thermal behaviour name
  const char* libraryALFENI =
      "src/libALFENI-generic.so";  // MFront ALFENI library
                                   // (behaviours/models)
  const char* libraryU3SI2 = "src/libU3SI2-generic.so";  // MFront U3SI2 library
                                                         // (behaviours/models)
  const char* solver_thermo =
      "HypreGMRES";  // linear solver used for the thermal problem
  const char* precond_thermo =
      "HypreBoomerAMG";  // preconditioner associated with the thermal solver
  const char* solver_meca =
      "HyprePCG";  // linear solver used for the mechanical problem
  const char* precond_meca =
      "HypreBoomerAMG";  // preconditioner associated with the mechanical solver
  int order = 1;         // finite element order
  bool parallel = true;  // enable parallel execution (MPI)
  bool debug = true;     // enable debug output/checks
  int refinement = 0;    // number of uniform mesh refinements
  int post_processing = 1;  // default value : disabled
  int verbosity_level = 0;  // default value : lower level

  // Physical properties
  double Ti = 293.15;           // initial temperature (K)
  double Te = 315.0;            // external/convection temperature (K)
  double source = 1e10;         // volumetric power source term
  double water_pressure = 0.0;  // imposed water pressure
  double duree = 1e5;           // total simulation duration
  int nbsteps = 1;              // number of time steps
  double h_conv = 5e4;          // thermal convection coefficient
};

struct GaussFieldStorage {
  std::shared_ptr<std::vector<double>>
      T_s0;  // temperature at the beginning of the time step
  std::shared_ptr<std::vector<double>>
      T_s1;  // temperature at the end of the time step

  std::shared_ptr<std::vector<double>>
      Pow_s0_sw;  // power density (swelling model) at the beginning of the time
                  // step
  std::shared_ptr<std::vector<double>>
      Pow_s1_sw;  // power density (swelling model) at the end of the time step

  std::shared_ptr<std::vector<double>>
      Pow_s0_mmc;  // power density (mechanical material) at the beginning of
                   // the time step
  std::shared_ptr<std::vector<double>>
      Pow_s1_mmc;  // power density (mechanical material) at the end of the time
                   // step

  std::shared_ptr<std::vector<double>>
      Pow_s0_th;  // power density (thermal material) at the beginning of the
                  // time step
  std::shared_ptr<std::vector<double>>
      Pow_s1_th;  // power density (thermal material) at the end of the time
                  // step
};

struct SetupPropertiesResult {
  std::vector<GaussFieldStorage>
      fields;  // Gauss point field storage, one entry per material
  std::shared_ptr<mfem_mgis::PointWiseModel>
      swelling_model;  // swelling model attached to material 1
};

/*!
 * \brief Configures materials, models, and field storages
 * \param[in,out] ctx: execution context, used for timing sections and error
 * handling
 * \param[in] p: test parameters (mesh, libraries, physical properties, ...)
 * \param[in,out] heat_transfer: non-linear heat transfer problem to configure
 * \param[in,out] mechanics: non-linear mechanical problem to configure
 * \param[in] power_history: function giving the power density as a function of
 * time
 * \return the field storages and the swelling model created during setup
 */
inline SetupPropertiesResult setup_properties(
    mgis::Context& ctx,
    const TestParameters& p,
    mfem_mgis::NonLinearEvolutionProblem& heat_transfer,
    mfem_mgis::NonLinearEvolutionProblem& mechanics,
    const std::function<double(double)>& power_history) {
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
    auto& m_th = heat_transfer.getMaterial(mat_id);
    auto& m_mc = mechanics.getMaterial(mat_id);

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

      auto& m_sw = sw_model->getMaterial();
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