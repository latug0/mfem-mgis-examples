#include <cmath>
#include <iostream>
#include <memory>

#include "MFEMMGIS/AbstractBoundaryCondition.hxx"
#include "MFEMMGIS/Material.hxx"
#include "MFEMMGIS/NonLinearEvolutionProblem.hxx"
#include "MFEMMGIS/Parameters.hxx"
#include "MFEMMGIS/Profiler.hxx"
#include "MFEMMGIS/UniformDirichletBoundaryCondition.hxx"
#include "mfem/fem/fe_coll.hpp"
#include "mfem/fem/nonlininteg.hpp"
#include "mfem/linalg/vector.hpp"
#include "mfem/mesh/mesh.hpp"

#include "../headers/RobinBC.hxx"

static constexpr int TAG_LEFT = 1;
static constexpr int TAG_RIGHT = 2;
static constexpr int TAG_MAT = 1;

int main(int argc, char **argv) {
  mfem_mgis::initialize(argc, argv);

  const double L = 1.0;
  const double lambda = 10.0;
  const double T0 = 500.0;
  const double T_inf = 300.0;
  const double h = 100.0;

  const double a = -h * (T0 - T_inf) / (lambda + h * L);
  const double TL_ref = T0 + a * L;
  std::cout << "Ref : T(L) = " << TL_ref << " K\n\n";

  int nx = 10;
  int order = 2;
  mfem::OptionsParser args(argc, argv);
  args.AddOption(&nx, "-nx", "--nx", "Elements selon x");
  args.AddOption(&order, "-o", "--order", "Ordre EF");
  args.Parse();
  if (!args.Good()) {
    args.PrintUsage(std::cout);
    return EXIT_FAILURE;
  }

  auto *m = new mfem::Mesh(mfem::Mesh::MakeCartesian3D(
      nx, 2, 2, mfem::Element::HEXAHEDRON, L, 1.0, 1.0));
  
  for (int i = 0; i < m->GetNBE(); ++i) {
    mfem::ElementTransformation *tr = m->GetBdrElementTransformation(i);
    mfem::IntegrationPoint ip;
    ip.Set3(0.5, 0.5, 0.0);
    tr->SetIntPoint(&ip);
    mfem::Vector center(3);
    tr->Transform(ip, center);
    if (center[0] < 1e-10)
      m->GetBdrElement(i)->SetAttribute(TAG_LEFT);
    else if (center[0] > L - 1e-10)
      m->GetBdrElement(i)->SetAttribute(TAG_RIGHT);
    else
      m->GetBdrElement(i)->SetAttribute(3);
  }
  m->SetAttributes();
  auto mesh = std::shared_ptr<mfem::Mesh>(m);

  auto fec = std::make_shared<mfem::H1_FECollection>(order, 3);
  auto fed = std::make_shared<mfem_mgis::FiniteElementDiscretization>(
      mesh, fec, mfem_mgis::size_type{1});

  mfem_mgis::NonLinearEvolutionProblem pb(
      fed, mgis::behaviour::Hypothesis::TRIDIMENSIONAL, {});

  pb.getUnknowns(mfem_mgis::bts) = T0;
  pb.getUnknowns(mfem_mgis::ets) = T0;

  pb.addBehaviourIntegrator("StationaryNonLinearHeatTransfer", TAG_MAT,
                            "src/libBehaviour.so",
                            "StationaryLinearHeatTransfer");

  auto &mat = pb.getMaterial(TAG_MAT);
  std::vector<mgis::real> T_field(mat.n, 0.0);
  mgis::behaviour::setExternalStateVariable(
      mat.s0, "Temperature", T_field,
      mgis::behaviour::MaterialStateManager::EXTERNAL_STORAGE,
      mgis::behaviour::MaterialStateManager::NOUPDATE);
  mgis::behaviour::setExternalStateVariable(
      mat.s1, "Temperature", T_field,
      mgis::behaviour::MaterialStateManager::EXTERNAL_STORAGE,
      mgis::behaviour::MaterialStateManager::NOUPDATE);

  pb.addBoundaryCondition(
      std::make_unique<mfem_mgis::UniformDirichletBoundaryCondition>(
          fed, TAG_LEFT, 0, [T0](const auto) noexcept { return T0; }));


  pb.addBoundaryCondition(
      std::make_unique<mfem_mgis::RobinBC>(fed, TAG_RIGHT, h, T_inf, nullptr));

  pb.setLinearSolver("UMFPackSolver", {});
  pb.setSolverParameters({{"VerbosityLevel", 1},
                          {"RelativeTolerance", 1e-10},
                          {"AbsoluteTolerance", 0},
                          {"MaximumNumberOfIterations", 100}});

  if (!pb.solve(0.0, 1.0))
    mfem_mgis::abort("Non-convergence");

  auto &fes = fed->getFiniteElementSpace<false>();
  mfem::GridFunction T_exact(&fes);
  mfem::FunctionCoefficient exact_coef(
      [&](const mfem::Vector &x) { return T0 + a * x[0]; });
  T_exact.ProjectCoefficient(exact_coef);

  pb.addPostProcessing("ParaviewExportResults",
                       {{"OutputFileName", "RobinTest"}});
  pb.executePostProcessings(0.0, 1.0);

  mfem::ParaViewDataCollection pv("RobinTest_exact", mesh.get());
  pv.RegisterField("T_exact", &T_exact);
  pv.SetCycle(0);
  pv.SetTime(0.0);
  pv.Save();

  const auto &T_sol = pb.getUnknowns(mfem_mgis::ets);
  double TL_num = 0.0;
  int n = 0;
  for (int i = 0; i < fes.GetNBE(); ++i) {
    if (mesh->GetBdrAttribute(i) != TAG_RIGHT)
      continue;
    mfem::Array<int> vdofs;
    fes.GetBdrElementVDofs(i, vdofs);
    for (int j = 0; j < vdofs.Size(); ++j) {
      TL_num += T_sol[vdofs[j]];
      ++n;
    }
  }
  TL_num /= n;

  std::cout << "Num : T(L) = " << TL_num << " K\n";
  std::cout << "Erreur relative     = " << std::abs(TL_num - TL_ref) / TL_ref << " \n";

  return EXIT_SUCCESS;
}