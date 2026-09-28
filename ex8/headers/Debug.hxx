#pragma once

#include <map>
#include <cmath>
#include <limits>
#include <string>
#include <vector>
#include <fstream>
#include <ostream>
#include <sstream>
#include <algorithm>

#include <mpi.h>

#include "mfem/fem/pgridfunc.hpp"
#include "MFEMMGIS/Config.hxx"
#include "MFEMMGIS/NonLinearEvolutionProblem.hxx"

#include "Setup.hxx"

//! \brief minimum, maximum and mean values of a field
struct FieldStatistics {
  double min = 0;
  double max = 0;
  double mean = 0;
};

/*!
 * \return the statistics of values distributed over all processes
 * \param[in] values: values of the current process
 * \param[in] comm: communicator
 */
inline FieldStatistics computeFieldStatistics(const std::vector<double>& values,
                                              MPI_Comm comm) {
  auto s = FieldStatistics{std::numeric_limits<double>::max(),
                           -std::numeric_limits<double>::max(), 0};
  for (const auto v : values) {
    s.min = std::min(s.min, v);
    s.max = std::max(s.max, v);
    s.mean += v;
  }
  auto n = static_cast<long long>(values.size());
  MPI_Allreduce(MPI_IN_PLACE, &s.min, 1, MPI_DOUBLE, MPI_MIN, comm);
  MPI_Allreduce(MPI_IN_PLACE, &s.max, 1, MPI_DOUBLE, MPI_MAX, comm);
  MPI_Allreduce(MPI_IN_PLACE, &s.mean, 1, MPI_DOUBLE, MPI_SUM, comm);
  MPI_Allreduce(MPI_IN_PLACE, &n, 1, MPI_LONG_LONG, MPI_SUM, comm);
  s.mean /= n;
  return s;
}  // end of computeFieldStatistics

/*!
 * \return the values of a nodal field at the nodes owned by the current
 * process, or the norm of the field at these nodes for a vector field
 * \param[in] f: nodal field
 */
inline std::vector<double> getOwnedNodalValues(const mfem::ParGridFunction& f) {
  const auto& fes = *(f.ParFESpace());
  auto values = std::vector<double>{};
  for (int i = 0; i != fes.GetNDofs(); ++i) {
    // a node shared between processes is counted by its owner only
    if (fes.GetLocalTDofNumber(fes.DofToVDof(i, 0)) < 0) {
      continue;
    }
    if (fes.GetVDim() == 1) {
      values.push_back(f(fes.DofToVDof(i, 0)));
    } else {
      auto n2 = 0.0;
      for (int d = 0; d != fes.GetVDim(); ++d) {
        const auto c = f(fes.DofToVDof(i, d));
        n2 += c * c;
      }
      values.push_back(std::sqrt(n2));
    }
  }
  return values;
}  // end of getOwnedNodalValues

/*!
 * \return the statistics of the temperature, of the norm of the
 * displacement, of the swelling and of the power density in the fuel at the
 * end of the simulation
 * \param[in] heat_transfer: heat transfer problem
 * \param[in] mechanics: mechanical problem
 * \param[in] setup: field storages and swelling model
 */
inline std::map<std::string, FieldStatistics> computePhysicsStatistics(
    mfem_mgis::NonLinearEvolutionProblem& heat_transfer,
    mfem_mgis::NonLinearEvolutionProblem& mechanics,
    const SetupPropertiesResult& setup) {
  auto stats = std::map<std::string, FieldStatistics>{};
  auto& thermal_fes = heat_transfer.getFiniteElementDiscretization()
                          .getFiniteElementSpace<true>();
  auto& mechanical_fes =
      mechanics.getFiniteElementDiscretization().getFiniteElementSpace<true>();
  const auto comm = thermal_fes.GetComm();
  // nodal fields, synchronized with the nodes shared between processes
  mfem::ParGridFunction T(&thermal_fes);
  T.SetFromTrueDofs(heat_transfer.getUnknowns(mfem_mgis::bts));
  stats["Temperature"] = computeFieldStatistics(getOwnedNodalValues(T), comm);
  mfem::ParGridFunction U(&mechanical_fes);
  U.SetFromTrueDofs(mechanics.getUnknowns(mfem_mgis::bts));
  stats["DisplacementNorm"] =
      computeFieldStatistics(getOwnedNodalValues(U), comm);
  // fields at the integration points of the fuel
  const auto& swelling =
      setup.swelling_model->getMaterial().s1.internal_state_variables;
  stats["Swelling"] = computeFieldStatistics(
      std::vector<double>(swelling.begin(), swelling.end()), comm);
  stats["PowerDensity"] =
      computeFieldStatistics(*(setup.fields[0].Pow_s1_sw), comm);
  return stats;
}  // end of computePhysicsStatistics

/*!
 * \brief print the statistics of the fields, one line per field, in the
 * format of the reference files
 * \param[in] os: output stream
 * \param[in] stats: statistics
 */
inline void printPhysicsStatistics(
    std::ostream& os, const std::map<std::string, FieldStatistics>& stats) {
  const auto precision = os.precision(15);
  os << "# field, minimum, maximum and mean values\n";
  for (const auto& [name, s] : stats) {
    os << name << ' ' << s.min << ' ' << s.max << ' ' << s.mean << '\n';
  }
  os.precision(precision);
}  // end of printPhysicsStatistics

/*!
 * \return if the swelling matches its exact value
 *
 * The swelling rate is proportional to the power density. The latter is
 * linear over each time step since the end of the power ramp is a time step
 * boundary, so the swelling computed by the `U3SI2_SolidSwelling` model with
 * the mean power density over each time step is exact.
 *
 * \param[in] s: statistics of the swelling
 * \param[in] p: parameters of the simulation
 */
inline bool checkSwelling(const FieldStatistics& s, const TestParameters& p) {
  // swelling per unit of energy released, see U3SI2_Swelling.mfront (6.2e-29
  // per fission, 200 MeV per fission)
  constexpr auto A = 6.2e-29 / (200 * 1.60218e-13);
  // energy released per unit of volume at the end of the simulation
  const auto t = p.duree;
  const auto E = (t <= p.t_ramp) ? p.source * t * t / (2 * p.t_ramp)
                                 : p.source * (t - p.t_ramp / 2);
  const auto S = A * E;
  if ((std::abs(s.min - S) > 1e-10 * S) || (std::abs(s.max - S) > 1e-10 * S)) {
    mfem_mgis::getErrorStream()
        << "the swelling does not match its exact value " << S << " (from "
        << s.min << " to " << s.max << ")\n";
    return false;
  }
  return true;
}  // end of checkSwelling

/*!
 * \return if the statistics of the fields match the reference values
 * \param[in] stats: statistics
 * \param[in] f: reference file, each line gives the name of a field followed
 * by its minimum, maximum and mean values, as printed by
 * `printPhysicsStatistics`
 */
inline bool checkPhysicsStatistics(
    const std::map<std::string, FieldStatistics>& stats, const std::string& f) {
  // relative tolerance, the values of a field are compared to its largest
  // absolute value
  constexpr auto eps = 1e-6;
  auto in = std::ifstream(f);
  if (!in) {
    mfem_mgis::getErrorStream() << "can't open file '" << f << "'\n";
    return false;
  }
  auto nfields = 0;
  auto line = std::string{};
  while (std::getline(in, line)) {
    if ((line.empty()) || (line[0] == '#')) {
      continue;
    }
    auto is = std::istringstream(line);
    auto name = std::string{};
    auto r = FieldStatistics{};
    if (!(is >> name >> r.min >> r.max >> r.mean)) {
      mfem_mgis::getErrorStream()
          << "invalid line '" << line << "' in '" << f << "'\n";
      return false;
    }
    const auto ps = stats.find(name);
    if (ps == stats.end()) {
      mfem_mgis::getErrorStream() << "unknown field '" << name << "'\n";
      return false;
    }
    const auto& s = ps->second;
    const auto scale = std::max(std::abs(r.min), std::abs(r.max));
    if ((std::abs(s.min - r.min) > eps * scale) ||
        (std::abs(s.max - r.max) > eps * scale) ||
        (std::abs(s.mean - r.mean) > eps * scale)) {
      mfem_mgis::getErrorStream()
          << "invalid statistics of field '" << name << "' (" << s.min << ' '
          << s.max << ' ' << s.mean << " vs " << r.min << ' ' << r.max << ' '
          << r.mean << ")\n";
      return false;
    }
    ++nfields;
  }
  if (nfields == 0) {
    mfem_mgis::getErrorStream() << "no reference values in '" << f << "'\n";
    return false;
  }
  return true;
}  // end of checkPhysicsStatistics
