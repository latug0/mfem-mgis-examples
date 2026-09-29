/*!
 * \file   CheckResultantForce.hxx
 * \brief  comparison of the resultant force on the upper boundary to
 * reference values
 */

#ifndef LIB_SSNA303_3D_CHECKRESULTANTFORCE_HXX
#define LIB_SSNA303_3D_CHECKRESULTANTFORCE_HXX

#include <cmath>
#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include <iostream>
#include "MFEMMGIS/Config.hxx"

//! \brief vertical component of the resultant force at a given time
struct VerticalForce {
  //! \brief time
  mfem_mgis::real t;
  //! \brief vertical component of the resultant force
  mfem_mgis::real fy;
};

/*!
 * \return the vertical component of the resultant force written by the
 * `ComputeResultantForceOnBoundary` post-processing at each time, or an empty
 * vector if the file can't be read
 * \param[in] f: file name
 */
inline std::vector<VerticalForce> readVerticalForce(const std::string& f) {
  auto forces = std::vector<VerticalForce>{};
  auto in = std::ifstream(f);
  auto line = std::string{};
  while (std::getline(in, line)) {
    if ((line.empty()) || (line[0] == '#')) {
      continue;
    }
    auto fx = mfem_mgis::real{};
    auto v = VerticalForce{};
    if (!(std::istringstream(line) >> v.t >> fx >> v.fy)) {
      std::cerr << "invalid line '" << line << "' in '" << f << "'\n";
      return {};
    }
    forces.push_back(v);
  }
  return forces;
}  // end of readVerticalForce

/*!
 * \return if the vertical component of the resultant force matches the
 * reference values. The computed times must be the first times of the
 * reference file, so that the beginning of a loading can be compared to the
 * reference values of the complete loading.
 * \param[in] f: file written by the `ComputeResultantForceOnBoundary`
 * post-processing
 * \param[in] r: reference file
 * \param[in] eps: relative tolerance. The default value is above the rounding
 * of the forces, which are written with 6 significant digits.
 */
inline bool checkVerticalForce(const std::string& f,
                               const std::string& r,
                               const mfem_mgis::real eps = 1e-4) {
  const auto values = readVerticalForce(f);
  const auto references = readVerticalForce(r);
  if (values.empty()) {
    std::cerr << "no value read in '" << f << "'\n";
    return false;
  }
  if (references.empty()) {
    std::cerr << "no value read in '" << r << "'\n";
    return false;
  }
  if (values.size() > references.size()) {
    std::cerr << "'" << f << "' has more values than '" << r << "'\n";
    return false;
  }
  for (std::size_t i = 0; i != values.size(); ++i) {
    // the times are written with 6 significant digits
    if (std::abs(values[i].t - references[i].t) >
        1e-6 * std::abs(references[i].t)) {
      std::cerr << "invalid time at time step " << i + 1 << " (" << values[i].t
                << " vs " << references[i].t << ")\n";
      return false;
    }
    if (std::abs(values[i].fy - references[i].fy) >
        eps * std::abs(references[i].fy)) {
      std::cerr << "invalid vertical force at time step " << i + 1 << " ("
                << values[i].fy << " vs " << references[i].fy << ")\n";
      return false;
    }
  }
  return true;
}  // end of checkVerticalForce

#endif /* LIB_SSNA303_3D_CHECKRESULTANTFORCE_HXX */
