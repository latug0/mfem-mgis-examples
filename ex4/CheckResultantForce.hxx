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

/*!
 * \return the vertical component of the resultant force written by the
 * `ComputeResultantForceOnBoundary` post-processing
 * \param[in] f: file name
 */
inline std::vector<mfem_mgis::real> readVerticalForce(const std::string& f) {
  auto fy = std::vector<mfem_mgis::real>{};
  auto in = std::ifstream(f);
  auto line = std::string{};
  while (std::getline(in, line)) {
    if ((line.empty()) || (line[0] == '#')) {
      continue;
    }
    auto t = mfem_mgis::real{};
    auto fx = mfem_mgis::real{};
    auto v = mfem_mgis::real{};
    std::istringstream(line) >> t >> fx >> v;
    fy.push_back(v);
  }
  return fy;
}  // end of readVerticalForce

/*!
 * \return if the vertical component of the resultant force matches the
 * reference values
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
  if ((references.empty()) || (values.size() != references.size())) {
    std::cerr << "'" << f << "' and '" << r
              << "' do not have the same number of values\n";
    return false;
  }
  for (std::size_t i = 0; i != values.size(); ++i) {
    if (std::abs(values[i] - references[i]) > eps * std::abs(references[i])) {
      std::cerr << "invalid vertical force at time step " << i + 1 << " ("
                << values[i] << " vs " << references[i] << ")\n";
      return false;
    }
  }
  return true;
}  // end of checkVerticalForce

#endif /* LIB_SSNA303_3D_CHECKRESULTANTFORCE_HXX */
