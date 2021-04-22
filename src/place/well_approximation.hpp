/**
 * @file well_approximation.h
 * @brief Approximate the well
 * @author Keren Zhu
 * @date 04/20/2021
 */

#pragma once

#include "global/global.h"
#include "util/box.hpp"

PROJECT_NAMESPACE_BEGIN

template<typename num_type>
struct BivariateGaussianParameters {
  BivariateGaussianParameters(num_type muX_, num_type sigmaX_, num_type muY_, num_type sigmaY_, num_type normalize_)
    : muX(muX_), sigmaX(sigmaX_), muY(muY_), sigmaY(sigmaY_), normalize(normalize_) {}
  num_type muX;
  num_type sigmaX;
  num_type muY;
  num_type sigmaY;
  num_type normalize; // mutliply this
};

template<typename num_type>
class BivariateGaussianWellApproximationGenerator {
  public:
    BivariateGaussianWellApproximationGenerator(
        const std::function<IndexType(void)> &getNumWells,
        const std::function<IndexType(IndexType)> &getNumRectsInWell,
        const std::function<Box<num_type>*(IndexType, IndexType)> &getRect)
      : _getNumWells(getNumWells), _getNumRectsInWell(getNumRectsInWell), _getRect(getRect)
    {}
    void generate(std::vector<BivariateGaussianParameters<num_type>> &vec);
  private:
    std::function<IndexType(void)> _getNumWells; ///< Get the number of wells
    std::function<IndexType(IndexType)> _getNumRectsInWell; ///< How many rectangles in a well
    std::function<Box<num_type>*(IndexType, IndexType)> _getRect; ///< rect = _getRect(wellIdx, rectInWellIdx)
};

template<typename num_type>
inline void BivariateGaussianWellApproximationGenerator<num_type>::generate(std::vector<BivariateGaussianParameters<num_type>> &vec) {
  vec.clear();
  IndexType numWell = _getNumWells();
  for (IndexType wellIdx = 0; wellIdx < numWell; ++wellIdx) {
    IndexType numRects = _getNumRectsInWell(wellIdx);
    // Note that theorically, the PDF value can exceed 1 given a set of very small disjointed wells
    // However, since well polygon is filtered by area in post-processing, the PDF is in generally controllable
    for (IndexType rectIdx = 0; rectIdx < numRects; ++rectIdx) {
      Box<num_type>*rect = _getRect(wellIdx, rectIdx);
      const num_type muX = static_cast<num_type>(rect->center().x());
      const num_type sigmaX = static_cast<num_type>(rect->xLen() / 2);
      const num_type muY = static_cast<num_type>(rect->center().y());
      const num_type sigmaY = static_cast<num_type>(rect->yLen() / 2);
      const num_type normalize =  (2 * 3.14 * sigmaX * sigmaY);
      vec.emplace_back(BivariateGaussianParameters<num_type>(muX, sigmaX, muY, sigmaY, normalize));
    }

  }
}

PROJECT_NAMESPACE_END
