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

constexpr RealType extendLookUpTable(RealType overlapAreaRatio) {
  if (overlapAreaRatio > 1.0) {
    return 5.0;
  }
  if (overlapAreaRatio > 0.5) {
    return 4.5;
  }
  if (overlapAreaRatio > 0.3) {
    return 4.2;
  }
  if (overlapAreaRatio > 0.2) {
    return 4.0;
  }
  if (overlapAreaRatio > 0.15) {
    return 3.9;
  }
  if (overlapAreaRatio > 0.1) {
    return 3.8;
  }
  if (overlapAreaRatio > 0.8) {
    return 3.7;
  }
  if (overlapAreaRatio > 0.5) {
    return 3.6;
  }
  if (overlapAreaRatio > 0.3) {
    return 3.4;
  }
  if (overlapAreaRatio > 0.2) {
    return 3.2;
  }
  if (overlapAreaRatio > 0.15) {
    return 3.0;
  }
  if (overlapAreaRatio > 0.1) {
    return 2.7;
  }
  if (overlapAreaRatio > 0.05) {
    return 2.4;
  }
  if (overlapAreaRatio > 0.04) {
    return 2.2;
  }
  if (overlapAreaRatio > 0.03) {
    return 2.0;
  }
  if (overlapAreaRatio > 0.02) {
    return 1.8;
  }
  if (overlapAreaRatio > 0.01) {
    return 1.6;
  }
  if (overlapAreaRatio > 0.005) {
    return 1.4;
  }
  return 1.2;
}

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
      const num_type sigmaX = static_cast<num_type>(rect->xLen()/ 2);
      const num_type muY = static_cast<num_type>(rect->center().y());
      const num_type sigmaY = static_cast<num_type>(rect->yLen()/ 2);
      const num_type normalize =  (2 * 3.14 * sigmaX * sigmaY);
      vec.emplace_back(BivariateGaussianParameters<num_type>(muX, sigmaX, muY, sigmaY, normalize));
    }

  }
}

PROJECT_NAMESPACE_END
