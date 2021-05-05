/**
 * @file different.h
 * @brief The numerical differentiable concepts and implementations
 * @author Keren Zhu
 * @date 02/25/2020
 */

#ifndef IDEAPLACE_DIFFERENT_H_
#define IDEAPLACE_DIFFERENT_H_

#include "db/Database.h"
#include "well_approximation.hpp"

#include <complex>

PROJECT_NAMESPACE_BEGIN

namespace diff {

enum class OpEnumType { hpwl, ovl, oob, asym, cosine, fence };
struct placement_differentiable_concept {};

template <typename ConceptType> struct is_placement_differentiable_concept {
  typedef std::false_type is_placement_differentiable_concept_type;
};

template <>
struct is_placement_differentiable_concept<placement_differentiable_concept> {
  typedef std::true_type is_placement_differentiable_concept_type;
};

template <typename DifType> struct placement_differentiable_traits {
  typedef DifType different_type;
  typedef typename different_type::numerical_type numerical_type;
  typedef typename different_type::coordinate_type coordinate_type;

  static numerical_type evaluate(const different_type &dif) {
    return dif.evaluate();
  }

  static void accumlateGradient(const different_type &dif) {
    dif.accumlateGradient();
  }
};

/// @namespace IDEAPLACE::op
/// @brief namespace for the operators
namespace op {
/// @brief abs smooth function with respect to 0
/// @return the smoothed result
template <typename NumType>
constexpr NumType logSumExp0(NumType var, NumType alpha) {
  return alpha * log(exp(var / alpha) + 1);
}

/// @brief calculate the log sum exp, used to smooth min max function
/// @return the calculated log sum exp
template <typename NumType>
constexpr NumType logSumExp(NumType var1, NumType var2, NumType alpha) {
  return alpha * log(exp(var1 / alpha) + exp(var2 / alpha));
}

/// @brief The partial of LSE(0, var) with respect to var
/// @return the gradient of logSumExp0
template <typename NumType>
constexpr NumType gradLogSumExp0(NumType var, NumType alpha) {
  return exp(var / alpha) / (exp(var / alpha) + 1);
}

namespace _conv_details {
template <typename, typename, bool> struct _conv_t {};
template <typename LhsType, typename _RhsType>
struct _conv_t<LhsType, _RhsType, true> {
  static constexpr LhsType _conv(_RhsType rhs) { return rhs; }
};
template <typename LhsType, typename _RhsType>
struct _conv_t<LhsType, _RhsType, false> {
  static constexpr LhsType _conv(_RhsType rhs) {
    return static_cast<LhsType>(rhs);
  }
};
} // namespace _conv_details

/// @brief convert rhs type to lhs type
template <typename LhsType, typename _RhsType>
constexpr LhsType conv(_RhsType rhs) {
  return _conv_details::_conv_t<
      LhsType, _RhsType, std::is_same<LhsType, _RhsType>::value>::_conv(rhs);
}

inline double realerf(double in) {
  return std::erf(in);
}


inline double realerf(std::complex<double> in) {
  return std::erf(std::real(in));
}
inline float realerf(std::complex<float> in) {
  return std::erf(std::real(in));
}
}; // namespace op

/// @brief LSE-smoothed HPWL
template <typename NumType, typename CoordType> struct LseHpwlDifferentiable {
  typedef NumType numerical_type;
  typedef CoordType coordinate_type;

  LseHpwlDifferentiable(const std::function<NumType(void)> &getAlphaFunc,
                        const std::function<NumType(void)> &getLambdaFunc) {
    _getAlphaFunc = getAlphaFunc;
    _getLambdaFunc = getLambdaFunc;
  }

  void setGetVarFunc(
      const std::function<CoordType(IndexType, Orient2DType)> &getVarFunc) {
    _getVarFunc = getVarFunc;
  }
  void setAccumulateGradFunc(
      const std::function<void(NumType, IndexType, Orient2DType)> &func) {
    _accumulateGradFunc = func;
  }
  void setGetAlphaFunc(const std::function<NumType(void)> &getAlphaFunc) {
    _getAlphaFunc = getAlphaFunc;
  }

  void setVirtualPin(const CoordType &x, const CoordType &y) {
    _validVirtualPin = 1;
    _virtualPinX = x;
    _virtualPinY = y;
  }
  void removeVirtualPin() { _validVirtualPin = 0; }
  void addVar(IndexType cellIdx, const CoordType &offsetX,
              const CoordType &offsetY) {
    _cells.emplace_back(cellIdx);
    _offsetX.emplace_back(offsetX);
    _offsetY.emplace_back(offsetY);
  }
  void setWeight(const NumType &weight) { _weight = weight; }
  bool validHpwl() const { return _cells.size() + _validVirtualPin > 1; }

  NumType evaluate() const {
    if (!validHpwl()) {
      return 0;
    }
    std::array<NumType, 4> max_val = {0, 0, 0, 0}; // xmax xin ymax ymin
    auto alpha = _getAlphaFunc();
    auto lambda = _getLambdaFunc();
    NumType *pMax = &max_val.front();
    for (IndexType pinIdx = 0; pinIdx < _cells.size(); ++pinIdx) {
      NumType x = op::conv<NumType>(
          _getVarFunc(_cells[pinIdx], Orient2DType::HORIZONTAL) +
          _offsetX[pinIdx]);
      NumType y = op::conv<NumType>(
          _getVarFunc(_cells[pinIdx], Orient2DType::VERTICAL) +
          _offsetY[pinIdx]);
      pMax[0] += exp(x / alpha);
      pMax[1] += exp(-x / alpha);
      pMax[2] += exp(y / alpha);
      pMax[3] += exp(-y / alpha);
    }
    if (_validVirtualPin == 1) {
      pMax[0] += exp(_virtualPinX / alpha);
      pMax[1] += exp(-_virtualPinX / alpha);
      pMax[2] += exp(_virtualPinY / alpha);
      pMax[3] += exp(-_virtualPinY / alpha);
    }
    NumType obj = 0;
    for (int i = 0; i < 4; ++i) {
      obj += log(pMax[i]);
    }
    return alpha * obj * _weight * lambda;
  }

  void accumlateGradient() const {
    if (!validHpwl()) {
      return;
    }
    std::array<NumType, 4> max_val = {0, 0, 0, 0}; // xmax xin ymax ymin
    auto alpha = _getAlphaFunc();
    auto lambda = _getLambdaFunc();
    NumType *pMax = &max_val.front();
    std::vector<std::array<NumType, 4>> exp_results(_cells.size());
    for (IndexType pinIdx = 0; pinIdx < _cells.size(); ++pinIdx) {
      NumType x = op::conv<NumType>(
          _getVarFunc(_cells[pinIdx], Orient2DType::HORIZONTAL) +
          _offsetX[pinIdx]);
      NumType y = op::conv<NumType>(
          _getVarFunc(_cells[pinIdx], Orient2DType::VERTICAL) +
          _offsetY[pinIdx]);
      exp_results[pinIdx][0] = exp(x / alpha);
      pMax[0] += exp_results[pinIdx][0];
      exp_results[pinIdx][1] = exp(-x / alpha);
      pMax[1] += exp_results[pinIdx][1];
      exp_results[pinIdx][2] = exp(y / alpha);
      pMax[2] += exp_results[pinIdx][2];
      exp_results[pinIdx][3] = exp(-y / alpha);
      pMax[3] += exp_results[pinIdx][3];
    }
    if (_validVirtualPin == 1) {
      pMax[0] += exp(_virtualPinX / alpha);
      pMax[1] += exp(-_virtualPinX / alpha);
      pMax[2] += exp(_virtualPinY / alpha);
      pMax[3] += exp(-_virtualPinY / alpha);
    }
    // avoid overflow
    for (IndexType i = 0; i < 4; ++i) {
      pMax[i] = std::max(pMax[i], op::conv<NumType>(1e-8));
    }
    for (IndexType pinIdx = 0; pinIdx < _cells.size(); ++pinIdx) {
      IndexType cellIdx = _cells[pinIdx];
      NumType xPartial = lambda * _weight;
      NumType yPartial = xPartial;
      xPartial *= (exp_results[pinIdx][0] / pMax[0]) -
                  (exp_results[pinIdx][1] / pMax[1]);
      yPartial *= (exp_results[pinIdx][2] / pMax[2]) -
                  (exp_results[pinIdx][3] / pMax[3]);
      _accumulateGradFunc(xPartial, cellIdx, Orient2DType::HORIZONTAL);
      _accumulateGradFunc(yPartial, cellIdx, Orient2DType::VERTICAL);
    }
  }


  IntType _validVirtualPin = 0;
  CoordType _virtualPinX = 0;
  CoordType _virtualPinY = 0;
  std::vector<IndexType> _cells;
  std::vector<CoordType> _offsetX;
  std::vector<CoordType> _offsetY;
  NumType _weight = 1;
  std::function<NumType(void)>
      _getAlphaFunc; ///< A function to get the current alpha
  std::function<NumType(void)>
      _getLambdaFunc; ///< A function to get the current lambda multiplier
  std::function<CoordType(IndexType cellIdx, Orient2DType orient)>
      _getVarFunc; ///< A function to get current variable value
  std::function<void(NumType, IndexType, Orient2DType)>
      _accumulateGradFunc; ///< A function to update partial
};

template <typename NumType, typename CoordType>
struct is_placement_differentiable_concept<
    LseHpwlDifferentiable<NumType, CoordType>> {
  typedef std::true_type is_placement_differentiable_concept_type;
};

// @brief pair-wise cell overlapping penalty
template <typename NumType, typename CoordType>
struct CellPairOverlapPenaltyDifferentiable {
  typedef NumType numerical_type;
  typedef CoordType coordinate_type;

  CellPairOverlapPenaltyDifferentiable() = default;

  CellPairOverlapPenaltyDifferentiable(
      IndexType cellIdxI, CoordType cellWidthI, CoordType cellHeightI,
      IndexType cellIdxJ, CoordType cellWidthJ, CoordType cellHeightJ,
      const std::function<NumType(void)> &getAlphaFunc,
      const std::function<NumType(void)> &getLambdaFunc) {
    _cellIdxI = cellIdxI;
    _cellWidthI = cellWidthI;
    _cellHeightI = cellHeightI;
    _cellIdxJ = cellIdxJ;
    _cellWidthJ = cellWidthJ;
    _cellHeightJ = cellHeightJ;
    _getAlphaFunc = getAlphaFunc;
    _getLambdaFunc = getLambdaFunc;
    _normal = op::conv<NumType>( _cellWidthI * _cellHeightI + _cellWidthJ * _cellHeightJ);
    _normal = 1.0 / std::sqrt(_normal);
  }

  void setGetVarFunc(
      const std::function<CoordType(IndexType, Orient2DType)> &getVarFunc) {
    _getVarFunc = getVarFunc;
  }
  void setAccumulateGradFunc(
      const std::function<void(NumType, IndexType, Orient2DType)> &func) {
    _accumulateGradFunc = func;
  }
  void setGetAlphaFunc(const std::function<NumType(void)> &getAlphaFunc) {
    _getAlphaFunc = getAlphaFunc;
  }

  NumType evaluate() const {
    const NumType xi =
        op::conv<NumType>(_getVarFunc(_cellIdxI, Orient2DType::HORIZONTAL));
    const NumType yi =
        op::conv<NumType>(_getVarFunc(_cellIdxI, Orient2DType::VERTICAL));
    const NumType wi = op::conv<NumType>(_cellWidthI);
    const NumType hi = op::conv<NumType>(_cellHeightI);
    const NumType wj = op::conv<NumType>(_cellWidthJ);
    const NumType hj = op::conv<NumType>(_cellHeightJ);
    NumType xj, yj;
    if (_considerOnlyOneCell) {
      xj = op::conv<NumType>(_cellJXLo);
      yj = op::conv<NumType>(_cellJYLo);
    } else {
      xj = op::conv<NumType>(_getVarFunc(_cellIdxJ, Orient2DType::HORIZONTAL));
      yj = op::conv<NumType>(_getVarFunc(_cellIdxJ, Orient2DType::VERTICAL));
    }
    const NumType alpha = _getAlphaFunc();
    const NumType lambda = _getLambdaFunc();

    const NumType ovl =
        pow(alpha, 2) *
        log(1 / (exp(-(hi + yi - yj) / alpha) + exp(-(hj - yi + yj) / alpha)) +
            1) *
        log(1 / (exp(-(wi + xi - xj) / alpha) + exp(-(wj - xi + xj) / alpha)) +
            1);
    return lambda * ovl * _weight ;
  }

  void accumlateGradient() const {
    /**
     * @brief syms xi xj wi wj alpha yi yj hi hj
     *
     *
     * syms xi xj wi wj alpha yi yj hi hj
     *
     * var1x = (xi + wi - xj);
     * var2x = (xj + wj - xi);
     * min_func_x = -alpha * log( exp(-var1x / alpha) + exp(-var2x / alpha));
     * max_func_x = alpha * log(exp(min_func_x/ alpha) +exp(0));
     *
     * var1y = (yi + hi - yj);
     * var2y = (yj + hj - yi);
     * min_func_y = -alpha * log( exp( - var1y / alpha) + exp( - var2y /
    alpha));
     * max_func_y = alpha * log(exp(min_func_y / alpha) + exp(0));
     *
     * ovl = max_func_x * max_func_y;
     *
     * dxi = diff(ovl, xi)
     * dxj = diff(ovl, xj)
     * dyi = diff(ovl, yi)
     * dyj = diff(ovl, yj)
     *
     *
    NumType dxi_gold = (alpha * alpha *log(1/(exp(-(hi + yi - yj)/alpha) +
    exp(-(hj - yi + yj)/alpha)) + 1)*(exp(-(wi + xi - xj)/alpha)/alpha -
    exp(-(wj - xi + xj)/alpha)/alpha))/((1/(exp(-(wi + xi - xj)/alpha) +
    exp(-(wj - xi + xj)/alpha)) + 1)*pow(exp(-(wi + xi - xj)/alpha) + exp(-(wj -
    xi + xj)/alpha), 2))
        ;
     NumType dxj_gold =
         -(alpha * alpha *log(1/(exp(-(hi + yi - yj)/alpha) + exp(-(hj - yi +
    yj)/alpha)) + 1)*(exp(-(wi + xi - xj)/alpha)/alpha - exp(-(wj - xi +
    xj)/alpha)/alpha))/((1/(exp(-(wi + xi - xj)/alpha) + exp(-(wj - xi +
    xj)/alpha)) + 1)*pow((exp(-(wi + xi - xj)/alpha) + exp(-(wj - xi +
    xj)/alpha)), 2))
         ;


    NumType dyi_gold =
        (alpha * alpha *log(1/(exp(-(wi + xi - xj)/alpha) + exp(-(wj - xi +
    xj)/alpha)) + 1)*(exp(-(hi + yi - yj)/alpha)/alpha - exp(-(hj - yi +
    yj)/alpha)/alpha))/((1/(exp(-(hi + yi - yj)/alpha) + exp(-(hj - yi +
    yj)/alpha)) + 1)*(pow(exp(-(hi + yi - yj)/alpha) + exp(-(hj - yi +
    yj)/alpha), 2)))
        ;

    NumType dyj_gold =
        -(alpha * alpha*log(1/(exp(-(wi + xi - xj)/alpha) + exp(-(wj - xi +
    xj)/alpha)) + 1)*(exp(-(hi + yi - yj)/alpha)/alpha - exp(-(hj - yi +
    yj)/alpha)/alpha))/((1/(exp(-(hi + yi - yj)/alpha) + exp(-(hj - yi +
    yj)/alpha)) + 1)*(pow((exp(-(hi + yi - yj)/alpha) + exp(-(hj - yi +
    yj)/alpha)), 2)))
     *
     */
    const NumType xi =
        op::conv<NumType>(_getVarFunc(_cellIdxI, Orient2DType::HORIZONTAL));
    const NumType yi =
        op::conv<NumType>(_getVarFunc(_cellIdxI, Orient2DType::VERTICAL));
    const NumType wi = op::conv<NumType>(_cellWidthI);
    const NumType hi = op::conv<NumType>(_cellHeightI);
    const NumType wj = op::conv<NumType>(_cellWidthJ);
    const NumType hj = op::conv<NumType>(_cellHeightJ);
    NumType xj, yj;
    if (_considerOnlyOneCell) {
      xj = op::conv<NumType>(_cellJXLo);
      yj = op::conv<NumType>(_cellJYLo);
    } else {
      xj = op::conv<NumType>(_getVarFunc(_cellIdxJ, Orient2DType::HORIZONTAL));
      yj = op::conv<NumType>(_getVarFunc(_cellIdxJ, Orient2DType::VERTICAL));
    }
    const NumType alpha = _getAlphaFunc();
    const NumType lambda = _getLambdaFunc();

    NumType dxi =
        (alpha * alpha *
         log(1 / (exp(-(hi + yi - yj) / alpha) + exp(-(hj - yi + yj) / alpha)) +
             1) *
         (exp(-(wi + xi - xj) / alpha) / alpha -
          exp(-(wj - xi + xj) / alpha) / alpha)) /
        ((1 / (exp(-(wi + xi - xj) / alpha) + exp(-(wj - xi + xj) / alpha)) +
          1) *
         pow(exp(-(wi + xi - xj) / alpha) + exp(-(wj - xi + xj) / alpha), 2));
    dxi *= (lambda);
    NumType dxj = -dxi;

    NumType dyi =
        (alpha * alpha *
         log(1 / (exp(-(wi + xi - xj) / alpha) + exp(-(wj - xi + xj) / alpha)) +
             1) *
         (exp(-(hi + yi - yj) / alpha) / alpha -
          exp(-(hj - yi + yj) / alpha) / alpha)) /
        ((1 / (exp(-(hi + yi - yj) / alpha) + exp(-(hj - yi + yj) / alpha)) +
          1) *
         (pow(exp(-(hi + yi - yj) / alpha) + exp(-(hj - yi + yj) / alpha), 2)));
    dyi *= (lambda);
    NumType dyj = -dyi;

    // accumulate the computed partials
    _accumulateGradFunc(dxi * _weight, _cellIdxI, Orient2DType::HORIZONTAL);
    _accumulateGradFunc(dyi * _weight, _cellIdxI, Orient2DType::VERTICAL);
    if (not _considerOnlyOneCell) {
      _accumulateGradFunc(dxj * _weight , _cellIdxJ, Orient2DType::HORIZONTAL);
      _accumulateGradFunc(dyj * _weight , _cellIdxJ, Orient2DType::VERTICAL);
    }
  }

  /// @brief only consider cell
  void configConsiderOnlyOneCell(CoordType xLo, CoordType yLo) {
    _considerOnlyOneCell = true;
    _cellJXLo = xLo;
    _cellJYLo = yLo;
  }
  
  void setWeight(NumType weight) { _weight  = weight; }

  void penalize() {
    _weight += scale;
  }

  IndexType _cellIdxI;
  CoordType _cellWidthI;
  CoordType _cellHeightI;
  IndexType _cellIdxJ;
  CoordType _cellWidthJ;
  CoordType _cellHeightJ;
  std::function<NumType(void)>
      _getAlphaFunc; ///< A function to get the current alpha
  std::function<NumType(void)>
      _getLambdaFunc; ///< A function to get the current lambda multiplier
  std::function<CoordType(IndexType cellIdx, Orient2DType orient)>
      _getVarFunc; ///< A function to get current variable value
  std::function<void(NumType, IndexType, Orient2DType)>
      _accumulateGradFunc; ///< A function to update partial
  bool _considerOnlyOneCell = false;
  CoordType _cellJXLo; ///< Valid only when _considerOnlyOneCell = true
  CoordType _cellJYLo; ///< Valid only when _consdierOnlyOneCell = true
  NumType _weight = 1.0;
  NumType _normal = 1.0; ///< Normalize between different cells
  static constexpr NumType scale = 1.1;  // Penalize this operator
};

template <typename op_type> struct place_overlap_trait {
  typedef typename op_type::coordinate_type coordinate_type;
  /// @brief calculate the overlap area of an operator
  static coordinate_type overlapArea(op_type &ovl) {
    const coordinate_type xi =
        ovl._getVarFunc(ovl._cellIdxI, Orient2DType::HORIZONTAL);
    const coordinate_type yi =
        ovl._getVarFunc(ovl._cellIdxI, Orient2DType::VERTICAL);
    const coordinate_type wi = ovl._cellWidthI;
    const coordinate_type hi = ovl._cellHeightI;
    const coordinate_type wj = ovl._cellWidthJ;
    const coordinate_type hj = ovl._cellHeightJ;
    coordinate_type xj, yj;
    if (ovl._considerOnlyOneCell) {
      xj = ovl._cellJXLo;
      yj = ovl._cellJYLo;
    } else {
      xj = ovl._getVarFunc(ovl._cellIdxJ, Orient2DType::HORIZONTAL);
      yj = ovl._getVarFunc(ovl._cellIdxJ, Orient2DType::VERTICAL);
    }

    const auto overlapX =
        std::max(std::min(xi + wi, xj + wj) - std::max(xi, xj),
                 op::conv<coordinate_type>(0.0));
    const auto overlapY =
        std::max(std::min(yi + hi, yj + hj) - std::max(yi, yj),
                 op::conv<coordinate_type>(0.0));
    return overlapX * overlapY;
  }
  static coordinate_type cellArea(op_type &ovl) {
    const coordinate_type wi = ovl._cellWidthI;
    const coordinate_type hi = ovl._cellHeightI;
    const coordinate_type wj = ovl._cellWidthJ;
    const coordinate_type hj = ovl._cellHeightJ;
    return wi * hi + wj * hj;
  }
  static coordinate_type smallCellArea(op_type &ovl) {
    const coordinate_type wi = ovl._cellWidthI;
    const coordinate_type hi = ovl._cellHeightI;
    const coordinate_type wj = ovl._cellWidthJ;
    const coordinate_type hj = ovl._cellHeightJ;
    return std::min(wi * hi, wj * hj);
  }
};

template <typename NumType, typename CoordType>
struct is_placement_differentiable_concept<
    CellPairOverlapPenaltyDifferentiable<NumType, CoordType>> {
  typedef std::true_type is_placement_differentiable_concept_type;
};

/// @brief the cell out of boundary penalty
template <typename NumType, typename CoordType>
struct CellOutOfBoundaryPenaltyDifferentiable {
  typedef NumType numerical_type;
  typedef CoordType coordinate_type;

  CellOutOfBoundaryPenaltyDifferentiable(
      IndexType cellIdx, CoordType cellWidth, CoordType cellHeight,
      Box<CoordType> *boundary,
      const std::function<NumType(void)> &getAlphaFunc,
      const std::function<NumType(void)> &getLambdaFunc) {
    _cellIdx = cellIdx;
    _cellWidth = cellWidth;
    _cellHeight = cellHeight;
    _boundary = boundary;
    _getAlphaFunc = getAlphaFunc;
    _getLambdaFunc = getLambdaFunc;
  }

  void setGetVarFunc(
      const std::function<CoordType(IndexType, Orient2DType)> &getVarFunc) {
    _getVarFunc = getVarFunc;
  }
  void setAccumulateGradFunc(
      const std::function<void(NumType, IndexType, Orient2DType)> &func) {
    _accumulateGradFunc = func;
  }
  void setGetAlphaFunc(const std::function<NumType(void)> &getAlphaFunc) {
    _getAlphaFunc = getAlphaFunc;
  }

  NumType evaluate() const {
    const NumType alpha = _getAlphaFunc();
    const NumType lambda = _getLambdaFunc();
    const CoordType xLo = _getVarFunc(_cellIdx, Orient2DType::HORIZONTAL);
    const CoordType yLo = _getVarFunc(_cellIdx, Orient2DType::VERTICAL);
    const CoordType xHi = xLo + _cellWidth;
    const CoordType yHi = yLo + _cellHeight;
    // Smooth abs xLo xHi
    NumType obXLo =
        op::logSumExp0(op::conv<NumType>(_boundary->xLo() - xLo), alpha);
    NumType obXHi =
        op::logSumExp0(op::conv<NumType>(xHi - _boundary->xHi()), alpha);
    // y
    NumType obYLo =
        op::logSumExp0(op::conv<NumType>(_boundary->yLo() - yLo), alpha);
    NumType obYHi =
        op::logSumExp0(op::conv<NumType>(yHi - _boundary->yHi()), alpha);
    return (obXLo + obXHi + obYLo + obYHi) * lambda;
  }

  void accumlateGradient() const {
    const NumType alpha = _getAlphaFunc();
    const NumType lambda = _getLambdaFunc();
    const CoordType xLo = _getVarFunc(_cellIdx, Orient2DType::HORIZONTAL);
    const CoordType yLo = _getVarFunc(_cellIdx, Orient2DType::VERTICAL);
    const CoordType xHi = xLo + _cellWidth;
    const CoordType yHi = yLo + _cellHeight;
    // max(lower - x/yLo, 0), max (x/yHi - upper, 0)
    NumType gradObX = -op::gradLogSumExp0( // negative comes from the derivative
        op::conv<NumType>(_boundary->xLo() - xLo), alpha);
    gradObX +=
        op::gradLogSumExp0(op::conv<NumType>(xHi - _boundary->xHi()), alpha);
    _accumulateGradFunc(gradObX * lambda, _cellIdx, Orient2DType::HORIZONTAL);
    // y
    NumType gradObY = -op::gradLogSumExp0( // negative comes from the derivative
        op::conv<NumType>(_boundary->yLo() - yLo), alpha);
    gradObY +=
        op::gradLogSumExp0(op::conv<NumType>(yHi - _boundary->yHi()), alpha);
    _accumulateGradFunc(gradObY * lambda, _cellIdx, Orient2DType::VERTICAL);
  }

  IndexType _cellIdx;
  CoordType _cellWidth;
  CoordType _cellHeight;
  Box<CoordType> *_boundary = nullptr;
  std::function<NumType(void)> _getAlphaFunc;
  std::function<NumType(void)>
      _getLambdaFunc; ///< A function to get the current lambda multiplier
  std::function<CoordType(IndexType cellIdx, Orient2DType orient)>
      _getVarFunc; ///< A function to get current variable value
  std::function<void(NumType, IndexType, Orient2DType)>
      _accumulateGradFunc; ///< A function to update partial

};

template <typename NumType, typename CoordType>
struct is_placement_differentiable_concept<
    CellOutOfBoundaryPenaltyDifferentiable<NumType, CoordType>> {
  typedef std::true_type is_placement_differentiable_concept_type;
};

template <typename op_type> struct place_out_of_boundary_trait {
  typedef typename op_type::coordinate_type coordinate_type;
  /// @brief calculate the out of boundary area of an operator
  static coordinate_type oobArea(op_type &oob) {
    const coordinate_type x =
        oob._getVarFunc(oob._cellIdx, Orient2DType::HORIZONTAL);
    const coordinate_type y =
        oob._getVarFunc(oob._cellIdx, Orient2DType::VERTICAL);
    const coordinate_type w = oob._cellWidth;
    const coordinate_type h = oob._cellHeight;
    const auto boxXHi = (*(oob._boundary)).xHi();
    const auto boxXLo = (*(oob._boundary)).xLo();
    const auto boxYHi = (*(oob._boundary)).yHi();
    const auto boxYLo = (*(oob._boundary)).yLo();

    const auto overlapX =
        std::max(std::min(x + w, boxXHi) - std::max(x, boxXLo),
                 op::conv<coordinate_type>(0.0));
    const auto overlapY =
        std::max(std::min(y + h, boxYHi) - std::max(y, boxYLo),
                 op::conv<coordinate_type>(0.0));
    return w * h - overlapX * overlapY;
  }
};

/// @brief Asymmetry penalty
template <typename NumType, typename CoordType> struct AsymmetryDifferentiable {
  typedef NumType numerical_type;
  typedef CoordType coordinate_type;

  AsymmetryDifferentiable(IndexType symGrpIdx,
                          const std::function<NumType(void)> &getLambdaFunc) {
    _symGrpIdx = symGrpIdx;
    _getLambdaFunc = getLambdaFunc;
  }

  void setGetVarFunc(
      const std::function<CoordType(IndexType, Orient2DType)> &getVarFunc) {
    _getVarFunc = getVarFunc;
  }
  void setAccumulateGradFunc(
      const std::function<void(NumType, IndexType, Orient2DType)> &func) {
    _accumulateGradFunc = func;
  }

  /// @brief add a symmetric pair. require the cell widths are the same
  void addSymPair(IndexType cellIdxI, IndexType cellIdxJ, CoordType width) {
    _pairCells.emplace_back(std::array<IndexType, 2>({cellIdxI, cellIdxJ}));
    _pairWidths.emplace_back(op::conv<NumType>(width));
  }
  void addSelfSym(IndexType cellIdx, CoordType width) {
    _selfSymCells.emplace_back(cellIdx);
    _selfSymWidths.emplace_back(op::conv<NumType>(width));
  }

  NumType evaluate() const {
    NumType lambda = _getLambdaFunc();
    NumType asym = 0;
    NumType symAxis =
        op::conv<NumType>(_getVarFunc(_symGrpIdx, Orient2DType::NONE));
    for (IndexType symPairIdx = 0; symPairIdx < _pairCells.size();
         ++symPairIdx) {
      const IndexType cellI = _pairCells[symPairIdx][0];
      const IndexType cellJ = _pairCells[symPairIdx][1];
      const NumType xi =
          op::conv<NumType>(_getVarFunc(cellI, Orient2DType::HORIZONTAL));
      const NumType yi =
          op::conv<NumType>(_getVarFunc(cellI, Orient2DType::VERTICAL));
      const NumType w = _pairWidths[symPairIdx];
      const NumType xj =
          op::conv<NumType>(_getVarFunc(cellJ, Orient2DType::HORIZONTAL));
      const NumType yj =
          op::conv<NumType>(_getVarFunc(cellJ, Orient2DType::VERTICAL));

      asym += pow(yi - yj, 2.0);
      asym += pow(xi + xj + w - 2 * symAxis, 2.0);
    }
    for (IndexType ssIdx = 0; ssIdx < _selfSymCells.size(); ++ssIdx) {
      NumType x = op::conv<NumType>(
          _getVarFunc(_selfSymCells[ssIdx], Orient2DType::HORIZONTAL));
      NumType w = _selfSymWidths[ssIdx];

      asym += pow(x + w / 2 - symAxis, 2.0);
    }
    return asym * lambda * _weight;
  }
  void accumlateGradient() const {
    NumType lambda = _getLambdaFunc();
    NumType symAxis =
        op::conv<NumType>(_getVarFunc(_symGrpIdx, Orient2DType::NONE));
    for (IndexType symPairIdx = 0; symPairIdx < _pairCells.size();
         ++symPairIdx) {
      const IndexType cellI = _pairCells[symPairIdx][0];
      const IndexType cellJ = _pairCells[symPairIdx][1];
      const NumType xi =
          op::conv<NumType>(_getVarFunc(cellI, Orient2DType::HORIZONTAL));
      const NumType yi =
          op::conv<NumType>(_getVarFunc(cellI, Orient2DType::VERTICAL));
      const NumType w = _pairWidths[symPairIdx];
      const NumType xj =
          op::conv<NumType>(_getVarFunc(cellJ, Orient2DType::HORIZONTAL));
      const NumType yj =
          op::conv<NumType>(_getVarFunc(cellJ, Orient2DType::VERTICAL));

      NumType partialX = 2.0 * (xi + xj + w - 2 * symAxis) * lambda * _weight;
      _accumulateGradFunc(partialX, cellI, Orient2DType::HORIZONTAL);
      _accumulateGradFunc(partialX, cellJ, Orient2DType::HORIZONTAL);
      _accumulateGradFunc(-2 * partialX, _symGrpIdx, Orient2DType::NONE);

      NumType partialYI = 2.0 * (yi - yj) * lambda * _weight;
      _accumulateGradFunc(partialYI, cellI, Orient2DType::VERTICAL);
      _accumulateGradFunc(-partialYI, cellJ, Orient2DType::VERTICAL);
    }
    for (IndexType ssIdx = 0; ssIdx < _selfSymCells.size(); ++ssIdx) {
      const NumType x = op::conv<NumType>(
          _getVarFunc(_selfSymCells[ssIdx], Orient2DType::HORIZONTAL));
      const NumType w = _selfSymWidths[ssIdx];

      NumType partial = 2.0 * (x + w / 2 - symAxis) * lambda;

      _accumulateGradFunc(partial * _weight, _selfSymCells[ssIdx],
                          Orient2DType::HORIZONTAL);
      _accumulateGradFunc(-partial * _weight, _symGrpIdx, Orient2DType::NONE);
    }
  }

  void setWeight(NumType weight) { _weight = weight; }

  void penalize() { _weight += _penalizeScale;  }

  IndexType _symGrpIdx;
  std::vector<std::array<IndexType, 2>> _pairCells;
  std::vector<NumType> _pairWidths;
  std::vector<IndexType> _selfSymCells;
  std::vector<NumType> _selfSymWidths;
  NumType _weight = 1.0;
  NumType _penalizeScale = 3.0;
  std::function<NumType(void)> _getLambdaFunc;
  std::function<CoordType(IndexType cellIdx, Orient2DType orient)>
      _getVarFunc; ///< A function to get current variable value
  std::function<void(NumType, IndexType, Orient2DType)>
      _accumulateGradFunc; ///< A function to update partial
};

template <typename NumType, typename CoordType>
struct is_placement_differentiable_concept<
    AsymmetryDifferentiable<NumType, CoordType>> {
  typedef std::true_type is_placement_differentiable_concept_type;
};

template <typename op_type> struct place_asym_trait {
  typedef typename op_type::coordinate_type coordinate_type;
  /// @brief calculate the asymmetry distance of an operator
  static coordinate_type asymDistance(op_type &asym) {
    coordinate_type dist = 0;
    auto symAxis = asym._getVarFunc(asym._symGrpIdx, Orient2DType::NONE);
    for (IndexType symPairIdx = 0; symPairIdx < asym._pairCells.size();
         ++symPairIdx) {
      const IndexType cellI = asym._pairCells[symPairIdx][0];
      const IndexType cellJ = asym._pairCells[symPairIdx][1];
      const auto xi = asym._getVarFunc(cellI, Orient2DType::HORIZONTAL);
      const auto yi = asym._getVarFunc(cellI, Orient2DType::VERTICAL);
      const auto w = asym._pairWidths[symPairIdx];
      const auto xj = asym._getVarFunc(cellJ, Orient2DType::HORIZONTAL);
      const auto yj = asym._getVarFunc(cellJ, Orient2DType::VERTICAL);

      dist += std::abs(xi + xj + w - 2 * symAxis);
      dist += std::abs(yi - yj);
    }
    for (IndexType ssIdx = 0; ssIdx < asym._selfSymCells.size(); ++ssIdx) {
      const auto x =
          asym._getVarFunc(asym._selfSymCells[ssIdx], Orient2DType::HORIZONTAL);
      const auto w = asym._selfSymWidths[ssIdx];

      dist += std::abs(x + w * 0.5 - symAxis);
    }
    return dist;
  }
  /// @brief calculate the normalized asymmetry distance of an operator
  static coordinate_type asymDistanceNormalized(op_type &asym) {
    return asymDistance(asym) /
           (asym._pairCells.size() + asym._selfSymCells.size());
  }
  static coordinate_type cellWidthSum(op_type &asym) {
    coordinate_type sum = 0;
    for (IndexType symPairIdx = 0; symPairIdx < asym._pairCells.size();
         ++symPairIdx) {
      sum +=  asym._pairWidths[symPairIdx];

    }
    for (IndexType ssIdx = 0; ssIdx < asym._selfSymCells.size(); ++ssIdx) {
      sum += asym._selfSymWidths[ssIdx];
    }
    return sum;
  }
};

/// @brief Each individual operator includes three cells and hence the two
/// vectors they compose
template <typename NumType, typename CoordType>
struct CosineDatapathDifferentiable {
  typedef NumType numerical_type;
  typedef CoordType coordinate_type;

  CosineDatapathDifferentiable(
      IndexType sCellIdx, const Point<CoordType> &sOffset, IndexType midCellIdx,
      const Point<CoordType> &midOffsetA, const Point<CoordType> &midOffsetB,
      IndexType tCellIdx, const Point<CoordType> &tOffset,
      const std::function<NumType(void)> &getLambdaFunc)
      : _sCellIdx(sCellIdx), _midCellIdx(midCellIdx), _tCellIdx(tCellIdx),
        _getLambdaFunc(getLambdaFunc) {
    _sOffset.setX(op::conv<NumType>(sOffset.x()));
    _sOffset.setY(op::conv<NumType>(sOffset.y()));
    _midOffsetA.setX(op::conv<NumType>(midOffsetA.x()));
    _midOffsetA.setY(op::conv<NumType>(midOffsetA.y()));
    _midOffsetB.setX(op::conv<NumType>(midOffsetB.x()));
    _midOffsetB.setY(op::conv<NumType>(midOffsetB.y()));
    _tOffset.setX(tOffset.x());
    _tOffset.setY(tOffset.y());
  }

  CosineDatapathDifferentiable(
      IndexType sCellIdx, const Point<CoordType> &sOffset, IndexType midCellIdx,
      const Point<CoordType> &midOffsetA, const Point<CoordType> &midOffsetB,
      const std::function<NumType(void)> &getLambdaFunc)
      : _sCellIdx(sCellIdx), _midCellIdx(midCellIdx),
        _getLambdaFunc(getLambdaFunc) {
    _sOffset.setX(op::conv<NumType>(sOffset.x()));
    _sOffset.setY(op::conv<NumType>(sOffset.y()));
    _midOffsetA.setX(op::conv<NumType>(midOffsetA.x()));
    _midOffsetA.setY(op::conv<NumType>(midOffsetA.y()));
    _midOffsetB.setX(op::conv<NumType>(midOffsetB.x()));
    _midOffsetB.setY(op::conv<NumType>(midOffsetB.y()));
    markTwoPin();
    _enable = false;
  }

  void setGetVarFunc(
      const std::function<CoordType(IndexType, Orient2DType)> &getVarFunc) {
    _getVarFunc = getVarFunc;
  }
  void setAccumulateGradFunc(
      const std::function<void(NumType, IndexType, Orient2DType)> &func) {
    _accumulateGradFunc = func;
  }

  BoolType isTwoPin() const { return _tCellIdx == INDEX_TYPE_MAX; }
  void markTwoPin() { _tCellIdx = INDEX_TYPE_MAX; }
  void setTwoPinEndingOffset(const Point<CoordType> &tOffset) {
    Assert(isTwoPin());
    _enable = true;
    _tOffset.setX(tOffset.x());
    _tOffset.setY(tOffset.y());
  }

  NumType evaluate() const;
  void accumlateGradient() const;

  void setWeight(NumType weight) { _weight = weight; }

  IndexType _sCellIdx = INDEX_TYPE_MAX;   ///< Source
  Point<NumType> _sOffset;                ///< The offset for x0
  IndexType _midCellIdx = INDEX_TYPE_MAX; ///< Middle
  Point<NumType> _midOffsetA;
  Point<NumType> _midOffsetB;
  IndexType _tCellIdx = INDEX_TYPE_MAX; ///< Target
  Point<NumType> _tOffset;
  std::function<NumType(void)>
      _getLambdaFunc; ///< A function to get the current lambda multiplier
  std::function<CoordType(IndexType cellIdx, Orient2DType orient)>
      _getVarFunc; ///< A function to get current variable value
  std::function<void(NumType, IndexType, Orient2DType)>
      _accumulateGradFunc; ///< A function to update partial
  NumType _weight = 1.0;
  bool _enable = true;
};

template <typename NumType, typename CoordType>
inline NumType
CosineDatapathDifferentiable<NumType, CoordType>::evaluate() const {
  if (not _enable) {
    return 0;
  }
  const NumType lambda = _getLambdaFunc();
  const NumType x1 =
      op::conv<NumType>(_getVarFunc(_sCellIdx, Orient2DType::HORIZONTAL));
  const NumType y1 =
      op::conv<NumType>(_getVarFunc(_sCellIdx, Orient2DType::VERTICAL));
  const NumType x2 =
      op::conv<NumType>(_getVarFunc(_midCellIdx, Orient2DType::HORIZONTAL));
  const NumType y2 =
      op::conv<NumType>(_getVarFunc(_midCellIdx, Orient2DType::VERTICAL));
  NumType x3 = 0;
  NumType y3 = 0;
  if (_tCellIdx != INDEX_TYPE_MAX) {
    x3 = op::conv<NumType>(_getVarFunc(_tCellIdx, Orient2DType::HORIZONTAL));
    y3 = op::conv<NumType>(_getVarFunc(_tCellIdx, Orient2DType::VERTICAL));
  }
  const NumType ox1 = _sOffset.x();
  const NumType oy1 = _sOffset.y();
  const NumType ox2a = _midOffsetA.x();
  const NumType oy2a = _midOffsetA.y();
  const NumType ox2b = _midOffsetB.x();
  const NumType oy2b = _midOffsetB.y();
  const NumType ox3 = _tOffset.x();
  const NumType oy3 = _tOffset.y();

  return (((ox1 - ox2a + x1 - x2) * (ox3 - ox2b - x2 + x3) +
           (oy1 - oy2a + y1 - y2) * (oy3 - oy2b - y2 + y3)) /
              (sqrt(pow(ox1 - ox2a + x1 - x2, 2.0) +
                    pow(oy1 - oy2a + y1 - y2, 2.0)) *
               sqrt(pow(ox3 - ox2b - x2 + x3, 2.0) +
                    pow(oy3 - oy2b - y2 + y3, 2.0))) +
          1) *
         lambda * _weight;
}

template <typename NumType, typename CoordType>
inline void
CosineDatapathDifferentiable<NumType, CoordType>::accumlateGradient() const {
  if (not _enable) {
    return;
  }
  const NumType lambda = _getLambdaFunc();
  const NumType x1 =
      op::conv<NumType>(_getVarFunc(_sCellIdx, Orient2DType::HORIZONTAL));
  const NumType y1 =
      op::conv<NumType>(_getVarFunc(_sCellIdx, Orient2DType::VERTICAL));
  const NumType x2 =
      op::conv<NumType>(_getVarFunc(_midCellIdx, Orient2DType::HORIZONTAL));
  const NumType y2 =
      op::conv<NumType>(_getVarFunc(_midCellIdx, Orient2DType::VERTICAL));
  NumType x3 = 0;
  NumType y3 = 0;
  if (_tCellIdx != INDEX_TYPE_MAX) {
    x3 = op::conv<NumType>(_getVarFunc(_tCellIdx, Orient2DType::HORIZONTAL));
    y3 = op::conv<NumType>(_getVarFunc(_tCellIdx, Orient2DType::VERTICAL));
  }
  const NumType ox1 = _sOffset.x();
  const NumType oy1 = _sOffset.y();
  const NumType ox2a = _midOffsetA.x();
  const NumType oy2a = _midOffsetA.y();
  const NumType ox2b = _midOffsetB.x();
  const NumType oy2b = _midOffsetB.y();
  const NumType ox3 = _tOffset.x();
  const NumType oy3 = _tOffset.y();

  NumType dx1 =
      (ox3 - ox2b - x2 + x3) / (sqrt(pow(ox1 - ox2a + x1 - x2, 2.0) +
                                     pow(oy1 - oy2a + y1 - y2, 2.0)) *
                                sqrt(pow(ox3 - ox2b - x2 + x3, 2.0) +
                                     pow(oy3 - oy2b - y2 + y3, 2.0))) -
      (((ox1 - ox2a + x1 - x2) * (ox3 - ox2b - x2 + x3) +
        (oy1 - oy2a + y1 - y2) * (oy3 - oy2b - y2 + y3)) *
       (2 * ox1 - 2 * ox2a + 2 * x1 - 2 * x2)) /
          (2 *
           pow(pow(ox1 - ox2a + x1 - x2, 2.0) + pow(oy1 - oy2a + y1 - y2, 2.0),
               1.5) *
           sqrt(pow(ox3 - ox2b - x2 + x3, 2.0) +
                pow(oy3 - oy2b - y2 + y3, 2.0)));

  dx1 *= lambda;

  NumType dy1 =
      (oy3 - oy2b - y2 + y3) / (sqrt(pow(ox1 - ox2a + x1 - x2, 2.0) +
                                     pow(oy1 - oy2a + y1 - y2, 2.0)) *
                                sqrt(pow(ox3 - ox2b - x2 + x3, 2.0) +
                                     pow(oy3 - oy2b - y2 + y3, 2.0))) -
      (((ox1 - ox2a + x1 - x2) * (ox3 - ox2b - x2 + x3) +
        (oy1 - oy2a + y1 - y2) * (oy3 - oy2b - y2 + y3)) *
       (2 * oy1 - 2 * oy2a + 2 * y1 - 2 * y2)) /
          (2 *
           pow(pow(ox1 - ox2a + x1 - x2, 2.0) + pow(oy1 - oy2a + y1 - y2, 2.0),
               1.5) *
           sqrt(pow(ox3 - ox2b - x2 + x3, 2.0) +
                pow(oy3 - oy2b - y2 + y3, 2.0)));
  dy1 *= lambda;

  NumType dx2 =
      (((ox1 - ox2a + x1 - x2) * (ox3 - ox2b - x2 + x3) +
        (oy1 - oy2a + y1 - y2) * (oy3 - oy2b - y2 + y3)) *
       (2 * ox1 - 2 * ox2a + 2 * x1 - 2 * x2)) /
          (2 *
           pow(pow(ox1 - ox2a + x1 - x2, 2.0) + pow(oy1 - oy2a + y1 - y2, 2.0),
               1.5) *
           sqrt(pow(ox3 - ox2b - x2 + x3, 2.0) +
                pow(oy3 - oy2b - y2 + y3, 2.0))) -
      (ox1 + ox3 - ox2a - ox2b + x1 - 2 * x2 + x3) /
          (sqrt(pow(ox1 - ox2a + x1 - x2, 2.0) +
                pow(oy1 - oy2a + y1 - y2, 2.0)) *
           sqrt(pow(ox3 - ox2b - x2 + x3, 2.0) +
                pow(oy3 - oy2b - y2 + y3, 2.0))) +
      (((ox1 - ox2a + x1 - x2) * (ox3 - ox2b - x2 + x3) +
        (oy1 - oy2a + y1 - y2) * (oy3 - oy2b - y2 + y3)) *
       (2 * ox3 - 2 * ox2b - 2 * x2 + 2 * x3)) /
          (2 *
           sqrt(pow(ox1 - ox2a + x1 - x2, 2.0) +
                pow(oy1 - oy2a + y1 - y2, 2.0)) *
           pow(pow(ox3 - ox2b - x2 + x3, 2.0) + pow(oy3 - oy2b - y2 + y3, 2.0),
               1.5));
  dx2 *= lambda;

  NumType dy2 =
      (((ox1 - ox2a + x1 - x2) * (ox3 - ox2b - x2 + x3) +
        (oy1 - oy2a + y1 - y2) * (oy3 - oy2b - y2 + y3)) *
       (2 * oy1 - 2 * oy2a + 2 * y1 - 2 * y2)) /
          (2 *
           pow(pow(ox1 - ox2a + x1 - x2, 2.0) + pow(oy1 - oy2a + y1 - y2, 2.0),
               1.5) *
           sqrt(pow(ox3 - ox2b - x2 + x3, 2.0) +
                pow(oy3 - oy2b - y2 + y3, 2.0))) -
      (oy1 + oy3 - oy2a - oy2b + y1 - 2 * y2 + y3) /
          (sqrt(pow(ox1 - ox2a + x1 - x2, 2.0) +
                pow(oy1 - oy2a + y1 - y2, 2.0)) *
           sqrt(pow(ox3 - ox2b - x2 + x3, 2.0) +
                pow(oy3 - oy2b - y2 + y3, 2.0))) +
      (((ox1 - ox2a + x1 - x2) * (ox3 - ox2b - x2 + x3) +
        (oy1 - oy2a + y1 - y2) * (oy3 - oy2b - y2 + y3)) *
       (2 * oy3 - 2 * oy2b - 2 * y2 + 2 * y3)) /
          (2 *
           sqrt(pow(ox1 - ox2a + x1 - x2, 2.0) +
                pow(oy1 - oy2a + y1 - y2, 2.0)) *
           pow(pow(ox3 - ox2b - x2 + x3, 2.0) + pow(oy3 - oy2b - y2 + y3, 2.0),
               1.5));
  dy2 *= lambda;

  NumType dx3 =
      (ox1 - ox2a + x1 - x2) / (sqrt(pow(ox1 - ox2a + x1 - x2, 2.0) +
                                     pow(oy1 - oy2a + y1 - y2, 2.0)) *
                                sqrt(pow(ox3 - ox2b - x2 + x3, 2.0) +
                                     pow(oy3 - oy2b - y2 + y3, 2.0))) -
      (((ox1 - ox2a + x1 - x2) * (ox3 - ox2b - x2 + x3) +
        (oy1 - oy2a + y1 - y2) * (oy3 - oy2b - y2 + y3)) *
       (2 * ox3 - 2 * ox2b - 2 * x2 + 2 * x3)) /
          (2 *
           sqrt(pow(ox1 - ox2a + x1 - x2, 2.0) +
                pow(oy1 - oy2a + y1 - y2, 2.0)) *
           pow(pow(ox3 - ox2b - x2 + x3, 2.0) + pow(oy3 - oy2b - y2 + y3, 2.0),
               1.5));
  dx3 *= lambda;

  NumType dy3 =
      (oy1 - oy2a + y1 - y2) / (sqrt(pow(ox1 - ox2a + x1 - x2, 2.0) +
                                     pow(oy1 - oy2a + y1 - y2, 2.0)) *
                                sqrt(pow(ox3 - ox2b - x2 + x3, 2.0) +
                                     pow(oy3 - oy2b - y2 + y3, 2.0))) -
      (((ox1 - ox2a + x1 - x2) * (ox3 - ox2b - x2 + x3) +
        (oy1 - oy2a + y1 - y2) * (oy3 - oy2b - y2 + y3)) *
       (2 * oy3 - 2 * oy2b - 2 * y2 + 2 * y3)) /
          (2 *
           sqrt(pow(ox1 - ox2a + x1 - x2, 2.0) +
                pow(oy1 - oy2a + y1 - y2, 2.0)) *
           pow(pow(ox3 - ox2b - x2 + x3, 2.0) + pow(oy3 - oy2b - y2 + y3, 2.0),
               1.5));
  dy3 *= lambda;

  _accumulateGradFunc(dx1 * _weight, _sCellIdx, Orient2DType::HORIZONTAL);
  _accumulateGradFunc(dy1 * _weight, _sCellIdx, Orient2DType::VERTICAL);
  _accumulateGradFunc(dx2 * _weight, _midCellIdx, Orient2DType::HORIZONTAL);
  _accumulateGradFunc(dy2 * _weight, _midCellIdx, Orient2DType::VERTICAL);
  _accumulateGradFunc(dx3 * _weight, _tCellIdx, Orient2DType::HORIZONTAL);
  _accumulateGradFunc(dy3 * _weight, _tCellIdx, Orient2DType::VERTICAL);
}

/// @brief LSE-smoothed HPWL
template <typename NumType, typename CoordType>
struct PowerVerQuadraticWireLengthDifferentiable {
  typedef NumType numerical_type;
  typedef CoordType coordinate_type;

  PowerVerQuadraticWireLengthDifferentiable(
      const std::function<NumType(void)> &getLambdaFunc) {
    _getLambdaFunc = getLambdaFunc;
  }

  void setGetVarFunc(
      const std::function<CoordType(IndexType, Orient2DType)> &getVarFunc) {
    _getVarFunc = getVarFunc;
  }
  void setAccumulateGradFunc(
      const std::function<void(NumType, IndexType, Orient2DType)> &func) {
    _accumulateGradFunc = func;
  }

  void setVirtualPin(const CoordType &x, const CoordType &y) {
    _validVirtualPin = 1;
    _virtualPinX = x;
    _virtualPinY = y;
  }
  void removeVirtualPin() { _validVirtualPin = 0; }
  void addVar(IndexType cellIdx, const CoordType &offsetX,
              const CoordType &offsetY) {
    _cells.emplace_back(cellIdx);
    _offsetX.emplace_back(offsetX);
    _offsetY.emplace_back(offsetY);
  }
  void setWeight(const NumType &weight) { _weight = weight; }
  bool validHpwl() const { return _cells.size() + _validVirtualPin > 1; }

  NumType evaluate() const {
    if (!validHpwl()) {
      return 0;
    }
    if (_validVirtualPin != 1) {
      return 0;
    }
    const NumType lambda = _getLambdaFunc();
    NumType obj = 0;
    for (IndexType pinIdx = 0; pinIdx < _cells.size(); ++pinIdx) {
      NumType y = op::conv<NumType>(
          _getVarFunc(_cells[pinIdx], Orient2DType::VERTICAL) +
          _offsetY[pinIdx]);
      obj += std::pow(y - _virtualPinY, 2.0);
    }
    return obj * _weight * lambda;
  }

  void accumlateGradient() const {
    if (!validHpwl()) {
      return;
    }
    if (_validVirtualPin != 1) {
      return;
    }
    const NumType lambda = _getLambdaFunc();
    for (IndexType pinIdx = 0; pinIdx < _cells.size(); ++pinIdx) {
      NumType y = op::conv<NumType>(
          _getVarFunc(_cells[pinIdx], Orient2DType::VERTICAL) +
          _offsetY[pinIdx]);
      const NumType yPartial = 2 * (y - _virtualPinY);
      _accumulateGradFunc(yPartial * _weight * lambda, _cells[pinIdx],
                          Orient2DType::VERTICAL);
    }
  }

  IntType _validVirtualPin = 0;
  CoordType _virtualPinX = 0;
  CoordType _virtualPinY = 0;
  std::vector<IndexType> _cells;
  std::vector<CoordType> _offsetX;
  std::vector<CoordType> _offsetY;
  NumType _weight = 1;
  std::function<NumType(void)>
      _getLambdaFunc; ///< A function to get the current lambda multiplier
  std::function<CoordType(IndexType cellIdx, Orient2DType orient)>
      _getVarFunc; ///< A function to get current variable value
  std::function<void(NumType, IndexType, Orient2DType)>
      _accumulateGradFunc; ///< A function to update partial
};

template <typename NumType, typename CoordType>
struct is_placement_differentiable_concept<
    PowerVerQuadraticWireLengthDifferentiable<NumType, CoordType>> {
  typedef std::true_type is_placement_differentiable_concept_type;
};

/// @brief Each individual operator includes three cells and hence the two
/// vectors they compose
template <typename NumType, typename CoordType>
struct CurrentFlowDifferentiable {
  typedef NumType numerical_type;
  typedef CoordType coordinate_type;

  CurrentFlowDifferentiable(IndexType sCellIdx, const CoordType sOffset,
                            IndexType tCellIdx, const CoordType tOffset,
                            const std::function<NumType(void)> &getLambdaFunc)
      : _sCellIdx(sCellIdx), _tCellIdx(tCellIdx),
        _getLambdaFunc(getLambdaFunc) {
    _sOffset = op::conv<NumType>(sOffset);
    _tOffset = op::conv<NumType>(tOffset);
  }

  void setGetVarFunc(
      const std::function<CoordType(IndexType, Orient2DType)> &getVarFunc) {
    _getVarFunc = getVarFunc;
  }
  void setAccumulateGradFunc(
      const std::function<void(NumType, IndexType, Orient2DType)> &func) {
    _accumulateGradFunc = func;
  }
  void setGetAlphaFunc(const std::function<NumType(void)> &getAlphaFunc) {
    _getAlphaFunc = getAlphaFunc;
  }

  NumType evaluate() const;
  void accumlateGradient() const;

  void setWeight(NumType weight) { _weight = weight; }

  IndexType _sCellIdx = INDEX_TYPE_MAX; ///< Source
  NumType _sOffset;                     ///< The offset for source y
  IndexType _tCellIdx = INDEX_TYPE_MAX; ///< Target
  NumType _tOffset;                     ///< The offset for target y
  std::function<NumType(void)>
      _getLambdaFunc; ///< A function to get the current lambda multiplier
  std::function<CoordType(IndexType cellIdx, Orient2DType orient)>
      _getVarFunc; ///< A function to get current variable value
  std::function<void(NumType, IndexType, Orient2DType)>
      _accumulateGradFunc; ///< A function to update partial
  NumType _weight = 1.0;
  std::function<NumType(void)> _getAlphaFunc;
};
template <typename NumType, typename CoordType>
inline NumType CurrentFlowDifferentiable<NumType, CoordType>::evaluate() const {
  const NumType lambda = _getLambdaFunc();
  const NumType alpha = _getAlphaFunc();
  const NumType y1 =
      op::conv<NumType>(_getVarFunc(_sCellIdx, Orient2DType::VERTICAL));
  const NumType y2 =
      op::conv<NumType>(_getVarFunc(_tCellIdx, Orient2DType::VERTICAL));
  const NumType oy1 = _sOffset;
  const NumType oy2 = _tOffset;

  return alpha * std::log(std::exp(-(oy1 - oy2 + y1 - y2) / alpha) + 1) *
         lambda * _weight;
}

template <typename NumType, typename CoordType>
inline void
CurrentFlowDifferentiable<NumType, CoordType>::accumlateGradient() const {
  const NumType lambda = _getLambdaFunc();
  const NumType alpha = _getAlphaFunc();
  const NumType y1 =
      op::conv<NumType>(_getVarFunc(_sCellIdx, Orient2DType::VERTICAL));
  const NumType y2 =
      op::conv<NumType>(_getVarFunc(_tCellIdx, Orient2DType::VERTICAL));
  const NumType oy1 = _sOffset;
  const NumType oy2 = _tOffset;

  NumType dy1 = -std::exp(-(oy1 - oy2 + y1 - y2) / alpha) /
                (std::exp(-(oy1 - oy2 + y1 - y2) / alpha) + 1);
  ;
  dy1 *= (lambda * _weight);

  NumType dy2 = std::exp(-(oy1 - oy2 + y1 - y2) / alpha) /
                (std::exp(-(oy1 - oy2 + y1 - y2) / alpha) + 1);
  dy2 *= (lambda * _weight);

  _accumulateGradFunc(dy1 * _weight, _sCellIdx, Orient2DType::VERTICAL);
  _accumulateGradFunc(dy2 * _weight, _tCellIdx, Orient2DType::VERTICAL);
}

/* Fence region */

/// @brief Model fence region penalty as recirocral of the overlapping area.
/// it first slices the polygon into rectangles, and
/// the cost is overlapping the sum of overlapping area of the cell and each
/// rectangle
template <typename NumType, typename CoordType>
struct FenceReciprocalOverlapSumBoxDifferentiable {
  typedef NumType numerical_type;
  typedef CoordType coordinate_type;

  FenceReciprocalOverlapSumBoxDifferentiable(
      IndexType cellIdx, const CoordType cellWidth, const CoordType cellHeight,
      const std::vector<Box<CoordType>> &boxes,
      const std::function<NumType(void)> &getAlphaFunc,
      const std::function<NumType(void)> &getLambdaFunc) {
    _cellIdx = cellIdx;
    _cellWidth = cellWidth;
    _cellHeight = cellHeight;
    _boxes = boxes;
    _getAlphaFunc = getAlphaFunc;
    _getLambdaFunc = getLambdaFunc;
  }

  void setGetVarFunc(
      const std::function<CoordType(IndexType, Orient2DType)> &getVarFunc) {
    _getVarFunc = getVarFunc;
  }
  void setAccumulateGradFunc(
      const std::function<void(NumType, IndexType, Orient2DType)> &func) {
    _accumulateGradFunc = func;
  }
  void setGetAlphaFunc(const std::function<NumType(void)> &getAlphaFunc) {
    _getAlphaFunc = getAlphaFunc;
  }

  NumType evaluate() const;
  void accumlateGradient() const;

  void setWeight(NumType weight) { _weight = weight; }

  IndexType _cellIdx;                 ///< The index of the cell
  CoordType _cellWidth;               ///< The width of the cell
  CoordType _cellHeight;              ///< The height of the cell
  std::vector<Box<CoordType>> _boxes; ///< The sliced polygon

  std::function<NumType(void)>
      _getLambdaFunc; ///< A function to get the current lambda multiplier
  std::function<CoordType(IndexType cellIdx, Orient2DType orient)>
      _getVarFunc; ///< A function to get current variable value
  std::function<void(NumType, IndexType, Orient2DType)>
      _accumulateGradFunc; ///< A function to update partial
  NumType _weight = 1.0;
  std::function<NumType(void)> _getAlphaFunc;
};

template <typename NumType, typename CoordType>
inline NumType
FenceReciprocalOverlapSumBoxDifferentiable<NumType, CoordType>::evaluate()
    const {
  const NumType lambda = _getLambdaFunc();
  const NumType alpha = _getAlphaFunc();
  const NumType xCell =
      op::conv<NumType>(_getVarFunc(_cellIdx, Orient2DType::HORIZONTAL));
  const NumType wCell = op::conv<NumType>(_cellWidth);
  const NumType yCell =
      op::conv<NumType>(_getVarFunc(_cellIdx, Orient2DType::VERTICAL));
  const NumType hCell = op::conv<NumType>(_cellHeight);

  NumType ovlArea = 0.0;

  for (const auto &fence : _boxes) {
    const NumType xFence = op::conv<NumType>(fence.xLo());
    const NumType wFence = op::conv<NumType>(fence.xLen());
    const NumType yFence = op::conv<NumType>(fence.yLo());
    const NumType hFence = op::conv<NumType>(fence.yLen());
    ovlArea += std::pow(alpha, 2.0) *
               std::log(1.0 / (std::exp(-(hCell + yCell - yFence) / alpha) +
                               std::exp(-(hFence - yCell + yFence) / alpha)) +
                        1.0) *
               std::log(1 / (std::exp(-(wCell + xCell - xFence) / alpha) +
                             std::exp(-(wFence - xCell + xFence) / alpha)) +
                        1.0);
  }
  return (1.0 / ovlArea) * _weight * lambda;
  ;
}

template <typename NumType, typename CoordType>
inline void FenceReciprocalOverlapSumBoxDifferentiable<
    NumType, CoordType>::accumlateGradient() const {
  const NumType lambda = _getLambdaFunc();
  const NumType alpha = _getAlphaFunc();
  const NumType xCell =
      op::conv<NumType>(_getVarFunc(_cellIdx, Orient2DType::HORIZONTAL));
  const NumType wCell = op::conv<NumType>(_cellWidth);
  const NumType yCell =
      op::conv<NumType>(_getVarFunc(_cellIdx, Orient2DType::VERTICAL));
  const NumType hCell = op::conv<NumType>(_cellHeight);

  NumType ovlArea = 0.0;
  NumType dOvlDX = 0.0;
  NumType dOvlDY = 0.0;

  for (const auto &fence : _boxes) {
    const NumType xFence = op::conv<NumType>(fence.xLo());
    const NumType wFence = op::conv<NumType>(fence.xLen());
    const NumType yFence = op::conv<NumType>(fence.yLo());
    const NumType hFence = op::conv<NumType>(fence.yLen());
    ovlArea += std::pow(alpha, 2.0) *
               std::log(1.0 / (std::exp(-(hCell + yCell - yFence) / alpha) +
                               std::exp(-(hFence - yCell + yFence) / alpha)) +
                        1.0) *
               std::log(1 / (std::exp(-(wCell + xCell - xFence) / alpha) +
                             std::exp(-(wFence - xCell + xFence) / alpha)) +
                        1.0);
    dOvlDX += (alpha * alpha *
               std::log(1.0 / (std::exp(-(hCell + yCell - yFence) / alpha) +
                               std::exp(-(hFence - yCell + yFence) / alpha)) +
                        1.0) *
               (std::exp(-(wCell + xCell - xFence) / alpha) / alpha -
                std::exp(-(wFence - xCell + xFence) / alpha) / alpha)) /
              ((1.0 / (std::exp(-(wCell + xCell - xFence) / alpha) +
                       std::exp(-(wFence - xCell + xFence) / alpha)) +
                1.0) *
               std::pow(std::exp(-(wCell + xCell - xFence) / alpha) +
                            std::exp(-(wFence - xCell + xFence) / alpha),
                        2.0));

    dOvlDY = (alpha * alpha *
              std::log(1.0 / (std::exp(-(wCell + xCell - xFence) / alpha) +
                              std::exp(-(wFence - xCell + xFence) / alpha)) +
                       1.0) *
              (std::exp(-(hCell + yCell - yFence) / alpha) / alpha -
               std::exp(-(hFence - yCell + yFence) / alpha) / alpha)) /
             ((1.0 / (std::exp(-(hCell + yCell - yFence) / alpha) +
                      std::exp(-(hFence - yCell + yFence) / alpha)) +
               1.0) *
              (std::pow(std::exp(-(hCell + yCell - yFence) / alpha) +
                            std::exp(-(hFence - yCell + yFence) / alpha),
                        2.0)));
  }

  const NumType dx = -(1.0 / std::pow(ovlArea, 2.0)) * 2 * ovlArea * dOvlDX;
  const NumType dy = -(1.0 / std::pow(ovlArea, 2.0)) * 2 * ovlArea * dOvlDY;

  _accumulateGradFunc(dx * _weight * lambda, _cellIdx,
                      Orient2DType::HORIZONTAL);
  _accumulateGradFunc(dy * _weight * lambda, _cellIdx, Orient2DType::VERTICAL);
}

template <typename op_type> struct place_fence_trait {
  typedef typename op_type::coordinate_type coordinate_type;
  static coordinate_type outFenceRegionArea(op_type &fence);
};

template <typename nlp_numerical_type, typename nlp_coordinate_type>
struct place_fence_trait<FenceReciprocalOverlapSumBoxDifferentiable<
    nlp_numerical_type, nlp_coordinate_type>> {
  typedef nlp_coordinate_type coordinate_type;
  typedef FenceReciprocalOverlapSumBoxDifferentiable<nlp_numerical_type,
                                                     nlp_coordinate_type>
      op_type;
  static coordinate_type outFenceRegionArea(op_type &fence) {
    const coordinate_type cellWidth = fence._cellWidth;
    const coordinate_type cellHeight = fence._cellHeight;
    const coordinate_type xCell =
        fence._getVarFunc(fence._cellIdx, Orient2DType::HORIZONTAL);
    const coordinate_type yCell =
        fence._getVarFunc(fence._cellIdx, Orient2DType::VERTICAL);
    coordinate_type overlapArea = 0.0;
    for (const auto &box : fence._boxes) {
      const auto xLo = box.xLo();
      const auto yLo = box.yLo();
      const auto xHi = box.xHi();
      const auto yHi = box.yHi();

      const auto overlapX =
          std::max(std::min(xCell + cellWidth, xHi) - std::max(xCell, xLo),
                   op::conv<coordinate_type>(0.0));
      const auto overlapY =
          std::max(std::min(yCell + cellHeight, yHi) - std::max(yCell, yLo),
                   op::conv<coordinate_type>(0.0));
      overlapArea += overlapX * overlapY;
    }
    return cellWidth * cellHeight -
           overlapArea; // Cell area - sum of overlapping area between cell and
                        // each splitted box
  }
};

struct FENCE_GAUSSIAN_INTEGRAL_COST {};
struct FENCE_SIGMOID_COST {};

namespace _fence_bivariate_gaussian_details {
  template<typename CostType> struct calc_trait{};

  template<>
    struct calc_trait<FENCE_SIGMOID_COST> {
      template<typename numerical_type>
        static numerical_type calcCost(
            numerical_type muX,
            numerical_type muY,
            numerical_type sigmaX,
            numerical_type sigmaY,
            numerical_type xLo,
            numerical_type yLo,
            numerical_type width,
            numerical_type height,
            numerical_type alpha,
            numerical_type 
            ) {
          using NumType = numerical_type;
          const NumType x = xLo + width / 2;
          const NumType y = yLo + height / 2;
          return .0/((std::exp(-alpha*(muX+sigmaX-x))+1.0)*(std::exp(-alpha*(-muX+sigmaX+x))+1.0)*(std::exp(-alpha*(muY+sigmaY-y))+1.0)*(std::exp(-alpha*(-muY+sigmaY+y))+1.0));;
        }
      template<typename numerical_type>
        static numerical_type calcDiffX(
            numerical_type muX,
            numerical_type muY,
            numerical_type sigmaX,
            numerical_type sigmaY,
            numerical_type xLo,
            numerical_type yLo,
            numerical_type width,
            numerical_type height,
            numerical_type alpha,
            numerical_type 
            ) {
          using NumType = numerical_type;
          const NumType x = xLo + width / 2;
          const NumType y = yLo + height / 2;
          return -(alpha*std::exp(-alpha*(muX+sigmaX-x))*1.0/pow(std::exp(-alpha*(muX+sigmaX-x))+1.0,2.0))/((std::exp(-alpha*(-muX+sigmaX+x))+1.0)*(std::exp(-alpha*(muY+sigmaY-y))+1.0)*(std::exp(-alpha*(-muY+sigmaY+y))+1.0))+(alpha*std::exp(-alpha*(-muX+sigmaX+x))*1.0/pow(std::exp(-alpha*(-muX+sigmaX+x))+1.0,2.0))/((std::exp(-alpha*(muX+sigmaX-x))+1.0)*(std::exp(-alpha*(muY+sigmaY-y))+1.0)*(std::exp(-alpha*(-muY+sigmaY+y))+1.0));
        }
      template<typename numerical_type>
        static numerical_type calcDiffY(
            numerical_type muX,
            numerical_type muY,
            numerical_type sigmaX,
            numerical_type sigmaY,
            numerical_type xLo,
            numerical_type yLo,
            numerical_type width,
            numerical_type height,
            numerical_type alpha,
            numerical_type 
            ) {
          using NumType = numerical_type;
          const NumType x = xLo + width / 2;
          const NumType y = yLo + height / 2;
          return -(alpha*std::exp(-alpha*(muY+sigmaY-y))*1.0/pow(std::exp(-alpha*(muY+sigmaY-y))+1.0,2.0))/((std::exp(-alpha*(muX+sigmaX-x))+1.0)*(std::exp(-alpha*(-muX+sigmaX+x))+1.0)*(std::exp(-alpha*(-muY+sigmaY+y))+1.0))+(alpha*std::exp(-alpha*(-muY+sigmaY+y))*1.0/pow(std::exp(-alpha*(-muY+sigmaY+y))+1.0,2.0))/((std::exp(-alpha*(muX+sigmaX-x))+1.0)*(std::exp(-alpha*(-muX+sigmaX+x))+1.0)*(std::exp(-alpha*(muY+sigmaY-y))+1.0));
        }
    };
  template<>
    struct calc_trait<FENCE_GAUSSIAN_INTEGRAL_COST> {
      template<typename numerical_type>
        static numerical_type calcCost(
            numerical_type muX,
            numerical_type muY,
            numerical_type sigmaX,
            numerical_type sigmaY,
            numerical_type xLo,
            numerical_type yLo,
            numerical_type width,
            numerical_type height,
            numerical_type ,
            numerical_type normalize
            ) {
          using NumType = numerical_type;
              std::complex<NumType> complexCost = (normalize*(op::realerf(sqrt(2.0)*sqrt(1.0/(sigmaX*sigmaX))*(-muX+width+xLo)*(1.0/2.0))+op::realerf(sqrt(2.0)*(muX-xLo)*sqrt(1.0/(sigmaX*sigmaX))*(1.0/2.0)))*1.0/sqrt(1.0/(sigmaX*sigmaX))*1.0/sqrt(std::complex<NumType>(-1.0/(sigmaY*sigmaY)))*(op::realerf(sqrt(2.0)*(1.0/(sigmaY*sigmaY)*(height+yLo)*sqrt(std::complex<NumType>(-1.0))-muY*1.0/(sigmaY*sigmaY)*sqrt(std::complex<NumType>(-1.0)))*1.0/sqrt(std::complex<NumType>(-1.0/(sigmaY*sigmaY)))*(1.0/2.0))+op::realerf(sqrt(2.0)*(muY*1.0/(sigmaY*sigmaY)*sqrt(std::complex<NumType>(-1.0))-1.0/(sigmaY*sigmaY)*yLo*sqrt(std::complex<NumType>(-1.0)))*1.0/sqrt(std::complex<NumType>(-1.0/(sigmaY*sigmaY)))*(1.0/2.0)))*2.5E-1*sqrt(std::complex<NumType>(-1.0)))/(sigmaX*sigmaY);
          complexCost /= (width * height);
          return std::real(complexCost);
        }
      template<typename numerical_type>
        static numerical_type calcDiffX(
            numerical_type muX,
            numerical_type muY,
            numerical_type sigmaX,
            numerical_type sigmaY,
            numerical_type xLo,
            numerical_type yLo,
            numerical_type width,
            numerical_type height,
            numerical_type ,
            numerical_type normalize
            ) {
          using NumType = numerical_type;
          std::complex<NumType> diff = (normalize*(sqrt(2.0)*1.0/sqrt(3.141592653589793)*exp(1.0/(sigmaX*sigmaX)*pow(muX-xLo,2.0)*(-1.0/2.0))*sqrt(1.0/(sigmaX*sigmaX))-sqrt(2.0)*1.0/sqrt(3.141592653589793)*exp(1.0/(sigmaX*sigmaX)*pow(-muX+width+xLo,2.0)*(-1.0/2.0))*sqrt(1.0/(sigmaX*sigmaX)))*1.0/sqrt(1.0/(sigmaX*sigmaX))*1.0/sqrt(std::complex<NumType>(-1.0/(sigmaY*sigmaY)))*(op::realerf(sqrt(2.0)*(1.0/(sigmaY*sigmaY)*(height+yLo)*sqrt(std::complex<NumType>(-1.0))-muY*1.0/(sigmaY*sigmaY)*sqrt(std::complex<NumType>(-1.0)))*1.0/sqrt(std::complex<NumType>(-1.0/(sigmaY*sigmaY)))*(1.0/2.0))+op::realerf(sqrt(2.0)*(muY*1.0/(sigmaY*sigmaY)*sqrt(std::complex<NumType>(-1.0))-1.0/(sigmaY*sigmaY)*yLo*sqrt(std::complex<NumType>(-1.0)))*1.0/sqrt(std::complex<NumType>(-1.0/(sigmaY*sigmaY)))*(1.0/2.0)))*-2.5E-1*sqrt(std::complex<NumType>(-1.0)))/(sigmaX*sigmaY);        
          diff /= (width * height);
          return std::real(diff);
        }
      template<typename numerical_type>
        static numerical_type calcDiffY(
            numerical_type muX,
            numerical_type muY,
            numerical_type sigmaX,
            numerical_type sigmaY,
            numerical_type xLo,
            numerical_type yLo,
            numerical_type width,
            numerical_type height,
            numerical_type ,
            numerical_type normalize
            ) {
          using NumType = numerical_type;
          std::complex<NumType> diff = (normalize*(sqrt(2.0)*1.0/(sigmaY*sigmaY)*1.0/sqrt(3.141592653589793)*exp((sigmaY*sigmaY)*pow(1.0/(sigmaY*sigmaY)*(height+yLo)*sqrt(std::complex<NumType>(-1.0))-muY*1.0/(sigmaY*sigmaY)*sqrt(std::complex<NumType>(-1.0)),2.0)*(1.0/2.0))*1.0/sqrt(std::complex<NumType>(-1.0/(sigmaY*sigmaY)))*sqrt(std::complex<NumType>(-1.0))-sqrt(std::complex<NumType>(2.0))*1.0/(sigmaY*sigmaY)*1.0/sqrt(3.141592653589793)*exp((sigmaY*sigmaY)*pow(muY*1.0/(sigmaY*sigmaY)*sqrt(std::complex<NumType>(-1.0))-1.0/(sigmaY*sigmaY)*yLo*sqrt(std::complex<NumType>(-1.0)),2.0)*(1.0/2.0))*1.0/sqrt(std::complex<NumType>(-1.0/(sigmaY*sigmaY)))*sqrt(std::complex<NumType>(-1.0)))*(op::realerf(sqrt(2.0)*sqrt(1.0/(sigmaX*sigmaX))*(-muX+width+xLo)*(1.0/2.0))+op::realerf(sqrt(2.0)*(muX-xLo)*sqrt(1.0/(sigmaX*sigmaX))*(1.0/2.0)))*1.0/sqrt(1.0/(sigmaX*sigmaX))*1.0/sqrt(std::complex<NumType>(-1.0/(sigmaY*sigmaY)))*2.5E-1*sqrt(std::complex<NumType>(-1.0)))/(sigmaX*sigmaY);
          diff /= (width * height);
          return std::real(diff);
        }
    };
} //namespace _fence_bivariate_gaussian_details

/// @brief Model fence region penalty as the sum of a set of Gaussian distribution.
/// it transformed the well bounding box into a Gaussian and the cost is the integral
// of the cost over cell area
template <typename NumType, typename CoordType, typename CostType=FENCE_SIGMOID_COST>
struct FenceBivariateGaussianDifferentiable {
  typedef NumType numerical_type;
  typedef CoordType coordinate_type;

  FenceBivariateGaussianDifferentiable(
      const std::vector<IndexType> &inFenceCellIdx,
      const std::vector<IndexType> &outFenceCellIdx,
      const std::vector<CoordType> &inFenceCellWidths,
      const std::vector<CoordType> &inFenceCellHeights,
      const std::vector<CoordType> &outFenceCellWidths,
      const std::vector<CoordType> &outFenceCellHeights,
      const std::function<NumType(void)> &getLambdaFunc) 
  : _inFenceCellIdx(inFenceCellIdx),_outFenceCellIdx(outFenceCellIdx),  
  _inFenceCellWidths(inFenceCellWidths), _inFenceCellHeights(inFenceCellHeights),
  _outFenceCellWidths(outFenceCellWidths), _outFenceCellHeights(outFenceCellHeights),
  _getLambdaFunc(getLambdaFunc) {
    _getLambdaFunc = getLambdaFunc;
    Assert(inFenceCellIdx.size() == inFenceCellWidths.size());
    Assert(inFenceCellIdx.size() == inFenceCellHeights.size());
    Assert(outFenceCellIdx.size() == outFenceCellWidths.size());
    Assert(outFenceCellIdx.size() == outFenceCellHeights.size());
    _getAlphaFunc = [](){ return 1.0; };
  }

  void setGetVarFunc(
      const std::function<CoordType(IndexType, Orient2DType)> &getVarFunc) {
    _getVarFunc = getVarFunc;
  }
  void setAccumulateGradFunc(
      const std::function<void(NumType, IndexType, Orient2DType)> &func) {
    _accumulateGradFunc = func;
  }
  void setGetAlphaFunc(const std::function<NumType(void)> &alphaFunc) {
    _getAlphaFunc = alphaFunc;
  }

  NumType evaluate() const;
  void accumlateGradient() const;


  std::vector<BivariateGaussianParameters<NumType>> _gaussianParameters;
  std::vector<IndexType> _inFenceCellIdx;
  std::vector<IndexType> _outFenceCellIdx;
  std::vector<CoordType> _inFenceCellWidths;
  std::vector<CoordType> _inFenceCellHeights;
  std::vector<CoordType> _outFenceCellWidths;
  std::vector<CoordType> _outFenceCellHeights;
  std::function<NumType(void)>
      _getLambdaFunc; ///< A function to get the current lambda multiplier
  std::function<CoordType(IndexType cellIdx, Orient2DType orient)>
      _getVarFunc; ///< A function to get current variable value
  std::function<void(NumType, IndexType, Orient2DType)>
      _accumulateGradFunc; ///< A function to update partial
  static constexpr NumType tol = 1.0; ///< In fence use tol - integral(PDF) as cost
  BoolType _considerOutFenceCells = false; ///< Wether calculate the cost for out of the fence cost
  NumType _weight = 1.0;
  std::function<NumType(void)> _getAlphaFunc;
};

template <typename NumType, typename CoordType, typename CostType>
inline NumType
FenceBivariateGaussianDifferentiable<NumType, CoordType, CostType>::evaluate()
    const {
      const NumType lambda = _getLambdaFunc();
      const NumType alpha = _getAlphaFunc();
      std::vector<NumType> costOut;
      costOut.resize(_inFenceCellIdx.size() + _outFenceCellIdx.size(), 0);
#pragma omp parallel for schedule(static)
      for (IndexType idx = 0; idx < _inFenceCellIdx.size(); ++idx) {
        IndexType cellIdx = _inFenceCellIdx[idx];
        const NumType xLo =
            op::conv<NumType>(_getVarFunc(cellIdx, Orient2DType::HORIZONTAL));
        const NumType width = op::conv<NumType>(_inFenceCellWidths[idx]);
        const NumType yLo =
            op::conv<NumType>(_getVarFunc(cellIdx, Orient2DType::VERTICAL));
        const NumType height = op::conv<NumType>(_inFenceCellHeights[idx]);
        for (IndexType gauIdx = 0; gauIdx < _gaussianParameters.size(); ++gauIdx) {
          const NumType muX = _gaussianParameters[gauIdx].muX;
          const NumType muY = _gaussianParameters[gauIdx].muY;
          const NumType sigmaX = _gaussianParameters[gauIdx].sigmaX * _gaussianParameters[gauIdx].extention;
          const NumType sigmaY = _gaussianParameters[gauIdx].sigmaY * _gaussianParameters[gauIdx].extention;
          const NumType normalize =_gaussianParameters[gauIdx].normalize;
          // Calculate cost
          const auto cost = _fence_bivariate_gaussian_details::calc_trait<CostType>::calcCost(
              muX, muY, sigmaX, sigmaY, xLo, yLo, width, height, alpha, normalize
              );
          costOut[idx] += lambda * (tol  - cost ) * _weight;
        }
      }
      if (not _considerOutFenceCells) {
        return std::accumulate(costOut.begin(), costOut.end(), 0.0);
      }
#pragma omp parallel for schedule(static)
      for (IndexType idx = 0; idx < _outFenceCellIdx.size(); ++idx) {
        IndexType cellIdx = _outFenceCellIdx[idx];
        const NumType xLo =
            op::conv<NumType>(_getVarFunc(cellIdx, Orient2DType::HORIZONTAL));
        const NumType width = op::conv<NumType>(_outFenceCellWidths[idx]);
        const NumType yLo =
            op::conv<NumType>(_getVarFunc(cellIdx, Orient2DType::VERTICAL));
        const NumType height = op::conv<NumType>(_outFenceCellHeights[idx]);
        for (IndexType gauIdx = 0; gauIdx < _gaussianParameters.size(); ++gauIdx) {
          const NumType muX = _gaussianParameters[gauIdx].muX;
          const NumType muY = _gaussianParameters[gauIdx].muY;
          const NumType sigmaX = _gaussianParameters[gauIdx].sigmaX;
          const NumType sigmaY = _gaussianParameters[gauIdx].sigmaY;
          const NumType normalize =_gaussianParameters[gauIdx].normalize;
          // Calculate cost
          const auto cost = _fence_bivariate_gaussian_details::calc_trait<CostType>::calcCost(
              muX, muY, sigmaX, sigmaY, xLo, yLo, width, height, alpha, normalize
              );

          costOut[idx + _inFenceCellIdx.size()] += lambda * cost * _weight;
        }
      }
      return std::accumulate(costOut.begin(), costOut.end(), 0.0);
  }

template <typename NumType, typename CoordType, typename CostType>
inline void
FenceBivariateGaussianDifferentiable<NumType, CoordType, CostType>::accumlateGradient()
     const {
      const NumType lambda = _getLambdaFunc();
      const NumType alpha = _getAlphaFunc();
#pragma omp parallel for schedule(static)
      for (IndexType idx = 0; idx < _inFenceCellIdx.size(); ++idx) {
        IndexType cellIdx = _inFenceCellIdx[idx];
        const NumType xLo =
            op::conv<NumType>(_getVarFunc(cellIdx, Orient2DType::HORIZONTAL));
        const NumType width = op::conv<NumType>(_inFenceCellWidths[idx]);
        const NumType yLo =
            op::conv<NumType>(_getVarFunc(cellIdx, Orient2DType::VERTICAL));
        const NumType height = op::conv<NumType>(_inFenceCellHeights[idx]);
        for (IndexType gauIdx = 0; gauIdx < _gaussianParameters.size(); ++gauIdx) {
          const NumType muX = _gaussianParameters[gauIdx].muX;
          const NumType muY = _gaussianParameters[gauIdx].muY;
          const NumType sigmaX = _gaussianParameters[gauIdx].sigmaX * _gaussianParameters[gauIdx].extention;
          const NumType sigmaY = _gaussianParameters[gauIdx].sigmaY * _gaussianParameters[gauIdx].extention;
          const NumType normalize = _gaussianParameters[gauIdx].normalize;
          // Calculate derivative
          const auto diffx = _fence_bivariate_gaussian_details::calc_trait<CostType>::calcDiffX(
              muX, muY, sigmaX, sigmaY, xLo, yLo, width, height, alpha, normalize
              );
          const auto diffy = _fence_bivariate_gaussian_details::calc_trait<CostType>::calcDiffY(
              muX, muY, sigmaX, sigmaY, xLo, yLo, width, height, alpha, normalize
              );
          _accumulateGradFunc( - lambda * diffx *_weight, cellIdx, Orient2DType::HORIZONTAL);
          _accumulateGradFunc( - lambda * diffy * _weight, cellIdx, Orient2DType::VERTICAL);
        }
      }
      if (not _considerOutFenceCells) {
        return;
      }
#pragma omp parallel for schedule(static)
      for (IndexType idx = 0; idx < _outFenceCellIdx.size(); ++idx) {
        IndexType cellIdx = _outFenceCellIdx[idx];
        const NumType xLo =
            op::conv<NumType>(_getVarFunc(cellIdx, Orient2DType::HORIZONTAL));
        const NumType width = op::conv<NumType>(_outFenceCellWidths[idx]);
        const NumType yLo =
            op::conv<NumType>(_getVarFunc(cellIdx, Orient2DType::VERTICAL));
        const NumType height = op::conv<NumType>(_outFenceCellHeights[idx]);
        for (IndexType gauIdx = 0; gauIdx < _gaussianParameters.size(); ++gauIdx) {
          const NumType muX = _gaussianParameters[gauIdx].muX;
          const NumType muY = _gaussianParameters[gauIdx].muY;
          const NumType sigmaX = _gaussianParameters[gauIdx].sigmaX;
          const NumType sigmaY = _gaussianParameters[gauIdx].sigmaY;
          const NumType normalize = _gaussianParameters[gauIdx].normalize;
          // Calculate cost
          const auto diffx = _fence_bivariate_gaussian_details::calc_trait<CostType>::calcDiffX(
              muX, muY, sigmaX, sigmaY, xLo, yLo, width, height, alpha, normalize
              );
          const auto diffy = _fence_bivariate_gaussian_details::calc_trait<CostType>::calcDiffY(
              muX, muY, sigmaX, sigmaY, xLo, yLo, width, height, alpha, normalize
              );
          _accumulateGradFunc(lambda * diffx *_weight, cellIdx, Orient2DType::HORIZONTAL);
          _accumulateGradFunc(lambda * diffy * _weight, cellIdx, Orient2DType::VERTICAL);
        }
      }
  }


/// @brief LSE-smoothed Area
template <typename NumType, typename CoordType> struct LseAreaDifferentiable {
  typedef NumType numerical_type;
  typedef CoordType coordinate_type;

  LseAreaDifferentiable(const std::function<NumType(void)> &getAlphaFunc,
                        const std::function<NumType(void)> &getLambdaFunc) {
    _getAlphaFunc = getAlphaFunc;
    _getLambdaFunc = getLambdaFunc;
  }

  void setGetVarFunc(
      const std::function<CoordType(IndexType, Orient2DType)> &getVarFunc) {
    _getVarFunc = getVarFunc;
  }
  void setAccumulateGradFunc(
      const std::function<void(NumType, IndexType, Orient2DType)> &func) {
    _accumulateGradFunc = func;
  }
  void setGetAlphaFunc(const std::function<NumType(void)> &getAlphaFunc) {
    _getAlphaFunc = getAlphaFunc;
  }

  void addVar(IndexType cellIdx, CoordType width, CoordType height) {
    _cells.emplace_back(cellIdx);
    _cellWidths.emplace_back(op::conv<NumType>(width));
    _cellHeights.emplace_back(op::conv<NumType>(height));
  }
  void setWeight(const NumType &weight) { _weight = weight; }

  NumType evaluate() const {
    const NumType lambda = _getLambdaFunc();
    const NumType alpha = _getAlphaFunc();
    std::array<NumType, 4> sumExpArray = {0.0, 0.0, 0.0, 0.0}; // sum exp(x/alpha), exp(-x/ alpha), exp(y/alpha), exp(-y/alpha)
    NumType *pSumExp = &sumExpArray.front();
    for (IndexType idx = 0; idx < _cells.size(); ++idx) {
      IndexType cellIdx = _cells[idx];
      const NumType x = op::conv<NumType>(
          _getVarFunc(cellIdx, Orient2DType::HORIZONTAL));
      const NumType y = op::conv<NumType>(
          _getVarFunc(cellIdx, Orient2DType::VERTICAL));
      const NumType xHi = x + _cellWidths[idx];
      const NumType yHi = y + _cellHeights[idx];
      pSumExp[0] += std::exp(x / alpha) + std::exp(xHi / alpha);
      pSumExp[1] += std::exp(-x / alpha) + std::exp(-yHi / alpha);
      pSumExp[2] += std::exp(y / alpha) + std::exp(yHi / alpha);
      pSumExp[3] += std::exp(-y / alpha) + std::exp(-yHi / alpha);
    }
    const NumType cost = (alpha * (std::log(pSumExp[0]) + std::log(pSumExp[1]))) * (alpha * (std::log(pSumExp[2]) + std::log(pSumExp[3]))) * lambda * _weight;
    return std::sqrt(cost);
  }

  void accumlateGradient() const {
    const NumType lambda = _getLambdaFunc();
    const NumType alpha = _getAlphaFunc();
    std::array<NumType, 4> sumExpArray = {0.0, 0.0, 0.0, 0.0}; // sum exp(x/alpha), exp(-x/ alpha), exp(y/alpha), exp(-y/alpha)
    NumType *pSumExp = &sumExpArray.front();
    for (IndexType idx = 0; idx < _cells.size(); ++idx) {
      IndexType cellIdx = _cells[idx];
      const NumType x = op::conv<NumType>(
          _getVarFunc(cellIdx, Orient2DType::HORIZONTAL));
      const NumType y = op::conv<NumType>(
          _getVarFunc(cellIdx, Orient2DType::VERTICAL));
      const NumType xHi = x + _cellWidths[idx];
      const NumType yHi = y + _cellHeights[idx];
      pSumExp[0] += std::exp(x / alpha) + std::exp(xHi / alpha);
      pSumExp[1] += std::exp(-x / alpha) + std::exp(-yHi / alpha);
      pSumExp[2] += std::exp(y / alpha) + std::exp(yHi / alpha);
      pSumExp[3] += std::exp(-y / alpha) + std::exp(-yHi / alpha);
    }
    std::array<NumType, 2> logSumExpArray = {0.0, 0.0}; // log exp sum x + log exp sum -x, log exp sum y + log exp sum -y
    logSumExpArray[0] = std::log(sumExpArray[0]) + std::log(sumExpArray[1]);
    logSumExpArray[1] = std::log(sumExpArray[2]) + std::log(sumExpArray[3]);
    const NumType costScale = 1 / std::sqrt((alpha * (std::log(pSumExp[0]) + std::log(pSumExp[1]))) * (alpha * (std::log(pSumExp[2]) + std::log(pSumExp[3]))) * lambda * _weight);
    for (IndexType idx = 0; idx < _cells.size(); ++idx) {
      IndexType cellIdx = _cells[idx];
      const NumType x = op::conv<NumType>(
          _getVarFunc(cellIdx, Orient2DType::HORIZONTAL));
      const NumType y = op::conv<NumType>(
          _getVarFunc(cellIdx, Orient2DType::VERTICAL));
      const NumType xHi = x + _cellWidths[idx];
      const NumType yHi = y + _cellHeights[idx];
      const NumType diffxLo = alpha * logSumExpArray[1] * ( std::exp(x / alpha ) /  sumExpArray[0] - std::exp(-x / alpha ) /  sumExpArray[1]);
      const NumType diffyLo = alpha * logSumExpArray[0] * ( std::exp(y / alpha ) /  sumExpArray[2] - std::exp(-y / alpha ) /  sumExpArray[3]);
      const NumType diffxHi = alpha * logSumExpArray[1] * ( std::exp(xHi / alpha ) /  sumExpArray[0] - std::exp(-xHi / alpha ) /  sumExpArray[1]);
      const NumType diffyHi = alpha * logSumExpArray[0] * ( std::exp(yHi / alpha ) /  sumExpArray[2] - std::exp(-yHi / alpha ) /  sumExpArray[3]);
      _accumulateGradFunc(costScale * lambda * (diffxLo + diffxHi) *_weight, cellIdx, Orient2DType::HORIZONTAL);
      _accumulateGradFunc(costScale * lambda * (diffyLo + diffyHi) * _weight, cellIdx, Orient2DType::VERTICAL);
    }
  }


  std::vector<IndexType> _cells;
  std::vector<NumType> _cellWidths;
  std::vector<NumType> _cellHeights;
  NumType _weight = 1;
  std::function<NumType(void)>
      _getAlphaFunc; ///< A function to get the current alpha
  std::function<NumType(void)>
      _getLambdaFunc; ///< A function to get the current lambda multiplier
  std::function<CoordType(IndexType cellIdx, Orient2DType orient)>
      _getVarFunc; ///< A function to get current variable value
  std::function<void(NumType, IndexType, Orient2DType)>
      _accumulateGradFunc; ///< A function to update partial
};
} // namespace diff

PROJECT_NAMESPACE_END

#endif // IDEAPLACE_DIFFERENT_H_
