/**
 * @file different.h
 * @brief The numerical differentiable concepts and implementations
 * @author Keren Zhu
 * @date 02/25/2020
 */

#ifndef IDEAPLACE_DIFFERENT_H_
#define IDEAPLACE_DIFFERENT_H_

#include "db/Database.h"

PROJECT_NAMESPACE_BEGIN

namespace diff {

enum class OpEnumType
{
    hpwl, ovl, oob, asym, cosine
};
struct placement_differentiable_concept {};

template <typename ConceptType>
struct is_placement_differentiable_concept
{
    typedef std::false_type is_placement_differentiable_concept_type;
};

template <>
struct is_placement_differentiable_concept<placement_differentiable_concept>
{
    typedef std::true_type  is_placement_differentiable_concept_type;
};

template <typename DifType>
struct placement_differentiable_traits
{
    typedef DifType different_type;
    typedef typename different_type::numerical_type numerical_type;
    typedef typename different_type::coordinate_type coordinate_type;

    static numerical_type evaluate( const different_type & dif, const std::function<coordinate_type(IndexType, Orient2DType)> &getVarFunc) 
    {
        return dif.evaluate(getVarFunc);
    }

    static void accumlateGradient(const different_type & dif, 
            const std::function<coordinate_type(IndexType cellIdx, Orient2DType orient)> &getVarFunc,
            const std::function<void(numerical_type, IndexType, Orient2DType)> &accumulateGradFunc)
    {
        dif.accumlateGradient(getVarFunc, accumulateGradFunc);
    }

};

/// @namespace IDEAPLACE::op
/// @brief namespace for the operators
namespace op
{
    /// @brief abs smooth function with respect to 0
    /// @return the smoothed result
    template<typename NumType>
    constexpr NumType logSumExp0(NumType var, NumType alpha)
    {
        return alpha * log( exp( var / alpha) + 1);
    }

    /// @brief calculate the log sum exp, used to smooth min max function
    /// @return the calculated log sum exp
    template<typename NumType>
    constexpr NumType logSumExp(NumType var1, NumType var2, NumType alpha)
    {
        return alpha * log(exp( var1 / alpha) + exp( var2 / alpha));
    }

    /// @brief The partial of LSE(0, var) with respect to var
    /// @return the gradient of logSumExp0
    template<typename NumType>
    constexpr NumType gradLogSumExp0(NumType var, NumType alpha)
    {
        return exp( var / alpha ) / ( exp( var / alpha) + 1);
    }

    namespace _conv_details
    {
        template<typename, typename, bool>
        struct _conv_t {};
        template<typename LhsType, typename _RhsType>
        struct _conv_t<LhsType, _RhsType, true>
        {
            static constexpr LhsType _conv(_RhsType rhs)
            {
                return rhs;
            }
        };
        template<typename LhsType, typename _RhsType>
        struct _conv_t<LhsType, _RhsType, false>
        {
            static constexpr LhsType _conv(_RhsType rhs)
            {
                return static_cast<LhsType>(rhs);
            }
        };
    } // namespace conv_details

    /// @brief convert rhs type to lhs type
    template<typename LhsType, typename _RhsType>
    constexpr LhsType conv(_RhsType rhs)
    {
        return _conv_details::_conv_t<LhsType, _RhsType, std::is_same<LhsType, _RhsType>::value>::_conv(rhs);
    }
};


/// @brief LSE-smoothed HPWL
template<typename NumType, typename CoordType>
struct LseHpwlDifferentiable
{
    typedef NumType numerical_type;
    typedef CoordType coordinate_type;

    LseHpwlDifferentiable(const std::function<NumType(void)> &getAlphaFunc, const std::function<NumType(void)> &getLambdaFunc) 
    { _getAlphaFunc = getAlphaFunc; _getLambdaFunc = getLambdaFunc; }

    void setVirtualPin(const CoordType &x, const CoordType &y) 
    { 
        _validVirtualPin = 1; 
        _virtualPinX = x;
        _virtualPinY = y;
    }
    void removeVirtualPin() { _validVirtualPin = 0; }
    void addVar(IndexType cellIdx, const CoordType &offsetX, const CoordType &offsetY)
    {
        _cells.emplace_back(cellIdx);
        _offsetX.emplace_back(offsetX);
        _offsetY.emplace_back(offsetY);
    }
    void setWeight(const NumType &weight) { _weight = weight; }
    bool validHpwl() const { return _cells.size() + _validVirtualPin > 1;}


    NumType evaluate( const std::function<CoordType(IndexType, Orient2DType)> &getVarFunc) const
    {
        if (! validHpwl())
        {
            return 0;
        }
        std::array<NumType, 4> max_val = { 0, 0, 0, 0}; // xmax xin ymax ymin
        auto alpha = _getAlphaFunc();
        auto lambda = _getLambdaFunc();
        NumType *pMax = &max_val.front(); 
        for (IndexType pinIdx = 0; pinIdx < _cells.size(); ++pinIdx)
        {
            NumType x = op::conv<NumType>(
                    getVarFunc(_cells[pinIdx], Orient2DType::HORIZONTAL) + _offsetX[pinIdx]
                    );
            NumType y = op::conv<NumType>(
                    getVarFunc(_cells[pinIdx], Orient2DType::VERTICAL) + _offsetY[pinIdx]
                    );
            pMax[0] += exp(x / alpha);
            pMax[1] += exp(- x / alpha);
            pMax[2] += exp(y / alpha);
            pMax[3] += exp(- y / alpha);
        }
        if (_validVirtualPin == 1)
        {
            pMax[0] += exp(_virtualPinX / alpha);
            pMax[1] += exp(- _virtualPinX / alpha);
            pMax[2] += exp(_virtualPinY / alpha);
            pMax[3] += exp(- _virtualPinY / alpha);
        }
        NumType obj = 0;
        for (int i = 0; i < 4; ++ i)
        {
            obj += log(pMax[i]);
        }
        return alpha * obj * _weight * lambda;
    }

    void accumlateGradient(
            const std::function<CoordType(IndexType cellIdx, Orient2DType orient)> &getVarFunc,
            const std::function<void(NumType, IndexType, Orient2DType)> &accumulateGradFunc) const
    {
        if (! validHpwl())
        {
            return;
        }
        std::array<NumType, 4> max_val = { 0, 0, 0, 0}; // xmax xin ymax ymin
        auto alpha = _getAlphaFunc();
        auto lambda = _getLambdaFunc();
        NumType *pMax = &max_val.front(); 
        std::vector<std::array<NumType, 4>> exp_results(_cells.size());
        for (IndexType pinIdx = 0; pinIdx < _cells.size(); ++pinIdx)
        {
            NumType x = op::conv<NumType>(
                    getVarFunc(_cells[pinIdx], Orient2DType::HORIZONTAL) + _offsetX[pinIdx]
                    );
            NumType y = op::conv<NumType>(
                    getVarFunc(_cells[pinIdx], Orient2DType::VERTICAL) + _offsetY[pinIdx]
                    );
            exp_results[pinIdx][0] = exp(x / alpha);
            pMax[0] += exp_results[pinIdx][0];
            exp_results[pinIdx][1] = exp(- x / alpha);
            pMax[1] += exp_results[pinIdx][1];
            exp_results[pinIdx][2] = exp( y / alpha);
            pMax[2] += exp_results[pinIdx][2];
            exp_results[pinIdx][3] = exp(- y / alpha);
            pMax[3] += exp_results[pinIdx][3];
        }
        if (_validVirtualPin == 1)
        {
            pMax[0] += exp(_virtualPinX / alpha);
            pMax[1] += exp(- _virtualPinX / alpha);
            pMax[2] += exp(_virtualPinY / alpha);
            pMax[3] += exp(- _virtualPinY / alpha);
        }
        for (IndexType pinIdx = 0; pinIdx < _cells.size(); ++pinIdx)
        {
            IndexType cellIdx = _cells[pinIdx];
            NumType xPartial = lambda * _weight;
            NumType yPartial = xPartial;
            xPartial *= (exp_results[pinIdx][0] / pMax[0]) - (exp_results[pinIdx][1] / pMax[1]);
            yPartial *= (exp_results[pinIdx][2] / pMax[2]) - (exp_results[pinIdx][3] / pMax[3]);
            accumulateGradFunc(xPartial, cellIdx, Orient2DType::HORIZONTAL);
            accumulateGradFunc(yPartial, cellIdx, Orient2DType::VERTICAL);
        }
    }

    IntType _validVirtualPin = 0;
    CoordType _virtualPinX = 0;
    CoordType _virtualPinY = 0;
    std::vector<IndexType> _cells;
    std::vector<CoordType> _offsetX;
    std::vector<CoordType> _offsetY;
    NumType _weight = 1;
    std::function<NumType(void)> _getAlphaFunc; ///< A function to get the current alpha
    std::function<NumType(void)> _getLambdaFunc; ///< A function to get the current lambda multiplier
};


template <typename NumType, typename CoordType>
struct is_placement_differentiable_concept<LseHpwlDifferentiable<NumType, CoordType>>
{
    typedef std::true_type  is_placement_differentiable_concept_type;
};

// @brief pair-wise cell overlapping penalty
template<typename NumType, typename CoordType>
struct CellPairOverlapPenaltyDifferentiable
{
    typedef NumType numerical_type;
    typedef CoordType coordinate_type;

    CellPairOverlapPenaltyDifferentiable(IndexType cellIdxI, CoordType cellWidthI, CoordType cellHeightI,
                           IndexType cellIdxJ, CoordType cellWidthJ, CoordType cellHeightJ,
                           const std::function<NumType(void)> &getAlphaFunc, const std::function<NumType(void)> &getLambdaFunc)
    {
        _cellIdxI = cellIdxI;
        _cellWidthI = cellWidthI;
        _cellHeightI = cellHeightI;
        _cellIdxJ = cellIdxJ;
        _cellWidthJ = cellWidthJ;
        _cellHeightJ = cellHeightJ;
        _getAlphaFunc = getAlphaFunc;
        _getLambdaFunc = getLambdaFunc;
    }

    NumType evaluate( const std::function<CoordType(IndexType, Orient2DType)> &getVarFunc) const
    {
        NumType alpha = _getAlphaFunc();
        NumType lambda  = _getLambdaFunc();
        CoordType xLoI = getVarFunc(_cellIdxI, Orient2DType::HORIZONTAL);
        CoordType xHiI = xLoI + _cellWidthI;
        CoordType yLoI = getVarFunc(_cellIdxI, Orient2DType::VERTICAL);
        CoordType yHiI = yLoI + _cellHeightI;
        CoordType xLoJ = getVarFunc(_cellIdxJ, Orient2DType::HORIZONTAL);
        CoordType xHiJ = xLoJ + _cellWidthJ;
        CoordType yLoJ = getVarFunc(_cellIdxJ, Orient2DType::VERTICAL);
        CoordType yHiJ = yLoJ + _cellHeightJ;
        // max (min(xHiI - xLoJ, xHiJ - xLoI), 0), vice versa
        // Notice that the calculation results for changing the order of i and j.
        // In the gradient, the results will be different
        CoordType var1X = xHiI - xLoJ;
        CoordType var2X = xHiJ - xLoI;
        NumType overlapX = op::logSumExp(
            op::conv<NumType>(var1X),
            op::conv<NumType>(var2X),
            - alpha
            );
        overlapX = op::logSumExp0(
                overlapX,
                alpha
                );
        // y
        CoordType var1Y = yHiI - yLoJ;
        CoordType var2Y = yHiJ - yLoI;
        NumType overlapY = op::logSumExp(
            op::conv<NumType>(var1Y),
            op::conv<NumType>(var2Y),
            - alpha
            );
        overlapY = op::logSumExp0(
                overlapY,
                alpha
                );
        return lambda * overlapX * overlapY;
    }

    void accumlateGradient(
            const std::function<CoordType(IndexType cellIdx, Orient2DType orient)> &getVarFunc,
            const std::function<void(NumType, IndexType, Orient2DType)> &accumulateGradFunc) const
    {
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
         * min_func_y = -alpha * log( exp( - var1y / alpha) + exp( - var2y / alpha));
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
        NumType dxi_gold = (alpha * alpha *log(1/(exp(-(hi + yi - yj)/alpha) + exp(-(hj - yi + yj)/alpha)) + 1)*(exp(-(wi + xi - xj)/alpha)/alpha - exp(-(wj - xi + xj)/alpha)/alpha))/((1/(exp(-(wi + xi - xj)/alpha) + exp(-(wj - xi + xj)/alpha)) + 1)*pow(exp(-(wi + xi - xj)/alpha) + exp(-(wj - xi + xj)/alpha), 2))
            ;
         NumType dxj_gold = 
             -(alpha * alpha *log(1/(exp(-(hi + yi - yj)/alpha) + exp(-(hj - yi + yj)/alpha)) + 1)*(exp(-(wi + xi - xj)/alpha)/alpha - exp(-(wj - xi + xj)/alpha)/alpha))/((1/(exp(-(wi + xi - xj)/alpha) + exp(-(wj - xi + xj)/alpha)) + 1)*pow((exp(-(wi + xi - xj)/alpha) + exp(-(wj - xi + xj)/alpha)), 2))
             ;


        NumType dyi_gold = 
            (alpha * alpha *log(1/(exp(-(wi + xi - xj)/alpha) + exp(-(wj - xi + xj)/alpha)) + 1)*(exp(-(hi + yi - yj)/alpha)/alpha - exp(-(hj - yi + yj)/alpha)/alpha))/((1/(exp(-(hi + yi - yj)/alpha) + exp(-(hj - yi + yj)/alpha)) + 1)*(pow(exp(-(hi + yi - yj)/alpha) + exp(-(hj - yi + yj)/alpha), 2)))
            ;

        NumType dyj_gold = 
            -(alpha * alpha*log(1/(exp(-(wi + xi - xj)/alpha) + exp(-(wj - xi + xj)/alpha)) + 1)*(exp(-(hi + yi - yj)/alpha)/alpha - exp(-(hj - yi + yj)/alpha)/alpha))/((1/(exp(-(hi + yi - yj)/alpha) + exp(-(hj - yi + yj)/alpha)) + 1)*(pow((exp(-(hi + yi - yj)/alpha) + exp(-(hj - yi + yj)/alpha)), 2)))
         *
         */
        NumType xi = op::conv<NumType>(getVarFunc(_cellIdxI, Orient2DType::HORIZONTAL));
        NumType yi = op::conv<NumType>(getVarFunc(_cellIdxI, Orient2DType::VERTICAL));
        NumType xj = op::conv<NumType>(getVarFunc(_cellIdxJ, Orient2DType::HORIZONTAL));
        NumType yj = op::conv<NumType>(getVarFunc(_cellIdxJ, Orient2DType::VERTICAL));
        NumType wi = op::conv<NumType>(_cellWidthI);
        NumType hi = op::conv<NumType>(_cellHeightI);
        NumType wj = op::conv<NumType>(_cellWidthJ);
        NumType hj = op::conv<NumType>(_cellHeightJ);
        NumType alpha = _getAlphaFunc();
        NumType lambda = _getLambdaFunc();

        NumType dxi = (alpha * alpha *log(1/(exp(-(hi + yi - yj)/alpha) + exp(-(hj - yi + yj)/alpha)) + 1)*(exp(-(wi + xi - xj)/alpha)/alpha - exp(-(wj - xi + xj)/alpha)/alpha))/((1/(exp(-(wi + xi - xj)/alpha) + exp(-(wj - xi + xj)/alpha)) + 1)*pow(exp(-(wi + xi - xj)/alpha) + exp(-(wj - xi + xj)/alpha), 2)) * lambda;
         NumType dxj = -dxi;

        NumType dyi = (alpha * alpha *log(1/(exp(-(wi + xi - xj)/alpha) + exp(-(wj - xi + xj)/alpha)) + 1)*(exp(-(hi + yi - yj)/alpha)/alpha - exp(-(hj - yi + yj)/alpha)/alpha))/((1/(exp(-(hi + yi - yj)/alpha) + exp(-(hj - yi + yj)/alpha)) + 1)*(pow(exp(-(hi + yi - yj)/alpha) + exp(-(hj - yi + yj)/alpha), 2))) * lambda;
        NumType dyj = -dyi;


        // accumulate the computed partials
        accumulateGradFunc(dxi , _cellIdxI, Orient2DType::HORIZONTAL);
        accumulateGradFunc(dxj , _cellIdxJ, Orient2DType::HORIZONTAL);
        accumulateGradFunc(dyi , _cellIdxI, Orient2DType::VERTICAL);
        accumulateGradFunc(dyj , _cellIdxJ, Orient2DType::VERTICAL);
    }

    IndexType _cellIdxI;
    CoordType _cellWidthI;
    CoordType _cellHeightI;
    IndexType _cellIdxJ;
    CoordType _cellWidthJ;
    CoordType _cellHeightJ;
    std::function<NumType(void)> _getAlphaFunc; ///< A function to get the current alpha
    std::function<NumType(void)> _getLambdaFunc; ///< A function to get the current lambda multiplier
};


template <typename NumType, typename CoordType>
struct is_placement_differentiable_concept<CellPairOverlapPenaltyDifferentiable<NumType, CoordType>>
{
    typedef std::true_type  is_placement_differentiable_concept_type;
};

/// @brief the cell out of boundary penalty
template<typename NumType, typename CoordType>
struct CellOutOfBoundaryPenaltyDifferentiable
{
    typedef NumType numerical_type;
    typedef CoordType coordinate_type;

    CellOutOfBoundaryPenaltyDifferentiable(IndexType cellIdx, CoordType cellWidth, CoordType cellHeight, Box<CoordType> *boundary,
            const std::function<NumType(void)> &getAlphaFunc, const std::function<NumType(void)> &getLambdaFunc)
    {
        _cellIdx = cellIdx;
        _cellWidth = cellWidth;
        _cellHeight = cellHeight;
        _boundary = boundary;
        _getAlphaFunc = getAlphaFunc;
        _getLambdaFunc = getLambdaFunc;
    }

    NumType evaluate( const std::function<CoordType(IndexType, Orient2DType)> &getVarFunc) const
    {
        NumType alpha  = _getAlphaFunc();
        NumType lambda = _getLambdaFunc();
        CoordType xLo = getVarFunc(_cellIdx, Orient2DType::HORIZONTAL);
        CoordType yLo = getVarFunc(_cellIdx, Orient2DType::VERTICAL);
        CoordType xHi = xLo + _cellWidth;
        CoordType yHi = yLo + _cellHeight;
        // Smooth abs xLo xHi
        NumType obXLo = op::logSumExp0(
                op::conv<NumType>(_boundary->xLo() - xLo),
                alpha);
        NumType obXHi = op::logSumExp0(
                op::conv<NumType>(xHi - _boundary->xHi()),
                alpha);
        // y
        NumType obYLo = op::logSumExp0(
                op::conv<NumType>(_boundary->yLo() - yLo),
                alpha
                );
        NumType obYHi = op::logSumExp0(
                op::conv<NumType>(yHi - _boundary->yHi()),
                alpha
                );
        return (obXLo + obXHi + obYLo + obYHi) * lambda;
    }

    void accumlateGradient(
            const std::function<CoordType(IndexType cellIdx, Orient2DType orient)> &getVarFunc,
            const std::function<void(NumType, IndexType, Orient2DType)> &accumulateGradFunc) const
    {
        NumType alpha  = _getAlphaFunc();
        NumType lambda = _getLambdaFunc();
        CoordType xLo = getVarFunc(_cellIdx, Orient2DType::HORIZONTAL);
        CoordType yLo = getVarFunc(_cellIdx, Orient2DType::VERTICAL);
        CoordType xHi = xLo + _cellWidth;
        CoordType yHi = yLo + _cellHeight;
        // max(lower - x/yLo, 0), max (x/yHi - upper, 0)
        NumType gradObX =
            - op::gradLogSumExp0( // negative comes from the derivative
                    op::conv<NumType>(_boundary->xLo() - xLo),
                    alpha
                    );
        gradObX +=
            op::gradLogSumExp0(
                    op::conv<NumType>(xHi - _boundary->xHi()),
                    alpha
                    );
        accumulateGradFunc(gradObX * lambda, _cellIdx, Orient2DType::HORIZONTAL);
        // y
        NumType gradObY =
            - op::gradLogSumExp0( // negative comes from the derivative
                    op::conv<NumType>(_boundary->yLo() - yLo),
                    alpha
                    );
        gradObY +=
            op::gradLogSumExp0( 
                    op::conv<NumType>(yHi - _boundary->yHi()),
                    alpha
                    );
        accumulateGradFunc(gradObY * lambda, _cellIdx, Orient2DType::VERTICAL);
    }


    IndexType _cellIdx;
    CoordType _cellWidth;
    CoordType _cellHeight;
    Box<CoordType> *_boundary = nullptr;
    std::function<NumType(void)> _getAlphaFunc;
    std::function<NumType(void)> _getLambdaFunc; ///< A function to get the current lambda multiplier
};


template <typename NumType, typename CoordType>
struct is_placement_differentiable_concept<CellOutOfBoundaryPenaltyDifferentiable<NumType, CoordType>>
{
    typedef std::true_type  is_placement_differentiable_concept_type;
};

/// @brief Asymmetry penalty
template<typename NumType, typename CoordType>
struct AsymmetryDifferentiable
{
    typedef NumType numerical_type;
    typedef CoordType coordinate_type;

    AsymmetryDifferentiable(IndexType symGrpIdx, const std::function<NumType(void)> &getLambdaFunc)
    {
        _symGrpIdx = symGrpIdx;
        _getLambdaFunc = getLambdaFunc;
    }
    /// @brief add a symmetric pair. require the cell widths are the same
    void addSymPair(IndexType cellIdxI, IndexType cellIdxJ, CoordType width)
    {
        _pairCells.emplace_back(std::array<IndexType, 2>({cellIdxI, cellIdxJ}));
        _pairWidths.emplace_back(op::conv<NumType>(width));
    }
    void addSelfSym(IndexType cellIdx, CoordType width)
    {
        _selfSymCells.emplace_back(cellIdx);
        _selfSymWidths.emplace_back(op::conv<NumType>(width));
    }

    NumType evaluate(const std::function<CoordType(IndexType, Orient2DType)> &getVarFunc) const
    {
        NumType lambda = _getLambdaFunc();
        NumType asym = 0;
        NumType symAxis = op::conv<NumType>(getVarFunc(_symGrpIdx, Orient2DType::NONE));
        for (IndexType symPairIdx = 0; symPairIdx < _pairCells.size(); ++symPairIdx)
        {
            IndexType cellI = _pairCells[symPairIdx][0];
            IndexType cellJ = _pairCells[symPairIdx][1];
            NumType xi = op::conv<NumType>(getVarFunc(cellI, Orient2DType::HORIZONTAL));
            NumType yi = op::conv<NumType>(getVarFunc(cellI, Orient2DType::VERTICAL));
            NumType w = _pairWidths[symPairIdx];
            NumType xj = op::conv<NumType>(getVarFunc(cellJ, Orient2DType::HORIZONTAL));
            NumType yj = op::conv<NumType>(getVarFunc(cellJ, Orient2DType::VERTICAL));

            asym += pow(yi - yj, 2.0);
            asym += pow(xi + xj + w - 2 *symAxis, 2.0);
        }
        for (IndexType ssIdx = 0; ssIdx < _selfSymCells.size(); ++ssIdx)
        {
            NumType x = op::conv<NumType>(getVarFunc(_selfSymCells[ssIdx], Orient2DType::HORIZONTAL));
            NumType w = _selfSymWidths[ssIdx];

            asym += pow( x + w / 2 - symAxis, 2.0);
        }
        return asym  * lambda;
    }
    void accumlateGradient(
            const std::function<CoordType(IndexType cellIdx, Orient2DType orient)> &getVarFunc,
            const std::function<void(NumType, IndexType, Orient2DType)> &accumulateGradFunc) const
    {
        NumType lambda = _getLambdaFunc();
        NumType symAxis = op::conv<NumType>(getVarFunc(_symGrpIdx, Orient2DType::NONE));
        for (IndexType symPairIdx = 0; symPairIdx < _pairCells.size(); ++symPairIdx)
        {
            IndexType cellI = _pairCells[symPairIdx][0];
            IndexType cellJ = _pairCells[symPairIdx][1];
            NumType xi = op::conv<NumType>(getVarFunc(cellI, Orient2DType::HORIZONTAL));
            NumType yi = op::conv<NumType>(getVarFunc(cellI, Orient2DType::VERTICAL));
            NumType w = _pairWidths[symPairIdx];
            NumType xj = op::conv<NumType>(getVarFunc(cellJ, Orient2DType::HORIZONTAL));
            NumType yj = op::conv<NumType>(getVarFunc(cellJ, Orient2DType::VERTICAL));

            NumType partialX = 2.0 * (xi + xj + w - 2 * symAxis) * lambda;
            accumulateGradFunc(partialX, cellI, Orient2DType::HORIZONTAL);
            accumulateGradFunc(partialX, cellJ, Orient2DType::HORIZONTAL);
            accumulateGradFunc( -2 *partialX, _symGrpIdx, Orient2DType::NONE);

            NumType partialYI = 2.0 * (yi - yj) * lambda;
            accumulateGradFunc(partialYI, cellI, Orient2DType::VERTICAL);
            accumulateGradFunc(- partialYI, cellJ, Orient2DType::VERTICAL);
        }
        for (IndexType ssIdx = 0; ssIdx < _selfSymCells.size(); ++ssIdx)
        {
            NumType x = op::conv<NumType>(getVarFunc(_selfSymCells[ssIdx], Orient2DType::HORIZONTAL));
            NumType w = _selfSymWidths[ssIdx];

            NumType partial = 2.0 * (x + w /2 - symAxis) * lambda;

            accumulateGradFunc(partial, _selfSymCells[ssIdx], Orient2DType::HORIZONTAL);
            accumulateGradFunc(partial, _symGrpIdx, Orient2DType::NONE);

        }
    }

    IndexType _symGrpIdx;
    std::vector<std::array<IndexType, 2>> _pairCells;
    std::vector<NumType> _pairWidths;
    std::vector<IndexType> _selfSymCells;
    std::vector<NumType> _selfSymWidths;
    std::function<NumType(void)> _getLambdaFunc;
};


template <typename NumType, typename CoordType>
struct is_placement_differentiable_concept<AsymmetryDifferentiable<NumType, CoordType>>
{
    typedef std::true_type  is_placement_differentiable_concept_type;
};


/// @brief Each individual operator includes three cells and hence the two vectors they compose
template<typename NumType, typename CoordType>
struct CosineDatapathDifferentiable
{
    typedef NumType numerical_type;
    typedef CoordType coordinate_type;

    CosineDatapathDifferentiable(
            IndexType sCellIdx, const XY<CoordType> &sOffset,
            IndexType midCellIdx, const XY<CoordType> &midOffsetA, const XY<CoordType> &midOffsetB,
            IndexType tCellIdx, const XY<CoordType> &tOffset,
            const std::function<NumType(void)> &getLambdaFunc)
        : _sCellIdx(sCellIdx),
          _midCellIdx(midCellIdx),
          _tCellIdx(tCellIdx), 
          _getLambdaFunc(getLambdaFunc)
        {
            _sOffset.setX(op::conv<NumType>(sOffset.x()));
            _sOffset.setY(op::conv<NumType>(sOffset.y()));
            _midOffsetA.setX(op::conv<NumType>(midOffsetA.x()));
            _midOffsetA.setY(op::conv<NumType>(midOffsetA.y()));
            _midOffsetB.setX(op::conv<NumType>(midOffsetB.x()));
            _midOffsetB.setY(op::conv<NumType>(midOffsetB.y()));
            _tOffset.setX(tOffset.x());
            _tOffset.setY(tOffset.y());
        }
    

    NumType evaluate(const std::function<CoordType(IndexType, Orient2DType)> &getVarFunc) const;
    void accumlateGradient(
            const std::function<CoordType(IndexType cellIdx, Orient2DType orient)> &getVarFunc,
            const std::function<void(NumType, IndexType, Orient2DType)> &accumulateGradFunc) const;

    IndexType _sCellIdx = INDEX_TYPE_MAX; ///< Source
    XY<NumType> _sOffset; ///< The offset for x0
    IndexType _midCellIdx = INDEX_TYPE_MAX; ///< Middle
    XY<NumType> _midOffsetA;
    XY<NumType> _midOffsetB;
    IndexType _tCellIdx = INDEX_TYPE_MAX; ///< Target
    XY<NumType> _tOffset;
    std::function<NumType(void)> _getLambdaFunc; ///< A function to get the current lambda multiplier
};

template<typename NumType, typename CoordType>
inline NumType CosineDatapathDifferentiable<NumType, CoordType>::evaluate(
        const std::function<CoordType(IndexType, Orient2DType)> &getVarFunc) const
{
    NumType lambda = _getLambdaFunc();
    NumType x1 = op::conv<NumType>(getVarFunc(_sCellIdx, Orient2DType::HORIZONTAL));
    NumType y1 = op::conv<NumType>(getVarFunc(_sCellIdx, Orient2DType::VERTICAL));
    NumType x2 = op::conv<NumType>(getVarFunc(_midCellIdx, Orient2DType::HORIZONTAL));
    NumType y2 = op::conv<NumType>(getVarFunc(_midCellIdx, Orient2DType::VERTICAL));
    NumType x3 = op::conv<NumType>(getVarFunc(_tCellIdx, Orient2DType::HORIZONTAL));
    NumType y3 = op::conv<NumType>(getVarFunc(_tCellIdx, Orient2DType::VERTICAL));
    NumType ox1 = _sOffset.x();
    NumType oy1 = _sOffset.y();
    NumType ox2a = _midOffsetA.x();
    NumType oy2a = _midOffsetA.y();
    NumType ox2b = _midOffsetB.x();
    NumType oy2b = _midOffsetB.y();
    NumType ox3 = _tOffset.x();
    NumType oy3 = _tOffset.y();

    return (((ox1 - ox2a + x1 - x2)*(ox3 - ox2b - x2 + x3) + (oy1 - oy2a + y1 - y2)*(oy3 - oy2b - y2 + y3))/(sqrt(pow(ox1 - ox2a + x1 - x2, 2.0) + pow(oy1 - oy2a + y1 - y2, 2.0))*(pow(ox3 - ox2b - x2 + x3, 2.0) + (oy3 - oy2b - y2 + y3)^2)^(1/2)) + 1) * lambda;
}

template<typename NumType, typename CoordType>
inline void CosineDatapathDifferentiable<NumType, CoordType>::accumlateGradient(
        const std::function<CoordType(IndexType cellIdx, Orient2DType orient)> &getVarFunc,
        const std::function<void(NumType, IndexType, Orient2DType)> &accumulateGradFunc) const
{
    NumType lambda = _getLambdaFunc();
    NumType x1 = op::conv<NumType>(getVarFunc(_sCellIdx, Orient2DType::HORIZONTAL));
    NumType y1 = op::conv<NumType>(getVarFunc(_sCellIdx, Orient2DType::VERTICAL));
    NumType x2 = op::conv<NumType>(getVarFunc(_midCellIdx, Orient2DType::HORIZONTAL));
    NumType y2 = op::conv<NumType>(getVarFunc(_midCellIdx, Orient2DType::VERTICAL));
    NumType x3 = op::conv<NumType>(getVarFunc(_tCellIdx, Orient2DType::HORIZONTAL));
    NumType y3 = op::conv<NumType>(getVarFunc(_tCellIdx, Orient2DType::VERTICAL));
    NumType ox1 = _sOffset.x();
    NumType oy1 = _sOffset.y();
    NumType ox2a = _midOffsetA.x();
    NumType oy2a = _midOffsetA.y();
    NumType ox2b = _midOffsetB.x();
    NumType oy2b = _midOffsetB.y();
    NumType ox3 = _tOffset.x();
    NumType oy3 = _tOffset.y();

    NumType dx1 = (ox3 - ox2b - x2 + x3)/(((ox1 - ox2a + x1 - x2)^2 + (oy1 - oy2a + y1 - y2)^2)^(1/2)*((ox3 - ox2b - x2 + x3)^2 + (oy3 - oy2b - y2 + y3)^2)^(1/2)) - (((ox1 - ox2a + x1 - x2)*(ox3 - ox2b - x2 + x3) + (oy1 - oy2a + y1 - y2)*(oy3 - oy2b - y2 + y3))*(2*ox1 - 2*ox2a + 2*x1 - 2*x2))/(2*((ox1 - ox2a + x1 - x2)^2 + (oy1 - oy2a + y1 - y2)^2)^(3/2)*((ox3 - ox2b - x2 + x3)^2 + (oy3 - oy2b - y2 + y3)^2)^(1/2));

    dx1 *= lambda;

    NumType dy1 = 
(oy3 - oy2b - y2 + y3)/(((ox1 - ox2a + x1 - x2)^2 + (oy1 - oy2a + y1 - y2)^2)^(1/2)*((ox3 - ox2b - x2 + x3)^2 + (oy3 - oy2b - y2 + y3)^2)^(1/2)) - (((ox1 - ox2a + x1 - x2)*(ox3 - ox2b - x2 + x3) + (oy1 - oy2a + y1 - y2)*(oy3 - oy2b - y2 + y3))*(2*oy1 - 2*oy2a + 2*y1 - 2*y2))/(2*((ox1 - ox2a + x1 - x2)^2 + (oy1 - oy2a + y1 - y2)^2)^(3/2)*((ox3 - ox2b - x2 + x3)^2 + (oy3 - oy2b - y2 + y3)^2)^(1/2))
;
    dy1 *= lambda;

    NumType dx2 = 
(((ox1 - ox2a + x1 - x2)*(ox3 - ox2b - x2 + x3) + (oy1 - oy2a + y1 - y2)*(oy3 - oy2b - y2 + y3))*(2*ox1 - 2*ox2a + 2*x1 - 2*x2))/(2*((ox1 - ox2a + x1 - x2)^2 + (oy1 - oy2a + y1 - y2)^2)^(3/2)*((ox3 - ox2b - x2 + x3)^2 + (oy3 - oy2b - y2 + y3)^2)^(1/2)) - (ox1 + ox3 - ox2a - ox2b + x1 - 2*x2 + x3)/(((ox1 - ox2a + x1 - x2)^2 + (oy1 - oy2a + y1 - y2)^2)^(1/2)*((ox3 - ox2b - x2 + x3)^2 + (oy3 - oy2b - y2 + y3)^2)^(1/2)) + (((ox1 - ox2a + x1 - x2)*(ox3 - ox2b - x2 + x3) + (oy1 - oy2a + y1 - y2)*(oy3 - oy2b - y2 + y3))*(2*ox3 - 2*ox2b - 2*x2 + 2*x3))/(2*((ox1 - ox2a + x1 - x2)^2 + (oy1 - oy2a + y1 - y2)^2)^(1/2)*((ox3 - ox2b - x2 + x3)^2 + (oy3 - oy2b - y2 + y3)^2)^(3/2))
;
    dx2 *= lambda;

    NumType dy2 = 
(((ox1 - ox2a + x1 - x2)*(ox3 - ox2b - x2 + x3) + (oy1 - oy2a + y1 - y2)*(oy3 - oy2b - y2 + y3))*(2*oy1 - 2*oy2a + 2*y1 - 2*y2))/(2*((ox1 - ox2a + x1 - x2)^2 + (oy1 - oy2a + y1 - y2)^2)^(3/2)*((ox3 - ox2b - x2 + x3)^2 + (oy3 - oy2b - y2 + y3)^2)^(1/2)) - (oy1 + oy3 - oy2a - oy2b + y1 - 2*y2 + y3)/(((ox1 - ox2a + x1 - x2)^2 + (oy1 - oy2a + y1 - y2)^2)^(1/2)*((ox3 - ox2b - x2 + x3)^2 + (oy3 - oy2b - y2 + y3)^2)^(1/2)) + (((ox1 - ox2a + x1 - x2)*(ox3 - ox2b - x2 + x3) + (oy1 - oy2a + y1 - y2)*(oy3 - oy2b - y2 + y3))*(2*oy3 - 2*oy2b - 2*y2 + 2*y3))/(2*((ox1 - ox2a + x1 - x2)^2 + (oy1 - oy2a + y1 - y2)^2)^(1/2)*((ox3 - ox2b - x2 + x3)^2 + (oy3 - oy2b - y2 + y3)^2)^(3/2));
    dy2 *= lambda;

    NumType dx3 = 
(ox1 - ox2a + x1 - x2)/(((ox1 - ox2a + x1 - x2)^2 + (oy1 - oy2a + y1 - y2)^2)^(1/2)*((ox3 - ox2b - x2 + x3)^2 + (oy3 - oy2b - y2 + y3)^2)^(1/2)) - (((ox1 - ox2a + x1 - x2)*(ox3 - ox2b - x2 + x3) + (oy1 - oy2a + y1 - y2)*(oy3 - oy2b - y2 + y3))*(2*ox3 - 2*ox2b - 2*x2 + 2*x3))/(2*((ox1 - ox2a + x1 - x2)^2 + (oy1 - oy2a + y1 - y2)^2)^(1/2)*((ox3 - ox2b - x2 + x3)^2 + (oy3 - oy2b - y2 + y3)^2)^(3/2))
;
    dx3 *= lambda;

    NumType dy3 = 
(oy1 - oy2a + y1 - y2)/(((ox1 - ox2a + x1 - x2)^2 + (oy1 - oy2a + y1 - y2)^2)^(1/2)*((ox3 - ox2b - x2 + x3)^2 + (oy3 - oy2b - y2 + y3)^2)^(1/2)) - (((ox1 - ox2a + x1 - x2)*(ox3 - ox2b - x2 + x3) + (oy1 - oy2a + y1 - y2)*(oy3 - oy2b - y2 + y3))*(2*oy3 - 2*oy2b - 2*y2 + 2*y3))/(2*((ox1 - ox2a + x1 - x2)^2 + (oy1 - oy2a + y1 - y2)^2)^(1/2)*((ox3 - ox2b - x2 + x3)^2 + (oy3 - oy2b - y2 + y3)^2)^(3/2))
;
    dy3 *= lambda;

    accumulateGradFunc(dx1, _sCellIdx, Orient2DType::HORIZONTAL);
    accumulateGradFunc(dy1, _sCellIdx, Orient2DType::VERTICAL);
    accumulateGradFunc(dx2, _midCellIdx, Orient2DType::HORIZONTAL);
    accumulateGradFunc(dy2, _midCellIdx, Orient2DType::VERTICAL);
    accumulateGradFunc(dx3, _tCellIdx, Orient2DType::HORIZONTAL);
    accumulateGradFunc(dy3, _tCellIdx, Orient2DType::VERTICAL);
}

} //namespace diff

PROJECT_NAMESPACE_END

#endif // IDEAPLACE_DIFFERENT_H_
