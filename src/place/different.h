/**
 * @file different.h
 * @brief The numerical differentiable concepts and implementations
 * @author Keren Zhu
 * @date 02/25/2020
 */

#ifndef IDEAPLACE_DIFFERENT_H_
#define IDEAPLACE_DIFFERENT_H_

#include "db/Database.h"
#include <math.h>

PROJECT_NAMESPACE_BEGIN

struct placement_differentiable_concept {};

template <typename ConceptType>
struct is_placement_differentiable_concept
{
    typedef NoType type;
};

template <>
struct is_placement_differentiable_concept<placement_differentiable_concept>
{
    typedef YesType type;
};

template <typename DifType>
struct placement_differentiable_traits
{
    typedef DifType different_type;
    typedef typename different_type::numerical_type numerical_type;
    typedef typename different_type::coordinate_type coordinate_type;

    static numerical_type evaluate( const different_type & dif, std::function<coordinate_type(IndexType, Orient2DType)> getVarFunc) 
    {
        return dif.evaluate(getVarFunc);
    }

    static void accumlateGradient(const different_type & dif, 
            std::function<coordinate_type(IndexType cellIdx, Orient2DType orient)> getVarFunc,
            std::function<void(numerical_type, IndexType, Orient2DType)> accumulateGradFunc)
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

    template<typename RhsType, typename LhsType>
    constexpr LhsType conv(RhsType rhs)
    {
        if (std::is_same<RhsType, LhsType>::value)
        {
            return rhs;
        }
        else
        {
            return static_cast<LhsType>(rhs);
        }
    }
};

// @brief pair-wise cell overlapping penalty
template<typename NumType, typename CoordType>
struct CellPairOverlapPenalty
{
    typedef NumType numerical_type;
    typedef CoordType coordinate_type;

    CellPairOverlapPenalty(IndexType cellIdxI, CoordType cellWidthI, CoordType cellHeightI,
                           IndexType cellIdxJ, CoordType cellWidthJ, CoordType cellHeightJ,
                           NumType *alpha, NumType *lambda)
    {
        _cellIdxI = cellIdxI;
        _cellWidthI = cellWidthI;
        _cellHeightI = cellHeightI;
        _cellIdxJ = cellIdxJ;
        _cellWidthJ = cellWidthJ;
        _cellHeightJ = cellHeightJ;
        _alpha = alpha;
        _lambda = lambda;
    }

    NumType evaluate( std::function<CoordType(IndexType, Orient2DType)> getVarFunc) const
    {
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
                op::conv<CoordType, NumType>(var1X),
                op::conv<CoordType, NumType>(var2X),
                - (*_alpha)
                );
            overlapX = op::logSumExp0(
                    overlapX,
                    (*_alpha)
                    );
            // y
            CoordType var1Y = yHiI - yLoJ;
            CoordType var2Y = yHiJ - yLoI;
            NumType overlapY = op::logSumExp(
                op::conv<CoordType, NumType>(var1Y),
                op::conv<CoordType, NumType>(var2Y),
                - (*_alpha)
                );
            overlapY = op::logSumExp0(
                    overlapY,
                    (*_alpha)
                    );
            return (*_lambda) * overlapX * overlapY;
    }

    void accumlateGradient(std::function<CoordType(IndexType cellIdx, Orient2DType orient)> getVarFunc,
            std::function<void(NumType, IndexType, Orient2DType)> accumulateGradFunc) const
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
        NumType xi = op::conv<CoordType, NumType>(getVarFunc(_cellIdxI, Orient2DType::HORIZONTAL));
        NumType yi = op::conv<CoordType, NumType>(getVarFunc(_cellIdxI, Orient2DType::VERTICAL));
        NumType xj = op::conv<CoordType, NumType>(getVarFunc(_cellIdxJ, Orient2DType::HORIZONTAL));
        NumType yj = op::conv<CoordType, NumType>(getVarFunc(_cellIdxJ, Orient2DType::VERTICAL));
        NumType wi = op::conv<CoordType, NumType>(_cellWidthI);
        NumType hi = op::conv<CoordType, NumType>(_cellHeightI);
        NumType wj = op::conv<CoordType, NumType>(_cellWidthJ);
        NumType hj = op::conv<CoordType, NumType>(_cellHeightJ);
        NumType alpha = (*_alpha);

        NumType dxi = (alpha * alpha *log(1/(exp(-(hi + yi - yj)/alpha) + exp(-(hj - yi + yj)/alpha)) + 1)*(exp(-(wi + xi - xj)/alpha)/alpha - exp(-(wj - xi + xj)/alpha)/alpha))/((1/(exp(-(wi + xi - xj)/alpha) + exp(-(wj - xi + xj)/alpha)) + 1)*pow(exp(-(wi + xi - xj)/alpha) + exp(-(wj - xi + xj)/alpha), 2))
            ;
         NumType dxj = -dxi;

        NumType dyi = (alpha * alpha *log(1/(exp(-(wi + xi - xj)/alpha) + exp(-(wj - xi + xj)/alpha)) + 1)*(exp(-(hi + yi - yj)/alpha)/alpha - exp(-(hj - yi + yj)/alpha)/alpha))/((1/(exp(-(hi + yi - yj)/alpha) + exp(-(hj - yi + yj)/alpha)) + 1)*(pow(exp(-(hi + yi - yj)/alpha) + exp(-(hj - yi + yj)/alpha), 2)));
        NumType dyj = -dyi;

        // accumulate the computed partials
        accumulateGradFunc(dxi * (*_lambda), _cellIdxI, Orient2DType::HORIZONTAL);
        accumulateGradFunc(dxj * (*_lambda), _cellIdxJ, Orient2DType::HORIZONTAL);
        accumulateGradFunc(dyi * (*_lambda), _cellIdxI, Orient2DType::VERTICAL);
        accumulateGradFunc(dyj * (*_lambda), _cellIdxJ, Orient2DType::VERTICAL);
    }

    IndexType _cellIdxI;
    CoordType _cellWidthI;
    CoordType _cellHeightI;
    IndexType _cellIdxJ;
    CoordType _cellWidthJ;
    CoordType _cellHeightJ;
    NumType *_alpha = nullptr;
    NumType *_lambda = nullptr;
};

/// @brief the cell out of boundary penalty
template<typename NumType, typename CoordType>
struct CellOutOfBoundaryPenalty
{
    typedef NumType numerical_type;
    typedef CoordType coordinate_type;

    CellOutOfBoundaryPenalty(IndexType cellIdx, CoordType cellWidth, CoordType cellHeight, Box<CoordType> *boundary, NumType *alpha, NumType *lambda)
    {
        _cellIdx = cellIdx;
        _cellWidth = cellWidth;
        _cellHeight = cellHeight;
        _boundary = boundary;
        _alpha = alpha;
        _lambda = lambda;
    }

    NumType evaluate( std::function<CoordType(IndexType, Orient2DType)> getVarFunc) const
    {
        CoordType xLo = getVarFunc(_cellIdx, Orient2DType::HORIZONTAL);
        CoordType yLo = getVarFunc(_cellIdx, Orient2DType::VERTICAL);
        CoordType xHi = xLo + _cellWidth;
        CoordType yHi = yLo + _cellHeight;
        // Smooth abs xLo xHi
        NumType obXLo = op::logSumExp0(
                op::conv<CoordType, NumType>(_boundary->xLo() - xLo),
                *_alpha);
        NumType obXHi = op::logSumExp0(
                op::conv<CoordType, NumType>(xHi - _boundary->xHi()),
                *_alpha);
        // y
        NumType obYLo = op::logSumExp0(
                op::conv<CoordType, NumType>(_boundary->yLo() - yLo),
                *_alpha
                );
        NumType obYHi = op::logSumExp0(
                op::conv<CoordType, NumType>(yHi - _boundary->yHi()),
                *_alpha
                );
        return (obXLo + obXHi + obYLo + obYHi) * (*_lambda);
    }

    void accumlateGradient(std::function<CoordType(IndexType cellIdx, Orient2DType orient)> getVarFunc,
            std::function<void(NumType, IndexType, Orient2DType)> accumulateGradFunc) const
    {
        CoordType xLo = getVarFunc(_cellIdx, Orient2DType::HORIZONTAL);
        CoordType yLo = getVarFunc(_cellIdx, Orient2DType::VERTICAL);
        CoordType xHi = xLo + _cellWidth;
        CoordType yHi = yLo + _cellHeight;
        // max(lower - x/yLo, 0), max (x/yHi - upper, 0)
        NumType gradObX =
            - op::gradLogSumExp0( // negative comes from the derivative
                    op::conv<CoordType, NumType>(_boundary->xLo() - xLo),
                    *_alpha
                    );
        gradObX +=
            op::gradLogSumExp0(
                    op::conv<CoordType, NumType>(xHi - _boundary->xHi()),
                    *_alpha
                    );
        accumulateGradFunc(gradObX * (*_lambda), _cellIdx, Orient2DType::HORIZONTAL);
        // y
        NumType gradObY =
            - op::gradLogSumExp0( // negative comes from the derivative
                    op::conv<CoordType, NumType>(_boundary->yLo() - yLo),
                    *_alpha
                    );
        gradObY +=
            op::gradLogSumExp0( 
                    op::conv<CoordType, NumType>(yHi - _boundary->yHi()),
                    *_alpha
                    );
        accumulateGradFunc(gradObY * (*_lambda), _cellIdx, Orient2DType::VERTICAL);
    }


    IndexType _cellIdx;
    CoordType _cellWidth;
    CoordType _cellHeight;
    Box<CoordType> *_boundary = nullptr;
    NumType *_alpha = nullptr;
    NumType *_lambda = nullptr;
};

PROJECT_NAMESPACE_END

#endif // IDEAPLACE_DIFFERENT_H_
