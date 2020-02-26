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

namespace op
{
    /// @brief abs smooth function with respect to 0
    /// @return the smoothed result
    template<typename NumType>
    constexpr NumType logSumExp0(NumType var, NumType alpha)
    {
        return alpha * log(exp(var / alpha) + 1);
    }

    /// @brief The gradient of some smooth function from Biying I cannot understand
    /// @return the gradient of logSumExp0
    template<typename NumType>
    constexpr NumType gradLogSumExp0(NumType var, NumType alpha)
    {
        return exp( var / alpha ) / ( exp(var / alpha) + 1);
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


/// @brief the cell out of boundary penalty
template<typename NumType, typename CoordType>
struct CellOutOfBoundaryPenalty
{
    typedef NumType numerical_type;
    typedef CoordType coordinate_type;

    CellOutOfBoundaryPenalty(IndexType cellIdx_, CoordType cellWidth_, CoordType cellHeight_, Box<CoordType> *boundary_, NumType *alpha_, NumType *lambda_)
    {
        cellIdx = cellIdx_;
        cellWidth = cellWidth_;
        cellHeight = cellHeight_;
        boundary = boundary_;
        alpha = alpha_;
        lambda = lambda_;
    }

    NumType evaluate( std::function<CoordType(IndexType, Orient2DType)> getVarFunc) const
    {
        CoordType xLo = getVarFunc(cellIdx, Orient2DType::HORIZONTAL);
        CoordType yLo = getVarFunc(cellIdx, Orient2DType::VERTICAL);
        CoordType xHi = xLo + cellWidth;
        CoordType yHi = yLo + cellHeight;
        // Smooth abs xLo xHi
        NumType obXLo = op::logSumExp0(
                op::conv<CoordType, NumType>(boundary->xLo() - xLo),
                *alpha);
        NumType obXHi = op::logSumExp0(
                op::conv<CoordType, NumType>(xHi - boundary->xHi()),
                *alpha);
        // y
        NumType obYLo = op::logSumExp0(
                op::conv<CoordType, NumType>(boundary->yLo() - yLo),
                *alpha
                );
        NumType obYHi = op::logSumExp0(
                op::conv<CoordType, NumType>(yHi - boundary->yHi()),
                *alpha
                );
        return (obXLo + obXHi + obYLo + obYHi) * (*lambda);
    }

    void accumlateGradient(std::function<CoordType(IndexType cellIdx, Orient2DType orient)> getVarFunc,
            std::function<void(NumType, IndexType, Orient2DType)> accumulateGradFunc) const
    {
        CoordType xLo = getVarFunc(cellIdx, Orient2DType::HORIZONTAL);
        CoordType yLo = getVarFunc(cellIdx, Orient2DType::VERTICAL);
        CoordType xHi = xLo + cellWidth;
        CoordType yHi = yLo + cellHeight;
        // max(lower - x/yLo, 0), max (x/yHi - upper, 0)
        NumType gradObX =
            - op::gradLogSumExp0( // negative comes from the derivative
                    op::conv<CoordType, NumType>(boundary->xLo() - xLo),
                    *alpha
                    );
        gradObX +=
            op::gradLogSumExp0(
                    op::conv<CoordType, NumType>(xHi - boundary->xHi()),
                    *alpha
                    );
        accumulateGradFunc(gradObX * (*lambda), cellIdx, Orient2DType::HORIZONTAL);
        // y
        NumType gradObY =
            - op::gradLogSumExp0( // negative comes from the derivative
                    op::conv<CoordType, NumType>(boundary->yLo() - yLo),
                    *alpha
                    );
        gradObY +=
            op::gradLogSumExp0( 
                    op::conv<CoordType, NumType>(yHi - boundary->yHi()),
                    *alpha
                    );
        accumulateGradFunc(gradObY * (*lambda), cellIdx, Orient2DType::VERTICAL);
    }


    IndexType cellIdx;
    CoordType cellWidth;
    CoordType cellHeight;
    Box<CoordType> *boundary = nullptr;
    NumType *alpha = nullptr;
    NumType *lambda = nullptr;
};

PROJECT_NAMESPACE_END

#endif // IDEAPLACE_DIFFERENT_H_
