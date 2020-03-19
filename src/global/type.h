/**
 * @file type.h
 * @brief Define some the types being used globally
 * @author Keren Zhu
 * @date 09/30/2019
 */

#ifndef IDEAPLACE_TYPE_H_
#define IDEAPLACE_TYPE_H_

#include <cstdint>
#include <string>
#include <sstream>
#include "namespace.h"
#include "util/lp_limbo.h"

PROJECT_NAMESPACE_BEGIN
// Built-in type aliases
using IndexType  = std::uint32_t;
using IntType    = std::int32_t;
using RealType   = double;
using Byte       = std::uint8_t;
using LocType    = std::int32_t; // Location/design unit // Location/design unit
 // Location/design unit
// Built-in type constants
constexpr IndexType INDEX_TYPE_MAX  = UINT32_MAX;
constexpr IntType INT_TYPE_MAX      = INT32_MAX;
constexpr IntType INT_TYPE_MIN      = INT32_MIN;
constexpr RealType REAL_TYPE_MAX    = 1e100;
constexpr RealType REAL_TYPE_MIN    = -1e100;
constexpr RealType REAL_TYPE_TOL    = 1e-6;
constexpr LocType LOC_TYPE_MAX      = INT32_MAX;
constexpr LocType LOC_TYPE_MIN      = INT32_MIN;

// Type aliases
//using CostTy     = double;
using CostType = RealType;
constexpr CostType COST_TYPE_INVALID = REAL_TYPE_MIN;
constexpr  CostType COST_TYPE_MAX = REAL_TYPE_MAX;


// Enums

///  @brief The type of symmetry
enum class SymType
{
    HORIZONTAL = 0,
    VERTICAL   = 1,
    NONE = 2
};

enum class Orient2DType
{
    HORIZONTAL = 0,
    VERTICAL   = 1,
    NONE = 2
};

enum class Direction2DType
{
    WEST = 0,
    EAST = 1,
    SOUTH = 2,
    NORTH  = 3,
    NONE = 4
};

struct YesType
{
    static const bool value = true;
};


struct NoType
{
    static const bool value = false;
};

// Select LP solver
namespace _lp {
    struct gurobi
    {
        static const int rank = 0;
#ifndef LP_NOT_USE_GUROBI
        static const bool if_enable = true;
        typedef lp::LimboLpGurobi LpModel;
        typedef lp::LimboLpGurobiTrait LpTrait;
#else
        static const bool if_enable = false;
#endif
    };
    struct lp_solve
    {
        static const int rank = 1;
#ifndef LP_NOT_USE_LPSOLVE
        static const bool if_enable = true;
        typedef lp::LimboLpsolve LpModel;
        typedef lp::LimboLpsolveTrait LpTrait;
#else
        static const bool if_enable = false;
#endif
    };

    struct undefined
    {
        static const bool if_enable = true;
        typedef _undefined LpModel;
        typedef linear_programming_trait<_undefined> LpTrait;
    };

    template<int rank>
    struct lp_rank
    {
        typedef undefined solver_type;
    };

    template<>
    struct lp_rank<0>
    {
        typedef gurobi solver_type;
    };

    template<>
    struct lp_rank<1>
    {
        typedef lp_solve solver_type;
    };

    template<int rank, bool>
    struct select_lp
    {
    };
    template<int rank>
    struct select_lp<rank, true>
    {
        typedef typename lp_rank<rank>::solver_type solver_type;
        typedef typename solver_type::LpModel LpModel;
        typedef typename solver_type::LpTrait LpTrait;
    };
    template<int rank>
    struct select_lp<rank, false>
    {
        typedef typename lp_rank<rank+1>::solver_type solver_type;
        typedef typename select_lp<rank+1, solver_type::if_enable>::LpModel LpModel;
        typedef typename select_lp<rank+1, solver_type::if_enable>::LpTrait LpTrait;
    };
}; //namspace _lp

namespace lp
{
    typedef typename _lp::select_lp<0, _lp::lp_rank<0>::solver_type::if_enable>::LpModel LpModel;
    typedef typename _lp::select_lp<0, _lp::lp_rank<0>::solver_type::if_enable>::LpTrait LpTrait;
};

PROJECT_NAMESPACE_END

#endif // AROUTER_TYPE_H_

