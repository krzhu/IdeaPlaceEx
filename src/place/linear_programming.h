/**
 * @file linear_programming.h
 * @brief The linear programming concept
 * @author Keren Zhu
 * @date 03/10/2020
 */

#ifndef IDEAPLACE_LINEAR_PROGRAMMING_H_
#define IDEAPLACE_LINEAR_PROGRAMMING_H_

#include "global/global.h"

PROJECT_NAMESPACE_BEGIN

template<typename solver_type>
struct linear_programming_trait
{
    typedef typename solver_type::variable_type variable_type; ///< The LP variable type
    typedef typename solver_type::value_type value_type; ///< The value type of LP problem. Usually variety of floats
    typedef typename solver_type::expr_type expr_type; ///< The type for expressions
    typedef typename solver_type::constr_type constr_type; ///< The type for constraints
    typedef typename solver_type::status_type status_type; ///< The LP solving status

    static variable_type addVar(solver_type &solver)
    {
        return solver.addVar();
    }
    static void addConstr(solver_type &solver, const constr_type &constr)
    {
        solver.addConstr(constr);
    }
    static void setVarLowerBound(solver_type &solver, const variable_type &var, const value_type &val)
    {
        solver.setVarLowerBound(var, val);
    }
    static void setVarUpperBound(solver_type &solver, const variable_type &var, const value_type &val)
    {
        solver.setVarUpperBound(var, val);
    }
    static void setVarInteger(solver_type &solver, const variable_type &var)
    {
        solver.setVarInteger(var);
    }
    static void setVarContinuous(solver_type &solver, const variable_type &var)
    {
        solver.setVarContinuous(var);
    }
    static void setObjectiveMaximize(solver_type &solver)
    {
        solver.setObjectiveMaximize();
    }
    static void setObjectiveMinimize(solver_type &solver)
    {
        solver.setObjectiveMinimize();
    }
    static void setObjective(solver_type &solver, const expr_type &expr)
    {
        solver.setObjective(expr);
    }
    static void solve(solver_type &solver)
    {
        solver.solve();
    }
    static status_type status(solver_type &solver)
    {
        return solver.status();
    }
    static bool isOptimal(solver_type &solver)
    {
        return status(solver).isOptimal();
    }
    static bool isSuboptimal(solver_type &solver)
    {
        return status(solver).isSuboptimal();
    }
    static bool isUnbounded(solver_type &solver)
    {
        return status(solver).isUnbounded();
    }
    static bool isInfeasible(solver_type &solver)
    {
        return status(solver).isInfeasible();
    }
    static value_type evaluateExpr(const solver_type &solver, const expr_type &expr)
    {
        return solver.evaluate(expr);
    }
    static value_type solution(const solver_type &solver, const variable_type &var)
    {
        return solver.solution(var);
    }
    static std::string statusStr(const solver_type &solver)
    {
        return status(solver).toStr();
    }
};

PROJECT_NAMESPACE_END

#endif // IDEAPLACE_LINEAR_PROGRAMMING_H_
