/**
 * @file lp_limbo_lpsolve.h
 * @brief The linear programming implementation using LIMBO::LPSolveApi
 * @author Keren Zhu
 * @date 03/10/2020
 */

#ifndef IDEAPLACE_LP_LIMBO_LPSOLVE_H_
#define IDEAPLACE_LP_LIMBO_LPSOLVE_H_

#include <limbo/solvers/api/LPSolveApi.h>
#include "place/linear_programming.h"

PROJECT_NAMESPACE_BEGIN

template<typename value_type=double>
class LimboLpsolve
{
    friend linear_programming_trait<LimboLpsolve<value_type>>;
    public:
        LimboLpsolve() {}
        ~LimboLpsolve() {}
        typedef typename limbo::solvers::LinearModel<value_type, value_type> model_type;
        typedef typename model_type::variable_type variable_type;
        typedef typename model_type::expression_type expr_type;
        typedef typename model_type::constraint_type constr_type;
        typedef limbo::solvers::SolverProperty status_type;
        typedef limbo::solvers::LPSolveLinearApi
            <typename model_type::coefficient_value_type, 
            typename model_type::variable_value_type>
                limbo_solver_type;
    private:
        model_type _model; ///< The LP problem model
        status_type _status; ///< The result status. e.g. OPTIMAL
        limbo::solvers::LPSolveParameters _params; ///< The parameters for LIMBO solver
};

template<typename value_type>
struct linear_programming_trait<LimboLpsolve<value_type>>
{
    typedef LimboLpsolve<value_type> solver_type;
    typedef typename solver_type::variable_type variable_type;
    typedef typename solver_type::expr_type expr_type;
    typedef typename solver_type::constr_type constr_type;
    typedef typename solver_type::status_type status_type;



    static variable_type addVar(solver_type &solver)
    {
        return solver._model.addVariable(0, std::numeric_limits<RealType>::max(),
                                                limbo::solvers::CONTINUOUS, 
                                                "x" + solver._model.numVariables());
    }
    static void addConstr(solver_type &solver, const constr_type &constr)
    {
        bool success = solver._model.addConstraint(constr, "CONSTR");
        AssertMsg(success, "Limbo lib LP solver add constraint failed\n");
    }
    static void setVarLowerBound(solver_type &solver, const variable_type &var, const value_type &val)
    {
        solver._model.updateVariableLowerBound(var, val);
    }
    static void setVarUpperBound(solver_type &solver, const variable_type &var, const value_type &val)
    {
        solver._model.updateVariableUpperBound(var, val);
    }
    static void setVarInteger(solver_type &solver, const variable_type &var)
    {
        solver._model.setVariableNumericType(var, limbo::solvers::INTEGER);
    }
    static void setVarContinuous(solver_type &solver, const variable_type &var)
    {
        solver._model.setVariableNumericType(var, limbo::solvers::CONTINUOUS);
    }
    static void setObjectiveMaximize(solver_type &solver)
    {
        solver._model.setOptimizeType(limbo::solvers::MAX);
    }
    static void setObjectiveMinimize(solver_type &solver)
    {
        solver._model.setOptimizeType(limbo::solvers::MIN);
    }
    static void setObjective(solver_type &solver, const expr_type &expr)
    {
        solver._model.setObjective(expr);
    }
    static void solve(solver_type &solver)
    {
        solver._params.setVerbose(2); // 2: SEVERE
        typename solver_type::limbo_solver_type sol(&solver._model);
        solver._status = sol(&solver._params);
    }
    static status_type status(solver_type &solver)
    {
        return solver._status;
    }
    static bool isOptimal(solver_type &solver)
    {
        return solver._status == limbo::solvers::OPTIMAL;
    }
    static bool isSuboptimal(solver_type &solver)
    {
        return solver._status == limbo::solvers::SUBOPTIMAL;
    }
    static bool isUnbounded(solver_type &solver)
    {
        return solver._status == limbo::solvers::UNBOUNDED;
    }
    static bool isInfeasible(solver_type &solver)
    {
        return solver._status == limbo::solvers::INFEASIBLE;
    }
    static value_type evaluateExpr(const solver_type &solver, const expr_type &expr)
    {
        return solver._model.evaluateExpression(expr, solver._model.variableSolutions());
    }
    static value_type solution(const solver_type &solver, const variable_type &var)
    {
        return solver._model.variableSolution(var);
    }
    static std::string statusStr(const solver_type &solver)
    {
        return limbo::solvers::toString(solver._status);
    }
};

PROJECT_NAMESPACE_END


#endif //IDEAPLACE_LP_LIMBO_LPSOLVE_H_
