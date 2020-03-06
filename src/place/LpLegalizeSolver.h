/**
 * @file LpLegalizeSolver.h
 * @brief The Lp solver for legalization
 * @author Keren Zhu
 * @date 03/06/2020
 */

#ifndef IDEAPLACE_LP_LEGALIZESOLVER_H_
#define IDEAPLACE_LP_LEGALIZESOLVER_H_

#include "db/Database.h"
#include <limbo/solvers/api/LPSolveApi.h>

PROJECT_NAMESPACE_BEGIN

class Constraints;

template<typename SolverType>
struct legalization_lp_traits
{
    static bool solve(SolverType &solver)
    {
        return solver.solve();
    }
    static void enableOptHpwl(SolverType &solver)
    {
        solver.enableOptHpwl();
    }
    static void enableOptArea(SolverType &solver)
    {
        solver.enableOptArea();
    }
    static void setHor(SolverType &solver)
    {
        solver.enbleHor();
    }
    static void setVer(SolverType &solver)
    {
        solver.setVer();
    }
    static void setMaxWidth(SolverType &solver, RealType wStar)
    {
        solver.setMaxWidth(wStar);
    }
    static void exportSolutions(SolverType &solver)
    {
        solver.exportSolutions();
    }
    static RealType evaluateObjective(SolverType &solver)
    {
        return solver.evaluateObjective();
    }

};
/// @brief The LP solver for legalization
class LpLegalizeSolver
{
    friend legalization_lp_traits<LpLegalizeSolver>;
    public:
        typedef limbo::solvers::LinearModel<RealType, RealType> LpModelType;
        typedef limbo::solvers::LPSolveLinearApi
            <LpModelType::coefficient_value_type, 
            LpModelType::variable_value_type>
                SolverType;
        explicit LpLegalizeSolver(Database &db, Constraints &constraints, bool isHor=true,
                IntType optHpwl=0, IntType optArea=1)
            : _db(db), _constrains(constraints), _isHor(isHor), _optHpwl(optHpwl), _optArea(optArea)
        {} //_solver = SolverType(&_ilpModel); }
        /// @brief solve the problem
        bool solve();
        // @brief dump out the solutions to the database
        void exportSolution();
        /// @brief evaluate the objective function and return the value
        RealType evaluateObj();
        /// @brief set the maximum width or height (_wStar)
        /// @param the maximum width or height in the hpwl optimization problem
        void setWStar(RealType wStar) { _wStar = wStar; }
    private:
        /// @brief solve the LP
        bool solveLp();
        /* Varibles functions */
        /// @brief add ILP variables
        void addIlpVars();
        /// @brief calculate the number of variables32
        IndexType numVars() const;
        /// @brief add location variables
        void addLocVars();
        /// @brief add wirelegth variables
        void addWirelengthVars();
        /// @brief add area variables
        void addAreaVars();
        /// @brief add sym group varibales
        void addSymVars();
        /* Obj functions */
        /// @brief set the objective function
        void configureObjFunc();
        /// @brief add wire length objective
        void addWirelengthObj();
        /// @brief add area objective
        void addAreaObj();
        /// @brief add relaxed symmetric penalty
        void addSymObj();
        /* Constraint functions */
        /// @brief add constraints
        void addIlpConstraints();
        /// @brief add area constraints
        void addBoundaryConstraints();
        /// @brief add topology constraints
        void addTopologyConstraints();
        /// @brief add symmetry constraints
        void addSymmetryConstraints();
        void addSymmetryConstraintsWithEqu();
        void addSymmetryConstraintsRex();
        /// @brief add hpwl constraints
        void addHpwlConstraints();
    private:
        /* Configurations - Inputs */
        Database &_db; ///< The database for the Ideaplace
        Constraints &_constrains; ///< The constraints edges to be honored
        bool _isHor = true; ///< Whether solving horizontal or vertical
        IntType _optHpwl = 0; ///< Whether optimizing HPWL in ILP problems
        IntType _optArea = 1; ///< Whether optimizing area in ILP problems
        /* Optimization supporting variables */
        LpModelType _ilpModel; ///< The ILP model
        LpModelType::expression_type _obj; ///< The objective function of the ILP model
        std::vector<LpModelType::variable_type> _locs; ///< The location variables of the ILP model
        std::vector<LpModelType::variable_type> _wlL; ///< The left wirelength variables of the ILP model
        std::vector<LpModelType::variable_type> _wlR; ///< The right wirelength variables of the ILP model
        LpModelType::variable_type _dim; ///< The variable for area optimization
        RealType _wStar = 0; ///< The optimal W found in legalization step
        std::vector<LpModelType::variable_type> _symLocs; ///< The variable for symmetric group axises
#ifdef MULTI_SYM_GROUP
        bool _isMultipleSymGrp = true;
#else
        bool _isMultipleSymGrp = false;
#endif 
        bool _relaxEqualityConstraint = false;
        //SolverType _solver; ///< Solver
        /*  Optimization Results */
        limbo::solvers::SolverProperty _optimStatus; ///< The resulting status
        limbo::solvers::LPSolveParameters _params;
        RealType _largeNum = 100000.0; ///< A large number

};


template<>
struct legalization_lp_traits<LpLegalizeSolver>
{
    static bool solve(LpLegalizeSolver &solver)
    {
        return solver.solve();
    }
    static void enableOptHpwl(LpLegalizeSolver &solver)
    {
        solver._optHpwl = 1;
        solver._optArea = 0;
    }
    static void enableOptArea(LpLegalizeSolver &solver)
    {
        solver._optHpwl = 0;
        solver._optArea = 1;
    }
    static void setHor(LpLegalizeSolver &solver)
    {
        solver._isHor = true;
    }
    static void setVer(LpLegalizeSolver &solver)
    {
        solver._isHor = false;
    }
    static void setMaxWidth(LpLegalizeSolver &solver, RealType wStar)
    {
        solver._wStar = wStar;
    }
    static void exportSolutions(LpLegalizeSolver &solver)
    {
        solver.exportSolution();
    }
    static RealType evaluateObjective(LpLegalizeSolver &solver)
    {
        return solver.evaluateObj();
    }

};

PROJECT_NAMESPACE_END

#endif //IDEAPLACE_LP_LEGALIZESOLVER_H_
