/**
 * @file NlpGPlacer.h
 * @brief The global placement solver with non-linear optimization
 * @author Keren Zhu
 * @date 03/29/2020
 */

#ifndef IDEAPLACE_NLPGPLACER_H_
#define IDEAPLACE_NLPGPLACER_H_

#include <Eigen/Dense>
#include "db/Database.h"
#include "place/different.h"

PROJECT_NAMESPACE_BEGIN

class NlpGPlacerBase;

namespace nlp{

}// namespace nlp

/// @brief non-linear programming-based analog global placement
class NlpGPlacerBase
{
    protected:
        using EigenMatrix = Eigen::Matrix<RealType, Eigen::Dynamic, Eigen::Dynamic>;
        using EigenVector = Eigen::Matrix<RealType, Eigen::Dynamic, 1>;
        using EigenXY = Eigen::Matrix<RealType, Eigen::Dynamic, 2, Eigen::ColMajor>;
        typedef RealType nlp_coordinate_type;
        typedef RealType nlp_numerical_type;
        typedef diff::LseHpwlDifferentiable<nlp_numerical_type, nlp_coordinate_type> nlp_hpwl_type;
        typedef diff::CellPairOverlapPenaltyDifferentiable<nlp_numerical_type, nlp_coordinate_type> nlp_ovl_type;
        typedef diff::CellOutOfBoundaryPenaltyDifferentiable<nlp_numerical_type, nlp_coordinate_type> nlp_oob_type;
        typedef diff::AsymmetryDifferentiable<nlp_numerical_type, nlp_coordinate_type> nlp_asym_type;
        typedef diff::CosineDatapathDifferentiable<nlp_numerical_type, nlp_coordinate_type> nlp_cos_type;
    protected:
        template<typename task_type>
        class Task
        {
            public: 
                explicit Task(task_type &task) : _task(task) {}
                explicit Task(task_type &&task) : _task(task) {}
                void run() { task_type::run(_task); }
                const task_type &taskData() const { return _task; } 
            private:
                task_type _task;
        };

        /// @brief Evaluating objective tasks
        class EvaObjTask
        {
            public:
                EvaObjTask(const std::function<nlp_numerical_type(void)> &func) : _evaFunc(func) {}
                EvaObjTask(const EvaObjTask & other) { _evaFunc = other._evaFunc; }
                static void run(EvaObjTask & task) { task._obj =  task._evaFunc(); }
                nlp_numerical_type obj() const { return _obj; }
            private:
                std::function<nlp_numerical_type(void)> _evaFunc; 
                nlp_numerical_type _obj;
        };
    public:
        explicit NlpGPlacerBase(Database &db) : _db(db) {}
        IntType solve();

    protected:
        /* Init functions */
        void initProblem();
        void initHyperParams();
        void initBoundaryParams();
        void initVariables();
        void initRandomPlacement();
        void initOperators();
        /* Output functions */
        void writeOut();
        /* construct tasks */
        void constructObjectiveCalculationTasks();
    protected:
        Database &_db; ///< The placement engine database
        /* NLP problem parameters */
        RealType _alpha; ///< Used in LSE approximation hyperparameter
        Box<RealType> _boundary; ///< The boundary constraint for the placement
        RealType _scale = 0.01; /// The scale ratio between float optimization kernel coordinate and placement database coordinate unit
        RealType _totalCellArea = 0; ///< The total cell area of the problem
        RealType _overlapThreshold = NLP_WN_CONJ_OVERLAP_THRESHOLD; ///< Threshold for whether increase penalty for overlapping penalty
        RealType _oobThreshold = NLP_WN_CONJ_OOB_THRESHOLD; ///< The threshold for wehther increasing the penalty for out of boundry
        RealType _asymThreshold = NLP_WN_CONJ_ASYM_THRESHOLD; ///< The threshold for whether increasing the penalty for asymmetry
        RealType _defaultSymAxis = 0.0; ///< The default symmetric axis
        /* Optimization internal results */
        RealType _curOvlRatio = 1.0; ///< The current overlapping ratio
        RealType _curOOBRatio = 1.0; ///< The current out of boundry ratio
        RealType _curAsymDist = 1.0; ///< The current asymmetric distance
        /* Optimization data */
        EigenXY _pl; ///< The placement solutions
        EigenVector _sym; ///< The symmetry axis variables
        /* Tasks */
        std::vector<Task<EvaObjTask>> _evaHpwlTasks; ///< The tasks for evaluating hpwl objectives
        std::vector<Task<EvaObjTask>> _evaOvlTasks; ///< The tasks for evaluating overlap objectives
        std::vector<Task<EvaObjTask>> _evaOobTasks; ///< The tasks for evaluating out of boundary objectives
        std::vector<Task<EvaObjTask>> _evaAsymTasks;  ///< The tasks for evaluating asymmetry objectives
        std::vector<Task<EvaObjTask>> _evaCosTasks;  ///< The tasks for evaluating signal path objectives
        /* Operators */
        std::vector<nlp_hpwl_type> _hpwlOps; ///< The HPWL cost 
        std::vector<nlp_ovl_type> _ovlOps; ///< The cell pair overlapping penalty operators
        std::vector<nlp_oob_type> _oobOps; ///< The cell out of boundary penalty operators 
        std::vector<nlp_asym_type> _asymOps; ///< The asymmetric penalty operators
        std::vector<nlp_cos_type> _cosOps;

};

PROJECT_NAMESPACE_END
#endif //IDEAPLACE_NLPGPLACER_H_
