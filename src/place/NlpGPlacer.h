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

namespace nlp {

    template<typename T>
    struct stop_condition_trait {};

    /// @brief stop condition with number of iterations
    struct stop_after_num_outer_iterations
    {
        IntType maxIter = 20;
        IntType curIter = 0;
    };
    
    template<>
    struct stop_condition_trait<stop_after_num_outer_iterations>
    {
        template<typename NlpType>
        static stop_after_num_outer_iterations construct(NlpType &) { return stop_after_num_outer_iterations(); }
        template<typename NlpType>
        static IntType stopPlaceCondition(stop_after_num_outer_iterations &stop, NlpType &)
        {
            if (stop.curIter >= stop.maxIter)
            {
                stop.curIter = 0;
                return 1;
            }
            ++stop.curIter;
            return 0;
        }
    };

    struct nlp_parameters
    {
        typedef stop_after_num_outer_iterations stop_condition_type;
    };

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

        /* parameters */
        typedef nlp::nlp_parameters::stop_condition_type stop_condition_type;
        typedef nlp::stop_condition_trait<stop_condition_type> stop_condition_trait;
        friend stop_condition_trait;
    protected:
        template<typename task_type>
        class Task
        {
            public: 
                explicit Task() {}
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

        /// @brief The tasks for only wrapping a function
        class FuncTask
        {
            public:
                FuncTask() { _func = [](){}; }
                FuncTask(const std::function<void(void)> &func) : _func(func) {}
                FuncTask(const FuncTask &other) { _func = other._func; }
                static void run(FuncTask &task) { task._func(); }
            private:
                std::function<void(void)> _func;
        };

        /// @brief The tasks for only wrapping a function returning integer
        class ConditionTask
        {
            public:
                ConditionTask() { _func = [](){ return 0; }; _cond = 0; }
                ConditionTask(const std::function<IntType(void)> &func) : _func(func) { _cond = 0;}
                ConditionTask(const ConditionTask &other) { _func = other._func; _cond = other._cond; }
                static void run(ConditionTask &task) { task._cond = task._func(); }
                IntType cond() const { return _cond; }
            private:
                std::function<IntType(void)> _func;
                IntType _cond = 0;
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
        void initOptimizationKernelMembers();
        /* Output functions */
        void writeOut();
        /* construct tasks */
        void constructTasks();
        // Obj-related
        void constructObjTasks();
        void constructObjectiveCalculationTasks();
        void constructSumObjTasks();
        void constructWrapObjTask();
        // Optimization kernel-related
        void constructOptimizationKernelTasks();
        void constructStopConditionTask();
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
        //RealType _curOvlRatio = 1.0; ///< The current overlapping ratio
        //RealType _curOOBRatio = 1.0; ///< The current out of boundry ratio
        //RealType _curAsymDist = 1.0; ///< The current asymmetric distance
        RealType _objHpwl = 0.0; ///< The current value for hpwl
        RealType _objOvl = 0.0; ///< The current value for overlapping penalty
        RealType _objOob = 0.0; ///< The current value for out of boundary penalty
        RealType _objAsym = 0.0; ///< The current value for asymmetry penalty
        RealType _objCos = 0.0; ///< The current value for the cosine signal path penalty
        RealType _obj = 0.0; ///< The current value for the total objective penalty
        /* NLP optimization kernel memebers */
        stop_condition_type _stopCondition;
        /* Optimization data */
        EigenXY _pl; ///< The placement solutions
        EigenVector _sym; ///< The symmetry axis variables
        /* Tasks */
        // Evaluating objectives
        std::vector<Task<EvaObjTask>> _evaHpwlTasks; ///< The tasks for evaluating hpwl objectives
        std::vector<Task<EvaObjTask>> _evaOvlTasks; ///< The tasks for evaluating overlap objectives
        std::vector<Task<EvaObjTask>> _evaOobTasks; ///< The tasks for evaluating out of boundary objectives
        std::vector<Task<EvaObjTask>> _evaAsymTasks;  ///< The tasks for evaluating asymmetry objectives
        std::vector<Task<EvaObjTask>> _evaCosTasks;  ///< The tasks for evaluating signal path objectives
        // Sum the objectives
        Task<FuncTask> _sumObjHpwlTask; ///< The task for summing hpwl objective
        Task<FuncTask> _sumObjOvlTask; ///< The task for summing the overlapping objective
        Task<FuncTask> _sumObjOobTask; ///< The task for summing the out of boundary objective
        Task<FuncTask> _sumObjAsymTask; ///< The task for summing the asymmetry objective
        Task<FuncTask> _sumObjCosTask; ///< The task for summing the cosine signal path objective
        Task<FuncTask> _sumObjAllTask; ///< The task for summing the different objectives together
        // Optimization kernel
        Task<ConditionTask> _checkStopConditionTask; ///< The task to check whether the optimization should stop
#ifdef DEBUG_SINGLE_THREAD_GP
        // Wrapper tasks for debugging
        Task<FuncTask> _wrapObjHpwlTask; ///< The task for wrap the objective 
        Task<FuncTask> _wrapObjOvlTask;
        Task<FuncTask> _wrapObjOobTask;
        Task<FuncTask> _wrapObjAsymTask;
        Task<FuncTask> _wrapObjCosTask;
        Task<FuncTask> _wrapObjAllTask;
#endif //DEBUG_SINGLE_THREAD_GP
        /* Operators */
        std::vector<nlp_hpwl_type> _hpwlOps; ///< The HPWL cost 
        std::vector<nlp_ovl_type> _ovlOps; ///< The cell pair overlapping penalty operators
        std::vector<nlp_oob_type> _oobOps; ///< The cell out of boundary penalty operators 
        std::vector<nlp_asym_type> _asymOps; ///< The asymmetric penalty operators
        std::vector<nlp_cos_type> _cosOps;
};

/// @brief first-order optimization
class NlpGPlacerFirstOrder : public NlpGPlacerBase
{
    public:
        NlpGPlacerFirstOrder(Database &db) : NlpGPlacerBase(db) {}
    protected:

        /// @brief The tasks for calculating the partials from a operator
        template<typename op_type>
        class CalculateOperatorPartialTask {};

        // @brief The tasks for transfer from CalculateOperatorPartialTask to a target matrix
        template<typename op_type>
        class UpdateGradientFromPartialTask {};

        template<>
        class CalculateOperatorPartialTask<nlp_hpwl_type>
        {
            friend UpdateGradientFromPartialTask<nlp_hpwl_type>;
            static constexpr IntType MAX_NUM_CELLS = IDEAPLACE_DEFAULT_MAX_NUM_CELLS;
            public:
                CalculateOperatorPartialTask() {  }
                CalculateOperatorPartialTask(const CalculateOperatorPartialTask<nlp_hpwl_type> &other) { _partials = other._partials; _op = other._op; }
                CalculateOperatorPartialTask(nlp_hpwl_type *op) 
                { 
                    _op = op; 
                    IndexType numPartials = op->_cells.size() * 2; // dx and dy for each cells
                    _partials.resize(numPartials, 2);
                    _inverseCellMap.resize(_op->_cells.size());
                    for (IndexType idx = 0; idx < numPartials; ++idx)
                    {
                        _partials(idx, 0) = 0.0;
                        _partials(idx, 1) = 0.0;
                    }
                    for (IndexType idx = 0; idx < op->_cells.size(); ++idx)
                    {
                        _cellMap[op->_cells[idx]] = idx;
                        _inverseCellMap[idx] = op->_cells[idx];
                    }
                    auto accumulateFunc = [&](nlp_numerical_type num, IndexType cellIdx, Orient2DType orient)
                    {
                        if (orient == Orient2DType::HORIZONTAL)
                        {
                            _partials(_cellMap[cellIdx], 0) += num;
                        }
                        else if (orient == Orient2DType::VERTICAL)
                        {
                            _partials(_cellMap[cellIdx], 1) += num;
                        }
                    };
                    op->setAccumulateGradFunc(accumulateFunc);
                }
                IndexType numCells() const { return _inverseCellMap.size(); }
                static void run(CalculateOperatorPartialTask &task) { diff::placement_differentiable_traits<nlp_hpwl_type>::accumlateGradient(*task._op); }
            private:
                EigenXY _partials;
                nlp_hpwl_type* _op = nullptr;
                std::array<IndexType, MAX_NUM_CELLS> _cellMap; ///< From db cell index to this class index
                std::vector<IndexType> _inverseCellMap;
        };

        template<>
        class UpdateGradientFromPartialTask<nlp_hpwl_type>
        {
            public:
                UpdateGradientFromPartialTask() {  }
                UpdateGradientFromPartialTask(const UpdateGradientFromPartialTask<nlp_hpwl_type> &other) { _target = other._target; _calcTask = other._calcTask; }
                UpdateGradientFromPartialTask(CalculateOperatorPartialTask<nlp_hpwl_type> *calcTask) { _calcTask = calcTask; }
                static void run(UpdateGradientFromPartialTask &task)
                {
                    for (IndexType idx = 0; idx < task._calcTask->numCells(); ++idx)
                    {
                        IndexType cellIdx = task._calcTask->_inverseCellMap[idx];
                        (*task._target)(cellIdx, 0) = task._calcTask->_partials(cellIdx, 0);
                    }
                }
            private:
                EigenXY *_target;
                CalculateOperatorPartialTask<nlp_hpwl_type> *_calcTask;
        };
    protected:
        EigenXY _grad; ///< The first order graident
};

PROJECT_NAMESPACE_END
#endif //IDEAPLACE_NLPGPLACER_H_
