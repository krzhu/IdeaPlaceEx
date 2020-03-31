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

namespace nlp {

    template<typename T>
    struct stop_condition_trait 
    {
        // static T construct(NlpType &)
        // static IntType stopPlaceCondition(T&, NlpType &)
    };

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

    struct nlp_types
    {
        typedef Eigen::Matrix<RealType, Eigen::Dynamic, Eigen::Dynamic> EigenMatrix;
        typedef Eigen::Matrix<RealType, Eigen::Dynamic, 1> EigenVector;
        typedef Eigen::Map<EigenVector> EigenMap;
        typedef RealType nlp_coordinate_type;
        typedef RealType nlp_numerical_type;
        typedef diff::LseHpwlDifferentiable<nlp_numerical_type, nlp_coordinate_type> nlp_hpwl_type;
        typedef diff::CellPairOverlapPenaltyDifferentiable<nlp_numerical_type, nlp_coordinate_type> nlp_ovl_type;
        typedef diff::CellOutOfBoundaryPenaltyDifferentiable<nlp_numerical_type, nlp_coordinate_type> nlp_oob_type;
        typedef diff::AsymmetryDifferentiable<nlp_numerical_type, nlp_coordinate_type> nlp_asym_type;
        typedef diff::CosineDatapathDifferentiable<nlp_numerical_type, nlp_coordinate_type> nlp_cos_type;

    };


}// namespace nlp

/// @brief namespace for nlp tasks
namespace nt
{
    template<typename task_type>
    class Task
    {
        public: 
            explicit Task() {}
            explicit Task(task_type &task) : _task(::klib::createSharedPtr(task)) { AssertMsg(false, "try point to a lvalue. Not tested\n"); }
            explicit Task(task_type &&task) : _task(std::make_shared<task_type>(std::move(task))) {}
            void run() { task_type::run(*_task); }
            const task_type &taskData() const { return *_task; } 
            task_type &taskData() { return *_task; }
            std::shared_ptr<task_type> taskDataPtr() { return _task; }
        private:
            std::shared_ptr<task_type> _task;
    };

    /// @brief Evaluating objective tasks
    class EvaObjTask
    {
        typedef nlp::nlp_types::nlp_numerical_type nlp_numerical_type;
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

    /// @brief The tasks for calculating the partials from a operator.
    /// @tparam the differentiable operator type
    template<typename op_type>
    class CalculateOperatorPartialTask;

    // @brief The tasks for transfer from CalculateOperatorPartialTask to a target matrix
    template<typename nlp_op_type>
    class UpdateGradientFromPartialTask
    {
        typedef typename nlp_op_type::numerical_type nlp_numerical_type;
        typedef typename nlp_op_type::coordinate_type nlp_coordiante_type;
        typedef nlp::nlp_types::EigenVector EigenVector;
        public:
            UpdateGradientFromPartialTask() = delete;
            UpdateGradientFromPartialTask(const UpdateGradientFromPartialTask<nlp_op_type> &other)  = delete;
            UpdateGradientFromPartialTask(UpdateGradientFromPartialTask<nlp_op_type> &other)  = delete;
            UpdateGradientFromPartialTask(UpdateGradientFromPartialTask<nlp_op_type> &&other)
                : _calcTask(std::move(other._calcTask)), _target(std::move(other._target)), _idxFunc(std::move(other._idxFunc))
            {
            }
            UpdateGradientFromPartialTask(std::shared_ptr<CalculateOperatorPartialTask<nlp_op_type>> calcTask, EigenVector *target,
                    const std::function<IndexType(IndexType, Orient2DType)> &idxFunc) 
            { 
                _calcTask = calcTask; 
                _target = target; 
                _idxFunc = idxFunc;  
            }
            static void run(UpdateGradientFromPartialTask &task)
            {
                for (IndexType idx = 0; idx < task._calcTask->numCells(); ++idx)
                {
                    IndexType cellIdx = task._calcTask->_inverseCellMap[idx];
                    (*task._target)(task._idxFunc(cellIdx, Orient2DType::HORIZONTAL)) = task._calcTask->_partialsX(idx);
                    (*task._target)(task._idxFunc(cellIdx, Orient2DType::VERTICAL)) = task._calcTask->_partialsY(idx);
                }
            }
        private:
            std::shared_ptr<CalculateOperatorPartialTask<nlp_op_type>> _calcTask;
            EigenVector *_target;
            std::function<IndexType(IndexType, Orient2DType)> _idxFunc; //< convert cell idx to eigen vector idx
    };

#ifdef MULTI_SYM_GROUP
    template<typename nlp_numerical_type, typename nlp_coordinate_type>
    class UpdateGradientFromPartialTask<diff::AsymmetryDifferentiable<nlp_numerical_type, nlp_coordinate_type>>
    {
        typedef diff::AsymmetryDifferentiable<nlp_numerical_type, nlp_coordinate_type> nlp_op_type;
        typedef nlp::nlp_types::EigenVector EigenVector;
        public:
            UpdateGradientFromPartialTask() {}
            UpdateGradientFromPartialTask(const UpdateGradientFromPartialTask<nlp_op_type> &other)  = delete;
            UpdateGradientFromPartialTask(UpdateGradientFromPartialTask<nlp_op_type> &other)  = delete;
            UpdateGradientFromPartialTask(UpdateGradientFromPartialTask<nlp_op_type> &&other)
                : _calcTask(std::move(other._calcTask)), _target(std::move(other._target)), _idxFunc(std::move(other._idxFunc))
            {
            }
            UpdateGradientFromPartialTask(std::shared_ptr<CalculateOperatorPartialTask<nlp_op_type>> calcTask, EigenVector *target,
                    const std::function<IndexType(IndexType, Orient2DType)> &idxFunc) 
            { 
                _calcTask = calcTask; 
                _target = target; 
                _idxFunc = idxFunc;  
            }
            static void run(UpdateGradientFromPartialTask &task)
            {
                for (IndexType idx = 0; idx < task._calcTask->numCells(); ++idx)
                {
                    IndexType cellIdx = task._calcTask->_inverseCellMap[idx];
                    (*task._target)(task._idxFunc(cellIdx, Orient2DType::HORIZONTAL)) = task._calcTask->_partialsX(idx);
                    (*task._target)(task._idxFunc(cellIdx, Orient2DType::VERTICAL)) = task._calcTask->_partialsY(idx);
                }
                (*task._target)(task._idxFunc(task._calcTask->_symGrpIdx, Orient2DType::HORIZONTAL)) = task._calcTask->_partialsX(idx);
            }
        private:
            std::shared_ptr<CalculateOperatorPartialTask<nlp_op_type>> _calcTask;
            EigenVector *_target;
            std::function<IndexType(IndexType, Orient2DType)> _idxFunc; //< convert cell idx to eigen vector idx
    };
#endif

    /// @brief trait template for building the _cellMap and _inverseCellMap from the different types opeartor. This template need partial specification.
    /// @tparam the differentiable operator type
    template<typename op_type>
    struct calc_operator_partial_build_cellmap_trait 
    {
        typedef op_type nlp_op_type;
        typedef CalculateOperatorPartialTask<nlp_op_type> calc_type;
        static void build(op_type *, calc_type *) {}
    };


    template<typename op_type>
    class CalculateOperatorPartialTask
    {
        typedef op_type nlp_op_type;
        typedef typename nlp_op_type::numerical_type nlp_numerical_type;
        typedef typename nlp_op_type::coordinate_type nlp_coordiante_type;
        typedef nlp::nlp_types::EigenVector EigenVector;
        friend UpdateGradientFromPartialTask<nlp_op_type>;
        friend calc_operator_partial_build_cellmap_trait<nlp_op_type>;
        static constexpr IntType MAX_NUM_CELLS = IDEAPLACE_DEFAULT_MAX_NUM_CELLS;
        public:
            CalculateOperatorPartialTask() = delete;
            CalculateOperatorPartialTask(CalculateOperatorPartialTask &other) = delete;
            CalculateOperatorPartialTask(CalculateOperatorPartialTask &&other)
                : _partialsX(std::move(other._partialsX)), _partialsY(std::move(other._partialsY)),
                _op(std::move(other._op)), _cellMap(std::move(other._cellMap)), _inverseCellMap(std::move(other._inverseCellMap)), 
                _numCells(std::move(other._numCells))
            {
                setAccumulateGradFunc();
            }
            CalculateOperatorPartialTask(nlp_op_type *op) 
            { 
                _op = op; 
                // Use this trait to speficify different number of cells for different operators
                calc_operator_partial_build_cellmap_trait<nlp_op_type>::build(*op, *this); 
                clear();
                setAccumulateGradFunc();
            }
            virtual void accumatePartial(nlp_numerical_type num, IndexType cellIdx, Orient2DType orient)
            {
                if (orient == Orient2DType::HORIZONTAL)
                {
                    _partialsX(_cellMap[cellIdx]) += num;
                }
                else if (orient == Orient2DType::VERTICAL)
                {
                    _partialsY(_cellMap[cellIdx]) += num;
                }
            }
            void setAccumulateGradFunc()
            {
                _op->setAccumulateGradFunc([&](nlp_numerical_type num, IndexType cellIdx, Orient2DType orient){ accumatePartial(num, cellIdx, orient);});
            }
            virtual void clear() 
            { 
                for (IndexType idx = 0; idx < numCells(); ++idx)
                {
                    _partialsX(idx) = 0.0;
                    _partialsY(idx) = 0.0;
                }
            }
            IndexType numCells() const { return _numCells; }
            static void run(CalculateOperatorPartialTask &task) 
            { 
                task.clear(); 
                diff::placement_differentiable_traits<nlp_op_type>::accumlateGradient(*(task._op)); 
            }
        protected:
            EigenVector _partialsX;
            EigenVector _partialsY;
            nlp_op_type* _op = nullptr;
            std::array<IndexType, MAX_NUM_CELLS> _cellMap; ///< From db cell index to this class index
            std::vector<IndexType> _inverseCellMap;
            IndexType _numCells;
    };

    template<typename nlp_numerical_type, typename nlp_coordinate_type>
    struct calc_operator_partial_build_cellmap_trait<diff::LseHpwlDifferentiable<nlp_numerical_type, nlp_coordinate_type>>
    {
        typedef diff::LseHpwlDifferentiable<nlp_numerical_type, nlp_coordinate_type> nlp_op_type;
        typedef CalculateOperatorPartialTask<nlp_op_type> calc_type;
        static void build(nlp_op_type &op, calc_type &calc)
        {
            calc._numCells = op._cells.size(); // dx and dy for each cells
            calc._partialsX.resize(calc._numCells);
            calc._partialsY.resize(calc._numCells);
            calc._inverseCellMap.resize(op._cells.size());
            for (IndexType idx = 0; idx < op._cells.size(); ++idx)
            {
                calc._cellMap[op._cells[idx]] = idx;
                calc._inverseCellMap[idx] = op._cells[idx];
            }
        }
    };

#ifdef MULTI_SYM_GROUP
    template<typename nlp_numerical_type, typename nlp_coordinate_type>
    class CalculateOperatorPartialTask<diff::AsymmetryDifferentiable<nlp_numerical_type, nlp_coordinate_type>> 
    : public CalculateOperatorPartialTask<diff::LseHpwlDifferentiable<nlp_numerical_type, nlp_coordinate_type>>
    {
        typedef diff::AsymmetryDifferentiable<nlp_numerical_type, nlp_coordinate_type> nlp_op_type;
        typedef nlp::nlp_types::EigenVector EigenVector;
        friend UpdateGradientFromPartialTask<nlp_op_type>;
        friend calc_operator_partial_build_cellmap_trait<nlp_op_type>;
        static constexpr IntType MAX_NUM_CELLS = IDEAPLACE_DEFAULT_MAX_NUM_CELLS;
        typedef CalculateOperatorPartialTask<diff::LseHpwlDifferentiable<nlp_numerical_type, nlp_coordinate_type>> base_type;
        public:
            virtual void accumatePartial(nlp_numerical_type num, IndexType cellIdx, Orient2DType orient)
            {
                if (orient == Orient2DType::HORIZONTAL)
                {
                    _partialsX(base_type::_cellMap[cellIdx]) += num;
                }
                else if (orient == Orient2DType::VERTICAL)
                {
                    _partialsY(base_type::_cellMap[cellIdx]) += num;
                }
                else
                {
                    _sym += num;
                }
            }
            void clear() override
            { 
                for (IndexType idx = 0; idx < base_type::numCells(); ++idx)
                {
                    base_type::_partialsX(idx) = 0.0;
                    base_type::_partialsY(idx) = 0.0;
                }
                _sym = 0.0;
            }
        protected:
            RealType _sym;
            IndexType _symGrpIdx;
    };
#endif
    template<typename nlp_numerical_type, typename nlp_coordinate_type>
    struct calc_operator_partial_build_cellmap_trait<diff::CellPairOverlapPenaltyDifferentiable<nlp_numerical_type, nlp_coordinate_type>>
    {
        typedef diff::CellPairOverlapPenaltyDifferentiable<nlp_numerical_type, nlp_coordinate_type> nlp_op_type;
        typedef CalculateOperatorPartialTask<nlp_op_type> calc_type;
        static void build(nlp_op_type &op, calc_type &calc)
        {
            calc._numCells = 2; // Always have exactly two cells
            calc._partialsX.resize(calc._numCells);
            calc._partialsY.resize(calc._numCells);
            calc._inverseCellMap.resize(2);
            calc._cellMap[op._cellIdxI] = 0;
            calc._inverseCellMap[0] = op._cellIdxI;
            calc._cellMap[op._cellIdxJ] = 1;
            calc._inverseCellMap[1] = op._cellIdxJ;
        }
    };

    template<typename nlp_numerical_type, typename nlp_coordinate_type>
    struct calc_operator_partial_build_cellmap_trait<diff::CellOutOfBoundaryPenaltyDifferentiable<nlp_numerical_type, nlp_coordinate_type>>
    {
        typedef diff::CellOutOfBoundaryPenaltyDifferentiable<nlp_numerical_type, nlp_coordinate_type> nlp_op_type;
        typedef CalculateOperatorPartialTask<nlp_op_type> calc_type;
        static void build(nlp_op_type &op, calc_type &calc)
        {
            calc._numCells = 1; // Always have exactly two cells
            calc._partialsX.resize(calc._numCells);
            calc._partialsY.resize(calc._numCells);
            calc._inverseCellMap.resize(1);
            calc._cellMap[op._cellIdx] = 0;
            calc._inverseCellMap[0] = op._cellIdx;
        }
    };


    template<typename nlp_numerical_type, typename nlp_coordinate_type>
    struct calc_operator_partial_build_cellmap_trait<diff::AsymmetryDifferentiable<nlp_numerical_type, nlp_coordinate_type>>
    {
        typedef diff::AsymmetryDifferentiable<nlp_numerical_type, nlp_coordinate_type> nlp_op_type;
        typedef CalculateOperatorPartialTask<nlp_op_type> calc_type;
        static void build(nlp_op_type &op, calc_type &calc)
        {
            calc._numCells = op._pairCells.size() * 2 + op._selfSymCells.size(); // Always have exactly two cells
            calc._partialsX.resize(calc._numCells);
            calc._partialsY.resize(calc._numCells);
            calc._inverseCellMap.resize(calc._numCells);
            IndexType idx = 0;
            for (IndexType pairIdx = 0; pairIdx < op._pairCells.size(); ++pairIdx)
            {
                calc._cellMap[op._pairCells[pairIdx][0]] = idx;
                calc._inverseCellMap[idx] = op._pairCells[pairIdx][0];
                ++idx;
                calc._cellMap[op._pairCells[pairIdx][1]] = idx;
                calc._inverseCellMap[idx] = op._pairCells[pairIdx][1];
                ++idx;
            }
            for (IndexType ssIdx = 0; ssIdx < op._selfSymCells.size(); ++ssIdx)
            {
                calc._cellMap[op._selfSymCells[ssIdx]] = idx;
                calc._inverseCellMap[idx] = op._selfSymCells[ssIdx];
                ++idx;
            }
#ifdef MULTI_SYM_GROUP
            calc._symGrpIdx = op._symGrpIdx;
#endif
        }
    };


}// namespace nt

/// @brief non-linear programming-based analog global placement
class NlpGPlacerBase
{
    public:
        typedef nlp::nlp_types::EigenMatrix EigenMatrix;
        typedef nlp::nlp_types::EigenVector EigenVector;
        typedef nlp::nlp_types::EigenMap EigenMap;
        typedef nlp::nlp_types::nlp_coordinate_type nlp_coordinate_type;
        typedef nlp::nlp_types::nlp_numerical_type nlp_numerical_type;
        typedef nlp::nlp_types::nlp_hpwl_type nlp_hpwl_type;
        typedef nlp::nlp_types::nlp_ovl_type nlp_ovl_type;
        typedef nlp::nlp_types::nlp_oob_type nlp_oob_type;
        typedef nlp::nlp_types::nlp_asym_type nlp_asym_type;
        typedef nlp::nlp_types::nlp_cos_type nlp_cos_type;

        /* parameters */
        typedef nlp::nlp_parameters::stop_condition_type stop_condition_type;
        typedef nlp::stop_condition_trait<stop_condition_type> stop_condition_trait;
        friend stop_condition_trait;
    
    public:
        explicit NlpGPlacerBase(Database &db) : _db(db) {}
        IntType solve();

    protected:
        /* Init functions */
        virtual void initProblem();
        void initHyperParams();
        void initBoundaryParams();
        void initVariables();
        void initRandomPlacement();
        void initOperators();
        void initOptimizationKernelMembers();
        /* Output functions */
        void writeOut();
        /* Util functions */
        IndexType plIdx(IndexType cellIdx, Orient2DType orient) 
        {
            if (orient == Orient2DType::HORIZONTAL)
            {
                return cellIdx;
            }
            else if (orient == Orient2DType::VERTICAL)
            {
                return cellIdx + _numCells;
            }
            else
            {
#ifdef MULTI_SYM_GROUP
                return cellIdx + 2 *  _numCells; // here cell index representing the idx of sym grp
#else
                return 2 * _numCells;
#endif
            }

        }
        /* construct tasks */
        virtual void constructTasks();
        // Obj-related
        void constructObjTasks();
        void constructObjectiveCalculationTasks();
        void constructSumObjTasks();
        void constructWrapObjTask();
        // Optimization kernel-related
        void constructOptimizationKernelTasks();
        void constructStopConditionTask();
        /* Optimization  kernel */
        virtual void optimize();
    protected:
        Database &_db; ///< The placement engine database
        /* NLP problem parameters */
        IndexType _numCells; ///< The number of cells
        RealType _alpha; ///< Used in LSE approximation hyperparameter
        Box<RealType> _boundary; ///< The boundary constraint for the placement
        RealType _scale = 0.01; /// The scale ratio between float optimization kernel coordinate and placement database coordinate unit
        RealType _totalCellArea = 0; ///< The total cell area of the problem
        RealType _overlapThreshold = NLP_WN_CONJ_OVERLAP_THRESHOLD; ///< Threshold for whether increase penalty for overlapping penalty
        RealType _oobThreshold = NLP_WN_CONJ_OOB_THRESHOLD; ///< The threshold for wehther increasing the penalty for out of boundry
        RealType _asymThreshold = NLP_WN_CONJ_ASYM_THRESHOLD; ///< The threshold for whether increasing the penalty for asymmetry
        RealType _defaultSymAxis = 0.0; ///< The default symmetric axis
        /* Optimization internal results */
        RealType _objHpwl = 0.0; ///< The current value for hpwl
        RealType _objOvl = 0.0; ///< The current value for overlapping penalty
        RealType _objOob = 0.0; ///< The current value for out of boundary penalty
        RealType _objAsym = 0.0; ///< The current value for asymmetry penalty
        RealType _objCos = 0.0; ///< The current value for the cosine signal path penalty
        RealType _obj = 0.0; ///< The current value for the total objective penalty
        /* NLP optimization kernel memebers */
        stop_condition_type _stopCondition;
        /* Optimization data */
        EigenVector _pl; ///< The placement solutions
        std::shared_ptr<EigenMap> _plx; ///< The placement solutions for x coodinates
        std::shared_ptr<EigenMap> _ply; ///< The placement solutions for y coordinates
        std::shared_ptr<EigenMap> _sym; ///< The symmetry axis variables
        /* Tasks */
        // Evaluating objectives
        std::vector<nt::Task<nt::EvaObjTask>> _evaHpwlTasks; ///< The tasks for evaluating hpwl objectives
        std::vector<nt::Task<nt::EvaObjTask>> _evaOvlTasks; ///< The tasks for evaluating overlap objectives
        std::vector<nt::Task<nt::EvaObjTask>> _evaOobTasks; ///< The tasks for evaluating out of boundary objectives
        std::vector<nt::Task<nt::EvaObjTask>> _evaAsymTasks;  ///< The tasks for evaluating asymmetry objectives
        std::vector<nt::Task<nt::EvaObjTask>> _evaCosTasks;  ///< The tasks for evaluating signal path objectives
        // Sum the objectives
        nt::Task<nt::FuncTask> _sumObjHpwlTask; ///< The task for summing hpwl objective
        nt::Task<nt::FuncTask> _sumObjOvlTask; ///< The task for summing the overlapping objective
        nt::Task<nt::FuncTask> _sumObjOobTask; ///< The task for summing the out of boundary objective
        nt::Task<nt::FuncTask> _sumObjAsymTask; ///< The task for summing the asymmetry objective
        nt::Task<nt::FuncTask> _sumObjCosTask; ///< The task for summing the cosine signal path objective
        nt::Task<nt::FuncTask> _sumObjAllTask; ///< The task for summing the different objectives together
        // Optimization kernel
        nt::Task<nt::ConditionTask> _checkStopConditionTask; ///< The task to check whether the optimization should stop
#ifdef DEBUG_SINGLE_THREAD_GP
        // Wrapper tasks for debugging
        nt::Task<nt::FuncTask> _wrapObjHpwlTask; ///< The task for wrap the objective 
        nt::Task<nt::FuncTask> _wrapObjOvlTask;
        nt::Task<nt::FuncTask> _wrapObjOobTask;
        nt::Task<nt::FuncTask> _wrapObjAsymTask;
        nt::Task<nt::FuncTask> _wrapObjCosTask;
        nt::Task<nt::FuncTask> _wrapObjAllTask;
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
        /* Init */
        virtual void initProblem() override;
        void initFirstOrderGrad();
        /* Construct tasks */
        void constructFirstOrderTasks();
        void constructCalcPartialsTasks();
        void constructUpdatePartialsTasks();
        /* optimization */
        virtual void optimize() override;

    protected:
        /* Optimization data */
        EigenVector _grad; ///< The first order graident
        EigenVector _gradHpwl; ///< The first order gradient of hpwl objective
        EigenVector _gradOvl; ///< The first order gradient  of overlapping objective
        EigenVector _gradOob; ///< The first order gradient of out of boundary objective
        EigenVector _gradAsym; ///< The first order gradient of asymmetry objective
        /* Tasks */
        virtual void constructTasks() override;
        // Calculate the partials
        std::vector<nt::Task<nt::CalculateOperatorPartialTask<nlp_hpwl_type>>> _calcHpwlPartialTasks;
        std::vector<nt::Task<nt::CalculateOperatorPartialTask<nlp_ovl_type>>> _calcOvlPartialTasks;
        std::vector<nt::Task<nt::CalculateOperatorPartialTask<nlp_oob_type>>> _calcOobPartialTasks;
        std::vector<nt::Task<nt::CalculateOperatorPartialTask<nlp_asym_type>>> _calcAsymPartialTasks;
        // Update the partials
        std::vector<nt::Task<nt::UpdateGradientFromPartialTask<nlp_hpwl_type>>> _updateHpwlPartialTasks;
        std::vector<nt::Task<nt::UpdateGradientFromPartialTask<nlp_ovl_type>>> _updateOvlPartialTasks;
        std::vector<nt::Task<nt::UpdateGradientFromPartialTask<nlp_oob_type>>> _updateOobPartialTasks;
        std::vector<nt::Task<nt::UpdateGradientFromPartialTask<nlp_asym_type>>> _updateAsymPartialTasks;
};

PROJECT_NAMESPACE_END
#endif //IDEAPLACE_NLPGPLACER_H_
