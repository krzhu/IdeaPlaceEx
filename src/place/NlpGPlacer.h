/**
 * @file NlpGPlacer.h
 * @brief The global placement solver with non-linear optimization
 * @author Keren Zhu
 * @date 03/29/2020
 */

#ifndef IDEAPLACE_NLPGPLACER_H_
#define IDEAPLACE_NLPGPLACER_H_

#include <Eigen/Dense>
#ifdef IDEAPLACE_TASKFLOR_FOR_GRAD_OBJ_
#include <taskflow/taskflow.hpp>
#endif // IDEAPLACE_TASKFLOR_FOR_GRAD_OBJ
#include "db/Database.h"
#include "place/different.h"
#include "place/nlp/nlpOuterOptm.hpp"
#include "place/nlp/nlpInitPlace.hpp"
#include "place/nlp/nlpTasks.hpp"
#include "place/nlp/nlpTypes.hpp"
#include "place/nlp/nlpOptmKernels.hpp"
#include "place/nlp/nlpFirstOrderKernel.hpp"
#include "place/nlp/conjugateGradientWnlib.hpp" // TODO: remove after no need

PROJECT_NAMESPACE_BEGIN

namespace nlp 
{
    /* The wrapper of settings */

    struct nlp_default_hyperparamters
    {
    };

    struct nlp_default_types
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
    
    struct nlp_default_zero_order_algorithms
    {
        typedef outer_stop_condition::stop_condition_list<
            outer_stop_condition::stop_after_violate_small,
            outer_stop_condition::stop_after_num_outer_iterations<200>
            > stop_condition_type;
        typedef init_place::init_random_placement_with_normal_distribution_near_center init_place_type;

        /* multipliers */
        typedef outer_multiplier::init::hard_code_init mult_init_type;
        typedef outer_multiplier::update::subgradient_normalized_by_init<nlp_default_types::nlp_numerical_type> mult_update_type;
        //typedef outer_multiplier::update::direct_subgradient mult_update_type;
        typedef outer_multiplier::mult_const_hpwl_cos_and_penalty_by_type<nlp_default_types::nlp_numerical_type, mult_init_type, mult_update_type> mult_type;
    };

    struct nlp_default_first_order_algorithms
    {
        typedef converge::converge_list<
                    converge::converge_grad_norm_by_init<nlp_default_types::nlp_numerical_type>,
                    converge::converge_criteria_max_iter<3000>
                        >
                converge_type;
        //typedef optm::first_order::naive_gradient_descent<converge_type> optm_type;
        typedef optm::first_order::adam<converge_type, nlp_default_types::nlp_numerical_type> optm_type;
        //typedef optm::first_order::conjugate_gradient_wnlib optm_type;
        
        /* multipliers */
        typedef outer_multiplier::init::init_by_matching_gradient_norm mult_init_type;
        //typedef outer_multiplier::update::subgradient_normalized_by_init<nlp_default_types::nlp_numerical_type> mult_update_type;
        typedef outer_multiplier::update::direct_subgradient mult_update_type;
        typedef outer_multiplier::mult_const_hpwl_cos_and_penalty_by_type<nlp_default_types::nlp_numerical_type, mult_init_type, mult_update_type> mult_type;
    };

    struct nlp_default_settings
    {
        typedef nlp_default_zero_order_algorithms nlp_zero_order_algorithms_type;
        typedef nlp_default_first_order_algorithms nlp_first_order_algorithms_type;
        typedef nlp_default_hyperparamters nlp_hyperparamters_type;
        typedef nlp_default_types nlp_types_type;
    };


}// namespace nlp

/// @brief non-linear programming-based analog global placement
template<typename nlp_settings>
class NlpGPlacerBase
{
    public:
        typedef typename nlp_settings::nlp_types_type nlp_types;
        typedef typename nlp_settings::nlp_zero_order_algorithms_type nlp_zero_order_algorithms;
        typedef typename nlp_settings::nlp_hyperparamters_type nlp_hyperparamters;

        typedef typename nlp_types::EigenMatrix EigenMatrix;
        typedef typename nlp_types::EigenVector EigenVector;
        typedef typename nlp_types::EigenMap EigenMap;
        typedef typename nlp_types::nlp_coordinate_type nlp_coordinate_type;
        typedef typename nlp_types::nlp_numerical_type nlp_numerical_type;
        typedef typename nlp_types::nlp_hpwl_type nlp_hpwl_type;
        typedef typename nlp_types::nlp_ovl_type nlp_ovl_type;
        typedef typename nlp_types::nlp_oob_type nlp_oob_type;
        typedef typename nlp_types::nlp_asym_type nlp_asym_type;
        typedef typename nlp_types::nlp_cos_type nlp_cos_type;


        /* algorithms */
        typedef typename nlp_zero_order_algorithms::stop_condition_type stop_condition_type;
        typedef nlp::outer_stop_condition::stop_condition_trait<stop_condition_type> stop_condition_trait;
        template<typename _T>
        friend struct nlp::outer_stop_condition::stop_condition_trait;
        typedef typename nlp_zero_order_algorithms::init_place_type init_placement_type;
        typedef nlp::init_place::init_place_trait<init_placement_type> init_place_trait;
        friend init_place_trait;
        typedef typename nlp_zero_order_algorithms::mult_init_type mult_init_type;
        typedef nlp::outer_multiplier::init::multiplier_init_trait<mult_init_type> mult_init_trait;
        friend mult_init_trait;
        typedef typename nlp_zero_order_algorithms::mult_update_type mult_update_type;
        typedef nlp::outer_multiplier::update::multiplier_update_trait<mult_update_type> mult_update_trait;
        friend mult_update_trait;
        typedef typename nlp_zero_order_algorithms::mult_type mult_type;
        typedef nlp::outer_multiplier::multiplier_trait<mult_type> mult_trait;
        friend mult_trait;
    
    public:
        explicit NlpGPlacerBase(Database &db) : _db(db) {}
        IntType solve();

    protected:
        /* calculating obj */
        void calcObj()
        {
            _wrapObjAllTask.run();
        }
        /* Init functions */
        virtual void initProblem();
        void initHyperParams();
        void initBoundaryParams();
        void initVariables();
        void initPlace();
        void initOperators();
        void initOptimizationKernelMembers();
        /* Output functions */
        void writeOut();
        /* Util functions */
        IndexType plIdx(IndexType cellIdx, Orient2DType orient);
        void alignToSym();
        /* construct tasks */
        virtual void constructTasks();
        // Obj-related
        void constructObjTasks();
        void constructObjectiveCalculationTasks();
        void constructSumObjTasks();
#ifdef DEBUG_SINGLE_THREAD_GP
        void constructWrapObjTask();
#endif
        /* Optimization  kernel */
        virtual void optimize();
        /* build the computational graph */
#ifdef IDEAPLACE_TASKFLOR_FOR_GRAD_OBJ_
        void regEvaHpwlTaskflow(tf::Taskflow & tfFlow);
        void regEvaOvlTaskflow(tf::Taskflow & tfFlow);
        void regEvaOobTaskflow(tf::Taskflow & tfFlow);
        void regEvaAsymTaskflow(tf::Taskflow & tfFlow);
        void regEvaCosTaskflow(tf::Taskflow & tfFlow);
        void regEvaAllObjTaskflow(tf::Taskflow & tfFlow);
#endif
        /* Debugging function */
#ifdef DEBUG_GR
#ifdef DEBUG_DRAW
        void drawCurrentLayout(const std::string &name);
#endif 
#endif
    protected:
        Database &_db; ///< The placement engine database
        /* NLP problem parameters */
        IndexType _numCells; ///< The number of cells
        RealType _alpha; ///< Used in LSE approximation hyperparameter
        Box<RealType> _boundary; ///< The boundary constraint for the placement
        nlp_coordinate_type _scale = 0.01; /// The scale ratio between float optimization kernel coordinate and placement database coordinate unit
        nlp_coordinate_type _totalCellArea = 0; ///< The total cell area of the problem
        nlp_coordinate_type _defaultSymAxis = 0.0; ///< The default symmetric axis
        IndexType _numVariables = 0; ///< The number of variables
        /* Optimization internal results */
        nlp_numerical_type _objHpwl = 0.0; ///< The current value for hpwl
        nlp_numerical_type _objOvl = 0.0; ///< The current value for overlapping penalty
        nlp_numerical_type _objOob = 0.0; ///< The current value for out of boundary penalty
        nlp_numerical_type _objAsym = 0.0; ///< The current value for asymmetry penalty
        nlp_numerical_type _objCos = 0.0; ///< The current value for the cosine signal path penalty
        nlp_numerical_type _obj = 0.0; ///< The current value for the total objective penalty
        /* NLP optimization kernel memebers */
        stop_condition_type _stopCondition;
        /* Optimization data */
        EigenVector _pl; ///< The placement solutions
        std::shared_ptr<EigenMap> _plx; ///< The placement solutions for x coodinates
        std::shared_ptr<EigenMap> _ply; ///< The placement solutions for y coordinates
        std::shared_ptr<EigenMap> _sym; ///< The symmetry axis variables
        /* Tasks */
        // Evaluating objectives
        std::vector<nt::Task<nt::EvaObjTask<nlp_numerical_type>>> _evaHpwlTasks; ///< The tasks for evaluating hpwl objectives
        std::vector<nt::Task<nt::EvaObjTask<nlp_numerical_type>>> _evaOvlTasks; ///< The tasks for evaluating overlap objectives
        std::vector<nt::Task<nt::EvaObjTask<nlp_numerical_type>>> _evaOobTasks; ///< The tasks for evaluating out of boundary objectives
        std::vector<nt::Task<nt::EvaObjTask<nlp_numerical_type>>> _evaAsymTasks;  ///< The tasks for evaluating asymmetry objectives
        std::vector<nt::Task<nt::EvaObjTask<nlp_numerical_type>>> _evaCosTasks;  ///< The tasks for evaluating signal path objectives
        // Sum the objectives
        nt::Task<nt::FuncTask> _sumObjHpwlTask; ///< The task for summing hpwl objective
        nt::Task<nt::FuncTask> _sumObjOvlTask; ///< The task for summing the overlapping objective
        nt::Task<nt::FuncTask> _sumObjOobTask; ///< The task for summing the out of boundary objective
        nt::Task<nt::FuncTask> _sumObjAsymTask; ///< The task for summing the asymmetry objective
        nt::Task<nt::FuncTask> _sumObjCosTask; ///< The task for summing the cosine signal path objective
        nt::Task<nt::FuncTask> _sumObjAllTask; ///< The task for summing the different objectives together
        // Wrapper tasks for debugging
        nt::Task<nt::FuncTask> _wrapObjHpwlTask; ///< The task for wrap the objective 
        nt::Task<nt::FuncTask> _wrapObjOvlTask;
        nt::Task<nt::FuncTask> _wrapObjOobTask;
        nt::Task<nt::FuncTask> _wrapObjAsymTask;
        nt::Task<nt::FuncTask> _wrapObjCosTask;
        nt::Task<nt::FuncTask> _wrapObjAllTask;
        /* Operators */
        std::vector<nlp_hpwl_type> _hpwlOps; ///< The HPWL cost 
        std::vector<nlp_ovl_type> _ovlOps; ///< The cell pair overlapping penalty operators
        std::vector<nlp_oob_type> _oobOps; ///< The cell out of boundary penalty operators 
        std::vector<nlp_asym_type> _asymOps; ///< The asymmetric penalty operators
        std::vector<nlp_cos_type> _cosOps;
#ifdef IDEAPLACE_TASKFLOR_FOR_GRAD_OBJ_
        /* taskflow */
        tf::Taskflow _taskflow; ///< The taskflow of cpp-taskflow
#endif
};

template<typename nlp_settings>
inline IndexType NlpGPlacerBase<nlp_settings>::plIdx(IndexType cellIdx, Orient2DType orient)
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

/// @brief first-order optimization
template<typename nlp_settings>
class NlpGPlacerFirstOrder : public NlpGPlacerBase<nlp_settings>
{
    public:
        typedef NlpGPlacerBase<nlp_settings> base_type;
        typedef typename base_type::EigenVector EigenVector;
        typedef typename base_type::nlp_hpwl_type nlp_hpwl_type;
        typedef typename base_type::nlp_ovl_type nlp_ovl_type;
        typedef typename base_type::nlp_oob_type nlp_oob_type;
        typedef typename base_type::nlp_asym_type nlp_asym_type;
        typedef typename base_type::nlp_cos_type nlp_cos_type;
        typedef typename nlp_settings::nlp_first_order_algorithms_type nlp_first_order_algorithms;
        typedef typename nlp_first_order_algorithms::converge_type converge_type;
        typedef typename nlp::converge::converge_criteria_trait<converge_type> converge_trait;
        typedef typename nlp_first_order_algorithms::optm_type optm_type;
        typedef typename nlp::optm::optm_trait<optm_type> optm_trait;
        
        friend converge_type;
        template<typename converge_criteria_type>
        friend struct nlp::converge::converge_criteria_trait;
        friend optm_type;
        friend optm_trait;

        typedef typename nlp_settings::nlp_first_order_algorithms_type::mult_init_type mult_init_type;
        typedef nlp::outer_multiplier::init::multiplier_init_trait<mult_init_type> mult_init_trait;
        friend mult_init_trait;
        typedef typename nlp_settings::nlp_first_order_algorithms_type::mult_update_type mult_update_type;
        typedef nlp::outer_multiplier::update::multiplier_update_trait<mult_update_type> mult_update_trait;
        friend mult_update_trait;
        typedef typename nlp_settings::nlp_first_order_algorithms_type::mult_type mult_type;
        typedef nlp::outer_multiplier::multiplier_trait<mult_type> mult_trait;
        friend mult_trait;

        NlpGPlacerFirstOrder(Database &db) : NlpGPlacerBase<nlp_settings>(db) {}
    protected:
        /* calculating gradient */
        void calcGrad()
        {
            _wrapCalcGradTask.run();
        }
        /* Init */
        virtual void initProblem() override;
        void initFirstOrderGrad();
        /* Construct tasks */
        virtual void constructTasks() override;
        void constructFirstOrderTasks();
        void constructCalcPartialsTasks();
        void constructUpdatePartialsTasks();
        void constructClearGradTasks();
        void constructSumGradTask();
        void constructWrapCalcGradTask();
        /* optimization */
        virtual void optimize() override;
        /* Build the computational graph */
#ifdef IDEAPLACE_TASKFLOR_FOR_GRAD_OBJ_
        void regCalcHpwlGradTaskFlow(tf::Taskflow &tfFlow);
        void regCalcOvlGradTaskFlow(tf::Taskflow &tfFlow);
        void regCalcOobGradTaskFlow(tf::Taskflow &tfFlow);
        void regCalcAsymGradTaskFlow(tf::Taskflow &tfFlow);
        void regCalcCosGradTaskFlow(tf::Taskflow &tfFlow);
        void regCalcAllGradTaskFlow(tf::Taskflow &tfFlow);
        void regCalcGradForLoop(tf::Taskflow &tfFlow);
#endif
        

    protected:
        /* Optimization data */
        EigenVector _grad; ///< The first order graident
        EigenVector _gradHpwl; ///< The first order gradient of hpwl objective
        EigenVector _gradOvl; ///< The first order gradient  of overlapping objective
        EigenVector _gradOob; ///< The first order gradient of out of boundary objective
        EigenVector _gradAsym; ///< The first order gradient of asymmetry objective
        EigenVector _gradCos; ///< The first order gradient of cosine signal path objective
        /* Tasks */
        // Calculate the partials
        std::vector<nt::Task<nt::CalculateOperatorPartialTask<nlp_hpwl_type, EigenVector>>> _calcHpwlPartialTasks;
        std::vector<nt::Task<nt::CalculateOperatorPartialTask<nlp_ovl_type,  EigenVector>>> _calcOvlPartialTasks;
        std::vector<nt::Task<nt::CalculateOperatorPartialTask<nlp_oob_type,  EigenVector>>> _calcOobPartialTasks;
        std::vector<nt::Task<nt::CalculateOperatorPartialTask<nlp_asym_type, EigenVector>>> _calcAsymPartialTasks;
        std::vector<nt::Task<nt::CalculateOperatorPartialTask<nlp_cos_type,  EigenVector>>> _calcCosPartialTasks;
        // Update the partials
        std::vector<nt::Task<nt::UpdateGradientFromPartialTask<nlp_hpwl_type, EigenVector>>> _updateHpwlPartialTasks;
        std::vector<nt::Task<nt::UpdateGradientFromPartialTask<nlp_ovl_type,  EigenVector>>> _updateOvlPartialTasks;
        std::vector<nt::Task<nt::UpdateGradientFromPartialTask<nlp_oob_type,  EigenVector>>> _updateOobPartialTasks;
        std::vector<nt::Task<nt::UpdateGradientFromPartialTask<nlp_asym_type, EigenVector>>> _updateAsymPartialTasks;
        std::vector<nt::Task<nt::UpdateGradientFromPartialTask<nlp_cos_type,  EigenVector>>> _updateCosPartialTasks;
        // Clear the gradient. Use to clear the _gradxxx records. Needs to call before updating the partials
        nt::Task<nt::FuncTask> _clearGradTask; //FIXME: not used right noe
        nt::Task<nt::FuncTask> _clearHpwlGradTask;
        nt::Task<nt::FuncTask> _clearOvlGradTask;
        nt::Task<nt::FuncTask> _clearOobGradTask;
        nt::Task<nt::FuncTask> _clearAsymGradTask;
        nt::Task<nt::FuncTask> _clearCosGradTask;
        // Sum the _grad from individual
        nt::Task<nt::FuncTask> _sumGradTask;
        nt::Task<nt::FuncTask> _sumHpwlGradTask;
        nt::Task<nt::FuncTask> _sumOvlGradTask;
        nt::Task<nt::FuncTask> _sumOobGradTask;
        nt::Task<nt::FuncTask> _sumAsymGradTask;
        nt::Task<nt::FuncTask> _sumCosGradTask;
        // all the grads has been calculated but have not updated
#ifdef IDEAPLACE_TASKFLOR_FOR_GRAD_OBJ_
        nt::Task<nt::EmptyTask> _endGradCalcTask;
        nt::Task<nt::EmptyTask> _beginGradCalcTask;
#endif
        nt::Task<nt::FuncTask> _wrapCalcGradTask; ///< For debugging: calculating the gradient and sum them
};

PROJECT_NAMESPACE_END
#endif //IDEAPLACE_NLPGPLACER_H_
