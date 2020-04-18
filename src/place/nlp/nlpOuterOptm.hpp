/**
 * @file nlpOuterOptm.hpp
 * @brief The non-lnear programming outer-problem optimization
 * @author Keren Zhu
 * @date 04/09/2020
 */

#pragma once

#include "global/global.h"
#include "nlpTypes.hpp"

PROJECT_NAMESPACE_BEGIN

namespace nlp 
{
    namespace outer_stop_condition
    {
        /* Stop outer problem condition */
        template<typename T>
        struct stop_condition_trait 
        {
            // static T construct(NlpType &)
            // static IntType stopPlaceCondition(T&, NlpType &)
        };

        /// @brief stop condition with number of iterations
        template<IntType MaxIter=10>
        struct stop_after_num_outer_iterations
        {
            static constexpr IntType maxIter = MaxIter;
            IntType curIter = 0;
        };
        
        template<IntType MaxIter>
        struct stop_condition_trait<stop_after_num_outer_iterations<MaxIter>>
        {
            template<typename NlpType>
            static stop_after_num_outer_iterations<MaxIter> construct(NlpType &) { return stop_after_num_outer_iterations<MaxIter>(); }
            template<typename NlpType>
            static IntType stopPlaceCondition(NlpType &, stop_after_num_outer_iterations<MaxIter> &stop)
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
    } // namespace outer_stop_condition

    namespace outer_multiplier
    {
        /// @brief the multiplier is categoried by types
        template<typename mult_type>
        struct is_mult_type_dependent_diff : std::false_type {};

        namespace init
        {
            template<typename T>
            struct multiplier_init_trait {};

            struct hard_code_init {};

            template<>
            struct multiplier_init_trait<hard_code_init>
            {
                template<typename nlp_type, typename mult_type>
                static void init(nlp_type &, mult_type& mult)
                {
                    mult._constMults.at(0) = 16; // hpwl
                    mult._constMults.at(1) = 16; // cos
                    mult._variedMults.at(0) = 1; // ovl
                    mult._variedMults.at(1) = 1; // oob
                    mult._variedMults.at(2) = 1; // asym
                }

            };

            /// @brief match by gradient norm
            struct init_by_matching_gradient_norm 
            {
                static constexpr RealType penaltyRatioToObj = 1.0; ///< The ratio of targeting penalty 
                static constexpr RealType small = 0.01;
            };

            
            template<>
            struct multiplier_init_trait<init_by_matching_gradient_norm>
            {
                typedef init_by_matching_gradient_norm init_type;
                template<typename nlp_type, typename mult_type, std::enable_if_t<nlp::is_first_order_diff<nlp_type>::value, void>* = nullptr>
                static void init(nlp_type &nlp, mult_type& mult)
                {
                    mult._constMults.at(0) = 1.0; // hpwl
                    const auto hpwlMult = mult._constMults.at(0);
                    const auto hpwlNorm = nlp._gradHpwl.norm();
                    const auto hpwlMultNorm = hpwlMult * hpwlNorm;
                    const auto hpwlMultNormPenaltyRatio = hpwlMultNorm * init_type::penaltyRatioToObj;
                    Assert(nlp._gradHpwl.norm() > REAL_TYPE_TOL);
                    const auto cosNorm = nlp._gradCos.norm();
                    const auto ovlNorm = nlp._gradOvl.norm();
                    const auto oobNorm = nlp._gradOob.norm();
                    const auto asymNorm = nlp._gradAsym.norm();
                    const auto maxPenaltyNorm = ovlNorm;
                    // Make a threshold on by referencing hpwl to determine whether one is small
                    const auto small  = init_type::small * hpwlNorm;
                    // match gradient norm for signal path
                    if (cosNorm > small)
                    {
                        mult._constMults.at(1) = hpwlMultNorm / cosNorm;
                    }
                    else
                    {
                        mult._constMults.at(1) = hpwlMult;
                    }
                    // overlap
                    if (ovlNorm > small)
                    {
                        mult._variedMults.at(0) = hpwlMultNormPenaltyRatio / ovlNorm;
                    }
                    else
                    {
                        mult._variedMults.at(0) = hpwlMultNormPenaltyRatio / maxPenaltyNorm;
                    }
                    // out of boundary
                    // Since we know oob is small at beginning, and it will take effect after a few iterations. Therefore it would be better to set it to resonable range first
                    //mult._variedMults.at(1) = hpwlMultNormPenaltyRatio;
                    if (oobNorm > small)
                    {
                        mult._variedMults.at(1) = hpwlMultNormPenaltyRatio / oobNorm;
                    }
                    else
                    {
                        mult._variedMults.at(1) = hpwlMultNormPenaltyRatio / maxPenaltyNorm;
                    }
                    // asym
                    if (asymNorm > small)
                    {
                        mult._variedMults.at(2) = hpwlMultNormPenaltyRatio / asymNorm;
                    }
                    else
                    {
                        mult._variedMults.at(2) = hpwlMultNormPenaltyRatio / maxPenaltyNorm;
                    }
                    DBG("init mult: hpwl %f cos %f \n",
                            mult._constMults[0], mult._constMults[1]);
                    DBG("init mult: ovl %f oob %f asym %f \n",
                            mult._variedMults[0], mult._variedMults[1], mult._variedMults[2]);
                }
            };
        }; // namespace init
        namespace update
        {
            template<typename T>
            struct multiplier_update_trait {};

            /// @brief increase the total amounts of penalty of by a constant
            struct shared_constant_increase_penalty
            {
                static constexpr RealType penalty  = 20;
            };

            template<>
            struct multiplier_update_trait<shared_constant_increase_penalty>
            {
                typedef shared_constant_increase_penalty update_type;

                template<typename nlp_type, typename mult_type, std::enable_if_t<is_mult_type_dependent_diff<mult_type>::value, void>* = nullptr>
                static update_type construct(nlp_type &, mult_type&) { return update_type(); }

                template<typename nlp_type, typename mult_type,  std::enable_if_t<is_mult_type_dependent_diff<mult_type>::value, void>* = nullptr>
                static void update(nlp_type &nlp, mult_type &mult, update_type &update)
                {
                    nlp._wrapObjAllTask.run();
                    const auto rawOvl = nlp._objOvl / mult._variedMults.at(0);
                    const auto rawOob = nlp._objOob / mult._variedMults.at(1);
                    const auto rawAsym = nlp._objAsym / mult._variedMults.at(2);
                    const auto fViolate = rawOvl + rawOob + rawAsym;
                    DBG("update mult: raw ovl %f oob %f asym %f total %f \n", rawOvl, rawOob, rawAsym, fViolate);
                    DBG("update mult:  before ovl %f oob %f asym %f \n",
                            mult._variedMults[0], mult._variedMults[1], mult._variedMults[2]);
                    mult._variedMults.at(0) += update.penalty * (rawOvl / fViolate);
                    mult._variedMults.at(1) += update.penalty * (rawOob / fViolate);
                    mult._variedMults.at(2) += update.penalty * (rawAsym / fViolate);
                    DBG("update mult: afterafter  ovl %f oob %f asym %f \n",
                            mult._variedMults[0], mult._variedMults[1], mult._variedMults[2]);
                }
            };

            /// @brief direct subgradient
            struct direct_subgradient
            {
                static constexpr RealType stepSize = 1;
            };

            template<>
            struct multiplier_update_trait<direct_subgradient>
            {
                typedef direct_subgradient update_type;

                template<typename nlp_type, typename mult_type, std::enable_if_t<is_mult_type_dependent_diff<mult_type>::value, void>* = nullptr>
                static update_type construct(nlp_type &, mult_type&) { return update_type(); }

                template<typename nlp_type, typename mult_type, std::enable_if_t<is_mult_type_dependent_diff<mult_type>::value, void>* = nullptr>
                static void init(nlp_type &, mult_type&, update_type &) { }

                template<typename nlp_type, typename mult_type,  std::enable_if_t<is_mult_type_dependent_diff<mult_type>::value, void>* = nullptr>
                static void update(nlp_type &nlp, mult_type &mult, update_type &update)
                {
                    nlp._wrapObjAllTask.run();
                    const auto rawOvl = nlp._objOvl / mult._variedMults.at(0);
                    const auto rawOob = nlp._objOob / mult._variedMults.at(1);
                    const auto rawAsym = nlp._objAsym / mult._variedMults.at(2);
                    DBG("update mult: raw ovl %f oob %f asym %f total %f \n", rawOvl, rawOob, rawAsym);
                    DBG("update mult:  before ovl %f oob %f asym %f \n",
                            mult._variedMults[0], mult._variedMults[1], mult._variedMults[2]);
                    mult._variedMults.at(0) += update.stepSize * (rawOvl );
                    mult._variedMults.at(1) += update.stepSize * (rawOob );
                    mult._variedMults.at(2) += update.stepSize * (rawAsym );
                    DBG("update mult: afterafter  ovl %f oob %f asym %f \n",
                            mult._variedMults[0], mult._variedMults[1], mult._variedMults[2]);
                }
            };

            /// @breif subgradient normalize by values in iter 0
            template<typename nlp_numerical_type>
            struct subgradient_normalized_by_init
            {
                static constexpr nlp_numerical_type stepSize = 10;
                std::vector<nlp_numerical_type> normalizeFactor;
            };

            template<typename nlp_numerical_type>
            struct multiplier_update_trait<subgradient_normalized_by_init<nlp_numerical_type>>
            {
                typedef subgradient_normalized_by_init<nlp_numerical_type> update_type;

                template<typename nlp_type, typename mult_type, std::enable_if_t<is_mult_type_dependent_diff<mult_type>::value, void>* = nullptr>
                static update_type construct(nlp_type &, mult_type&) { return update_type(); }

                template<typename nlp_type, typename mult_type, std::enable_if_t<is_mult_type_dependent_diff<mult_type>::value, void>* = nullptr>
                static void init(nlp_type &nlp, mult_type &mult, update_type &update) 
                { 
                    update.normalizeFactor.resize(3);
                    update.normalizeFactor.at(0) = mult._variedMults.at(0) / nlp._objOvl;
                    update.normalizeFactor.at(1) = 1;// mult._variedMults.at(1) / nlp._objOob;
                    update.normalizeFactor.at(2) = mult._variedMults.at(2) / nlp._objAsym;
                    //update.normalizeFactor.at(0) = 1 / nlp._objOvl;
                    //update.normalizeFactor.at(1) = 1;// / nlp._objOob;
                    //update.normalizeFactor.at(2) = 1 / nlp._objAsym;
                }

                template<typename nlp_type, typename mult_type,  std::enable_if_t<is_mult_type_dependent_diff<mult_type>::value, void>* = nullptr>
                static void update(nlp_type &nlp, mult_type &mult, update_type &update)
                {
                    nlp._wrapObjAllTask.run();
                    const auto rawOvl = nlp._objOvl / mult._variedMults.at(0);
                    const auto rawOob = nlp._objOob / mult._variedMults.at(1);
                    const auto rawAsym = nlp._objAsym / mult._variedMults.at(2);
                    const auto normalizedOvl = rawOvl * update.normalizeFactor.at(0);
                    const auto normalizedOob = rawOob * update.normalizeFactor.at(1);
                    const auto normalizedAsym = rawAsym  * update.normalizeFactor.at(2);
                    DBG("update mult: raw ovl %f oob %f asym %f total %f \n", normalizedOvl, normalizedOob, normalizedAsym);
                    DBG("update mult:  before ovl %f oob %f asym %f \n",
                            mult._variedMults[0], mult._variedMults[1], mult._variedMults[2]);
                    mult._variedMults.at(0) += update.stepSize * (normalizedOvl );
                    mult._variedMults.at(1) += update.stepSize * (normalizedOob );
                    mult._variedMults.at(2) += update.stepSize * (normalizedAsym );
                    DBG("update mult: afterafter  ovl %f oob %f asym %f \n",
                            mult._variedMults[0], mult._variedMults[1], mult._variedMults[2]);
                }
            };


        } // namespace update
        template<typename T> 
        struct multiplier_trait {};

        template<typename nlp_numerical_type, typename init_type, typename update_type>
        struct mult_const_hpwl_cos_and_penalty_by_type
        {
            typedef init_type mult_init_type;
            typedef update_type mult_update_type;
            std::vector<nlp_numerical_type> _constMults; ///< constant mults
            std::vector<nlp_numerical_type> _variedMults; ///< varied penalty multipliers
            update_type update;
        };

        template<typename nlp_numerical_type, typename init_type, typename update_type>
        struct is_mult_type_dependent_diff<mult_const_hpwl_cos_and_penalty_by_type<nlp_numerical_type, init_type, update_type>> : std::true_type {};

        template<typename nlp_numerical_type, typename init_type, typename update_type>
        struct multiplier_trait<mult_const_hpwl_cos_and_penalty_by_type<nlp_numerical_type, init_type, update_type>>
        {
            typedef mult_const_hpwl_cos_and_penalty_by_type<nlp_numerical_type, init_type, update_type> mult_type;

            template<typename nlp_type>
            static mult_type construct(nlp_type &nlp)
            {
                mult_type mult;
                mult._constMults.resize(2, 0.0);
                mult._variedMults.resize(3, 0.0);
                mult.update = update::multiplier_update_trait<update_type>::construct(nlp, mult);
                return mult;
            }

            template<typename nlp_type>
            static void init(nlp_type &nlp, mult_type &mult)
            {
                init::multiplier_init_trait<init_type>::init(nlp, mult);
                update::multiplier_update_trait<update_type>::init(nlp, mult, mult.update);
                for (auto &op : nlp._hpwlOps) { op._getLambdaFunc = [&](){ return mult._constMults[0]; }; }
                for (auto &op : nlp._cosOps) { op._getLambdaFunc = [&](){ return mult._constMults[1]; }; }
                for (auto &op : nlp._ovlOps) { op._getLambdaFunc = [&](){ return mult._variedMults[0]; }; }
                for (auto &op : nlp._oobOps) { op._getLambdaFunc = [&](){ return mult._variedMults[1]; }; }
                for (auto &op : nlp._asymOps) { op._getLambdaFunc = [&](){ return mult._variedMults[2]; }; }
            }

            template<typename nlp_type>
            static void update(nlp_type &nlp, mult_type &mult)
            {
                update::multiplier_update_trait<update_type>::update(nlp, mult, mult.update);
            }

        };
    } // namespace outer_multiplier
} //namespace nlp
PROJECT_NAMESPACE_END
