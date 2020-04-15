/**
 * @file nlpOuterOptm.hpp
 * @brief The non-lnear programming outer-problem optimization
 * @author Keren Zhu
 * @date 04/09/2020
 */

#pragma once

#include "global/global.h"

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
                    const auto rawOvl = nlp._objOvl / mult._variedMults.at(0);
                    const auto rawOob = nlp._objOob / mult._variedMults.at(1);
                    const auto rawAsym = nlp._objAsym / mult._variedMults.at(2);
                    const auto fViolate = rawOvl + rawOob + rawAsym;
                    mult._variedMults.at(0) += update.penalty * (rawOvl / fViolate);
                    mult._variedMults.at(1) += update.penalty * (rawOob / fViolate);
                    mult._variedMults.at(2) += update.penalty * (rawAsym / fViolate);
                    DBG("update mult: ovl %f oob %f asym %f \n",
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
