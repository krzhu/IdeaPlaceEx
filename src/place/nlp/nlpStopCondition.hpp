/**
 * @file nlpStopCondition.hpp
 * @brief The non-lnear programming check whether stop outer optimization algorithms for global plaement
 * @author Keren Zhu
 * @date 04/02/2020
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
        struct stop_after_num_outer_iterations
        {
            static constexpr IntType maxIter = 20;
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
    } // namespace outer_stop_condition
} //namespace nlp

PROJECT_NAMESPACE_END

