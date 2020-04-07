/**
 * @file nlpOptmKernel.hpp
 * @brief The non-lnear programming basic optimization kernel for global plaement
 * @author Keren Zhu
 * @date 04/05/2020
 */

#pragma once

#include "global/global.h"
#include "nlpTypes.hpp"

PROJECT_NAMESPACE_BEGIN

namespace nlp
{
    namespace converge
    {
        template<typename converge_type>
        struct converge_criteria_trait {};

        /// @brief used for place holder
        struct converge_criteria_empty {};

        template<>
        struct converge_criteria_trait<converge_criteria_empty>
        {
            typedef converge_criteria_empty converge_type;
            template<typename nlp_type, typename optm_type>
            static constexpr BoolType stopCriteria(nlp_type &, optm_type &, converge_type&) { return false;}
            template<typename converge_type>
            static constexpr void clear(converge_type&) {}
        };


        template<IntType max_iter=20>
        struct converge_criteria_max_iter
        {
            IntType _iter = 0;
        };

        template<IntType max_iter>
        struct converge_criteria_trait<converge_criteria_max_iter<max_iter>>
        {
            typedef converge_criteria_max_iter<max_iter> converge_type;
            static void clear(converge_type &c)
            {
                c._iter = 0;
            }
            template<typename nlp_type, typename optm_type>
            static BoolType stopCriteria(nlp_type &, optm_type &, converge_type &c)
            {
                if (c._iter >= max_iter)
                {
                    clear(c);
                    return true;
                }
                else
                {
                    ++c._iter;
                    return false;
                }
            }
        };

        template<typename nlp_numerical_type=RealType>
        struct converge_criteria_stop_improve_less_in_last_two_step
        {
            static constexpr RealType _stopThreshold = 0.0001;
            nlp_numerical_type _fLastStep = - 1.0;
            nlp_numerical_type _fLastLastStep = - 1.0;
        };
        template<typename nlp_numerical_type>
        struct converge_criteria_trait<converge_criteria_stop_improve_less_in_last_two_step<nlp_numerical_type>>
        {
            typedef converge_criteria_stop_improve_less_in_last_two_step<nlp_numerical_type> converge_type;
            static void clear(converge_type &c)
            {
                c._fLastStep = - 1.0;
                c._fLastLastStep = - 1.0;
            }
            template<typename nlp_type, typename optm_type>
            static BoolType stopCriteria(nlp_type &n, optm_type &, converge_type &c)
            {
                if (c._fLastLastStep >= 0.0)
                {
                    if ((c._fLastLastStep - n._obj) / c._fLastLastStep < c._stopThreshold)
                    {
                        //clear(c);
                        //return true;
                    }
                }
                c._fLastLastStep = c._fLastStep;
                c._fLastStep = n._obj;
                return false;
            }
        };

        template<typename converge_type_1=converge_criteria_empty,
            typename converge_type_2=converge_criteria_empty,
            typename converge_type_3=converge_criteria_empty,
            typename converge_type_4=converge_criteria_empty,
            typename converge_type_5=converge_criteria_empty>
        struct converge_criteria_list
        {
            typedef converge_type_1 c1_type;
            typedef converge_type_2 c2_type;
            typedef converge_type_3 c3_type;
            typedef converge_type_4 c4_type;
            typedef converge_type_5 c5_type;
            converge_type_1 _c1;
            converge_type_2 _c2;
            converge_type_3 _c3;
            converge_type_4 _c4;
            converge_type_5 _c5;
        };
        template<typename converge_type_1,
            typename converge_type_2,
            typename converge_type_3,
            typename converge_type_4,
            typename converge_type_5>
        struct converge_criteria_trait<converge_criteria_list<converge_type_1, converge_type_2, converge_type_3, converge_type_4, converge_type_5>>
        {
            typedef converge_criteria_list<converge_type_1, converge_type_2, converge_type_3, converge_type_4, converge_type_5> converge_type;
            typedef typename converge_type::c1_type c1_type;
            typedef typename converge_type::c2_type c2_type;
            typedef typename converge_type::c3_type c3_type;
            typedef typename converge_type::c4_type c4_type;
            typedef typename converge_type::c5_type c5_type;
            static void clear(converge_type &c)
            {
                converge_criteria_trait<c1_type>::clear(c._c1);
                converge_criteria_trait<c2_type>::clear(c._c2);
                converge_criteria_trait<c3_type>::clear(c._c3);
                converge_criteria_trait<c4_type>::clear(c._c4);
                converge_criteria_trait<c5_type>::clear(c._c5);
            }
            template<typename nlp_type, typename optm_type>
            static constexpr BoolType stopCriteria(nlp_type &n, optm_type &o, converge_type&c) 
            { 
                BoolType stop = false;
                stop = converge_criteria_trait<c1_type>::stopCriteria(n, o, c._c1);
                if (!stop) { stop = converge_criteria_trait<c2_type>::stopCriteria(n, o, c._c2); }
                if (!stop) { stop = converge_criteria_trait<c3_type>::stopCriteria(n, o, c._c3); }
                if (!stop) { stop = converge_criteria_trait<c4_type>::stopCriteria(n, o, c._c4); }
                if (!stop) { stop = converge_criteria_trait<c5_type>::stopCriteria(n, o, c._c5); }
                if (stop) { clear(c); } // might redundantly clear the criteria. But it's okay for now.
                return stop;
            }
        };
    } // namespace converge
} // namespace nlp

PROJECT_NAMESPACE_END
