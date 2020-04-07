/**
 * @file nlpFirstOrderKernel.hpp
 * @brief The non-lnear programming first order optimization kernel for global plaement
 * @author Keren Zhu
 * @date 04/05/2020
 */

#pragma once

#include "global/global.h"
#include "nlpOptmKernels.hpp"

PROJECT_NAMESPACE_BEGIN

namespace nlp
{
    namespace optm
    {
        template<typename optm_type>
        struct optm_trait {};
        namespace first_order
        {
            /// @brief Naive gradient descent. It will stop if maxIter reaches or the improvement 
            template<typename converge_criteria_type>
            struct naive_gradient_descent
            {
                typedef converge_criteria_type converge_type;
                static constexpr RealType _stepSize = 0.01;
                converge_criteria_type _converge;
                friend nlp::converge::converge_criteria_trait<converge_criteria_type>;
                friend nlp::optm::optm_trait<naive_gradient_descent<converge_criteria_type>>;
            };
        } // namspace first_order
        template<typename converge_criteria_type>
        struct optm_trait<first_order::naive_gradient_descent<converge_criteria_type>>
        {
            typedef first_order::naive_gradient_descent<converge_criteria_type> optm_type;
            typedef typename optm_type::converge_type converge_type;
            typedef nlp::converge::converge_criteria_trait<converge_type> converge_trait;
            template<typename nlp_type, std::enable_if_t<nlp::is_first_order_diff<nlp_type>::value, void>* = nullptr>
            static void optimize(nlp_type &n, optm_type &o)
            {
                converge_trait::clear(o._converge);
                do 
                {
                    n.calcGrad();
                    n._pl -= o._stepSize * n._grad;
                    n.calcObj();
                } while (!converge_trait::stopCriteria(n, o, o._converge));
            }
        };
    } // namespace optm
} // namespace nlp

PROJECT_NAMESPACE_END
