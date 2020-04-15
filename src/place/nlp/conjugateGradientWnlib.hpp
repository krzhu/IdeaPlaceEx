
#pragma once

#include "global/global.h"
#include "nlpOptmKernels.hpp"
#include "nlpFirstOrderKernel.hpp"

PROJECT_NAMESPACE_BEGIN

namespace nlp
{
    namespace optm
    {
        namespace first_order
        {
            namespace _wnlib
            {
                template<typename nlp_numerical_type>
                struct wnlibWrapper
                {
                    wnlibWrapper(std::function<nlp_numerical_type(void)> evaObj, std::function<void(void)> evaGrad, nlp_numerical_type *pl, nlp_numerical_type *grad, IndexType numVariables)
                        : _evaObj(evaObj), _evaGrad(evaGrad),  _numVariables(numVariables)
                    {
                        _pl = pl;
                        _grad = grad;
                    }
                    void optimize();
                    nlp_numerical_type objFunc(nlp_numerical_type *sol)
                    {
                        for (IndexType i = 0; i < _numVariables; ++i)
                        {
                            _pl[i] = sol[i];
                        }
                        return _evaObj();
                    }
                    void gradFunc(nlp_numerical_type *sol,  nlp_numerical_type *grad)
                    {
                        for (IndexType i = 0; i < _numVariables; ++i)
                        {
                            _pl[i] = sol[i];
                        }
                        _evaGrad();
                        for (IndexType i = 0; i < _numVariables; ++i)
                        {
                            grad[i] = _grad[i];
                        }
                    }
                    std::function<nlp_numerical_type(void)> _evaObj;
                    std::function<void(void)> _evaGrad;
                    nlp_numerical_type *_pl;
                    nlp_numerical_type *_grad;
                    IndexType _numVariables;
                    double _lowerBoundary = -100.0;
                };
                extern wnlibWrapper<double>  *wrapperwrapperPtr;
            } // namespace _wnlib
            /// @brief conjugate gradient method using wn lib
            struct conjugate_gradient_wnlib 
            {
                template<typename nlp_type>
                void optimize(nlp_type &n)
                {
                    auto evaObj = [&](){ n._wrapObjAllTask.run(); return n._obj; };
                    auto evaGrad = [&](){ n._wrapCalcGradTask.run(); };
                    _wnlib::wnlibWrapper<typename nlp_type::nlp_numerical_type> wrapper(evaObj, evaGrad, n._pl.data(), n._grad.data(), n._numVariables);
                    _wnlib::wrapperwrapperPtr = &wrapper;
                    wrapper.optimize();
                }
            };
        } // namspace first_order
        template<>
        struct optm_trait<first_order::conjugate_gradient_wnlib>
        {
            typedef first_order::conjugate_gradient_wnlib optm_type;
            template<typename nlp_type, std::enable_if_t<nlp::is_first_order_diff<nlp_type>::value, void>* = nullptr>
            static void optimize(nlp_type &n, optm_type &o)
            {
                o.optimize(n);
            };
        };
    } // namespace optm
} // namespace nlp

PROJECT_NAMESPACE_END
