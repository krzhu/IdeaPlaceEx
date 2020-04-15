/**
 * @file conjugateGradientWnlib.hpp
 * @brief Implementation of wnlib interface for conjugate gradient
 * @author Keren Zhu
 * @date 04/15/2020
 */

#include "conjugateGradientWnlib.hpp"
#include "wnconj.h"

PROJECT_NAMESPACE_BEGIN

namespace nlp { namespace optm { namespace first_order { namespace _wnlib {

wnlibWrapper<double>  *wrapperwrapperPtr;
RealType objFuncWrapper(RealType *vec)
{
    return wrapperwrapperPtr->objFunc(vec);
}

void gradFuncWrapper(RealType *vec1, RealType *vec2)
{
    return wrapperwrapperPtr->gradFunc(vec1, vec2);
}

template<typename nlp_numerical_type>
void wnlibWrapper<nlp_numerical_type>::optimize()
{
    int code;
    wn_conj_gradient_method(&code, &_lowerBoundary, _pl, _numVariables, objFuncWrapper, gradFuncWrapper, 5000);
}

template struct wnlibWrapper<double>;
}}}}//namespace nlp { namespace optm { namespace first_order { namespace _wnlib {
PROJECT_NAMESPACE_END
