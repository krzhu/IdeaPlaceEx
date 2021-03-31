/**
 * @file conjugateGradientWnlib.hpp
 * @brief Implementation of wnlib interface for conjugate gradient
 * @author Keren Zhu
 * @date 04/15/2020
 */

#include "conjugateGradientWnlib.hpp"
#include "wnconj.h"

PROJECT_NAMESPACE_BEGIN

namespace nlp {
namespace optm {
namespace first_order {
namespace _wnlib {

wnlibWrapper<double> *wrapperwrapperPtr;
RealType objFuncWrapper(RealType *vec) {
  return wrapperwrapperPtr->objFunc(vec);
}

void gradFuncWrapper(RealType *vec1, RealType *vec2) {
  wrapperwrapperPtr->gradFunc(vec1, vec2);
}

template <typename nlp_numerical_type>
void wnlibWrapper<nlp_numerical_type>::optimize() {
  int code;
  double valMin;
  RealType *sol;
  sol = (RealType *)malloc(sizeof(RealType) * _numVariables);
  for (IndexType i = 0; i < _numVariables; ++i) {
    sol[i] = 1.0;
  }

  wn_conj_gradient_method(&code, &valMin, sol, _numVariables, objFuncWrapper,
                          gradFuncWrapper, 3000);
  // wn_conj_direction_method(&code, &valMin, _pl, initial_coord_x0s,
  // _numVariables, objFuncWrapper, 1000);
  for (IndexType i = 0; i < _numVariables; ++i) {
    _pl[i] = sol[i];
  }
  free(sol);
}

template struct wnlibWrapper<double>;
} // namespace _wnlib
} // namespace first_order
} // namespace optm
} // namespace nlp
PROJECT_NAMESPACE_END
