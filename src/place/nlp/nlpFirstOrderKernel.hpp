/**
 * @file nlpFirstOrderKernel.hpp
 * @brief The non-lnear programming first order optimization kernel for global
 * plaement
 * @author Keren Zhu
 * @date 04/05/2020
 */

#pragma once

#include "global/global.h"
#include "nlpOptmKernels.hpp"

PROJECT_NAMESPACE_BEGIN

namespace nlp {
namespace optm {
template <typename optm_type> struct optm_trait {};
namespace first_order {
/// @brief Naive gradient descent. It will stop if maxIter reaches or the
/// improvement
template <typename converge_criteria_type, typename nlp_numerical_type> struct naive_gradient_descent {
  typedef converge_criteria_type converge_type;
  static constexpr nlp_numerical_type _stepSize = 0.0003;
  converge_criteria_type _converge;
};
/// @brief adam
template <typename converge_criteria_type, typename nlp_numerical_type>
struct adam {
  typedef converge_criteria_type converge_type;
  converge_criteria_type _converge;
  static constexpr nlp_numerical_type alpha = 0.0003;
  static constexpr nlp_numerical_type beta1 = 0.9;
  static constexpr nlp_numerical_type beta2 = 0.999;
  static constexpr nlp_numerical_type epsilon = 1e-8;

  static constexpr nlp_numerical_type naiveGradientDescentStepSize = 0.001;
};
/// @brief nesterov accelerated gradient
template <typename converge_criteria_type, typename nlp_numerical_type>
struct nesterov {
  typedef converge_criteria_type converge_type;
  converge_criteria_type _converge;
  static constexpr nlp_numerical_type eta = 0.003;
};
} // namespace first_order
template <typename converge_criteria_type, typename nlp_numerical_type>
struct optm_trait<first_order::naive_gradient_descent<converge_criteria_type, nlp_numerical_type>> {
  typedef first_order::naive_gradient_descent<converge_criteria_type, nlp_numerical_type> optm_type;
  typedef typename optm_type::converge_type converge_type;
  typedef nlp::converge::converge_criteria_trait<converge_type> converge_trait;
  template <typename nlp_type,
            std::enable_if_t<nlp::is_first_order_diff<nlp_type>::value, void>
                * = nullptr>
  static void optimize(nlp_type &n, optm_type &o) {
    converge_trait::clear(o._converge);
    IndexType iter = 0;
    do {
      n.calcGrad();
      n._pl -= optm_type::_stepSize * n._grad;
      ++iter;
    } while (!converge_trait::stopCriteria(n, o, o._converge));
#ifdef DEBUG_GR
    DBG("naive gradient decesent: %f hpwl %f cos %f ovl %f oob %f asym %f \n",
        n._obj, n._objHpwl, n._objCos, n._objOvl, n._objOob, n._objAsym);
    DBG("gradient norm %f \n", n._grad.norm());
    DBG("converge at iter %d \n", iter);
#endif
  }
};
template <typename converge_criteria_type, typename nlp_numerical_type>
struct optm_trait<
    first_order::adam<converge_criteria_type, nlp_numerical_type>> {
  typedef first_order::adam<converge_criteria_type, nlp_numerical_type>
      optm_type;
  typedef typename optm_type::converge_type converge_type;
  typedef nlp::converge::converge_criteria_trait<converge_type> converge_trait;
  template <typename nlp_type,
            std::enable_if_t<nlp::is_first_order_diff<nlp_type>::value, void>
                * = nullptr>
  static void optimize(nlp_type &n, optm_type &o) {
    converge_trait::clear(o._converge);
    const IndexType numVars = n._numVariables;
    typename nlp_type::EigenVector m, v;
    m.resize(numVars);
    m.setZero();
    v.resize(numVars);
    v.setZero();
    IndexType iter = 0;
    double initNorm = 1;
    do {
      n.calcGrad();


      m = o.beta1 * m + (1 - o.beta1) * n._grad;
      v = o.beta2 * v + (1 - o.beta2) * n._grad.cwiseProduct(n._grad);
      auto mt = m / (1 - pow(o.beta1, iter+1));
      auto vt = v / (1 - pow(o.beta2, iter+1));
      auto bot = vt.array().sqrt() + o.epsilon;
      n._pl = n._pl - o.alpha * (mt.array() / bot).matrix();
#if 0
      for (int i = 0; i < n._pl.size(); ++i) {
        if (std::isnan(n._pl(i))) {
          std::cout<<"\n\n\nNAN!!!\n\n\n";
          //std::cout<<"hpwl "<< n._gradHpwl<<std::endl;
          //std::cout<<"ovl "<< n._gradOvl<<std::endl;
          //std::cout<<"asyn "<< n._gradAsym<<std::endl;
          //std::cout<<"fence "<< n._gradFence<<std::endl;
          //std::cout<<"grad" <<n._grad<<std::endl;
          std::cout<< "pl"<<n._pl<<std::endl;
          assert(false);
        }
      }
#endif

      // n.calcObj();
      if (iter == 0) {
        initNorm = n._grad.norm();
         DBG(" init norm %f \n", initNorm);
      }
      ++iter;
      // DBG("adam: %f hpwl %f cos %f ovl %f oob %f asym %f \n", n._obj,
      // n._objHpwl, n._objCos, n._objOvl, n._objOob, n._objAsym);
    } while (!converge_trait::stopCriteria(n, o, o._converge));
    n.calcObj();
//#ifdef DEBUG_GR
    DBG("adam: %f hpwl %f cos %f ovl %f oob %f asym %f fence %f \n", n._obj, n._objHpwl,
        n._objCos, n._objOvl, n._objOob, n._objAsym, n._objFence);
    DBG("gradient norm %f, ratio %f \n", n._grad.norm(), n._grad.norm() / initNorm);
    DBG("converge at iter %d \n", iter);
    DBG("fence norm %f \n", n._gradFence.norm());
//#endif
  }
};

template <typename converge_criteria_type, typename nlp_numerical_type>
struct optm_trait<
    first_order::nesterov<converge_criteria_type, nlp_numerical_type>> {
  typedef first_order::nesterov<converge_criteria_type, nlp_numerical_type>
      optm_type;
  typedef typename optm_type::converge_type converge_type;
  typedef nlp::converge::converge_criteria_trait<converge_type> converge_trait;
  template <typename nlp_type,
            std::enable_if_t<nlp::is_first_order_diff<nlp_type>::value, void>
                * = nullptr>
  static void optimize(nlp_type &n, optm_type &o) {
    converge_trait::clear(o._converge);
    const IndexType numVars = n._numVariables;
    IndexType iter = 0;
    nlp_numerical_type lambdaPrev = 0;
    nlp_numerical_type lambdaCurr = 1;
    nlp_numerical_type gamma = 1;
    typename nlp_type::EigenVector yCurr, yPrev;
    yCurr.resize(numVars);
    yCurr.setZero();
    yPrev.resize(numVars);
    yPrev = n._pl;
    do {
      ++iter;
      n.calcGrad();
      yCurr = n._pl - o.eta * n._grad;
      n._pl = (1 - gamma) * yCurr + gamma * yPrev;

      yPrev = yCurr;

      const auto lambdaTemp = lambdaCurr;
      lambdaCurr = (1 + std::sqrt(1 + 4 * lambdaPrev * lambdaPrev)) * 0.5;
      lambdaPrev = lambdaTemp;

      gamma = (1 - lambdaPrev) / lambdaCurr;
      gamma = std::min(gamma, 0.999999999);
      gamma = std::max(gamma, 1e-8);
    } while (!converge_trait::stopCriteria(n, o, o._converge));
    n.calcObj();
#ifdef DEBUG_GR
    DBG("nesterov: %f hpwl %f cos %f ovl %f oob %f asym %f \n", n._obj,
        n._objHpwl, n._objCos, n._objOvl, n._objOob, n._objAsym);
    DBG("gradient norm %f \n", n._grad.norm());
    DBG("converge at iter %d \n", iter);
#endif
  }
};

} // namespace optm
} // namespace nlp

PROJECT_NAMESPACE_END
