/**
 * @file nlpOuterOptm.hpp
 * @brief The non-lnear programming outer-problem optimization
 * @author Keren Zhu
 * @date 04/09/2020
 */

#pragma once

#include "stop_condition.hpp"
#include "multiplier.hpp"

PROJECT_NAMESPACE_BEGIN

namespace nlp {

/// @brief alpha
namespace alpha {
namespace update {
template <typename T> struct alpha_update_trait {};

} // namespace update

template <typename T> struct alpha_trait {};

template <typename nlp_numerical_type> struct alpha_hpwl_ovl_oob {
  std::vector<nlp_numerical_type> _alpha;
};

template <typename nlp_numerical_type>
struct alpha_trait<alpha_hpwl_ovl_oob<nlp_numerical_type>> {
  typedef alpha_hpwl_ovl_oob<nlp_numerical_type> alpha_type;

  template <typename nlp_type> static alpha_type construct(nlp_type &) {
    alpha_type alpha;
    alpha._alpha.resize(5, 1.0);
    return alpha;
  }

  template <typename nlp_type>
  static void init(nlp_type &nlp, alpha_type &alpha) {
    for (auto &op : nlp._hpwlOps) {
      op.setGetAlphaFunc([&]() { return alpha._alpha[0]; });
    }
    for (auto &op : nlp._ovlOps) {
      op.setGetAlphaFunc([&]() { return alpha._alpha[1]; });
    }
    for (auto &op : nlp._oobOps) {
      op.setGetAlphaFunc([&]() { return alpha._alpha[2]; });
    }
    for (auto &op : nlp._crfOps) {
      op.setGetAlphaFunc([&]() { return alpha._alpha[3]; });
    }
    for (auto &op : nlp._fenceOps) {
      op.setGetAlphaFunc([&]() { return alpha._alpha[4]; });
    }
  }
  static std::function<nlp_numerical_type(void)>
  fenceGetAlphaFunc(alpha_type &alpha) {
    return [&]() { return alpha._alpha[4]; };
  }
  static std::function<nlp_numerical_type(void)>
  ovlGetAlphaFunc(alpha_type &alpha) {
    return [&]() { return alpha._alpha[1]; };
  }
};

namespace update {
/// @breif update the alpha that mapping objective function to alpha, from [0,
/// init_obj] -> [min, max]
/// @tparam the index of which alpha to update
template <typename nlp_numerical_type, IndexType alphaIdx>
struct exponential_by_obj {
  static constexpr nlp_numerical_type alphaMax = 1.5;
  static constexpr nlp_numerical_type alphaMin = 0.3;
  static constexpr nlp_numerical_type alphaMin_minus_one = alphaMin - 1;
  static constexpr nlp_numerical_type log_alphaMax_minus_alphaMin_plus_1 =
      std::log(alphaMax - alphaMin_minus_one);
  nlp_numerical_type theConstant =
      0.0; ///< log(alpha_max - alpha_min + 1) / init
};

template <typename nlp_numerical_type, IndexType alphaIdx>
struct alpha_update_trait<exponential_by_obj<nlp_numerical_type, alphaIdx>> {
  typedef exponential_by_obj<nlp_numerical_type, alphaIdx> update_type;

  template <typename nlp_type>
  static constexpr typename nlp_type::nlp_numerical_type obj(nlp_type &nlp) {
    switch (alphaIdx) {
    case 0:
      return nlp._objHpwlRaw;
      break;
    case 1:
      return nlp._objOvlRaw;
      break;
    case 2:
      return nlp._objOobRaw;
      break;
    case 3:
      return nlp._objCrfRaw;
      break;
    default:
      return nlp._objFenceRaw;
      break;
    }
  }

  template <typename nlp_type>
  static constexpr update_type
  construct(nlp_type &, alpha_hpwl_ovl_oob<nlp_numerical_type> &) {
    return update_type();
  }

  template <typename nlp_type>
  static constexpr void init(nlp_type &nlp,
                             alpha_hpwl_ovl_oob<nlp_numerical_type> &alpha,
                             update_type &update) {
    if (obj(nlp) < REAL_TYPE_TOL) {
      update.theConstant = -1;
      alpha._alpha[alphaIdx] = update.alphaMax;
      return;
    }
    update.theConstant = update.log_alphaMax_minus_alphaMin_plus_1 / obj(nlp);
    alpha._alpha[alphaIdx] = update.alphaMax;
  }

  template <typename nlp_type>
  static constexpr void update(nlp_type &nlp,
                               alpha_hpwl_ovl_oob<nlp_numerical_type> &alpha,
                               update_type &update) {
    if (update.theConstant < REAL_TYPE_TOL) {
      return;
    }
    if (obj(nlp) < REAL_TYPE_TOL) {
      alpha._alpha[alphaIdx] = update.alphaMin;
      return;
    }
    alpha._alpha[alphaIdx] =
        std::exp(update.theConstant * obj(nlp)) + update.alphaMin - 1;
#ifdef DEBUG_GR
    DBG("new alpha idx %d %f \n", alphaIdx, alpha._alpha[alphaIdx]);
    DBG("obj %f , the const %f \n", obj(nlp), update.theConstant);
#endif
  }
};

/// @breif update the alpha that mapping objective function to alpha, from [0,
/// init_obj] -> [min, max]
/// @tparam the index of which alpha to update
template <typename nlp_numerical_type, IndexType alphaIdx>
struct reciprocal_by_obj {
  // alpha = a / (x - k * obj_init) + b
  static constexpr nlp_numerical_type alphaMax = 2.0;
  static constexpr nlp_numerical_type alphaMin = 0.2;
  static constexpr nlp_numerical_type k = 100; ///< k > 1.0
  nlp_numerical_type a =
      -1.0; ///< (k ^2 - k) * obj_init * (alphaMax - alphaMin), should > 0
  nlp_numerical_type b =
      1.0; // -k * alphaMax + k * alphaMin + alphaMax; < 0 for most k
  nlp_numerical_type kObjInit = -1.0;
};

template <typename nlp_numerical_type, IndexType alphaIdx>
struct alpha_update_trait<reciprocal_by_obj<nlp_numerical_type, alphaIdx>> {
  typedef reciprocal_by_obj<nlp_numerical_type, alphaIdx> update_type;

  template <typename nlp_type>
  static constexpr typename nlp_type::nlp_numerical_type obj(nlp_type &nlp) {
    switch (alphaIdx) {
    case 0:
      return nlp._objHpwlRaw;
      break;
    case 1:
      return nlp._objOvlRaw;
      break;
    case 2:
      return nlp._objOobRaw;
      break;
    case 3:
      return nlp._objCrfRaw;
      break;
    default:
      return nlp._objFenceRaw;
      break;
    }
  }

  template <typename nlp_type>
  static constexpr update_type
  construct(nlp_type &, alpha_hpwl_ovl_oob<nlp_numerical_type> &) {
    return update_type();
  }

  template <typename nlp_type>
  static constexpr void init(nlp_type &nlp,
                             alpha_hpwl_ovl_oob<nlp_numerical_type> &alpha,
                             update_type &update) {
    if (obj(nlp) < REAL_TYPE_TOL) {
      alpha._alpha[alphaIdx] = update.alphaMax;
      return;
    }
    auto objInit = obj(nlp);
    alpha._alpha[alphaIdx] = update.alphaMax;
    update.a = (update.k * update.k - update.k) * objInit *
               (update.alphaMax - update.alphaMin);
    update.b = -update.k * update.alphaMax + update.k * update.alphaMin +
               update.alphaMax;
    update.kObjInit = update.k * objInit;
  }

  template <typename nlp_type>
  static constexpr void update(nlp_type &nlp,
                               alpha_hpwl_ovl_oob<nlp_numerical_type> &alpha,
                               update_type &update) {
    if (update.kObjInit < REAL_TYPE_TOL) {
      init(nlp, alpha, update);
      return;
    }
    if (obj(nlp) < REAL_TYPE_TOL) {
      alpha._alpha[alphaIdx] = update.alphaMin;
      return;
    }
    if (update.kObjInit < update.k * obj(nlp)) {
      init(nlp, alpha, update);
      return;
    }
    alpha._alpha[alphaIdx] =
        -update.a / (obj(nlp) - update.kObjInit) + update.b;
#ifdef DEBUG_GR
    DBG("update alpha:: new alpha idx %d %f \n", alphaIdx,
        alpha._alpha[alphaIdx]);
    DBG("update alpha:: a %f obj %f kObjInit %f b %f \n", update.a, obj(nlp),
        update.kObjInit, update.b);
#endif
  }
};

/// @breif update the alpha that mapping objective function to alpha, from [0,
/// init_obj] -> [min, max]
/// @tparam the index of which alpha to update
template <typename nlp_numerical_type, IndexType alphaIdx>
struct linear_by_obj {
  // alpha = a / (x - k * obj_init) + b
  static constexpr nlp_numerical_type alphaMax = 2.0;
  static constexpr nlp_numerical_type alphaMin = 0.3;
  nlp_numerical_type objInit = -1.0;
};

template <typename nlp_numerical_type, IndexType alphaIdx>
struct alpha_update_trait<linear_by_obj<nlp_numerical_type, alphaIdx>> {
  typedef linear_by_obj<nlp_numerical_type, alphaIdx> update_type;

  template <typename nlp_type>
  static constexpr typename nlp_type::nlp_numerical_type obj(nlp_type &nlp) {
    switch (alphaIdx) {
    case 0:
      return nlp._objHpwlRaw;
      break;
    case 1:
      return nlp._objOvlRaw;
      break;
    case 2:
      return nlp._objOobRaw;
      break;
    case 3:
      return nlp._objCrfRaw;
      break;
    default:
      return nlp._objFenceRaw;
      break;
    }
  }

  template <typename nlp_type>
  static constexpr update_type
  construct(nlp_type &, alpha_hpwl_ovl_oob<nlp_numerical_type> &) {
    return update_type();
  }

  template <typename nlp_type>
  static constexpr void init(nlp_type &nlp,
                             alpha_hpwl_ovl_oob<nlp_numerical_type> &alpha,
                             update_type &update) {
    if (obj(nlp) < REAL_TYPE_TOL) {
      alpha._alpha[alphaIdx] = update.alphaMax;
      return;
    }
    auto objInit = obj(nlp);
    alpha._alpha[alphaIdx] = update.alphaMax;
    update.objInit = objInit;
  }

  template <typename nlp_type>
  static constexpr void update(nlp_type &nlp,
                               alpha_hpwl_ovl_oob<nlp_numerical_type> &alpha,
                               update_type &update) {
    auto objCurr = obj(nlp);
    if (update.objInit < objCurr) {
      init(nlp, alpha, update);
      return;
    }
    if (update.objInit < REAL_TYPE_TOL) {
      init(nlp, alpha, update);
      return;
    }
    if (objCurr < REAL_TYPE_TOL) {
      alpha._alpha[alphaIdx] = update.alphaMin;
      return;
    }
    alpha._alpha[alphaIdx] =
        ((update.alphaMax - update.alphaMin) / update.objInit) * objCurr +
        update.alphaMin;
#ifdef DEBUG_GR
    DBG("update alpha:: new alpha idx %d %f \n", alphaIdx,
        alpha._alpha[alphaIdx]);
#endif
  }
};

/// @brief a convenient wrapper for combining different types of stop_condition
/// condition. the list in the template will be check one by one and return
/// converge if any of them say so
template <typename alpha_update_type, typename... others>
struct alpha_update_list {
  typedef alpha_update_list<others...> base_type;
  alpha_update_type _update;
  alpha_update_list<others...> _list;
};

template <typename alpha_update_type>
struct alpha_update_list<alpha_update_type> {
  alpha_update_type _update;
};

template <typename alpha_update_type, typename... others>
struct alpha_update_trait<alpha_update_list<alpha_update_type, others...>> {
  typedef alpha_update_list<alpha_update_type, others...> list_type;
  typedef typename alpha_update_list<alpha_update_type, others...>::base_type
      base_type;

  template <typename nlp_type, typename alpha_type>
  static list_type construct(nlp_type &n, alpha_type &a) {
    list_type list;
    list._update =
        std::move(alpha_update_trait<alpha_update_type>::construct(n, a));
    list._list = std::move(alpha_update_trait<base_type>::construct(n, a));
    return list;
  }

  template <typename nlp_type, typename alpha_type>
  static void init(nlp_type &n, alpha_type &a, list_type &c) {
    alpha_update_trait<alpha_update_type>::init(n, a, c._update);
    alpha_update_trait<base_type>::init(n, a, c._list);
  }

  template <typename nlp_type, typename alpha_type>
  static void update(nlp_type &n, alpha_type &alpha, list_type &c) {
    alpha_update_trait<alpha_update_type>::update(n, alpha, c._update);
    alpha_update_trait<base_type>::update(n, alpha, c._list);
  }
};

template <typename alpha_update_type>
struct alpha_update_trait<alpha_update_list<alpha_update_type>> {
  typedef alpha_update_list<alpha_update_type> list_type;

  template <typename nlp_type, typename alpha_type>
  static void init(nlp_type &n, alpha_type &a, list_type &c) {
    alpha_update_trait<alpha_update_type>::init(n, a, c._update);
  }

  template <typename nlp_type, typename alpha_type>
  static list_type construct(nlp_type &n, alpha_type &alpha) {
    list_type list;
    list._update =
        std::move(alpha_update_trait<alpha_update_type>::construct(n, alpha));
    return list;
  }

  template <typename nlp_type, typename alpha_type>
  static void update(nlp_type &n, alpha_type &alpha, list_type &c) {
    alpha_update_trait<alpha_update_type>::update(n, alpha, c._update);
  }
};

} // namespace update
} // namespace alpha
} // namespace nlp
PROJECT_NAMESPACE_END
