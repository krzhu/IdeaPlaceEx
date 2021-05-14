/**
 * @file stop_condition.hpp
 * @brief The non-lnear programming outer-problem optimization
 * @author Keren Zhu
 * @date 04/09/2020
 */

#pragma once

#include "global/global.h"
#include "nlpTypes.hpp"
#include "place/different.h"

PROJECT_NAMESPACE_BEGIN

namespace nlp {


namespace outer_multiplier {
/// @brief the multiplier is categoried by types
template <typename mult_type>
struct is_mult_type_dependent_diff : std::false_type {};

namespace init {
template <typename T> struct multiplier_init_trait {};

struct hard_code_init {};

template <> struct multiplier_init_trait<hard_code_init> {
  template <typename nlp_type, typename mult_type>
  static void init(nlp_type &, mult_type &mult) {
    mult._constMults.at(0) = 16; // hpwl
    mult._constMults.at(1) = 1;  // cos
    mult._constMults.at(2) = 16; // power
    mult._variedMults.at(0) = 1; // ovl
    mult._variedMults.at(1) = 1; // oob
    mult._variedMults.at(2) = 1; // asym
    mult._variedMults.at(3) = 1; // crf
    mult._variedMults.at(4) = 1; // Fence
  }
};

/// @brief match by gradient norm
struct init_by_matching_gradient_norm {
  static constexpr RealType penaltyRatioToObj =
      0.02; ///< The ratio of targeting penalty
  static constexpr RealType small = 0.01;
};

template <> struct multiplier_init_trait<init_by_matching_gradient_norm> {
  typedef init_by_matching_gradient_norm init_type;
  template <typename nlp_type, typename mult_type,
            std::enable_if_t<nlp::is_first_order_diff<nlp_type>::value, void>
                * = nullptr>
  static void init(nlp_type &nlp, mult_type &mult) {
    RealType totalHpwlWeights = 0.0;
    RealType totalCosWeights = 0.0;
    RealType totalPowerWlWeights = 0.0;
    RealType totalCrfWeights = 0.0;
    RealType totalFenceWeights = 0.0;
    for (const auto &op : nlp._hpwlOps) {
      totalHpwlWeights += op._weight;
    }
    for (const auto &op : nlp._cosOps) {
      totalCosWeights += op._weight;
    }
    for (const auto &op : nlp._powerWlOps) {
      totalPowerWlWeights += op._weight;
    }
    for (const auto &op : nlp._crfOps) {
      totalCrfWeights += op._weight;
    }
    for (const auto &op : nlp._fenceOps) {
      totalFenceWeights += op._weight;
    }
    mult._constMults.at(0) = 1.0; // hpwl
    const auto hpwlMult = mult._constMults.at(0);
    const auto hpwlNorm = nlp._gradHpwl.norm();
    const auto hpwlMultNorm = hpwlMult * hpwlNorm;
    const auto hpwlMultNormPenaltyRatio =
        hpwlMultNorm * init_type::penaltyRatioToObj;
    const auto cosNorm = nlp._gradCos.norm();
    const auto powerWlNorm = nlp._gradPowerWl.norm();
    const auto ovlNorm = nlp._gradOvl.norm();
    const auto oobNorm = nlp._gradOob.norm();
    const auto asymNorm = nlp._gradAsym.norm();
    const auto crfNorm = nlp._gradCrf.norm();
    const auto fenceNorm = nlp._gradFence.norm();
    const auto maxPenaltyNorm = ovlNorm;
    // Make a threshold on by referencing hpwl to determine whether one is small
    const auto small = init_type::small * hpwlNorm;

    // Fix corner case that may happen when the placement is very small
    if (hpwlNorm < REAL_TYPE_TOL) {
      mult._constMults.resize(3, 1);
      mult._variedMults.resize(5, 1);
      WRN("Ideaplace: NLP global placement: init multipliers: wire length  "
          "gradient norm is very small %f!, ",
          hpwlNorm);
      return;
    }
    // match gradient norm for signal path
    if (cosNorm > small) {
      mult._constMults.at(1) =
          hpwlMultNorm * totalCosWeights / totalHpwlWeights / cosNorm;
    } else {
      mult._constMults.at(1) = hpwlMult * 30;
    }
    // match gradient norm for signal path
    if (powerWlNorm > small) {
      mult._constMults.at(2) =
          hpwlMultNorm * totalCosWeights / totalHpwlWeights / powerWlNorm;
    } else {
      mult._constMults.at(2) = hpwlMult * 30;
    }
    // overlap
    if (ovlNorm > small) {
      mult._variedMults.at(0) = hpwlMultNormPenaltyRatio / ovlNorm;
    } else {
      mult._variedMults.at(0) = hpwlMultNormPenaltyRatio / maxPenaltyNorm;
    }
    // out of boundary
    // Since we know oob is small at beginning, and it will take effect after a
    // few iterations. Therefore it would be better to set it to resonable range
    // first
    // mult._variedMults.at(1) = hpwlMultNormPenaltyRatio;
    if (oobNorm > small) {
      mult._variedMults.at(1) = hpwlMultNormPenaltyRatio / oobNorm;
    } else {
      mult._variedMults.at(1) = hpwlMultNormPenaltyRatio / maxPenaltyNorm;
    }
    // asym
    if (asymNorm > small) {
      mult._variedMults.at(2) = hpwlMultNormPenaltyRatio / asymNorm;
    } else {
      mult._variedMults.at(2) = hpwlMultNormPenaltyRatio / maxPenaltyNorm;
    }
    // crf
    if (crfNorm > small) {
      mult._variedMults.at(3) = hpwlMultNormPenaltyRatio * totalCrfWeights /
                                totalHpwlWeights / crfNorm;
    } else {
      mult._variedMults.at(3) = hpwlMultNormPenaltyRatio / maxPenaltyNorm;
    }
    // fence
    if (fenceNorm > small) {
      mult._variedMults.at(4) = hpwlMultNormPenaltyRatio * totalFenceWeights /
                                totalHpwlWeights / fenceNorm;
    } else {
      mult._variedMults.at(4) = hpwlMultNormPenaltyRatio / maxPenaltyNorm;
    }
#ifdef DEBUG_GR
    // crf
    DBG("init mult: hpwl %f cos %f power wl %f \n", mult._constMults[0],
        mult._constMults[1], mult._constMults[2]);
    DBG("init mult: ovl %f oob %f asym %f current flow %f fence %f\n",
        mult._variedMults[0], mult._variedMults[1], mult._variedMults[2],
        mult._variedMults[3], mult._variedMults[4]);
#endif
  }
};
}; // namespace init
namespace update {
template <typename T> struct multiplier_update_trait {};

struct no_update {};
template <> struct multiplier_update_trait<no_update> {
  typedef no_update update_type;

  template <typename nlp_type, typename mult_type,
            std::enable_if_t<is_mult_type_dependent_diff<mult_type>::value,
                             void> * = nullptr>
  static update_type construct(nlp_type &, mult_type &) {
    return update_type();
  }

  template <typename nlp_type, typename mult_type,
            std::enable_if_t<is_mult_type_dependent_diff<mult_type>::value,
                             void> * = nullptr>
  static void init(nlp_type &, mult_type &, update_type &) {}

  template <typename nlp_type, typename mult_type,
            std::enable_if_t<is_mult_type_dependent_diff<mult_type>::value,
                             void> * = nullptr>
  static void update(nlp_type &, mult_type &, update_type &) {}
};

/// @brief increase the total amounts of penalty of by a constant
struct shared_constant_increase_penalty {
  static constexpr RealType penalty = 20;
};

template <> struct multiplier_update_trait<shared_constant_increase_penalty> {
  typedef shared_constant_increase_penalty update_type;

  template <typename nlp_type, typename mult_type,
            std::enable_if_t<is_mult_type_dependent_diff<mult_type>::value,
                             void> * = nullptr>
  static update_type construct(nlp_type &, mult_type &) {
    return update_type();
  }

  template <typename nlp_type, typename mult_type,
            std::enable_if_t<is_mult_type_dependent_diff<mult_type>::value,
                             void> * = nullptr>
  static void update(nlp_type &nlp, mult_type &mult, update_type &update) {
    nlp._wrapObjAllTask.run();
    const auto rawOvl = nlp._objOvl / mult._variedMults.at(0);
    const auto rawOob = nlp._objOob / mult._variedMults.at(1);
    const auto rawAsym = nlp._objAsym / mult._variedMults.at(2);
    const auto rawFence = nlp._objFence / mult._variedMults.at(4);
    const auto fViolate = rawOvl + rawOob + rawAsym + rawFence;
#ifdef DEBUG_GR
    DBG("update mult: raw ovl %f oob %f asym %f total %f \n", rawOvl, rawOob,
        rawAsym, fViolate);
    DBG("update mult:  before ovl %f oob %f asym %f \n", mult._variedMults[0],
        mult._variedMults[1], mult._variedMults[2]);
#endif
    mult._variedMults.at(0) += update.penalty * (rawOvl / fViolate);
    mult._variedMults.at(1) += update.penalty * (rawOob / fViolate);
    mult._variedMults.at(2) += update.penalty * (rawAsym / fViolate);
    mult._variedMults.at(4) += update.penalty * (rawFence / fViolate);
#ifdef DEBUG_GR
    DBG("update mult: afterafter  ovl %f oob %f asym %f \n",
        mult._variedMults[0], mult._variedMults[1], mult._variedMults[2]);
#endif
  }
};

/// @brief direct subgradient
struct direct_subgradient {
  static constexpr RealType stepSize = 1e-3;


};

template <> struct multiplier_update_trait<direct_subgradient> {
  typedef direct_subgradient update_type;

  template <typename nlp_type, typename mult_type,
            std::enable_if_t<is_mult_type_dependent_diff<mult_type>::value,
                             void> * = nullptr>
  static update_type construct(nlp_type &, mult_type &) {
    return update_type();
  }

  template <typename nlp_type, typename mult_type,
            std::enable_if_t<is_mult_type_dependent_diff<mult_type>::value,
                             void> * = nullptr>
  static void init(nlp_type &, mult_type &, update_type &) {}

  template <typename nlp_type, typename mult_type,
            std::enable_if_t<is_mult_type_dependent_diff<mult_type>::value,
                             void> * = nullptr>
  static void update(nlp_type &nlp, mult_type &mult, update_type &update) {
    nlp._wrapObjAllTask.run();
    const auto rawOvl = nlp._objOvl / mult._variedMults.at(0);
    const auto rawOob = nlp._objOob / mult._variedMults.at(1);
    const auto rawAsym = nlp._objAsym / mult._variedMults.at(2);
    const auto rawCos = nlp._objCos / mult._constMults.at(1);
    const auto rawCrf = nlp._objCrf / mult._variedMults.at(3);
    const auto rawFence = nlp._objFence / mult._variedMults.at(4);
#ifdef DEBUG_GR
    DBG("update mult: raw ovl %f oob %f asym %f cos %f powerWl %f current "
        "flow\n",
        rawOvl, rawOob, rawAsym, rawCos, nlp._objPowerWl, rawCrf);
    DBG("update mult:  before ovl %f oob %f asym %f current flow %f cos %f \n",
        mult._variedMults[0], mult._variedMults[1], mult._variedMults[2],
        mult._variedMults[3], mult._constMults[1]);
#endif
    mult._variedMults.at(0) += update.stepSize * (rawOvl);
    mult._variedMults.at(1) += update.stepSize * (rawOob);
    mult._variedMults.at(2) += update.stepSize * (rawAsym);
    mult._variedMults.at(3) += update.stepSize * (rawCrf);
    mult._variedMults.at(4) += update.stepSize * (rawFence);
    mult._constMults.at(1) += update.stepSize * (rawCos);
#ifdef DEBUG_GR
    DBG("update mult: afterafter  ovl %f oob %f asym %f cos %f current flow %f "
        "Fence %f \n",
        mult._variedMults[0], mult._variedMults[1], mult._variedMults[2],
        mult._constMults[1], mult._variedMults[3], mult._variedMults[4]);
#endif
  }
};

template <typename nlp_numerical_type> struct match_grad_const_multipliers {
  nlp_numerical_type totalHpwlWeights = 0.0;
  nlp_numerical_type totalCosWeights = 0.0;
  nlp_numerical_type totalPowerWlWeights = 0.0;
  static constexpr nlp_numerical_type maxMult = 500;
  bool _recordedInit =
      false; ///< Whether the init multipliers have been recorded
  nlp_numerical_type ratio =
      1.0; ///< The ratio of matched part vs constant part
  static constexpr nlp_numerical_type ratioDecayRate =
      0.3; ///< The decay factor of "ratio"

  nlp_numerical_type hpwlInitMult = 1.0;
  nlp_numerical_type cosInitMult = 1.0;
  nlp_numerical_type powerWlInitMult = 1.0;
  void decay() { ratio *= ratioDecayRate; }
};

template <typename nlp_numerical_type>
struct multiplier_update_trait<
    match_grad_const_multipliers<nlp_numerical_type>> {
  typedef match_grad_const_multipliers<nlp_numerical_type> update_type;
  template <typename nlp_type, typename mult_type,
            std::enable_if_t<is_mult_type_dependent_diff<mult_type>::value,
                             void> * = nullptr>
  static update_type construct(nlp_type &, mult_type &) {
    return update_type();
  }

  template <typename nlp_type, typename mult_type,
            std::enable_if_t<is_mult_type_dependent_diff<mult_type>::value,
                             void> * = nullptr>
  static void init(nlp_type &n, mult_type &, update_type &u) {
    for (const auto &op : n._hpwlOps) {
      u.totalHpwlWeights += op._weight;
    }
    for (const auto &op : n._cosOps) {
      u.totalCosWeights += op._weight;
    }
    for (const auto &op : n._powerWlOps) {
      u.totalPowerWlWeights += op._weight;
    }
  }

  template <typename nlp_type, typename mult_type,
            std::enable_if_t<is_mult_type_dependent_diff<mult_type>::value,
                             void> * = nullptr>
  static void update(nlp_type &nlp, mult_type &mult, update_type &u) {
    nlp._wrapObjAllTask.run();
    // Record the init multipliers
    if (not u._recordedInit) {
      u._recordedInit = true;
      u.hpwlInitMult = mult._constMults.at(0);
      u.cosInitMult =
          mult._constMults.at(0) * u.totalCosWeights / u.totalHpwlWeights;
      u.powerWlInitMult =
          mult._constMults.at(0) * u.totalPowerWlWeights / u.totalHpwlWeights;
    }
    const auto hpwlMult = mult._constMults.at(0);
    const auto hpwlNorm = nlp._gradHpwl.norm();
    const auto hpwlMultNorm = hpwlMult * hpwlNorm / u.totalHpwlWeights;
    const auto cosNorm = nlp._gradCos.norm() / mult._constMults.at(1);
    const auto powerWlNorm = nlp._gradPowerWl.norm() / mult._constMults.at(2);
    // Make a threshold on by referencing hpwl to determine whether one is small
    const auto small = 0.001 * hpwlNorm;

    // Fix corner case that may happen when the placement is very small
    if (hpwlNorm < REAL_TYPE_TOL) {
      return;
    }
    nlp_numerical_type cosMatch;
    // match gradient norm for signal path
    if (cosNorm > small) {
      cosMatch = std::min(hpwlMultNorm * u.totalCosWeights / cosNorm,
                          update_type::maxMult);
    } else {
      cosMatch = hpwlMult;
    }
    // match gradient norm for signal path
    nlp_numerical_type powerWlMatch;
    if (powerWlNorm > small) {
      powerWlMatch =
          std::min(hpwlMultNorm * u.totalPowerWlWeights / powerWlNorm,
                   update_type::maxMult);
    } else {
      powerWlMatch = hpwlMult;
    }
    mult._constMults.at(1) = u.ratio * cosMatch + (1 - u.ratio) * u.cosInitMult;
    mult._constMults.at(2) =
        u.ratio * powerWlMatch + (1 - u.ratio) * u.powerWlInitMult;
    u.decay();
#ifdef DEBUG_GR
    DBG("match_grad_const_multipliers: multipliers hpwl %f cos %f power wl %f "
        "\n",
        mult._constMults.at(0), mult._constMults.at(1), mult._constMults.at(2));
#endif
  }
};

/// @breif subgradient normalize by values in iter 0
template <typename nlp_numerical_type> struct subgradient_normalized_by_init {
  static constexpr nlp_numerical_type stepSize = 10;
  std::vector<nlp_numerical_type> normalizeFactor;
};

template <typename nlp_numerical_type>
struct multiplier_update_trait<
    subgradient_normalized_by_init<nlp_numerical_type>> {
  typedef subgradient_normalized_by_init<nlp_numerical_type> update_type;

  template <typename nlp_type, typename mult_type,
            std::enable_if_t<is_mult_type_dependent_diff<mult_type>::value,
                             void> * = nullptr>
  static update_type construct(nlp_type &, mult_type &) {
    return update_type();
  }

  template <typename nlp_type, typename mult_type,
            std::enable_if_t<is_mult_type_dependent_diff<mult_type>::value,
                             void> * = nullptr>
  static void init(nlp_type &nlp, mult_type &mult, update_type &update) {
    update.normalizeFactor.resize(3);
    update.normalizeFactor.at(0) = mult._variedMults.at(0) / nlp._objOvl;
    update.normalizeFactor.at(1) = 1; // mult._variedMults.at(1) / nlp._objOob;
    update.normalizeFactor.at(2) = mult._variedMults.at(2) / nlp._objAsym;
    update.normalizeFactor.at(3) = mult._variedMults.at(3) / nlp._objCrf;
    update.normalizeFactor.at(4) = mult._variedMults.at(4) / nlp._objFence;
    // update.normalizeFactor.at(0) = 1 / nlp._objOvl;
    // update.normalizeFactor.at(1) = 1;// / nlp._objOob;
    // update.normalizeFactor.at(2) = 1 / nlp._objAsym;
  }

  template <typename nlp_type, typename mult_type,
            std::enable_if_t<is_mult_type_dependent_diff<mult_type>::value,
                             void> * = nullptr>
  static void update(nlp_type &nlp, mult_type &mult, update_type &update) {
    nlp._wrapObjAllTask.run();
    const auto rawOvl = nlp._objOvl / mult._variedMults.at(0);
    const auto rawOob = nlp._objOob / mult._variedMults.at(1);
    const auto rawAsym = nlp._objAsym / mult._variedMults.at(2);
    const auto rawCrf = nlp._objCrf / mult._variedMults.at(3);
    const auto rawFence = nlp._objFence / mult._variedMults.at(4);
    const auto normalizedOvl = rawOvl * update.normalizeFactor.at(0);
    const auto normalizedOob = rawOob * update.normalizeFactor.at(1);
    const auto normalizedAsym = rawAsym * update.normalizeFactor.at(2);
    const auto normalizedCrf = rawCrf * update.normalizeFactor.at(3);
    const auto normalizedFence = rawFence * update.normalizeFactor.at(4);
#ifdef DEBUG_GR
    DBG("update mult: raw ovl %f oob %f asym %f total %f \n", normalizedOvl,
        normalizedOob, normalizedAsym);
    DBG("update mult:  before ovl %f oob %f asym %f \n", mult._variedMults[0],
        mult._variedMults[1], mult._variedMults[2]);
#endif
    mult._variedMults.at(0) += update.stepSize * (normalizedOvl);
    mult._variedMults.at(1) += update.stepSize * (normalizedOob);
    mult._variedMults.at(2) += update.stepSize * (normalizedAsym);
    mult._variedMults.at(3) += update.stepSize * (normalizedCrf);
    mult._variedMults.at(4) += update.stepSize * (normalizedFence);
#ifdef DEBUG_GR
    DBG("update mult: afterafter  ovl %f oob %f asym %f \n",
        mult._variedMults[0], mult._variedMults[1], mult._variedMults[2]);
#endif
  }
};

} // namespace update
template <typename T> struct multiplier_trait {};

template <typename nlp_numerical_type, typename init_type, typename update_type>
struct mult_const_hpwl_cos_and_penalty_by_type {
  typedef init_type mult_init_type;
  typedef update_type mult_update_type;
  std::vector<nlp_numerical_type> _constMults;  ///< constant mults
  std::vector<nlp_numerical_type> _variedMults; ///< varied penalty multipliers
  update_type update;
};

template <typename nlp_numerical_type, typename init_type, typename update_type>
struct is_mult_type_dependent_diff<mult_const_hpwl_cos_and_penalty_by_type<
    nlp_numerical_type, init_type, update_type>> : std::true_type {};

template <typename nlp_numerical_type, typename init_type, typename update_type>
struct multiplier_trait<mult_const_hpwl_cos_and_penalty_by_type<
    nlp_numerical_type, init_type, update_type>> {
  typedef mult_const_hpwl_cos_and_penalty_by_type<nlp_numerical_type, init_type,
                                                  update_type>
      mult_type;

  template <typename nlp_type> static mult_type construct(nlp_type &nlp) {
    mult_type mult;
    mult._constMults.resize(3, 0.0);
    mult._variedMults.resize(5, 0.0);
    mult.update =
        update::multiplier_update_trait<update_type>::construct(nlp, mult);
    return mult;
  }

  template <typename nlp_type>
  static void init(nlp_type &nlp, mult_type &mult) {
    init::multiplier_init_trait<init_type>::init(nlp, mult);
    update::multiplier_update_trait<update_type>::init(nlp, mult, mult.update);
    for (auto &op : nlp._hpwlOps) {
      op._getLambdaFunc = [&]() { return mult._constMults[0]; };
    }
    for (auto &op : nlp._cosOps) {
      op._getLambdaFunc = [&]() { return mult._constMults[1]; };
    }
    for (auto &op : nlp._ovlOps) {
      op._getLambdaFunc = [&]() { return mult._variedMults[0]; };
    }
    for (auto &op : nlp._oobOps) {
      op._getLambdaFunc = [&]() { return mult._variedMults[1]; };
    }
    for (auto &op : nlp._asymOps) {
      op._getLambdaFunc = [&]() { return mult._variedMults[2]; };
    }
    for (auto &op : nlp._powerWlOps) {
      op._getLambdaFunc = [&]() { return mult._constMults[2]; };
    }
    for (auto &op : nlp._crfOps) {
      op._getLambdaFunc = [&]() { return mult._variedMults[3]; };
    }
    for (auto &op : nlp._fenceOps) {
      op._getLambdaFunc = [&]() { return mult._constMults[0]; }; // SAME as HPWL
    }
    for (auto &op : nlp._areaOps) {
      op._getLambdaFunc = [&]() { return mult._constMults[0]; }; // Area use the same lambda as wirelength
    }
  }

  template <typename nlp_type>
  static void update(nlp_type &nlp, mult_type &mult) {
    update::multiplier_update_trait<update_type>::update(nlp, mult,
                                                         mult.update);
  }

  template <typename nlp_type>
  static void recordRaw(nlp_type &nlp, mult_type &mult) {
    nlp._objHpwlRaw = nlp._objHpwl / mult._constMults[0];
    nlp._objCosRaw = nlp._objCos / mult._constMults[1];
    nlp._objOvlRaw = nlp._objOvl / mult._variedMults[0];
    nlp._objOobRaw = nlp._objOob / mult._variedMults[1];
    nlp._objAsymRaw = nlp._objAsym / mult._variedMults[2];
    nlp._objCrfRaw = nlp._objCrf / mult._variedMults[3];
    nlp._objFenceRaw = nlp._objFence / mult._variedMults[4];
  }

  static std::function<nlp_numerical_type(void)>
  fenceGetLambdaFunc(mult_type &mult) {
    return [&]() { return mult._constMults[0]; };
  }

  static std::function<nlp_numerical_type(void)>
  ovlGetLambdaFunc(mult_type &mult) {
    return [&]() { return mult._variedMults[0]; };
  }
};
} // namespace outer_multiplier
} // namespace nlp

PROJECT_NAMESPACE_END
