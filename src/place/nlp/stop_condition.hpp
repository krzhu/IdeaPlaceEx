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

namespace outer_stop_condition {
/* Stop outer problem condition */
template <typename T> struct stop_condition_trait {
  // static T construct(NlpType &)
  // static IntType stopPlaceCondition(T&, NlpType &)
};

/// @brief stop condition with number of iterations
template <IntType MaxIter = 10> struct stop_after_num_outer_iterations {
  static constexpr IntType maxIter = MaxIter;
  IntType curIter = 0;
};

template <IntType MaxIter>
struct stop_condition_trait<stop_after_num_outer_iterations<MaxIter>> {
  template <typename NlpType>
  static stop_after_num_outer_iterations<MaxIter> construct(NlpType &) {
    return stop_after_num_outer_iterations<MaxIter>();
  }

  template <typename NlpType>
  static void init(NlpType &, stop_after_num_outer_iterations<MaxIter> &stop) {
    stop.curIter = 0;
  }

  static void clear(stop_after_num_outer_iterations<MaxIter> &stop) {
    stop.curIter = 0;
  }

  template <typename NlpType>
  static BoolType
  stopPlaceCondition(NlpType &,
                     stop_after_num_outer_iterations<MaxIter> &stop) {
    if (stop.curIter >= stop.maxIter) {
      stop.curIter = 0;
      return 1;
    }
    ++stop.curIter;
    return 0;
  }
};

/// @brief stop condition that is enable only when the placer is in fast mode
/// @tparam the slave stop condition. The slave is checked when in fast mode
template <typename stop_slave_type> struct stop_enable_if_fast_mode {
  bool activate = false;
  stop_slave_type slave;
};

template <typename stop_slave_type>
struct stop_condition_trait<stop_enable_if_fast_mode<stop_slave_type>> {
  typedef stop_enable_if_fast_mode<stop_slave_type> stop_type;
  template <typename NlpType> static stop_type construct(NlpType &nlp) {
    stop_type stop;
    stop.slave = stop_condition_trait<stop_slave_type>::construct(nlp);
    stop.activate = nlp._db.parameters().isFastMode();
    return std::move(stop);
  }

  template <typename NlpType> static void init(NlpType &n, stop_type &stop) {
    if (stop.activate) {
      stop_condition_trait<stop_slave_type>::init(n, stop.slave);
    }
  }

  static void clear(stop_type &stop) {
    if (stop.activate) {
      stop_condition_trait<stop_slave_type>::clear(stop.slave);
    }
  }

  template <typename NlpType>
  static BoolType stopPlaceCondition(NlpType &n, stop_type &stop) {
    if (not stop.activate)
      return 0;
    return stop_condition_trait<stop_slave_type>::stopPlaceCondition(
        n, stop.slave);
  }
};
/// @brief stop after violating is small enough
struct stop_after_violate_small {
  static constexpr RealType overlapRatio =
      0.03; ///< with respect to total cell area
  static constexpr RealType outOfBoundaryRatio =
      0.1; ///< with respect to boundary
  static constexpr RealType asymRatio =
      0.02; ///< with respect to sqrt(total cell area)
  static constexpr RealType outFenceRatio =
      1; ///< With respect to total cell area needs well
  IntType curIter = 0;
  static constexpr IntType minIter = 0;
};

template <> struct stop_condition_trait<stop_after_violate_small> {
  typedef stop_after_violate_small stop_type;

  template <typename NlpType> static stop_type construct(NlpType &) {
    return stop_type();
  }

  template <typename NlpType> static void init(NlpType &, stop_type &s) {
    s.curIter = 0;
  }

  static void clear(stop_type &s) { s.curIter = 0; }

  template <typename NlpType>
  static BoolType stopPlaceCondition(NlpType &n, stop_type &stop) {
    using CoordType = typename NlpType::nlp_coordinate_type;


    ++stop.curIter;
    if (stop.curIter < stop.minIter) {
      return false;
    }

    // check whether overlapping is small than threshold
    CoordType ovlArea = 0;
    const CoordType ovlThreshold = stop.overlapRatio * n._totalCellArea;
    for (auto &op : n._ovlOps) {
      ovlArea += diff::place_overlap_trait<
          typename NlpType::nlp_ovl_type>::overlapArea(op);
    }
    if (ovlArea > ovlThreshold) {
      return false;
    }
    // Check whether out of boundary is smaller than threshold
    CoordType oobArea = 0;
    const CoordType oobThreshold = stop.outOfBoundaryRatio * n._boundary.area();
    for (auto &op : n._oobOps) {
      oobArea += diff::place_out_of_boundary_trait<
          typename NlpType::nlp_oob_type>::oobArea(op);
      if (oobArea > oobThreshold) {
#ifdef DEBUG_GR
        DBG("fail on oob \n");
#endif
        return false;
      }
    }
    // Check whether asymmetry distance is smaller than threshold
    CoordType asymDist = 0;
    const CoordType asymThreshold =
        stop.asymRatio * std::sqrt(n._totalCellArea);
    for (auto &op : n._asymOps) {
      asymDist += diff::place_asym_trait<
          typename NlpType::nlp_asym_type>::asymDistanceNormalized(op);
      if (asymDist > asymThreshold) {
#ifdef DEBUG_GR
        DBG("fail on asym \n");
#endif
        return false;
      }
    }
#if 0
    // Check whether out of fence region area is smaller than threshold
    CoordType outFenceArea = 0.0;
    CoordType inWellCellArea = 0.0;
    for (auto &op : n._fenceOps) {
      outFenceArea += diff::place_fence_trait<
          typename NlpType::nlp_fence_type>::outFenceRegionArea(op);
      inWellCellArea += op._cellWidth * op._cellHeight;
    }
    if (outFenceArea / inWellCellArea > stop.outFenceRatio) {
#ifdef DEBUG_GR
      DBG("fail on fence \n");
#endif
      return false;
    }
#endif

#ifdef DEBUG_GR
    DBG("ovl area %f target %f \n oob area %f target %f \n asym dist %f target "
        "%f \n",
        ovlArea, ovlThreshold, oobArea, oobThreshold, asymDist, asymThreshold);
#endif
    return true;
  }
};
/// @brief a convenient wrapper for combining different types of stop_condition
/// condition. the list in the template will be check one by one and return
/// converge if any of them say so
template <typename stop_condition_type, typename... others>
struct stop_condition_list {
  typedef stop_condition_list<others...> base_type;
  stop_condition_type _stop;
  stop_condition_list<others...> _list;
};

template <typename stop_condition_type>
struct stop_condition_list<stop_condition_type> {
  stop_condition_type _stop;
};

template <typename stop_condition_type, typename... others>
struct stop_condition_trait<
    stop_condition_list<stop_condition_type, others...>> {
  typedef stop_condition_list<stop_condition_type, others...> list_type;
  typedef
      typename stop_condition_list<stop_condition_type, others...>::base_type
          base_type;

  static void clear(list_type &c) {
    stop_condition_trait<stop_condition_type>::clear(c._stop);
    stop_condition_trait<base_type>::clear(c._list);
  }

  template <typename nlp_type> static list_type construct(nlp_type &n) {
    list_type list;
    list._stop =
        std::move(stop_condition_trait<stop_condition_type>::construct(n));
    list._list = std::move(stop_condition_trait<base_type>::construct(n));
    return list;
  }

  template <typename nlp_type>
  static BoolType stopPlaceCondition(nlp_type &n, list_type &c) {
    BoolType stop = false;
    if (stop_condition_trait<stop_condition_type>::stopPlaceCondition(
            n, c._stop)) {
      stop = true;
    }
    if (stop_condition_trait<base_type>::stopPlaceCondition(n, c._list)) {
      stop = true;
    }
    if (stop) {
      stop_condition_trait<stop_condition_type>::clear(c._stop);
    }
    return stop;
  }
};

template <typename stop_condition_type>
struct stop_condition_trait<stop_condition_list<stop_condition_type>> {
  typedef stop_condition_list<stop_condition_type> list_type;
  static void clear(stop_condition_list<stop_condition_type> &c) {
    stop_condition_trait<stop_condition_type>::clear(c._stop);
  }

  template <typename nlp_type> static list_type construct(nlp_type &n) {
    list_type list;
    list._stop =
        std::move(stop_condition_trait<stop_condition_type>::construct(n));
    return list;
  }

  template <typename nlp_type>
  static BoolType
  stopPlaceCondition(nlp_type &n, stop_condition_list<stop_condition_type> &c) {
    if (stop_condition_trait<stop_condition_type>::stopPlaceCondition(
            n, c._stop)) {
      stop_condition_trait<stop_condition_type>::clear(c._stop);
      return true;
    }
    return false;
  }
};

} // namespace outer_stop_condition
} // namespace nlp

PROJECT_NAMESPACE_END
