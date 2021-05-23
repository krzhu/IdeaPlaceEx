/**
 * @file NlpGPlacer.h
 * @brief The global placement solver with non-linear optimization
 * @author Keren Zhu
 * @date 03/29/2020
 */

#ifndef IDEAPLACE_NLPGPLACER_H_
#define IDEAPLACE_NLPGPLACER_H_

#include <chrono>
#include <Eigen/Dense>
#ifdef IDEAPLACE_TASKFLOR_FOR_GRAD_OBJ_
#include <taskflow/taskflow.hpp>
#endif // IDEAPLACE_TASKFLOR_FOR_GRAD_OBJ
#include "db/Database.h"
#include "pinassign/VirtualPinAssigner.h"
#include "place/different.h"
#include "place/differentSecondOrder.hpp"
#include "place/nlp/nlpFirstOrderKernel.hpp"
#include "place/nlp/nlpInitPlace.hpp"
#include "place/nlp/nlpOptmKernels.hpp"
#include "place/nlp/nlpOuterOptm.hpp"
#include "place/nlp/nlpSecondOrderKernels.hpp"
#include "place/nlp/nlpTasks.hpp"
#include "place/nlp/nlpTypes.hpp"
#include "util/Polygon2Rect.h"
PROJECT_NAMESPACE_BEGIN

namespace nlp {
/* The wrapper of settings */

struct nlp_default_hyperparamters {};

struct nlp_default_types {
  typedef RealType nlp_coordinate_type;
  typedef RealType nlp_numerical_type;
  typedef Eigen::Matrix<nlp_numerical_type, Eigen::Dynamic, Eigen::Dynamic>
      EigenMatrix;
  typedef Eigen::Matrix<nlp_numerical_type, Eigen::Dynamic, 1> EigenVector;
  typedef Eigen::Map<EigenVector> EigenMap;
  typedef Eigen::DiagonalMatrix<nlp_numerical_type, Eigen::Dynamic>
      EigenDiagonalMatrix;
  typedef diff::LseHpwlDifferentiable<nlp_numerical_type, nlp_coordinate_type>
      nlp_hpwl_type;
  typedef diff::CellPairOverlapPenaltyDifferentiable<nlp_numerical_type,
                                                     nlp_coordinate_type>
      nlp_ovl_type;
  typedef diff::CellOutOfBoundaryPenaltyDifferentiable<nlp_numerical_type,
                                                       nlp_coordinate_type>
      nlp_oob_type;
  typedef diff::AsymmetryDifferentiable<nlp_numerical_type, nlp_coordinate_type>
      nlp_asym_type;
  typedef diff::CosineDatapathDifferentiable<nlp_numerical_type,
                                             nlp_coordinate_type>
      nlp_cos_type;
  typedef diff::PowerVerQuadraticWireLengthDifferentiable<nlp_numerical_type,
                                                          nlp_coordinate_type>
      nlp_power_wl_type;
  typedef diff::CurrentFlowDifferentiable<nlp_numerical_type,
                                          nlp_coordinate_type>
      nlp_crf_type;
  typedef diff::FenceBivariateGaussianDifferentiable<nlp_numerical_type,
                                                           nlp_coordinate_type>
      nlp_fence_type;
  typedef diff::LseAreaDifferentiable<nlp_numerical_type, nlp_coordinate_type> 
    nlp_area_type;
};

struct nlp_default_zero_order_algorithms {
  typedef outer_stop_condition::stop_condition_list<
      outer_stop_condition::stop_after_violate_small>
      //outer_stop_condition::stop_after_num_outer_iterations<300>,
      //outer_stop_condition::stop_enable_if_fast_mode<
      //    outer_stop_condition::stop_after_num_outer_iterations<50>>>
      stop_condition_type;
  typedef init_place::init_random_placement_with_normal_distribution_near_center
      init_place_type;

  /* multipliers */
  typedef outer_multiplier::init::hard_code_init mult_init_type;
  typedef outer_multiplier::update::subgradient_normalized_by_init<
      nlp_default_types::nlp_numerical_type>
      mult_update_type;
  // typedef outer_multiplier::update::direct_subgradient mult_update_type;
  typedef outer_multiplier::mult_const_hpwl_cos_and_penalty_by_type<
      nlp_default_types::nlp_numerical_type, mult_init_type, mult_update_type>
      mult_type;
};

struct nlp_default_first_order_algorithms {
  typedef converge::converge_list<
      converge::converge_grad_norm_by_init<
          nlp_default_types::nlp_numerical_type>,
      converge::converge_criteria_max_iter<30000>,
      converge::converge_criteria_enable_if_fast_mode<
          converge::converge_criteria_max_iter<200>>,
      converge::converge_criteria_enable_if_exceed_iter<
          10, converge::converge_criteria_max_iter<10000>>,
      converge::converge_criteria_stop_when_large_variable_changes<nlp_default_types::nlp_numerical_type, nlp_default_types::EigenVector>>
      converge_type;
   //typedef optm::first_order::naive_gradient_descent<converge_type, nlp_default_types::nlp_numerical_type> optm_type;
  typedef optm::first_order::adam<converge_type, nlp_default_types::nlp_numerical_type> optm_type;
  // typedef optm::first_order::nesterov<converge_type,
  // nlp_default_types::nlp_numerical_type> optm_type;

  /* multipliers */
  typedef outer_multiplier::init::init_by_matching_gradient_norm mult_init_type;
  // typedef
  // outer_multiplier::update::subgradient_normalized_by_init<nlp_default_types::nlp_numerical_type>
  // mult_update_type;
  typedef outer_multiplier::update::direct_subgradient mult_update_type;
  typedef outer_multiplier::mult_const_hpwl_cos_and_penalty_by_type<
      nlp_default_types::nlp_numerical_type, mult_init_type, mult_update_type>
      mult_type;

  typedef outer_multiplier::update::match_grad_const_multipliers<
      nlp_default_types::nlp_numerical_type>
      mult_adjust_type;

  /* alpha */
  typedef alpha::alpha_hpwl_ovl_oob<nlp_default_types::nlp_numerical_type>
      alpha_type;
  typedef alpha::update::alpha_update_list<
      alpha::update::reciprocal_by_obj<nlp_default_types::nlp_numerical_type,
                                       1>,
      alpha::update::reciprocal_by_obj<nlp_default_types::nlp_numerical_type,
                                       2>,
      alpha::update::reciprocal_by_obj<nlp_default_types::nlp_numerical_type,
                                       3>,
      alpha::update::fence_update<nlp_default_types::nlp_numerical_type>
     >
      alpha_update_type;
};

template <typename nlp_types> struct nlp_default_second_order_settings {
  typedef diff::jacobi_hessian_approx_trait<typename nlp_types::nlp_hpwl_type>
      hpwl_hessian_trait;
  typedef diff::jacobi_hessian_approx_trait<typename nlp_types::nlp_ovl_type>
      ovl_hessian_trait;
  typedef diff::jacobi_hessian_approx_trait<typename nlp_types::nlp_oob_type>
      oob_hessian_trait;
  typedef diff::jacobi_hessian_approx_trait<typename nlp_types::nlp_asym_type>
      asym_hessian_trait;
  typedef diff::jacobi_hessian_approx_trait<typename nlp_types::nlp_cos_type>
      cos_hessian_trait;
  typedef diff::jacobi_hessian_approx_trait<
      typename nlp_types::nlp_power_wl_type>
      power_wl_hessian_trait;
};

struct nlp_default_second_order_algorithms {
  typedef converge::converge_list<converge::converge_grad_norm_by_init<
                                      nlp_default_types::nlp_numerical_type>,
                                  converge::converge_criteria_max_iter<3000>>
      converge_type;
  // typedef optm::second_order::naive_gradient_descent<converge_type>
  // optm_type;
  typedef optm::second_order::adam<converge_type,
                                   nlp_default_types::nlp_numerical_type>
      optm_type;
  // typedef optm::second_order::nesterov<converge_type,
  // nlp_default_types::nlp_numerical_type> optm_type;

  /* multipliers */
  typedef outer_multiplier::init::init_by_matching_gradient_norm mult_init_type;
  // typedef
  // outer_multiplier::update::subgradient_normalized_by_init<nlp_default_types::nlp_numerical_type>
  // mult_update_type;
  typedef outer_multiplier::update::direct_subgradient mult_update_type;
  typedef outer_multiplier::mult_const_hpwl_cos_and_penalty_by_type<
      nlp_default_types::nlp_numerical_type, mult_init_type, mult_update_type>
      mult_type;
  typedef outer_multiplier::update::match_grad_const_multipliers<
      nlp_default_types::nlp_numerical_type>
      mult_adjust_type;

  /* alpha */
  typedef alpha::alpha_hpwl_ovl_oob<nlp_default_types::nlp_numerical_type>
      alpha_type;
  typedef alpha::update::alpha_update_list<
      alpha::update::reciprocal_by_obj<nlp_default_types::nlp_numerical_type,
                                       1>,
      alpha::update::reciprocal_by_obj<nlp_default_types::nlp_numerical_type,
                                       2>>
      alpha_update_type;
};

struct nlp_default_settings {
  typedef nlp_default_zero_order_algorithms nlp_zero_order_algorithms_type;
  typedef nlp_default_first_order_algorithms nlp_first_order_algorithms_type;
  typedef nlp_default_hyperparamters nlp_hyperparamters_type;
  typedef nlp_default_types nlp_types_type;
  typedef nlp_default_second_order_settings<nlp_types_type>
      nlp_second_order_setting_type;
  typedef nlp_default_second_order_algorithms nlp_second_order_algorithms_type;
};

} // namespace nlp

namespace _nlp_details {

template <typename fence_type> struct construct_fence_type_trait {};

template <typename nlp_numerical_type, typename nlp_coordinate_type>
struct construct_fence_type_trait<
    diff::FenceReciprocalOverlapSumBoxDifferentiable<nlp_numerical_type,
                                                     nlp_coordinate_type>> {
  static diff::FenceReciprocalOverlapSumBoxDifferentiable<nlp_numerical_type,
                                                          nlp_coordinate_type>
  constructFenceOperator(
      IndexType cellIdx, nlp_coordinate_type scale, const Cell &cell,
      const Well &well,
      const std::function<nlp_numerical_type(void)> &getAlphaFunc,
      const std::function<nlp_numerical_type(void)> &getLambdaFunc) {
    std::vector<Box<nlp_coordinate_type>> boxes; // Splited polygon
    std::vector<Box<LocType>> boxesUnScaled;
    if (not klib::convertPolygon2Rects(well.shape().outer(), boxesUnScaled)) {
      ERR("NlpGPlacer:: cannot split well polygon! \n");
      Assert(false);
    }
    for (const auto &box : boxesUnScaled) {
      boxes.emplace_back(
          Box<nlp_coordinate_type>(box.xLo() * scale, box.yLo() * scale,
                                   box.xHi() * scale, box.yHi() * scale));
    }
    return diff::FenceReciprocalOverlapSumBoxDifferentiable<
        nlp_numerical_type, nlp_coordinate_type>(
        cellIdx, cell.cellBBox().xLen() * scale, cell.cellBBox().yLen() * scale,
        boxes, getAlphaFunc, getLambdaFunc);
  }
  template<typename nlp_type>
    static void construct_operators(nlp_type &nlp,
      const std::function<nlp_numerical_type(void)> &getAlphaFunc,
      const std::function<nlp_numerical_type(void)> &getLambdaFunc,
      const std::function<nlp_coordinate_type(IndexType, Orient2DType)> getVarFunc
      ) {
      for (const auto &well : nlp._db.vWells()) {
        for (IndexType cellIdx : well.sCellIds()) {
          nlp._fenceOps.emplace_back(
              _nlp_details::construct_fence_type_trait<
                  typename nlp_type::nlp_fence_type>::constructFenceOperator(cellIdx, nlp._scale,
                                                          nlp._db.cell(cellIdx), well,
                                                          getAlphaFunc,
                                                          getLambdaFunc));
          nlp._fenceOps.back().setGetVarFunc(getVarFunc);
          nlp._fenceOps.back().setWeight(nlp._db.parameters().defaultWellWeight());
        }
      }
    }
};

template <typename nlp_numerical_type, typename nlp_coordinate_type>
struct construct_fence_type_trait<
    diff::FenceBivariateGaussianDifferentiable<nlp_numerical_type,
                                                     nlp_coordinate_type>> {
  typedef diff::FenceBivariateGaussianDifferentiable<nlp_numerical_type,
                                                          nlp_coordinate_type> 
                                                            op_type;
  static void generateDistributionParameters(std::vector<BivariateGaussianParameters<nlp_numerical_type>> & result, const Database &db, nlp_numerical_type scale, IndexType wellType) {
#if 0
    std::vector<std::vector<Box<nlp_coordinate_type>>> boxes; // Splited polygon
    std::vector<std::vector<Box<LocType>>> boxesUnScaled;
    boxes.resize(db.vWells().size());
    boxesUnScaled.resize(db.vWells().size());
    for (IndexType wellIdx = 0; wellIdx < db.vWells().size(); ++wellIdx) {
      if (not klib::convertPolygon2Rects(db.well(wellIdx).shape().outer(), boxesUnScaled[wellIdx])) {
        ERR("NlpGPlacer:: cannot split well polygon! \n");
        Assert(false);
      }
      for (const auto &box : boxesUnScaled[wellIdx]) {
        boxes[wellIdx].emplace_back(
            Box<nlp_coordinate_type>(box.xLo() * scale, box.yLo() * scale,

                                     box.xHi() * scale, box.yHi() * scale));
      }
    }
#else
    std::vector<std::vector<Box<nlp_coordinate_type>>> boxes;
    boxes.resize(db.vWells().size());
    for (IndexType wellIdx = 0; wellIdx < db.vWells().size(); ++wellIdx) {
      if (db.well(wellIdx).wellType() != wellType) {
        continue;
      }
      auto box = db.well(wellIdx).boundingBox();
      boxes[wellIdx].emplace_back(
            Box<nlp_coordinate_type>(box.xLo() * scale, box.yLo() * scale,
                                     box.xHi() * scale, box.yHi() * scale));
    }
#endif
    auto getNumWells = [&]() { return db.vWells().size(); };
    auto getNumRects = [&](IndexType wellIdx) { return boxes.at(wellIdx).size(); };
    auto getRect = [&](IndexType wellIdx, IndexType rectIdx) {  return &boxes.at(wellIdx).at(rectIdx); };
    BivariateGaussianWellApproximationGenerator<nlp_numerical_type> gen(getNumWells, getNumRects, getRect);
    gen.generate(result);
  }

  template<typename nlp_type>
    static void construct_operators(nlp_type &nlp,
      const std::function<nlp_numerical_type(void)> &,
      const std::function<nlp_numerical_type(void)> &getLambdaFunc,
      const std::function<nlp_coordinate_type(IndexType, Orient2DType)> getVarFunc
      ) {
      const Database &db = nlp._db;
      std::vector<IndexType> inFenceCellIdx;
      std::vector<nlp_coordinate_type> inWidths;
      std::vector<nlp_coordinate_type> inHeights;
      std::vector<IndexType> outFenceCellIdx;
      std::vector<nlp_coordinate_type> outWidths;
      std::vector<nlp_coordinate_type> outHeights;
      for (IndexType wellType = 0;  wellType < db.numWellTypes(); ++wellType) {
        for (IndexType cellIdx = 0; cellIdx < db.numCells(); ++cellIdx) {
          const auto &cell = db.cell(cellIdx);
          nlp_coordinate_type width = cell.cellBBox().width() * nlp._scale;
          nlp_coordinate_type height = cell.cellBBox().height() * nlp._scale;
          if (db.cell(cellIdx).wellType() == wellType) {
            inFenceCellIdx.emplace_back(cellIdx);
            inWidths.emplace_back(width);
            inHeights.emplace_back(height);
          }
          else {
            outFenceCellIdx.emplace_back(cellIdx);
            outWidths.emplace_back(width);
            outHeights.emplace_back(height);
          }
        }
        nlp._fenceOps.emplace_back(op_type(inFenceCellIdx, outFenceCellIdx,
              inWidths, inHeights,
              outWidths, outHeights,
               getLambdaFunc));
        nlp._fenceOps.back().setGetVarFunc(getVarFunc);
        nlp._fenceOps.back()._weight = db.parameters().defaultWellWeight();
        generateDistributionParameters(nlp._fenceOps.back()._gaussianParameters, db, nlp._scale, wellType);
      }
    }
};


template<typename fence_type> struct reinit_well_trait {};

template <typename nlp_numerical_type, typename nlp_coordinate_type>
struct reinit_well_trait<
    diff::FenceReciprocalOverlapSumBoxDifferentiable<nlp_numerical_type,
    nlp_coordinate_type>> {

    template<typename nlp_type> 
      static void reinit_well_operators(nlp_type &nlp) {
        nlp._fenceOps.clear();
        nlp._ovlOps.resize(nlp._numCellOvlOps);
        auto getVarFunc = [&](IndexType cellIdx, Orient2DType orient) {
#ifdef MULTI_SYM_GROUP
          return _pl(plIdx(cellIdx, orient));
#else
          if (orient == Orient2DType::NONE) {
            return nlp._defaultSymAxis;
          }
          return nlp._pl(nlp.plIdx(cellIdx, orient));
#endif
        };

        auto getFenceAlphaFunc = nlp_type::alpha_trait::fenceGetAlphaFunc(*nlp._alpha);
        auto getFenceLambdaFunc = nlp_type::mult_trait::fenceGetLambdaFunc(*nlp._multiplier);
        auto getOvlAlphaFunc = nlp_type::alpha_trait::ovlGetAlphaFunc(*nlp._alpha);
        auto getOvlLambdaFunc = nlp_type::mult_trait::ovlGetLambdaFunc(*nlp._multiplier);

        if (not nlp._hasInitFirstOrderOuterProblem) {
          // Before init first order outer problem
          // The alpha, lambda have not been constructed
          // Just use naive values
          getFenceAlphaFunc = [&]() { return 1.0; };
          getFenceLambdaFunc = [&]() { return  1.0; };
          getOvlAlphaFunc = [&](){ return 1.0; };
          getOvlLambdaFunc = [&]() { return  1.0; };

        }

        // Pair-wise well-to-cell overlapping
        if (nlp._useWellCellOvl) {
          for (IndexType wellIdx = 0; wellIdx < nlp._db.vWells().size();
               ++wellIdx) {
            const auto &well = nlp._db.vWells().at(wellIdx);
            std::vector<Box<typename nlp_type::base_type::nlp_coordinate_type>>
                boxes; // Splited polygon
            std::vector<Box<LocType>> boxesUnScaled;
            if (not klib::convertPolygon2Rects(well.shape().outer(), boxesUnScaled)) {
              ERR("NlpGPlacer:: cannot split well polygon! \n");
              Assert(false);
            }
            for (const auto &box : boxesUnScaled) {
              boxes.emplace_back(
                  Box<typename nlp_type::base_type::nlp_coordinate_type>(
                  (box.xLo()  - box.xLen() / 2)* nlp._scale, (box.yLo() - box.yLen() / 2)* nlp._scale,
                  (box.xHi()  + box.xLen() / 2)* nlp._scale, (box.yHi() + box.yLen() / 2)* nlp._scale));
            }
            for (IndexType cellIdx = 0; cellIdx < nlp._db.numCells();
                 ++cellIdx) {
              if (well.sCellIds().find(cellIdx) == well.sCellIds().end()) {
                const auto cellBBox = nlp._db.cell(cellIdx).cellBBox();
                for (const auto &wellBox : boxes) {
                  nlp._ovlOps.emplace_back(
                      typename nlp_type::nlp_ovl_type(cellIdx, cellBBox.xLen() * nlp._scale,
                                   cellBBox.yLen() * nlp._scale,
                                   666666, // Doesn't matter
                                   wellBox.xLen(), wellBox.yLen(), getOvlAlphaFunc,
                                   getOvlLambdaFunc));
                  nlp._ovlOps.back().configConsiderOnlyOneCell(wellBox.xLo(),
                                                                      wellBox.yLo());
                  nlp._ovlOps.back().setGetVarFunc(getVarFunc);
                  nlp._ovlOps.back().setWeight(nlp._db.parameters().defaultWellWeight());
                }
              }
            }
          }
        }
        for (const auto &well : nlp._db.vWells()) {
          for (IndexType cellIdx : well.sCellIds()) {
            nlp._fenceOps.emplace_back(
                _nlp_details::construct_fence_type_trait<typename nlp_type::nlp_fence_type>::
                    constructFenceOperator(cellIdx, nlp._scale,
                                           nlp._db.cell(cellIdx), well,
                                           getFenceAlphaFunc, getFenceLambdaFunc));
            nlp._fenceOps.back().setGetVarFunc(getVarFunc);
            nlp._fenceOps.back().setWeight(
                nlp._db.parameters().defaultWellWeight());
          }
        }
        nlp._fenceOps.back()._considerOutFenceCells = nlp._useWellCellOvl;

        // Re-construct tasks
        nlp.clearTasks();
        nlp.constructTasks();
      }
};
template <typename nlp_numerical_type, typename nlp_coordinate_type>
struct reinit_well_trait<
    diff::FenceBivariateGaussianDifferentiable<nlp_numerical_type,
    nlp_coordinate_type>> {
      typedef diff::FenceBivariateGaussianDifferentiable<nlp_numerical_type,
    nlp_coordinate_type> op_type;

    template<typename nlp_type> 
      static void reinit_well_operators(nlp_type &nlp) {
        for (IndexType wellType = 0; wellType < nlp._db.numWellTypes(); ++wellType) {
          construct_fence_type_trait<op_type>::generateDistributionParameters(nlp._fenceOps.back()._gaussianParameters, nlp._db, nlp._scale, wellType);
          nlp._fenceOps.back()._considerOutFenceCells = nlp._useWellCellOvl;
        }
      }
};
} // namespace _nlp_details

/// @brief non-linear programming-based analog global placement
template <typename nlp_settings> class NlpGPlacerBase {
public:
  typedef typename nlp_settings::nlp_types_type nlp_types;
  typedef typename nlp_settings::nlp_zero_order_algorithms_type
      nlp_zero_order_algorithms;
  typedef typename nlp_settings::nlp_hyperparamters_type nlp_hyperparamters;

  typedef typename nlp_types::EigenMatrix EigenMatrix;
  typedef typename nlp_types::EigenVector EigenVector;
  typedef typename nlp_types::EigenMap EigenMap;
  typedef typename nlp_types::nlp_coordinate_type nlp_coordinate_type;
  typedef typename nlp_types::nlp_numerical_type nlp_numerical_type;
  typedef typename nlp_types::nlp_hpwl_type nlp_hpwl_type;
  typedef typename nlp_types::nlp_ovl_type nlp_ovl_type;
  typedef typename nlp_types::nlp_oob_type nlp_oob_type;
  typedef typename nlp_types::nlp_asym_type nlp_asym_type;
  typedef typename nlp_types::nlp_cos_type nlp_cos_type;
  typedef typename nlp_types::nlp_power_wl_type nlp_power_wl_type;
  typedef typename nlp_types::nlp_crf_type nlp_crf_type;
  typedef typename nlp_types::nlp_fence_type nlp_fence_type;
  typedef typename nlp_types::nlp_area_type nlp_area_type;

  /* algorithms */
  typedef typename nlp_zero_order_algorithms::stop_condition_type
      stop_condition_type;
  typedef nlp::outer_stop_condition::stop_condition_trait<stop_condition_type>
      stop_condition_trait;
  template <typename _T>
  friend struct nlp::outer_stop_condition::stop_condition_trait;
  typedef
      typename nlp_zero_order_algorithms::init_place_type init_placement_type;
  typedef nlp::init_place::init_place_trait<init_placement_type>
      init_place_trait;
  friend init_place_trait;
  typedef typename nlp_zero_order_algorithms::mult_init_type mult_init_type;
  typedef nlp::outer_multiplier::init::multiplier_init_trait<mult_init_type>
      mult_init_trait;
  friend mult_init_trait;
  typedef typename nlp_zero_order_algorithms::mult_update_type mult_update_type;
  typedef nlp::outer_multiplier::update::multiplier_update_trait<
      mult_update_type>
      mult_update_trait;
  friend mult_update_trait;
  typedef typename nlp_zero_order_algorithms::mult_type mult_type;
  typedef nlp::outer_multiplier::multiplier_trait<mult_type> mult_trait;
  friend mult_trait;


  /* Implementation details */
  friend _nlp_details::construct_fence_type_trait<nlp_fence_type>;

public:
  explicit NlpGPlacerBase(Database &db) : _db(db) {}
  IntType solve();

  /* Input/Output functions */
  /// @brief write the GPlacer locations to database
  void writeLocs();
  /// @brief read the cell locations from database
  void readLocs();
  /* Well-related initialization */
  virtual void reinitWellOperators();
  /// @brief get the overlapping area ratio
  RealType overlapAreaRatio() {
    nlp_coordinate_type overlapArea = 0.0;
    for (auto &op : _ovlOps) {
      overlapArea += diff::place_overlap_trait<
          nlp_ovl_type>::overlapArea(op);
      }
    return overlapArea / _totalCellArea;
  }

  void debug() {
    double xmin =99999;
    double xmax =-99999;
    double ymin = 99999;
    double ymax = - 999999;
    DBG("PRINT Stats \n");
    for (IndexType cellIdx = 0; cellIdx <  _db.numCells(); ++cellIdx) {
      DBG("Cell %d x %f y %f \n", cellIdx, _pl(plIdx(cellIdx, Orient2DType::HORIZONTAL)), _pl(plIdx(cellIdx, Orient2DType::VERTICAL)));
      xmin = std::min(xmin, _pl(plIdx(cellIdx, Orient2DType::HORIZONTAL)));
      xmax = std::max(xmax, _pl(plIdx(cellIdx, Orient2DType::HORIZONTAL)));
      ymin = std::min(ymin, _pl(plIdx(cellIdx, Orient2DType::VERTICAL)));
      ymax = std::max(ymax, _pl(plIdx(cellIdx, Orient2DType::VERTICAL)));
    }
    for (IndexType wellIdx = 0; wellIdx < _db.vWells().size(); ++wellIdx) {
      DBG("Well %d \n", wellIdx);
      _db.well(wellIdx).printInfo();
    }
    DBG("_scale %f \n", _scale);
    DBG("xmin %f xmax %f, ymin %f ymax %f \n", xmin, xmax, ymin, ymax);
    std::string normal_distr = "./debug/gaussian.txt";
    std::string cost_in = "./debug/cost_in.txt";
    std::string grad_in = "./debug/grad_in.txt";
    std::ofstream fGaussian, fCostIn, fGradIn;
    fGaussian.open(normal_distr);
    fCostIn.open(cost_in);
    fGradIn.open(grad_in);
    std::vector<BivariateGaussianParameters<nlp_numerical_type>> &paras = _fenceOps.back()._gaussianParameters;
    for (IndexType gauIdx  = 0; gauIdx < paras.size(); ++gauIdx) {
      const auto muX = paras[gauIdx].muX;
      const auto muY = paras[gauIdx].muY;
      const auto sigmaX = paras[gauIdx].sigmaX;
      const auto sigmaY = paras[gauIdx].sigmaY;
      const auto normalize = paras[gauIdx].normalize;
      fGaussian<< muX << " "<< sigmaX << " "<< muY << " "<< sigmaY <<  " "  << normalize << "\n";
    }
    paras.clear();
    paras.emplace_back(BivariateGaussianParameters<nlp_numerical_type>(6.5, 2.5, 4.0, 2.2, 1.0));
    paras.emplace_back(BivariateGaussianParameters<nlp_numerical_type>(-6.5, 3.0, -7.5, 2.0, 1.0));
    fGaussian.close();
    xmax = 10;
    xmin = -10;
    ymax = 10;
    ymin = -10;
    double width = (xmax - xmin) / 100;
    double height = (ymax - ymin) / 100;
    double size = std::max(width, height);
    for (int xStep = 0; xStep < 100; ++xStep) {
      double xLo = xmin + size * xStep;
      for (int yStep = 0; yStep < 100; ++yStep) {
        double yLo = ymin + size * yStep;
        double cost = 0;
        double xDiff = 0;
        double yDiff = 0;
        for (IndexType gauIdx  = 0; gauIdx < paras.size(); ++gauIdx) {
          const auto muX = paras[gauIdx].muX;
          const auto muY = paras[gauIdx].muY;
          const auto sigmaX = paras[gauIdx].sigmaX;
          const auto sigmaY = paras[gauIdx].sigmaY;
          cost += diff::_fence_bivariate_gaussian_details::calc_trait<diff::FENCE_SIGMOID_COST>::calcCost(
              muX, muY, sigmaX, sigmaY, xLo, yLo, width, height, 0.2, 1.0
              );

          const auto diffx = diff::_fence_bivariate_gaussian_details::calc_trait<diff::FENCE_SIGMOID_COST>::calcDiffX(
              muX, muY, sigmaX, sigmaY, xLo, yLo, width, height, 0.2, 1.0
              );
          const auto diffy = diff::_fence_bivariate_gaussian_details::calc_trait<diff::FENCE_SIGMOID_COST>::calcDiffY(
              muX, muY, sigmaX, sigmaY, xLo, yLo, width, height, 0.2, 1.0
              );
          xDiff += diffx;
          yDiff += diffy;
        }
        fCostIn << xStep << " " << yStep << " "<< xLo <<" "<< yLo <<" "<< cost << "\n";
        if (xStep % 10 == 0 and yStep % 10 == 0) {
          fGradIn << xStep << " " << yStep << " "<< xLo <<" "<< yLo <<" "<< xDiff << " "<< yDiff << "\n";
        }
      }
    }
    fCostIn.close();
    fGradIn.close();
  }
protected:

  /* construct tasks */
  virtual void constructTasks();
  virtual void clearTasks();
  void assignIoPins();
  /* calculating obj */
  void calcObj() { _wrapObjAllTask.run(); }
  /* Init functions */
  virtual void initProblem();
  void initPlace();
  void initOperators();
  void initHyperParams();
  void initBoundaryParams();
  void initVariables();
  void initOptimizationKernelMembers();
  /* Util functions */
  IndexType plIdx(IndexType cellIdx, Orient2DType orient);
  void alignToSym();
  // Obj-related
  void constructObjTasks();
  void constructObjectiveCalculationTasks();
  void constructSumObjTasks();
#ifdef DEBUG_SINGLE_THREAD_GP
  void constructWrapObjTask();
#endif
  /* Optimization  kernel */
  virtual void optimize();
  /* Debugging function */
#ifdef DEBUG_GR
#ifdef DEBUG_DRAW
  void drawCurrentLayout(const std::string &name);
#endif
#endif
protected:
  Database &_db; ///< The placement engine database
  /* NLP problem parameters */
  IndexType _numCells; ///< The number of cells
  Box<nlp_coordinate_type>
      _boundary; ///< The boundary constraint for the placement
  nlp_coordinate_type _scale =
      0.01; /// The scale ratio between float optimization kernel coordinate and
            /// placement database coordinate unit
  nlp_coordinate_type _totalCellArea =
      0; ///< The total cell area of the problem
  nlp_coordinate_type _defaultSymAxis = 0.0; ///< The default symmetric axis
  IndexType _numVariables = 0;               ///< The number of variables
  /* Optimization internal results */
  nlp_numerical_type _objHpwl = 0.0; ///< The current value for hpwl
  nlp_numerical_type _objOvl =
      0.0; ///< The current value for overlapping penalty
  nlp_numerical_type _objOob =
      0.0; ///< The current value for out of boundary penalty
  nlp_numerical_type _objAsym =
      0.0; ///< The current value for asymmetry penalty
  nlp_numerical_type _objCos =
      0.0; ///< The current value for the cosine signal path penalty
  nlp_numerical_type _objPowerWl = 0.0; ///< power wire length
  nlp_numerical_type _objCrf = 0.0;     ///< Current flow
  nlp_numerical_type _objFence = 0.0;   ///< Fence region
  nlp_numerical_type _objArea = 0.0; ///< Area objective
  nlp_numerical_type _obj =
      0.0; ///< The current value for the total objective penalty
  nlp_numerical_type _objHpwlRaw = 0.0; ///< The current value for hpwl
  nlp_numerical_type _objOvlRaw =
      0.0; ///< The current value for overlapping penalty
  nlp_numerical_type _objOobRaw =
      0.0; ///< The current value for out of boundary penalty
  nlp_numerical_type _objAsymRaw =
      0.0; ///< The current value for asymmetry penalty
  nlp_numerical_type _objCosRaw =
      0.0; ///< The current value for the cosine signal path penalty
  nlp_numerical_type _objPowrWlRaw = 0.0; ///< Power wire length
  nlp_numerical_type _objCrfRaw = 0.0;    ///< Current flow
  nlp_numerical_type _objFenceRaw = 0.0;  ///< Fence region
  nlp_numerical_type _objAreaRaw = 0.0; ///< Area
  /* NLP optimization kernel memebers */
  stop_condition_type _stopCondition;
  /* Optimization data */
  EigenVector _pl; ///< The placement solutions
  /* Tasks */
  // Evaluating objectives
  std::vector<nt::Task<nt::EvaObjTask<nlp_numerical_type>>> _evaHpwlTasks; ///< The tasks for evaluating hpwl objectives
  std::vector<nt::Task<nt::EvaObjTask<nlp_numerical_type>>> _evaOvlTasks; ///< The tasks for evaluating overlap objectives
  std::vector<nt::Task<nt::EvaObjTask<nlp_numerical_type>>> _evaOobTasks; ///< The tasks for evaluating out of boundary objectives
  std::vector<nt::Task<nt::EvaObjTask<nlp_numerical_type>>> _evaAsymTasks; ///< The tasks for evaluating asymmetry objectives
  std::vector<nt::Task<nt::EvaObjTask<nlp_numerical_type>>> _evaCosTasks; ///< The tasks for evaluating signal path objectives
  std::vector<nt::Task<nt::EvaObjTask<nlp_numerical_type>>> _evaPowerWlTasks; ///< The tasks for evaluating power wirelength objectives
  std::vector<nt::Task<nt::EvaObjTask<nlp_numerical_type>>> _evaCrfTasks; ///< The tasks for evaluating current flow objectives
  std::vector<nt::Task<nt::EvaObjTask<nlp_numerical_type>>> _evaFenceTasks; ///< The tasks for evaluating fence region objectives
  std::vector<nt::Task<nt::EvaObjTask<nlp_numerical_type>>> _evaAreaTasks; ///< The tasks for evaluating fence region objectives
  // Sum the objectives
  nt::Task<nt::FuncTask> _sumObjHpwlTask; ///< The task for summing hpwl objective
  nt::Task<nt::FuncTask> _sumObjOvlTask; ///< The task for summing the overlapping objective
  nt::Task<nt::FuncTask> _sumObjOobTask; ///< The task for summing the out of boundary objective
  nt::Task<nt::FuncTask> _sumObjAsymTask; ///< The task for summing the asymmetry objective
  nt::Task<nt::FuncTask> _sumObjCosTask; ///< The task for summing the cosine signal path objective
  nt::Task<nt::FuncTask> _sumObjPowerWlTask; ///< The task for summing the
                                             ///< cosine signal path objective
  nt::Task<nt::FuncTask> _sumObjCrfTask; ///< The task for summing the current flow objective
  nt::Task<nt::FuncTask> _sumObjFenceTask; ///< The task for summing the fence region objecitve
  nt::Task<nt::FuncTask> _sumObjAreaTask; ///< The task for summing the fence region objecitve
  nt::Task<nt::FuncTask> _sumObjAllTask; ///< The task for summing the different
                                         ///< objectives together
  // Wrapper tasks for debugging
  nt::Task<nt::FuncTask>
      _wrapObjHpwlTask; ///< The task for wrap the objective for wirelength
  nt::Task<nt::FuncTask> _wrapObjOvlTask;
  nt::Task<nt::FuncTask> _wrapObjOobTask;
  nt::Task<nt::FuncTask> _wrapObjAsymTask;
  nt::Task<nt::FuncTask> _wrapObjCosTask;
  nt::Task<nt::FuncTask> _wrapObjPowerWlTask;
  nt::Task<nt::FuncTask> _wrapObjCrfTask;   ///< The wrapper for caculating the
                                            ///< current flow objective
  nt::Task<nt::FuncTask> _wrapObjFenceTask; ///< The wrapper for calculating the
                                            ///< fence region objective
  nt::Task<nt::FuncTask> _wrapObjAreaTask; ///< The wrapper for calculating the
                                            ///< fence region objective
  nt::Task<nt::FuncTask> _wrapObjAllTask;
  /* Operators */
  std::vector<nlp_hpwl_type> _hpwlOps; ///< The HPWL cost
  std::vector<nlp_ovl_type>
      _ovlOps;              ///< The cell pair overlapping penalty operators
  IndexType _numCellOvlOps; ///< The number of ovl ops that are overlapping
                            ///< between cells
  std::vector<nlp_oob_type>
      _oobOps; ///< The cell out of boundary penalty operators
  std::vector<nlp_asym_type> _asymOps; ///< The asymmetric penalty operators
  std::vector<nlp_cos_type> _cosOps;   ///< The signal flow operators
  std::vector<nlp_power_wl_type> _powerWlOps;
  std::vector<nlp_crf_type> _crfOps;     ///< The current flow operators
  std::vector<nlp_fence_type> _fenceOps; ///< The fence region operators
  std::vector<nlp_area_type> _areaOps; ///< The area objective operators
};

template <typename nlp_settings>
inline IndexType NlpGPlacerBase<nlp_settings>::plIdx(IndexType cellIdx,
                                                     Orient2DType orient) {
  if (orient == Orient2DType::HORIZONTAL) {
    return cellIdx;
  } else if (orient == Orient2DType::VERTICAL) {
    return cellIdx + _numCells;
  } else {
#ifdef MULTI_SYM_GROUP
    return cellIdx +
           2 * _numCells; // here cell index representing the idx of sym grp
#else
    return 2 * _numCells;
#endif
  }
}

/// @brief first-order optimization
template <typename nlp_settings>
class NlpGPlacerFirstOrder : public NlpGPlacerBase<nlp_settings> {
public:
  typedef NlpGPlacerBase<nlp_settings> base_type;
  typedef typename base_type::EigenVector EigenVector;
  typedef typename base_type::nlp_hpwl_type nlp_hpwl_type;
  typedef typename base_type::nlp_ovl_type nlp_ovl_type;
  typedef typename base_type::nlp_oob_type nlp_oob_type;
  typedef typename base_type::nlp_asym_type nlp_asym_type;
  typedef typename base_type::nlp_cos_type nlp_cos_type;
  typedef typename base_type::nlp_power_wl_type nlp_power_wl_type;
  typedef typename base_type::nlp_crf_type nlp_crf_type;
  typedef typename base_type::nlp_fence_type nlp_fence_type;
  typedef typename base_type::nlp_area_type nlp_area_type;

  typedef typename nlp_settings::nlp_first_order_algorithms_type
      nlp_first_order_algorithms;
  typedef typename nlp_first_order_algorithms::converge_type converge_type;
  typedef typename nlp::converge::converge_criteria_trait<converge_type>
      converge_trait;
  typedef typename nlp_first_order_algorithms::optm_type optm_type;
  typedef typename nlp::optm::optm_trait<optm_type> optm_trait;

  friend converge_type;
  template <typename converge_criteria_type>
  friend struct nlp::converge::converge_criteria_trait;
  friend optm_type;
  friend optm_trait;

  typedef typename nlp_settings::nlp_first_order_algorithms_type::mult_init_type
      mult_init_type;
  typedef nlp::outer_multiplier::init::multiplier_init_trait<mult_init_type>
      mult_init_trait;
  friend mult_init_trait;
  typedef
      typename nlp_settings::nlp_first_order_algorithms_type::mult_update_type
          mult_update_type;
  typedef nlp::outer_multiplier::update::multiplier_update_trait<
      mult_update_type>
      mult_update_trait;
  friend mult_update_trait;
  typedef typename nlp_settings::nlp_first_order_algorithms_type::mult_type
      mult_type;
  typedef nlp::outer_multiplier::multiplier_trait<mult_type> mult_trait;
  friend mult_trait;

  typedef
      typename nlp_settings::nlp_first_order_algorithms_type::mult_adjust_type
          mult_adjust_type;
  typedef nlp::outer_multiplier::update::multiplier_update_trait<
      mult_adjust_type>
      mult_adjust_trait;
  friend mult_adjust_trait;

  /* updating alpha parameters */
  typedef typename nlp_settings::nlp_first_order_algorithms_type::alpha_type
      alpha_type;
  typedef nlp::alpha::alpha_trait<alpha_type> alpha_trait;
  template <typename T> friend struct nlp::alpha::alpha_trait;
  typedef
      typename nlp_settings::nlp_first_order_algorithms_type::alpha_update_type
          alpha_update_type;
  typedef nlp::alpha::update::alpha_update_trait<alpha_update_type>
      alpha_update_trait;
  template <typename T> friend struct nlp::alpha::update::alpha_update_trait;

  /* Implementation details */
  friend _nlp_details::reinit_well_trait<nlp_fence_type>;
  friend _nlp_details::construct_fence_type_trait<nlp_fence_type>;

  NlpGPlacerFirstOrder(Database &db) : NlpGPlacerBase<nlp_settings>(db) {}
  void writeoutCsv() {
    std::string ver1f = "ver1.csv";
    std::string ver2f = "ver2.csv";
    std::ofstream ver1(ver1f.c_str());
    std::ofstream ver2(ver2f.c_str());
    ver1 << "x y val\n";
    ver2 << "x y val\n";

    for (RealType x = -8; x < 8; x += (16.0 / 300)) {
      for (RealType y = -8; y < 8; y += (16.0 / 300)) {
        this->_pl(this->plIdx(0, Orient2DType::HORIZONTAL)) = x;
        this->_pl(this->plIdx(0, Orient2DType::VERTICAL)) = y;
        this->_wrapObjAllTask.run();
        auto obj = this->_obj;
        ver1 << x << " " << y << " " << obj << "\n";
      }
    }
    for (auto &op : this->_hpwlOps) {
      op._getLambdaFunc = [&]() { return 1.0; };
    }
    for (auto &op : this->_cosOps) {
      op._getLambdaFunc = [&]() { return 1.0; };
    }
    for (auto &op : this->_ovlOps) {
      op._getLambdaFunc = [&]() { return 1.0; };
    }
    for (auto &op : this->_oobOps) {
      op._getLambdaFunc = [&]() { return 1.0; };
    }
    for (auto &op : this->_asymOps) {
      op._getLambdaFunc = [&]() { return 1.0; };
    }
    for (auto &op : this->_powerWlOps) {
      op._getLambdaFunc = [&]() { return 1.0; };
    }
    for (auto &op : this->_crfOps) {
      op._getLambdaFunc = [&]() { return 1.0; };
    }
    for (RealType x = -8; x < 8; x += (16.0 / 300)) {
      for (RealType y = -8; y < 8; y += (16.0 / 300)) {
        this->_pl(this->plIdx(0, Orient2DType::HORIZONTAL)) = x;
        this->_pl(this->plIdx(0, Orient2DType::VERTICAL)) = y;
        this->_wrapObjAllTask.run();
        auto obj = this->_obj;
        ver2 << x << " " << y << " " << obj << "\n";
      }
    }
  }
  /* Well related */
  virtual void reinitWellOperators() override;
  /* Parameters */
  void openUseWellCellOvl() { _useWellCellOvl = true; }
  void closeUseWellCellOvl() { _useWellCellOvl = false; }
  /* Init */
  virtual void prepareWellAwarePlace();
  /* Core Optimization*/
  /// @brief Finish one iteration of optimization and update the problem
  void stepOptmIter();
  /// @brief query whether meet stop condition
  /// @return true: meet stop condition. false: does not meet
  BoolType meetStopCondition() {
    return base_type::stop_condition_trait::stopPlaceCondition(
      *this, this->_stopCondition);
  }
  void cleanupMode();
protected:
  /* Construct tasks */
  virtual void constructTasks() override;
  virtual void clearTasks() override;
  /* calculating gradient */
  void calcGrad() { _wrapCalcGradTask.run(); }
  /* Init */
  virtual void initProblem() override;
  void initFirstOrderGrad();
  void constructFirstOrderTasks();
  void constructCalcPartialsTasks();
  void constructUpdatePartialsTasks();
  void constructClearGradTasks();
  void constructSumGradTask();
  void constructWrapCalcGradTask();
  /* optimization */
  virtual void optimize() override;
  void printGradNorm() {
      DBG("grad hpwl %f, ovl %f, oob %f, asym %f, cos %f, pwl %f, crf %f, fence %f \n",
            _gradHpwl.norm(), _gradOvl.norm(), _gradOob.norm(), _gradAsym.norm(), _gradCos.norm(),
            _gradPowerWl.norm(), _gradCrf.norm(), _gradFence.norm());
      auto &grad = _gradOvl;
      double sumgrad = 0;
      std::cout<<"Grad X: ";
      for (IndexType cellIdx = 0; cellIdx < this->_db.numCells();++cellIdx) {
        sumgrad +=grad(this->plIdx(cellIdx, Orient2DType::HORIZONTAL));
        std::cout<<std::to_string(cellIdx)<< ": "<<grad(this->plIdx(cellIdx, Orient2DType::HORIZONTAL));
        std::cout<<". ";
        if (cellIdx %5 == 0) {
          std::cout<<"\n";
        }
      }
      std::cout<<"\n\nSum: "<< sumgrad<<"\n\n";
      std::cout<<"Grad Y: ";
      sumgrad = 0;
      for (IndexType cellIdx = 0; cellIdx < this->_db.numCells();++cellIdx) {
        sumgrad +=grad(this->plIdx(cellIdx, Orient2DType::VERTICAL));
        //std::cout<<std::to_string(cellIdx)<< ": "<<grad(this->plIdx(cellIdx, Orient2DType::VERTICAL));
        //std::cout<<". ";
        //if (cellIdx %5 == 0) {
        //  std::cout<<"\n";
        //}
      }
      //std::cout<<"\n\nSum: "<<sumgrad<<"\n\n"<<std::endl;
  }


private:
  void initFirstOrderOuterProblem() {
    if (_hasInitFirstOrderOuterProblem) { return; }
    _hasInitFirstOrderOuterProblem = true;
    // setting up the multipliers
    this->_wrapObjAllTask.run();
    _wrapCalcGradTask.run();


    _multiplier = std::make_shared<mult_type>(mult_trait::construct(*this));
    mult_trait::init(*this, *_multiplier);
    mult_trait::recordRaw(*this, *_multiplier);

    _multAdjuster = std::make_shared<mult_adjust_type>(
        mult_adjust_trait::construct(*this, *_multiplier));
    mult_adjust_trait::init(*this, *_multiplier, *_multAdjuster);

    _alpha = std::make_shared<alpha_type>(alpha_trait::construct(*this));
    alpha_trait::init(*this, *_alpha);
    _alphaUpdate = std::make_shared<alpha_update_type>(
        alpha_update_trait::construct(*this, *_alpha));
    alpha_update_trait::init(*this, *_alpha, *_alphaUpdate);

  }
protected:
  /* Optimization data */
  EigenVector _grad;     ///< The first order graident
  EigenVector _gradHpwl; ///< The first order gradient of hpwl objective
  EigenVector _gradOvl;  ///< The first order gradient  of overlapping objective
  EigenVector _gradOob; ///< The first order gradient of out of boundary objective
  EigenVector _gradAsym; ///< The first order gradient of asymmetry objective
  EigenVector _gradCos; ///< The first order gradient of cosine signal path objective
  EigenVector _gradPowerWl;
  EigenVector _gradCrf;
  EigenVector _gradFence; ///< The gradient of fence region objective
  EigenVector _gradArea; ///< The gradient of area objective
  /* Tasks */
  // Calculate the partials
  std::vector<
      nt::Task<nt::CalculateOperatorPartialTask<nlp_hpwl_type, EigenVector>>>
      _calcHpwlPartialTasks;
  std::vector<
      nt::Task<nt::CalculateOperatorPartialTask<nlp_ovl_type, EigenVector>>>
      _calcOvlPartialTasks;
  std::vector<
      nt::Task<nt::CalculateOperatorPartialTask<nlp_oob_type, EigenVector>>>
      _calcOobPartialTasks;
  std::vector<
      nt::Task<nt::CalculateOperatorPartialTask<nlp_asym_type, EigenVector>>>
      _calcAsymPartialTasks;
  std::vector<
      nt::Task<nt::CalculateOperatorPartialTask<nlp_cos_type, EigenVector>>>
      _calcCosPartialTasks;
  std::vector<nt::Task<
      nt::CalculateOperatorPartialTask<nlp_power_wl_type, EigenVector>>>
      _calcPowerWlPartialTasks;
  std::vector<
      nt::Task<nt::CalculateOperatorPartialTask<nlp_crf_type, EigenVector>>>
      _calcCrfPartialTasks;
  std::vector<
      nt::Task<nt::CalculateOperatorPartialTask<nlp_fence_type, EigenVector>>>
      _calcFencePartialTasks;
  std::vector<
      nt::Task<nt::CalculateOperatorPartialTask<nlp_area_type, EigenVector>>>
      _calcAreaPartialTasks;
  // Update the partials
  std::vector<
      nt::Task<nt::UpdateGradientFromPartialTask<nlp_hpwl_type, EigenVector>>>
      _updateHpwlPartialTasks;
  std::vector<
      nt::Task<nt::UpdateGradientFromPartialTask<nlp_ovl_type, EigenVector>>>
      _updateOvlPartialTasks;
  std::vector<
      nt::Task<nt::UpdateGradientFromPartialTask<nlp_oob_type, EigenVector>>>
      _updateOobPartialTasks;
  std::vector<
      nt::Task<nt::UpdateGradientFromPartialTask<nlp_asym_type, EigenVector>>>
      _updateAsymPartialTasks;
  std::vector<
      nt::Task<nt::UpdateGradientFromPartialTask<nlp_cos_type, EigenVector>>>
      _updateCosPartialTasks;
  std::vector<nt::Task<
      nt::UpdateGradientFromPartialTask<nlp_power_wl_type, EigenVector>>>
      _updatePowerWlPartialTasks;
  std::vector<
      nt::Task<nt::UpdateGradientFromPartialTask<nlp_crf_type, EigenVector>>>
      _updateCrfPartialTasks;
  std::vector<
      nt::Task<nt::UpdateGradientFromPartialTask<nlp_fence_type, EigenVector>>>
      _updateFencePartialTasks;
  std::vector<
      nt::Task<nt::UpdateGradientFromPartialTask<nlp_area_type, EigenVector>>>
      _updateAreaPartialTasks;
  // Clear the gradient. Use to clear the _gradxxx records. Needs to call before
  // updating the partials
  nt::Task<nt::FuncTask> _clearGradTask; // FIXME: not used right noe
  nt::Task<nt::FuncTask> _clearHpwlGradTask;
  nt::Task<nt::FuncTask> _clearOvlGradTask;
  nt::Task<nt::FuncTask> _clearOobGradTask;
  nt::Task<nt::FuncTask> _clearAsymGradTask;
  nt::Task<nt::FuncTask> _clearCosGradTask;
  nt::Task<nt::FuncTask> _clearPowerWlGradTask;
  nt::Task<nt::FuncTask> _clearCrfGradTask;
  nt::Task<nt::FuncTask> _clearFenceGradTask;
  nt::Task<nt::FuncTask> _clearAreaGradTask;
  // Sum the _grad from individual
  nt::Task<nt::FuncTask> _sumGradTask;
  nt::Task<nt::FuncTask> _sumHpwlGradTask;
  nt::Task<nt::FuncTask> _sumOvlGradTask;
  nt::Task<nt::FuncTask> _sumOobGradTask;
  nt::Task<nt::FuncTask> _sumAsymGradTask;
  nt::Task<nt::FuncTask> _sumCosGradTask;
  nt::Task<nt::FuncTask> _sumPowerWlTaskGradTask;
  nt::Task<nt::FuncTask> _sumCrfGradTask;
  nt::Task<nt::FuncTask> _sumFenceGradTask;
  nt::Task<nt::FuncTask> _sumAreaGradTask;
  // all the grads has been calculated but have not updated
  nt::Task<nt::FuncTask> _wrapCalcGradTask; ///<  calculating the gradient and sum them
  /* optim and problem */
  optm_type _optm;
  std::shared_ptr<mult_type> _multiplier;
  std::shared_ptr<mult_adjust_type> _multAdjuster;
  std::shared_ptr<alpha_type> _alpha;
  std::shared_ptr<alpha_update_type> _alphaUpdate;
  /* Parameters */
  BoolType _useWellCellOvl = false; ///<
  bool _useSimpleOptm = false;
  /* Misc. */
  BoolType _hasInitFirstOrderOuterProblem = false; ///< Whether has init multiplier, alpha, etc.
  IndexType _iter = 0;
  BoolType _convergeWithLargeVariableChange = false; ///< The last termination of optimization is because of large location changes
  std::vector<IndexType> _maskOutCells; ///< Cells that should not update values
};

//// @brief some helper function for NlpGPlacerSecondOrder
namespace _nlp_second_order_details {
template <typename nlp_settings, BoolType is_diagonal>
struct is_diagonal_select {};

template <typename nlp_settings> struct is_diagonal_select<nlp_settings, true> {
  typedef typename nlp_settings::nlp_types_type::EigenMatrix matrix_type;
  static void resize(matrix_type &matrix, IntType size) {
    matrix.resize(size, size);
  }

  static decltype(auto) inverse(matrix_type &matrix) {
    return matrix.diagonal().cwiseInverse().asDiagonal();
  }
};

template <typename nlp_settings>
struct is_diagonal_select<nlp_settings, false> {
  typedef typename nlp_settings::nlp_types_type::EigenMatrix matrix_type;
  static void resize(matrix_type &matrix, IntType size) {
    matrix.resize(size, size);
  }
  static decltype(auto) inverse(matrix_type &matrix) {
    Assert(false);
    return matrix.diagonal().cwiseInverse();
  }
};

template <typename hessian_target_type> struct update_hessian {
  hessian_target_type &target;
};
}; // namespace _nlp_second_order_details

/// @brief first-order optimization
template <typename nlp_settings>
class NlpGPlacerSecondOrder : public NlpGPlacerFirstOrder<nlp_settings> {
public:
  typedef typename NlpGPlacerFirstOrder<nlp_settings>::base_type base_type;
  typedef NlpGPlacerFirstOrder<nlp_settings> first_order_type;

  typedef typename first_order_type::nlp_numerical_type nlp_numerical_type;
  typedef typename first_order_type::nlp_coordinate_type nlp_coordinate_type;

  typedef typename nlp_settings::nlp_second_order_setting_type
      second_order_setting_type;

  typedef typename first_order_type::EigenMatrix EigenMatrix;

  typedef typename first_order_type::nlp_hpwl_type nlp_hpwl_type;
  typedef typename first_order_type::nlp_ovl_type nlp_ovl_type;
  typedef typename first_order_type::nlp_oob_type nlp_oob_type;
  typedef typename first_order_type::nlp_asym_type nlp_asym_type;
  typedef typename first_order_type::nlp_cos_type nlp_cos_type;
  typedef typename first_order_type::nlp_power_wl_type nlp_power_wl_type;

  typedef
      typename second_order_setting_type::hpwl_hessian_trait hpwl_hessian_trait;
  typedef
      typename second_order_setting_type::ovl_hessian_trait ovl_hessian_trait;
  typedef
      typename second_order_setting_type::oob_hessian_trait oob_hessian_trait;
  typedef
      typename second_order_setting_type::asym_hessian_trait asym_hessian_trait;
  typedef
      typename second_order_setting_type::cos_hessian_trait cos_hessian_trait;
  typedef typename second_order_setting_type::power_wl_hessian_trait
      power_wl_hessian_trait;

  /* figure out the types for storing the hessian */
  // Determine whether the operators are return a diagonal hessian
  constexpr static BoolType isHpwlHessianDiagonal =
      diff::is_diagnol_matrix<hpwl_hessian_trait>::value;
  constexpr static BoolType isOvlHessianDiagonal =
      diff::is_diagnol_matrix<ovl_hessian_trait>::value;
  constexpr static BoolType isOobHessianDiagonal =
      diff::is_diagnol_matrix<oob_hessian_trait>::value;
  constexpr static BoolType isAsymHessianDiagonal =
      diff::is_diagnol_matrix<asym_hessian_trait>::value;
  constexpr static BoolType isCosHessianDiagonal =
      diff::is_diagnol_matrix<cos_hessian_trait>::value;
  constexpr static BoolType isPowerWlHessianDiagonal =
      diff::is_diagnol_matrix<power_wl_hessian_trait>::value;
  constexpr static BoolType isHessianDiagonal =
      isHpwlHessianDiagonal and isOvlHessianDiagonal and
      isOobHessianDiagonal and isAsymHessianDiagonal and
      isCosHessianDiagonal and isPowerWlHessianDiagonal;

  // define the supporting trait
  typedef _nlp_second_order_details::is_diagonal_select<nlp_settings,
                                                        isHpwlHessianDiagonal>
      hpwl_hessian_diagonal_selector;
  friend hpwl_hessian_diagonal_selector;
  typedef _nlp_second_order_details::is_diagonal_select<nlp_settings,
                                                        isOvlHessianDiagonal>
      ovl_hessian_diagonal_selector;
  friend ovl_hessian_diagonal_selector;
  typedef _nlp_second_order_details::is_diagonal_select<nlp_settings,
                                                        isOobHessianDiagonal>
      oob_hessian_diagonal_selector;
  friend oob_hessian_diagonal_selector;
  typedef _nlp_second_order_details::is_diagonal_select<nlp_settings,
                                                        isAsymHessianDiagonal>
      asym_hessian_diagonal_selector;
  friend asym_hessian_diagonal_selector;
  typedef _nlp_second_order_details::is_diagonal_select<nlp_settings,
                                                        isCosHessianDiagonal>
      cos_hessian_diagonal_selector;
  friend cos_hessian_diagonal_selector;
  typedef _nlp_second_order_details::is_diagonal_select<
      nlp_settings, isPowerWlHessianDiagonal>
      power_wl_hessian_diagonal_selector;
  friend cos_hessian_diagonal_selector;
  typedef _nlp_second_order_details::is_diagonal_select<nlp_settings,
                                                        isHessianDiagonal>
      hessian_diagonal_selector;
  friend hessian_diagonal_selector;

  typedef
      typename hpwl_hessian_diagonal_selector::matrix_type hpwl_hessian_matrix;
  typedef
      typename ovl_hessian_diagonal_selector::matrix_type ovl_hessian_matrix;
  typedef
      typename oob_hessian_diagonal_selector::matrix_type oob_hessian_matrix;
  typedef
      typename asym_hessian_diagonal_selector::matrix_type asym_hessian_matrix;
  typedef
      typename cos_hessian_diagonal_selector::matrix_type cos_hessian_matrix;
  typedef typename power_wl_hessian_diagonal_selector::matrix_type
      power_wl_hessian_matrix;
  typedef typename hessian_diagonal_selector::matrix_type hessian_matrix;

  /* define the algorithms */
  typedef typename nlp_settings::nlp_second_order_algorithms_type
      nlp_second_order_algorithms;
  typedef typename nlp_second_order_algorithms::converge_type converge_type;
  typedef typename nlp::converge::converge_criteria_trait<converge_type>
      converge_trait;
  typedef typename nlp_second_order_algorithms::optm_type optm_type;
  typedef typename nlp::optm::optm_trait<optm_type> optm_trait;

  friend converge_type;
  template <typename converge_criteria_type>
  friend struct nlp::converge::converge_criteria_trait;
  friend optm_type;
  friend optm_trait;

  typedef
      typename nlp_settings::nlp_second_order_algorithms_type::mult_init_type
          mult_init_type;
  typedef nlp::outer_multiplier::init::multiplier_init_trait<mult_init_type>
      mult_init_trait;
  friend mult_init_trait;
  typedef
      typename nlp_settings::nlp_second_order_algorithms_type::mult_update_type
          mult_update_type;
  typedef nlp::outer_multiplier::update::multiplier_update_trait<
      mult_update_type>
      mult_update_trait;
  friend mult_update_trait;
  typedef typename nlp_settings::nlp_second_order_algorithms_type::mult_type
      mult_type;
  typedef nlp::outer_multiplier::multiplier_trait<mult_type> mult_trait;
  friend mult_trait;

  typedef
      typename nlp_settings::nlp_second_order_algorithms_type::mult_adjust_type
          mult_adjust_type;
  typedef nlp::outer_multiplier::update::multiplier_update_trait<
      mult_adjust_type>
      mult_adjust_trait;
  friend mult_adjust_trait;

  /* updating alpha parameters */
  typedef typename nlp_settings::nlp_second_order_algorithms_type::alpha_type
      alpha_type;
  typedef nlp::alpha::alpha_trait<alpha_type> alpha_trait;
  template <typename T> friend struct nlp::alpha::alpha_trait;
  typedef
      typename nlp_settings::nlp_second_order_algorithms_type::alpha_update_type
          alpha_update_type;
  typedef nlp::alpha::update::alpha_update_trait<alpha_update_type>
      alpha_update_trait;
  template <typename T> friend struct nlp::alpha::update::alpha_update_trait;

  static constexpr nlp_numerical_type hessianMinBound = 0.01;
  static constexpr nlp_numerical_type hessianMaxBound = 10;

  NlpGPlacerSecondOrder(Database &db)
      : NlpGPlacerFirstOrder<nlp_settings>(db) {}

public:
  decltype(auto) inverseHessian() {
    return hessian_diagonal_selector::inverse(_hessian);
  }
  void calcHessian() {
    _clearHessian();
    _calcAllHessians();
    _updateAllHessian();
    clipHessian();
  }

protected:
  virtual void initProblem() override {
    first_order_type::initProblem();
    initSecondOrder();
  }
  void initSecondOrder();
  /* Construct tasks */
  virtual void optimize() override {
    WATCH_QUICK_START();
    // setting up the multipliers
    this->assignIoPins();
    this->_wrapObjAllTask.run();
    this->_wrapCalcGradTask.run();
    calcHessian();

    optm_type optm;
    mult_type multiplier = mult_trait::construct(*this);
    mult_trait::init(*this, multiplier);
    mult_trait::recordRaw(*this, multiplier);

    mult_adjust_type multAdjuster =
        mult_adjust_trait::construct(*this, multiplier);
    mult_adjust_trait::init(*this, multiplier, multAdjuster);

    alpha_type alpha = alpha_trait::construct(*this);
    alpha_trait::init(*this, alpha);
    alpha_update_type alphaUpdate = alpha_update_trait::construct(*this, alpha);
    alpha_update_trait::init(*this, alpha, alphaUpdate);
    std::cout << "nlp address " << this << std::endl;

    IntType iter = 0;
    do {
      std::string debugGdsFilename = "./debug/";
      debugGdsFilename += "gp_iter_" + std::to_string(iter) + ".gds";
      optm_trait::optimize(*this, optm);
      mult_trait::update(*this, multiplier);
      mult_trait::recordRaw(*this, multiplier);
      mult_adjust_trait::update(*this, multiplier, multAdjuster);

      alpha_update_trait::update(*this, alpha, alphaUpdate);
      this->assignIoPins();
      DBG("obj %f hpwl %f ovl %f oob %f asym %f cos %f \n", this->_obj,
          this->_objHpwl, this->_objOvl, this->_objOob, this->_objAsym,
          this->_objCos);
      ++iter;
    } while (not base_type::stop_condition_trait::stopPlaceCondition(
        *this, this->_stopCondition));
    auto end = WATCH_QUICK_END();
    // std::cout<<"grad"<<"\n"<< _grad <<std::endl;
    std::cout << "time " << end / 1000 << " ms" << std::endl;
    this->writeLocs();
  }

private:
  void constructCalcHessianTasks() {
    using hpwl =
        nt::CalculateOperatorHessianTask<nlp_hpwl_type, hpwl_hessian_trait,
                                         EigenMatrix, hpwl_hessian_matrix>;
    using ovl =
        nt::CalculateOperatorHessianTask<nlp_ovl_type, ovl_hessian_trait,
                                         EigenMatrix, ovl_hessian_matrix>;
    using oob =
        nt::CalculateOperatorHessianTask<nlp_oob_type, oob_hessian_trait,
                                         EigenMatrix, oob_hessian_matrix>;
    using asym =
        nt::CalculateOperatorHessianTask<nlp_asym_type, asym_hessian_trait,
                                         EigenMatrix, asym_hessian_matrix>;
    using cos =
        nt::CalculateOperatorHessianTask<nlp_cos_type, cos_hessian_trait,
                                         EigenMatrix, cos_hessian_matrix>;
    using pwl =
        nt::CalculateOperatorHessianTask<nlp_power_wl_type,
                                         power_wl_hessian_trait, EigenMatrix,
                                         power_wl_hessian_matrix>;
    auto getIdxFunc = [&](IndexType cellIdx, Orient2DType orient) {
      return this->plIdx(cellIdx, orient);
    }; // wrapper the convert cell idx to pl idx
    for (IndexType i = 0; i < this->_hpwlOps.size(); ++i) {
      _calcHpwlHessianTasks.emplace_back(
          hpwl(&this->_hpwlOps[i], &_hessianHpwl, getIdxFunc));
    }
    for (auto &op : this->_ovlOps) {
      _calcOvlHessianTasks.emplace_back(ovl(&op, &_hessianOvl, getIdxFunc));
    }
    for (auto &op : this->_oobOps) {
      _calcOobHessianTasks.emplace_back(oob(&op, &_hessianOob, getIdxFunc));
    }
    for (auto &op : this->_asymOps) {
      _calcAsymHessianTasks.emplace_back(asym(&op, &_hessianAsym, getIdxFunc));
    }
    for (auto &op : this->_cosOps) {
      _calcCosHessianTasks.emplace_back(cos(&op, &_hessianCos, getIdxFunc));
    }
    for (auto &op : this->_powerWlOps) {
      _calcPowerWlHessianTasks.emplace_back(
          pwl(&op, &_hessianPowerWl, getIdxFunc));
    }
  }
  void _clearHessian() {
    _hessian.setZero();
    _hessianHpwl.setZero();
    _hessianOvl.setZero();
    _hessianOob.setZero();
    _hessianAsym.setZero();
    _hessianCos.setZero();
    _hessianPowerWl.setZero();
  }
  void _calcAllHessians() {
#pragma omp parallel for schedule(static)
    for (IndexType i = 0; i < _calcHpwlHessianTasks.size(); ++i) {
      _calcHpwlHessianTasks[i].calc();
    }
#pragma omp parallel for schedule(static)
    for (IndexType i = 0; i < _calcOvlHessianTasks.size(); ++i) {
      _calcOvlHessianTasks[i].calc();
    }
#pragma omp parallel for schedule(static)
    for (IndexType i = 0; i < _calcOobHessianTasks.size(); ++i) {
      _calcOobHessianTasks[i].calc();
    }
#pragma omp parallel for schedule(static)
    for (IndexType i = 0; i < _calcAsymHessianTasks.size(); ++i) {
      _calcAsymHessianTasks[i].calc();
    }
#pragma omp parallel for schedule(static)
    for (IndexType i = 0; i < _calcCosHessianTasks.size(); ++i) {
      _calcCosHessianTasks[i].calc();
    }
#pragma omp parallel for schedule(static)
    for (IndexType i = 0; i < _calcPowerWlHessianTasks.size(); ++i) {
      _calcPowerWlHessianTasks[i].calc();
    }
  }
  void _updateAllHessian() {
#pragma omp parallel for schedule(static)
    for (IndexType i = 0; i < 6; ++i) {
      if (i == 0) {
        for (auto &calc : _calcHpwlHessianTasks) {
          calc.update();
        }
      } else if (i == 1) {
        for (auto &calc : _calcOvlHessianTasks) {
          calc.update();
        }
      } else if (i == 2) {
        for (auto &calc : _calcOobHessianTasks) {
          calc.update();
        }
      } else if (i == 3) {
        for (auto &calc : _calcAsymHessianTasks) {
          calc.update();
        }
      } else if (i == 4) {
        for (auto &calc : _calcCosHessianTasks) {
          calc.update();
        }
      } else {
        for (auto &calc : _calcPowerWlHessianTasks) {
          calc.update();
        }
      }
    }
    _hessian = _hessianHpwl + _hessianOvl + _hessianOob + _hessianAsym +
               _hessianCos + _hessianPowerWl;
  }

  void clipHessian() {
    _hessian = _hessian.cwiseMin(hessianMinBound).cwiseMax(hessianMaxBound);
  }
  virtual void constructTasks() override {
    first_order_type::constructTasks();

    this->constructCalcHessianTasks();
  }

protected:
  hessian_matrix _hessian;          ///< The hessian for the objective function
  hpwl_hessian_matrix _hessianHpwl; ///< The hessian for the hpwl function
  ovl_hessian_matrix _hessianOvl; ///< The hessian for the overlapping function
  oob_hessian_matrix
      _hessianOob; ///< The hessian for the out of boundary function
  asym_hessian_matrix _hessianAsym; ///< The hessian for the asymmetry function
  cos_hessian_matrix _hessianCos; ///< The hessian for the signal path function
  power_wl_hessian_matrix _hessianPowerWl;
  /* Tasks */
  std::vector<nt::CalculateOperatorHessianTask<
      nlp_hpwl_type, hpwl_hessian_trait, EigenMatrix, hpwl_hessian_matrix>>
      _calcHpwlHessianTasks; ///< calculate and update the hessian
  std::vector<nt::CalculateOperatorHessianTask<nlp_ovl_type, ovl_hessian_trait,
                                               EigenMatrix, ovl_hessian_matrix>>
      _calcOvlHessianTasks; ///< calculate and update the hessian
  std::vector<nt::CalculateOperatorHessianTask<nlp_oob_type, oob_hessian_trait,
                                               EigenMatrix, oob_hessian_matrix>>
      _calcOobHessianTasks; ///< calculate and update the hessian
  std::vector<nt::CalculateOperatorHessianTask<
      nlp_asym_type, asym_hessian_trait, EigenMatrix, asym_hessian_matrix>>
      _calcAsymHessianTasks; ///< calculate and update the hessian
  std::vector<nt::CalculateOperatorHessianTask<nlp_cos_type, cos_hessian_trait,
                                               EigenMatrix, cos_hessian_matrix>>
      _calcCosHessianTasks; ///< calculate and update the hessian
  std::vector<nt::CalculateOperatorHessianTask<
      nlp_power_wl_type, power_wl_hessian_trait, EigenMatrix,
      power_wl_hessian_matrix>>
      _calcPowerWlHessianTasks; ///< calculate and update the hessian
};

template <typename nlp_settings>
inline void NlpGPlacerSecondOrder<nlp_settings>::initSecondOrder() {
  const IntType size = this->_numVariables;
  hpwl_hessian_diagonal_selector::resize(_hessianHpwl, size);
  ovl_hessian_diagonal_selector::resize(_hessianOvl, size);
  oob_hessian_diagonal_selector::resize(_hessianOob, size);
  asym_hessian_diagonal_selector::resize(_hessianAsym, size);
  cos_hessian_diagonal_selector::resize(_hessianCos, size);
  power_wl_hessian_diagonal_selector::resize(_hessianPowerWl, size);
}

PROJECT_NAMESPACE_END
#endif // IDEAPLACE_NLPGPLACER_H_
