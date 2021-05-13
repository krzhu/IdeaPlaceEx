/**
 * @file Tech.h
 * @brief The placement technology-related data structure
 * @author Keren Zhu
 * @date 09/30/2019
 */

#ifndef IDEAPLACE_TECH_H_
#define IDEAPLACE_TECH_H_

#include "global/global.h"
#include "util/Vector2D.h"
#include "util/box.hpp"
#include <unordered_map>

PROJECT_NAMESPACE_BEGIN

class WpeSpacingRule {
  public:
    explicit WpeSpacingRule() = default;
    /// @brief Add vertical spacing rule. If device width > width, it need spacing
    void addVerticalSpacingRule(LocType width, LocType spacing) {
      _fingerWidth.emplace_back(width);
      _verSpacing.emplace_back(spacing);
      if (_fingerWidth.size() >= 2) {
        if (_fingerWidth.at(_fingerWidth.size() - 2) > _fingerWidth.back()) {
          AssertMsg(false, "Check WPE spacing setting. Channel width is not in ascending order.");
        }
        if (_verSpacing.at(_verSpacing.size() - 2) > _verSpacing.back()) {
          AssertMsg(false, "Check WPE spacing setting. Longer width results in lower spacing. ");
        }
      }
    }
    /// @brief Add horizontal spacing rule. If device length > length, it need spacing
    void addHorizontalSpacingRule(LocType length, LocType spacing) {
      _fingerLength.emplace_back(length);
      _horSpacing.emplace_back(spacing);
      if (_fingerLength.size() >= 2) {
        if (_fingerLength.at(_fingerLength.size() - 2) > _fingerLength.back()) {
          AssertMsg(false, "Check WPE spacing setting. Channel length is not in ascending order.");
        }
        if (_horSpacing.at(_horSpacing.size() - 2) > _horSpacing.back()) {
          AssertMsg(false, "Check WPE spacing setting. Longer length results in lower spacing. ");
        }
      }
    }
    /// @brief get the required vertical spacing. 
    /// @param finger channel width of the device
    LocType verticalSpacing(LocType width) const {
      AssertMsg(width >= 0, "Check device width. Maybe it wasn't set correctely?");
      for (IndexType idx = 0; idx < _fingerWidth.size(); ++idx) {
        if (_fingerWidth.at(idx) > width) {
          AssertMsg(idx > 0, "Check WPE spacing setting. Include width=0 case.");
          return _verSpacing.at(idx - 1);
        }
      }
      return _verSpacing.back();
    }
    /// @brief get the required horizontal spacing. 
    /// @param finger channel length of the device
    LocType horizontalSpacing(LocType length) const {
      AssertMsg(length >= 0, "Check device length. Maybe it wasn't set correctely?");
      for (IndexType idx = 0; idx < _fingerLength.size(); ++idx) {
        if (_fingerLength.at(idx) > length) {
          AssertMsg(idx > 0, "Check WPE spacing setting. Include length=0 case.");
          return _horSpacing.at(idx - 1);
        }
      }
      return _horSpacing.back();
    }
  private:
    std::vector<LocType> _fingerWidth; ///< Finger width determine vertical spacing
    std::vector<LocType> _verSpacing; 
    std::vector<LocType> _fingerLength;
    std::vector<LocType> _horSpacing;
};

/// @brief VDD contact templates
class VddContactRule {
  public:
    explicit VddContactRule() = default;
    /// @brief Set the required spacing rule (from device to contact)
    void setRequiredSpacing(LocType requiredSpacing) { _requiredSpacing = requiredSpacing; }
    /// @brief Get the required spacing rule (from device to contact)
    LocType requiredSpacing() const { return _requiredSpacing; }
    /// @brief Add a contact shape candidate
    void addContactTemplate(const Box<LocType> &templ, IntType weight) { _contactTemplates.emplace_back(templ); _contactWeights.emplace_back(weight); }
    /// @brief Add a contact shape candidate
    void addContactTemplate(LocType xLo, LocType yLo, LocType xHi, LocType yHi, IntType weight) { addContactTemplate(Box<LocType>(xLo, yLo, xHi, yHi), weight); }
    /// @brief Get the number of templates available
    IndexType numContactTemplates() const { return _contactTemplates.size(); }
    /// @brief Get a template
    const Box<LocType> & contractTemplate(IndexType idx) const { return _contactTemplates.at(idx); }
    /// @brief Get the weight of a template
    IntType contactWeight(IndexType idx) const { return _contactWeights.at(idx); }
  private:
    std::vector<Box<LocType>> _contactTemplates;
    std::vector<IntType> _contactWeights; ///< The weight for each contact. Higher the preferred
    LocType _requiredSpacing = 0;
};

class Tech {
public:
  /// @brief default constructor
  explicit Tech() = default;
  /*------------------------------*/
  /* Construct the tech           */
  /*------------------------------*/
  /// @brief add a gds layer
  /// @param a gds layer ID
  /// @return the placement layer index
  IndexType addGdsLayer(IndexType gdsLayer) {
    IndexType layerIdx = _gdsLayerIdxArray.size();
    _gdsLayerIdxArray.emplace_back(gdsLayer);
    _layerIdxMap[gdsLayer] = layerIdx;
    return layerIdx;
  }
  /// @brief initialize the rule data structures
  void initRuleDataStructure();
  /// @brief set the width rule for a layer
  /// @param first: the layer index
  /// @param second: the width rule for the layer
  void setWidthRule(IndexType layerIdx, LocType width) {
    _widthRule.at(layerIdx) = width;
  }
  /// @brief set the area rule for a layer
  /// @param first: the layer index
  /// @param second: the area rule for the layer
  void setAreaRule(IndexType layerIdx, LocType area) {
    _areaRule.at(layerIdx) = area;
  }
  /// @brief set the spacing rule for a intra-layer
  /// @param first: the layer
  void setSpacingRule(IndexType layerIdx, LocType spacing) {
    _spacingRule.at(layerIdx, layerIdx) = spacing;
  }
  /// @brief set the spacing for inter-layer case
  /// @param first: the first layer
  /// @param second: the second layer
  /// @param third: the spacing required
  /// For now, just make two layers interchangable
  void setSpacingRule(IndexType layerIdxFirst, IndexType layerIdxSecond,
                      LocType spacing) {
    _spacingRule.at(layerIdxFirst, layerIdxSecond) = spacing;
    _spacingRule.at(layerIdxSecond, layerIdxFirst) = spacing;
  }
  /// @brief set the dbu/ database unit
  /// @param the databse unit
  void setDbu(LocType dbu) { _dbu = dbu; }
  /// @brief Set the n-well layer index
  /// @param The layer index of the N-well layer
  void setNwellLayerIdx(IndexType nwellLayerIdx) { _nwellLayerIdx = nwellLayerIdx; } 
  /// @brief Set the required vdd contact spacing
  void setVddContactRequiredSpacing(LocType require) { _contactRule.setRequiredSpacing(require); }
  /// @brief Add a VDD contact template
  void addVddContactTemplate(LocType xLo, LocType yLo, LocType xHi, LocType yHi, IntType weight) { _contactRule.addContactTemplate(xLo, yLo, xHi, yHi, weight); }
  /// @brief Add ad VDD contact template
  void addVddContactTemplate(const Box<LocType> &box, IntType weight) { _contactRule.addContactTemplate(box, weight); }
  /// @brief Add WPE vertical spacing rule. If device width > width, it need spacing
  void addWpeVerticalSpacingRule(LocType width, LocType spacing) {
    _wpe.addVerticalSpacingRule(width, spacing);
  }
  /// @brief Add WPE horizontal spacing rule. If device length > length, it need spacing
  void addWpeHorizontalSpacingRule(LocType length, LocType spacing) {
    _wpe.addHorizontalSpacingRule(length, spacing);
  }
  /*------------------------------*/
  /* Query the tech               */
  /*------------------------------*/
  /// @brief get the layer index from GDS techlayer
  /// @param the GDS techlayer ID
  /// @return the corresponding layer index in the placement
  IndexType gdsLayerToLayerIdx(IndexType gdsLayer) const {
    return _layerIdxMap.at(gdsLayer);
  }
  /// @brief get whether there is width rule in the given layer
  /// @param the layer index
  /// @return where there is width rule defined in the given layer
  bool hasWidthRule(IndexType layerIdx) const {
    return _widthRule.at(layerIdx) != -1;
  }
  /// @brief get the width rule of the given layer
  /// @param the layer index
  /// @return the width rule
  LocType widthRule(IndexType layerIdx) const {
    Assert(this->hasWidthRule(layerIdx));
    return _widthRule.at(layerIdx);
  }
  /// @brief get whether there is area rule defined in the given layer
  /// @param the layer index
  /// @return where there is area rule defined in the given layer
  bool hasAreaRule(IndexType layerIdx) const {
    return _areaRule.at(layerIdx) != -1;
  }
  /// @brief get the area rule of the given layer
  /// @param the layer index
  /// @return the area rule
  LocType areaRule(IndexType layerIdx) const {
    Assert(this->hasAreaRule(layerIdx));
    return _areaRule.at(layerIdx);
  }
  /// @brief get the spacing rule of intra-layer
  /// @param the layer index
  /// @return the intra-layer spacing rule
  LocType spacingRule(IndexType layerIdx) const {
    return _spacingRule.at(layerIdx, layerIdx);
  }
  /// @brief get the inter-spacing rule
  /// @param first: the first layer index
  /// @param second: the second layer index
  /// @return the spacing rule defined
  LocType spacingRule(IndexType layerIdxFirst, IndexType layerIdxSecond) const {
    return _spacingRule.at(layerIdxFirst, layerIdxSecond);
  }
  /// @brief get the database unit
  /// @return the database unit
  LocType dbu() const { return _dbu; }
  /// @brief get the number of layers
  /// @return the number of layers
  IndexType numLayers() const { return _gdsLayerIdxArray.size(); }
  /// @brief Get the required spacing for VDD contact
  LocType vddContactSpacing() const { return _contactRule.requiredSpacing(); }
  /// @brief Get number of VDD contact templates available
  IndexType numVddContactTemplates() const { return _contactRule.numContactTemplates(); }
  /// @brief Get a vdd contact template
  const Box<LocType> & vddContactTemplate(IndexType idx) const { return _contactRule.contractTemplate(idx); }
  /// @brief Get a vdd contact weight
  IntType vddContactWeight(IndexType idx) const { return _contactRule.contactWeight(idx); }
  /// @brief get the required vertical spacing. 
  /// @param finger channel width of the device
  LocType wpeVerticalSpacing(LocType width) const {
    return _wpe.verticalSpacing(width);
  }
  /// @brief get the required horizontal spacing. 
  /// @param finger channel length of the device
  LocType wpeHorizontalSpacing(LocType length) const {
    return _wpe.horizontalSpacing(length);
  }
  /*------------------------------*/
  /* Getters                      */
  /*------------------------------*/
  /// @brief get the map from gds techlayer to IDEAPLACE layer
  /// @return the map from gds techlayer to IDEAPLACE layer
  const std::unordered_map<IndexType, IndexType> &layerIdxMap() const {
    return _layerIdxMap;
  }
  /// @brief get the map from gds techlayer to IDEAPLACE layer
  /// @return the map from gds techlayer to IDEAPLACE layer
  std::unordered_map<IndexType, IndexType> &layerIdxMap() {
    return _layerIdxMap;
  }
  /// @brief Get if the nwell layer index is set
  bool isNwellLayerSet() const { return _nwellLayerIdx != INDEX_TYPE_MAX; }
  /// @brief Get the n-well layer index
  IndexType nwellLayerIdx() const { return _nwellLayerIdx; }

private:
  LocType _dbu = 1000; ///< 1 um = _dbu database units
  std::unordered_map<IndexType, IndexType>
      _layerIdxMap; ///< A map from gds layer to IDEAPLACE layer.
                    ///< _layerIdxMap[techlayer] = layer index
  std::vector<IndexType>
      _gdsLayerIdxArray; ///< The tech layer in GDS. _gdsLayerIdx[idx of layer]
                         ///< = techlayer in GDS
  std::vector<LocType> _widthRule; ///< The width rule of the layer
  std::vector<LocType> _areaRule;  ///< The area rule of the layer
  Vector2D<LocType>
      _spacingRule; ///< The intra/inter-layer spacing rules.
                    ///< _spacingRule[layer1][layer2] = spacing rules btween
                    ///< shapes of layer 1 and layer 2.
  IndexType _nwellLayerIdx = INDEX_TYPE_MAX;
  VddContactRule _contactRule; ///< The rules for VDD contact
  WpeSpacingRule _wpe; ///< The rule for WPE spacing
};

inline void Tech::initRuleDataStructure() {
  IndexType numLayers = _gdsLayerIdxArray.size();
  _widthRule.resize(numLayers, -1);
  _areaRule.resize(numLayers, -1);
  _spacingRule.resize(numLayers, numLayers, 0);
}

PROJECT_NAMESPACE_END

#endif /// IDEAPLACE_TECH_H_
