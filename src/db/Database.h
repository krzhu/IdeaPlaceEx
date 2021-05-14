/**
 * @file Database.h
 * @brief The placement database data structure
 * @author Keren Zhu
 * @date 10/02/2019
 */

#ifndef IDEAPLACE_DATABASE_H_
#define IDEAPLACE_DATABASE_H_

#include "Cell.h"
#include "Constraints.h"
#include "Net.h"
#include "Parameters.h"
#include "Pin.h"
#include "Tech.h"

#include <numeric>

PROJECT_NAMESPACE_BEGIN

/// @class IDEAPLACE::Database
/// @brief the database class of the placement engine
class Database {
public:
  /// @brief default database
  explicit Database() = default;
  /*------------------------------*/
  /* Initialization               */
  /*------------------------------*/
  /// @brief initializing cells. Allocate the correct number of layers
  /// @return if successful
  bool initCells();
  /// @brief initializing a cell
  bool initCell(IndexType cellIdx);
  /*------------------------------*/
  /* Getters                      */
  /*------------------------------*/
  /// @brief get the technology-wrapper of the placement
  /// @return the technology
  const Tech &tech() const { return _tech; }
  /// @brief get the technology-wrapper of the placement
  /// @return the technology
  Tech &tech() { return _tech; }
  /// @brief get the placement parameter wrapper
  /// @return the placement parameter wrapper
  const Parameters &parameters() const { return _para; }
  /// @brief get the placement parameter wrapper
  /// @return the placement parameter wrapper
  Parameters &parameters() { return _para; }
  /*------------------------------*/
  /* Vector operations            */
  /*------------------------------*/
  /// @brief allocate a new symmetric group
  /// @return the index of the symmetric group
  IndexType allocateSymGrp() {
    _symGroups.emplace_back(SymGroup());
    return _symGroups.size() - 1;
  }
  /// @brief get the number of symmetric groups
  /// @return the number of the symmetric groups
  IndexType numSymGroups() const { return _symGroups.size(); }
  /// @brief get a symmetric group
  /// @param the index of the symmetric group
  const SymGroup &symGroup(IndexType idx) const { return _symGroups.at(idx); }
  /// @brief get a symmetric group
  /// @param the index of the symmetric group
  SymGroup &symGroup(IndexType idx) { return _symGroups.at(idx); }
  /// @brief get the vector of symmetric groups
  /// @return the vector of symmetric groups
  const std::vector<SymGroup> &vSymGrpArray() const { return _symGroups; }
  /// @brief get the number of cells
  /// @return the number of cells
  IndexType numCells() const { return _cellArray.size(); }
  /// @brief get one cell from the array
  /// @param the index of the cell
  /// @return the cell of the index
  const Cell &cell(IndexType cellIdx) const { return _cellArray.at(cellIdx); }
  /// @brief get one cell from the array
  /// @param the index of the cell
  /// @return the cell of the index
  Cell &cell(IndexType cellIdx) { return _cellArray.at(cellIdx); }
  /// @brief allocate a new cell in the array
  /// @return the index of the new cell
  IndexType allocateCell() {
    _cellArray.emplace_back(Cell());
    return _cellArray.size() - 1;
  }
  /// @brief get the cell array
  /// @return the cell array
  const std::vector<Cell> &vCellArray() const { return _cellArray; }
  /// @brief get the cell array
  /// @return the cell array
  std::vector<Cell> &vCellArray() { return _cellArray; }
  /// @brief set name of cell
  /// @param first: cellIdx
  /// @param second: cellName
  void setCellName(IndexType cellIdx, const std::string name) {
    _cellArray.at(cellIdx).setName(name);
  }
  /// @brief get the number of nets
  /// @param the index of the nets
  /// @return the net of the index
  IndexType numNets() const { return _netArray.size(); }
  /// @brief get one net from the array
  /// @param the index of the net
  /// @return the net of the index
  const Net &net(IndexType netIdx) const { return _netArray.at(netIdx); }
  /// @brief get one net from the array
  /// @param the index of the net
  /// @return the net of the index
  Net &net(IndexType netIdx) { return _netArray.at(netIdx); }
  /// @brief allocate a new net in the array
  /// @return index of the new net
  IndexType allocateNet() {
    _netArray.emplace_back(Net());
    return _netArray.size() - 1;
  }
  const std::vector<Net> &nets() const { return _netArray; }
  std::vector<Net> &nets() { return _netArray; }
  /// @brief get the number of pins
  /// @return the number of pins
  IndexType numPins() const { return _pinArray.size(); }
  /// @brief get one pin from the array
  /// @param the index of the pin
  /// @return the pin of the index
  const Pin &pin(IndexType pinIdx) const { return _pinArray.at(pinIdx); }
  /// @brief get one pin from the array
  /// @param the index of the pin
  /// @return the pin of the index
  Pin &pin(IndexType pinIdx) { return _pinArray.at(pinIdx); }
  /// @brief allocate a new pin in the array
  /// @return the index of the pin
  IndexType allocatePin() {
    _pinArray.emplace_back(Pin());
    return _pinArray.size() - 1;
  }
  const std::vector<Pin> &vPinArray() const { return _pinArray; }
  std::vector<Pin> &vPinArray() { return _pinArray; }
  /// @brief get the vector of proximity groups
  /// @return the vector of proximity groups
  const std::vector<ProximityGroup> &proximityGrps() const {
    return _proximityGrps;
  }
  /// @brief get the vector of proximity groups
  /// @return the vector of proximity groups
  std::vector<ProximityGroup> &proximityGrps() { return _proximityGrps; }
  IndexType allocateProximityGroup() {
    _proximityGrps.emplace_back(ProximityGroup());
    return _proximityGrps.size() - 1;
  }
  const ProximityGroup &proximityGrp(IndexType grpIdx) const {
    return _proximityGrps.at(grpIdx);
  }
  ProximityGroup &proximityGrp(IndexType grpIdx) {
    return _proximityGrps.at(grpIdx);
  }

  const std::vector<SignalPath> &vSignalPaths() const { return _signalPaths; }
  std::vector<SignalPath> &vSignalPaths() { return _signalPaths; }
  IndexType allocateSignalPath() {
    _signalPaths.emplace_back(SignalPath());
    return _signalPaths.size() - 1;
  }
  SignalPath &signalPath(IndexType idx) { return _signalPaths.at(idx); }
  const SignalPath &signalPath(IndexType idx) const {
    return _signalPaths.at(idx);
  }
  /// @brief split the signal path to seperate the sym pairs
  void splitSignalPathsBySymPairs();
  /// @brief add a new pin to a signal path
  /// @param first: the index of signal path
  /// @param second: the name of cell
  /// @param third: the name of pin
  void addPinToSignalPath(IndexType pathIdx, const std::string &cellName,
                          const std::string &pinName) {
    auto &sig = this->signalPath(pathIdx);
    IndexType pinIdx = INDEX_TYPE_MAX;
    IndexType cellIdx = INDEX_TYPE_MAX;
    for (IndexType idx = 0; idx < this->numCells(); ++idx) {
      const auto &cell = this->cell(idx);
      if (cell.name() == cellName) {
        cellIdx = idx;
        for (IndexType jdx = 0; jdx < cell.numPinIdx(); ++jdx) {
          const auto &pin = this->pin(cell.pinIdx(jdx));
          if (pin.name() == pinName) {
            pinIdx = cell.pinIdx(jdx);
            break;
          }
        }
        break;
      }
    }
    AssertMsg(cellIdx != INDEX_TYPE_MAX,
              "Unknown pin! cellname: %s, pinname %s \n", cellName.c_str(),
              pinName.c_str());
    if (pinIdx == INDEX_TYPE_MAX) {
      WRN("%s unknown pin name. cell name %s pin name %s \n", __FUNCTION__,
          cellName.c_str(), pinName.c_str());
    }
    sig.addPinIdx(pinIdx);
  }

  // Well
  const std::vector<Well> &vWells() const { return _wellArray; }
  Well &well(const IndexType i) { return _wellArray.at(i); }
  const Well &well(const IndexType i) const { return _wellArray.at(i); }

  IndexType allocateWell() {
    Well w(_wellArray.size());
    _wellArray.emplace_back(w);
    return w.idx();
  }
  bool removeWell(const Well &w) {
    return _wellArray.erase(
               std::remove(_wellArray.begin(), _wellArray.end(), w),
               _wellArray.end()) != _wellArray.end();
  }
  void clearWells() { 
    _wellArray.clear(); 
    for (auto &cell : _cellArray) {
      cell.clearWell();
    }
    _wellRectSize.clear();
  }
  void assignCellToWellAndRemoveUnusedWell();
  /// @brief split the well polygon into rectangles
  void splitWells(IntType mode = 0) {
    _wellRectSize.clear();
    for (auto &well : _wellArray) {
      well.splitIntoRects(mode);
      _wellRectSize.emplace_back(well.rects().size());
    }
  }
  IndexType numWellRects() const {
    return std::accumulate(_wellRectSize.begin(), _wellRectSize.end(), 0);
  }
  /// @return well index, rect index in the well
  std::pair<IndexType, IndexType> getWellRectIdx(IndexType rectIdx) const {
    for (IndexType wellIdx = 0; wellIdx < _wellRectSize.size(); ++wellIdx) {
      if (rectIdx < _wellRectSize[wellIdx]) {
        return std::make_pair(wellIdx, rectIdx);
      }
      rectIdx -= _wellRectSize[wellIdx];
    }
    Assert(false);
    return std::make_pair(INDEX_TYPE_MAX, INDEX_TYPE_MAX);
  }
  IndexType numRectInWell(IndexType wellIdx) const {
    return _wellRectSize.at(wellIdx);
  }

  /*------------------------------*/
  /* Query function wrappers      */
  /*------------------------------*/
  /// @brief get the pin offset with regarded to the ll corner of cell
  Point<LocType> pinOffsetToCell(IndexType pinIdx) const {
    const auto &pin = this->pin(pinIdx);
    IndexType cellIdx = pin.cellIdx();
    const auto &cell = this->cell(cellIdx);
    // Get the cell location from the input arguments
    if (_para.ifUseRealPinLoc()) {
      return pin.midLoc() - cell.cellBBox().ll();
    } else {
      // Use the cell mid location
      return cell.cellBBox().center();
    }
  }
  /*------------------------------*/
  /* Technology-dependent         */
  /*------------------------------*/
  /// @brief get the spacing requirement between two cells
  /// @param index for cell 1
  /// @param index for cell 2
  /// @return A box object representing the spacing requirements. xyLo for cell
  /// 1 is left and bottm. xyHi for ... 2...
  Box<LocType> cellSpacing(IndexType cellIdx1, IndexType cellIdx2) const;
  /// @brief calculate the spacing between cells
  void calculateCellSpacings();
  /// @brief Get the reqired WPE spacing
  /// @param cell index
  /// @return Pair. .first=horizontal spacing .second=vertical spacing
  std::pair<LocType, LocType> wpeSpacing(IndexType cellIdx) {
    const LocType ver = _tech.wpeVerticalSpacing(cell(cellIdx).fingerChannelWidth());
    const LocType hor = _tech.wpeHorizontalSpacing(cell(cellIdx).fingerChannelLength());
    return std::make_pair(hor, ver);
  }
  /*------------------------------*/
  /* Supporting functions         */
  /*------------------------------*/
  /// @brief calcuate the total cell area
  /// @return the total cell area
  RealType calculateTotalCellArea() const;
  /// @brief calculate and return the HPWL
  /// @return HPWL
  LocType hpwl() const;
  LocType hpwlWithVitualPins() const;
  void expandCellToGridSize(LocType gridSize) {
    for (auto &cell : _cellArray) {
      cell.forceExtendToGrid(gridSize);
    }
  }
  RealType area() const;
  bool checkSym();
  LocType findSymAxis();
  /// @brief Move the layout so that everything is located >= 0
  void offsetLayout();
  /// @brief individual well mode, only for illustration
  void individualWell() {
    for (IndexType cellIdx = 0; cellIdx < numCells(); ++cellIdx) {
      if (not cell(cellIdx).needWell()) {
        continue;
      }
      cell(cellIdx).calculateCellBBox();
      auto wpe = wpeSpacing(cellIdx);
      Box<LocType> bbox = cell(cellIdx).cellBBox();
      bbox.expandX(std::max(wpe.first, 500 + 340 *2)); // FIXME: typical guard ring well extension. Only for experimental comparsion 
      bbox.expandY(std::max(wpe.second, 500 + 340 *2)); // 500 + 340*2 align with other experimental setting (guard ring + required spacing)
      cell(cellIdx).unionBBox(_tech.nwellLayerIdx(), bbox);
      cell(cellIdx).calculateCellBBox();
    }
  }
  /*------------------------------*/
  /* Debug functions              */
  /*------------------------------*/
#ifdef DEBUG_DRAW
  /// @brief draw the cell blocks
  /// @param the filename for saving
  void drawCellBlocks(const std::string &name);
#endif // DEBUG_DRAW

private:
  std::vector<Cell> _cellArray;     ///< The cells of the placement problem
  std::vector<Net> _netArray;       ///< The nets of the placement problem
  std::vector<Pin> _pinArray;       ///< The pins of the placement problem
  std::vector<SymGroup> _symGroups; ///< The symmetric groups
  std::vector<ProximityGroup>
      _proximityGrps;                   ///< The proximity group constraints
  std::vector<SignalPath> _signalPaths; ///< The signal/current paths
  std::vector<Well> _wellArray;         ///< The wells for PMOS
  std::vector<IndexType> _wellRectSize; ///< Record the number of rectangles in each well

  Tech _tech;       ///< The tech information
  Parameters _para; ///< The parameters for the placement engine
};

inline RealType Database::calculateTotalCellArea() const {
  RealType area = 0;
  for (const auto &cell : _cellArray) {
    area += cell.cellBBox().area();
  }
  return area;
}

inline void Database::calculateCellSpacings() {
  for (IndexType cellIdx  = 0; cellIdx < numCells(); ++cellIdx) {
    // Each cell only record the spacing for the cells after in index, size numCells() - cellIdx - 1
    // The last one is against N-well
    // In total: numCells() - cellIdx
    _cellArray.at(cellIdx).spacingToCells().resize(numCells() - cellIdx);
  }
  for (IndexType cellIdx1 = 0; cellIdx1 < numCells(); ++cellIdx1) {
    auto &cell1 = _cellArray.at(cellIdx1);
    for (IndexType cellIdx2 = cellIdx1 + 1; cellIdx2 < numCells(); ++cellIdx2) {
      auto &cell2 = _cellArray.at(cellIdx2);
      auto &spacings = cell1.spacingToCells().at(cellIdx2 - cellIdx1 - 1);
      spacings = Box<LocType>(0, 0, 0, 0);
      for (IndexType layer1 = 0; layer1 < _tech.numLayers(); ++layer1) {
        for (IndexType layer2 = 0; layer2 < _tech.numLayers(); ++layer2) {
          // Check spacing for layer1 in cell1 against layer2 in cell2
          bool cell1HasLayer = cell1.layerHasShape(layer1);
          bool cell2HasLayer = cell2.layerHasShape(layer2);
          if (!(cell1HasLayer && cell2HasLayer)) {
            continue;
          }
          LocType spacingRule = _tech.spacingRule(layer1, layer2);
          // Calculate the layer shape to the cell boundry
          LocType cell1ToBoundXLo =
              cell1.bbox(layer1).xLo() - cell1.cellBBox().xLo();
          LocType cell1ToBoundXHi =
              cell1.cellBBox().xHi() - cell1.bbox(layer1).xHi();
          LocType cell1ToBoundYLo =
              cell1.bbox(layer1).yLo() - cell1.cellBBox().yLo();
          LocType cell1ToBoundYHi =
              cell1.cellBBox().yHi() - cell1.bbox(layer1).yHi();
          LocType cell2ToBoundXLo =
              cell2.bbox(layer2).xLo() - cell2.cellBBox().xLo();
          LocType cell2ToBoundXHi =
              cell2.cellBBox().xHi() - cell2.bbox(layer2).xHi();
          LocType cell2ToBoundYLo =
              cell2.bbox(layer2).yLo() - cell2.cellBBox().yLo();
          LocType cell2ToBoundYHi =
              cell2.cellBBox().yHi() - cell2.bbox(layer2).yHi();
          // cell 1 is left
          LocType xLo = -cell1ToBoundXHi + spacingRule - cell2ToBoundXLo;
          // cell 1 is lower
          LocType yLo = -cell1ToBoundYHi + spacingRule - cell2ToBoundYLo;
          // cell 1 is right
          LocType xHi = -cell1ToBoundXLo + spacingRule - cell2ToBoundXHi;
          // cell 1 is higher
          LocType yHi = -cell1ToBoundYLo + spacingRule - cell2ToBoundYHi;
          // Update the cell-wise spacing
          spacings.setXLo(std::max(spacings.xLo(), xLo));
          spacings.setYLo(std::max(spacings.yLo(), yLo));
          spacings.setXHi(std::max(spacings.xHi(), xHi));
          spacings.setYHi(std::max(spacings.yHi(), yHi));
        }
      }
    }
  }
  // Finally, calculate the spacing for each cell against nwell
  for (IndexType cellIdx = 0; cellIdx < numCells(); ++cellIdx) {
    auto &cell1 = this->cell(cellIdx);
    auto &spacings = cell1.spacingToCells().back(); // Record in the last one
    spacings = Box<LocType>(0, 0, 0, 0);
    if (not _tech.isNwellLayerSet()) {
      continue;
    }
    for (IndexType layerIdx = 0; layerIdx < _tech.numLayers(); ++layerIdx) {
      bool cell1HasLayer = cell1.layerHasShape(layerIdx);
      if (!(cell1HasLayer)) {
        continue;
      }
      LocType spacingRule = _tech.spacingRule(layerIdx, _tech.nwellLayerIdx());
      // Calculate the layer shape to the cell boundry
      LocType cell1ToBoundXLo =
          cell1.bbox(layerIdx).xLo() - cell1.cellBBox().xLo();
      LocType cell1ToBoundXHi =
          cell1.cellBBox().xHi() - cell1.bbox(layerIdx).xHi();
      LocType cell1ToBoundYLo =
          cell1.bbox(layerIdx).yLo() - cell1.cellBBox().yLo();
      LocType cell1ToBoundYHi =
          cell1.cellBBox().yHi() - cell1.bbox(layerIdx).yHi();
      // cell 1 is left
      LocType xLo = -cell1ToBoundXHi + spacingRule;
      // cell 1 is lower
      LocType yLo = -cell1ToBoundYHi + spacingRule;
      // cell 1 is right
      LocType xHi = -cell1ToBoundXLo + spacingRule;
      // cell 1 is higher
      LocType yHi = -cell1ToBoundYLo + spacingRule;
      // Update the cell-wise spacing
      spacings.setXLo(std::max(spacings.xLo(), xLo));
      spacings.setYLo(std::max(spacings.yLo(), yLo));
      spacings.setXHi(std::max(spacings.xHi(), xHi));
      spacings.setYHi(std::max(spacings.yHi(), yHi));
    }
  }
}

inline Box<LocType> Database::cellSpacing(IndexType cellIdx1,
                                          IndexType cellIdx2) const {
  if (cellIdx1 > cellIdx2) {
    std::swap(cellIdx1, cellIdx2);
  }
  if (cellIdx1 >= numCells() and cellIdx2 >= numCells()) {
    // They are both referring to NW
    if (not _tech.isNwellLayerSet()) {
      return Box<LocType>(0, 0, 0, 0);
    }
    // If they are in the same well
    const auto rectIdx1 = getWellRectIdx(cellIdx1 - numCells());
    const auto rectIdx2 = getWellRectIdx(cellIdx2 - numCells());
    if (rectIdx1.first == rectIdx2.first) {
      return Box<LocType>(0, 0, 0, 0);
    }
    // If they are in the different wells
    
    LocType spacing = _tech.spacingRule(_tech.nwellLayerIdx());
    return Box<LocType>(spacing, spacing, spacing, spacing);
  }
  else if (cellIdx2 >= numCells()) {
    // idx1 refer to a cell, idx2 refer to nwell
    return cell(cellIdx1).spacingToCells().back();
  }
  else {
    return cell(cellIdx1).spacingToCells().at(cellIdx2 - cellIdx1 - 1);
  }
}

PROJECT_NAMESPACE_END

#endif /// IDEAPLACE_DATABASE_H_
