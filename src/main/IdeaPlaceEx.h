/**
 * @file IdeaPlaceEx.h
 * @brief wrapper of Everything
 * @author Keren Zhu
 * @date 10/02/2019
 */

#ifndef IDEAPLACE_IDEAPLACEEX_H_
#define IDEAPLACE_IDEAPLACEEX_H_

#include "db/Database.h"
/* Solver */
#include "place/CGLegalizer.h"
#include "place/NlpGPlacer.h"
#include "place/ProximityMgr.h"

PROJECT_NAMESPACE_BEGIN

/// @class IDEAPLACE::IdeaPlaceEx
/// @brief the main wrapper for the placement engine
class IdeaPlaceEx {
public:
  using GPlacer = NlpGPlacerFirstOrder<nlp::nlp_default_settings>;
  /// @brief default constructor
  explicit IdeaPlaceEx() = default;
  /// @brief the file-based input
  /// @param the system arguments
  /// @return if the parsing is successful
  bool parseFileBased(int argc, char **argv);
  /// @brief run the placement algorithm
  /// @return whether the placement is successful
  LocType solve(LocType gridSize = -1);
  /// @brief prepare input before gp. Should not be used together with solve()
  void prepareGp();
  /// @brief the file-based output
  /// @param the system arguments
  /// @return if the writing is successful
  bool outputFileBased(int argc, char **argv);
  /// @brief end the placement, should not called together with solve()
  LocType endPlace();
  /// @brief align the placement to grid
  /// @param grid stepsize
  LocType alignToGrid(LocType gridStepSize);
  /// @brief initialize an global placer
  /// @return A global placer object
  GPlacer initGlobalPlacer();
  CGLegalizer initLegalizer();
  /*------------------------------*/
  /* File-based input interface   */
  /*------------------------------*/
  /// @brief read in tech simple file
  /// @param the filename for the tech simple file
  void readTechSimpleFile(const std::string &techFileSimple);
  /// @brief read pin file
  /// @param the filename for the pin file
  void readPinFile(const std::string &pinFile);
  /// @brief read connection file
  /// @param the connection file
  void readConnectionFile(const std::string &connectionFile);
  /// @brief read net weight file
  /// @param net weight file
  void readNetWgtFile(const std::string &netWgtFile);
  /// @brief read sym file
  /// @param sym file
  void readSymFile(const std::string &symFile);
  /// @brief read a gds for a cell
  /// @param first: the filename for gds layout file
  /// @param second: the cell index. if INDEX_TYPE_MAX, then will try to match
  /// the gds cellname
  void readGdsLayout(const std::string &gdsFile,
                     IndexType cellIdx = INDEX_TYPE_MAX);
  /// @brief read the symnet file
  /// @param the filename for the symnet file
  void readSymNetFile(const std::string &symnetFile);
  /// @brief read the sigpath file
  /// @param the filename for the sigpath file
  void readSigpathFile(const std::string &sigpathFile);
  /*------------------------------*/
  /* paramters                    */
  /*------------------------------*/
  void setNumThreads(IndexType numThreads);
  void setIoPinBoundaryExtension(LocType ext) {
    _db.parameters().setVirtualBoundaryExtension(ext);
  }
  void setIoPinInterval(LocType interval) {
    _db.parameters().setVirtualPinInterval(interval);
  }
  /// @brief set the parameter to use the real pin location
  void openRealPinLocation() { _db.parameters().openRealPinLocation(); }
  /// @brief set the parameter to use the mid of cell as pin location
  void closeRealPinLocation() { _db.parameters().closeRealPinLocation(); }
  /// @brief open fast mode. Notice that the constraints might be violated
  void openFastMode() { _db.parameters().openFastMode(); }
  /// @brief close fast mode.
  void closeFastMode() { _db.parameters().closeFastMode(); }
  /// @brief set the boundary constraint.
  void setBoundaryConstraint(LocType xLo, LocType yLo, LocType xHi,
                             LocType yHi) {
    _db.parameters().setBoundaryConstraint(xLo, yLo, xHi, yHi);
  }
  void setGridStep(LocType gridStep) {
    if (gridStep > 0) {
      _db.parameters().setGridStep(gridStep);
    }
  }
  /*------------------------------*/
  /* tech input interface         */
  /*------------------------------*/
  /// @brief add gds layer
  /// @param gds layer index
  /// @return the placement layer index
  IndexType addGdsLayer(IndexType gdsLayer) {
    return _db.tech().addGdsLayer(gdsLayer);
  }
  /// @brief finish adding gds layers
  void finishAddingGdsLayer() { _db.tech().initRuleDataStructure(); }
  /// @brief Set a same layer spacing rules
  void setIntraLayerSpacing(IndexType layerIdx, LocType spacing) {
    _db.tech().setSpacingRule(layerIdx, spacing);
  }
  /// @brief Set a inter-layer spacing rules
  /// @param first: first layer. The order of layers does not matter.
  /// @param second: second layer. The order of layers does not matter.
  /// @param third: spacing requirement
  void setInterLayerSpacing(IndexType firstLayer, IndexType secondLayer, LocType spacing) {
    _db.tech().setSpacingRule(firstLayer, secondLayer, spacing);
  }
  /// @brief Set the N-well layer index
  void setNwellLayerIdx(IndexType layerIdx) {
    _db.tech().setNwellLayerIdx(layerIdx);
  }

  /*------------------------------*/
  /* Standard input interface     */
  /*------------------------------*/
  /// @brief add a new cell
  /// @return the index for that cell
  IndexType allocateCell() {
    IndexType cellIdx = _db.allocateCell();
    _db.initCell(cellIdx);
    return cellIdx;
  }
  /// @brief set this cell as it needs well
  /// @param the index of the cell
  void flagCellAsNeedWell(IndexType cellIdx) { _db.cell(cellIdx).setNeedWell(true); }
  /// @brief Set the channel width of a cell
  void setCellFingerChannelWidth(IndexType cellIdx, LocType width) { _db.cell(cellIdx).setFingerChannelWidth(width); }
  /// @brief Set the channel length of a cell
  void setCellFingerChannelLength(IndexType cellIdx, LocType length) { _db.cell(cellIdx).setFingerChannelLength(length); }
  /// @brief add a new pin
  /// @param the index of cell this pin belonging to
  /// @return the index of the pin allocated
  IndexType allocatePin(IndexType cellIdx) {
    auto pinIdx = _db.allocatePin();
    _db.pin(pinIdx).setCellIdx(cellIdx);
    _db.cell(cellIdx).addPin(pinIdx);
    return pinIdx;
  }
  /// @brief set the name of a pin
  /// @param the index of the pin
  /// @param the name of the pin
  void setPinName(IndexType pinIdx, const std::string &name) {
    _db.pin(pinIdx).setName(name);
  }
  /// @brief add pin shape for a pin
  /// @param the index of the pin
  /// @param xLo
  /// @param yLo
  /// @param xHi
  /// @param yHi
  void addPinShape(IndexType pinIdx, LocType xLo, LocType yLo, LocType xHi,
                   LocType yHi) {
    _db.pin(pinIdx).shape().unionBox(Box<LocType>(xLo, yLo, xHi, yHi));
  }
  /// @brief allocate a new net
  /// @return the index for the net
  IndexType allocateNet() { return _db.allocateNet(); }
  /// @brief set the name of a net
  /// @param first: the index of the net
  /// @param second: the net name
  void setNetName(IndexType netIdx, const std::string &netName) {
    _db.net(netIdx).setName(netName);
  }
  /// @brief add pin to net
  /// @param first: the pin index
  /// @param second: the net index
  void addPinToNet(IndexType pinIdx, IndexType netIdx) {
    _db.net(netIdx).addPin(pinIdx);
    _db.pin(pinIdx).addNetIdx(netIdx);
  }
  /// @brief set the external net bounding box
  /// @param first: index of net
  /// @param second: xLo
  /// @param third: yLo
  /// @param fourth: xHi
  /// @param fifth: yHi
  void setNetExternalBBox(IndexType netIdx, LocType xLo, LocType yLo,
                          LocType xHi, LocType yHi) {
    _db.net(netIdx).setExternalBBox(Box<LocType>(xLo, yLo, xHi, yHi));
  }
  /// @brief set the net weight
  /// @param the index of the net
  /// @param the weight for the net
  void setNetWgt(IndexType netIdx, IntType weight) {
    _db.net(netIdx).setWeight(weight);
  }

  /// @brief assign to left
  /// @param the net index
  void fpIoPinAssignLeft(IndexType netIdx) {
    _db.net(netIdx).fpIoPinAssignLeft();
  }
  /// @brief assign to right
  /// @param the net index
  void fpIoPinAssignRight(IndexType netIdx) {
    _db.net(netIdx).fpIoPinAssignRight();
  }
  /// @brief clear io pin assignment status from floorplan
  /// @param the net index
  void clearFpIoPinAssign(IndexType netIdx) {
    _db.net(netIdx).clearFpIoPinAssign();
  }
  /// @brief allocate a symmetric group
  /// @return the index of the symmetric group
  IndexType allocateSymGrp() { return _db.allocateSymGrp(); }
  /// @brief add a symmetric pair to a symmetric group
  /// @param first: symmetric group index
  /// @param second: the index of the first cell. The order of the cell indices
  /// does not matter.
  /// @param third: the index of the second cell. The order of the cell indices
  /// does not matter.
  void addSymPair(IndexType symGrpIdx, IndexType cellIdx1, IndexType cellIdx2) {
    _db.symGroup(symGrpIdx).addSymPair(cellIdx1, cellIdx2);
  }
  /// @brief add a self symmetric cell to a symmetric group
  /// @param first: the symmetric group index
  /// @param second: the cell index
  void addSelfSym(IndexType symGrpIdx, IndexType selfSymCellIdx) {
    _db.symGroup(symGrpIdx).addSelfSym(selfSymCellIdx);
  }
  /// @brief add cell layout shape
  /// @param the index of the cell
  /// @param the index of the placer layer
  /// @param xLo
  /// @param yLo
  /// @param xHi
  /// @param yHi
  void addCellShape(IndexType cellIdx, IndexType placerLayer, LocType xLo,
                    LocType yLo, LocType xHi, LocType yHi) {
    _db.cell(cellIdx).unionBBox(placerLayer, Box<LocType>(xLo, yLo, xHi, yHi));
  }
  void setCellFlip(IndexType cellIdx) { _db.cell(cellIdx).flip(); }
  /// @brief allocate a new proximity group
  /// @return the index of the proximity group
  IndexType allocateProximityGroup() { return _db.allocateProximityGroup(); }
  /// @brief add an cell to a prixmityGroup
  /// @param first: the index of the cell
  /// @param second: the index of the proximity group
  void addCellToProximityGroup(IndexType cellIdx, IndexType proximityGrpIdx) {
    _db.proximityGrp(proximityGrpIdx).addCell(cellIdx);
  }
  /// @brief set a pair of symmetric nets
  /// @param first and second: two net indices. The order does not matter
  void addSymNetPair(IndexType netIdx1, IndexType netIdx2) {
    _db.net(netIdx1).setSymNet(netIdx2, true);
    _db.net(netIdx2).setSymNet(netIdx1, false);
  }
  /// @brief allocate a signal path
  /// @return the index of the signal path
  IndexType allocateSignalPath() { return _db.allocateSignalPath(); }
  /// @brief mark a signal path as power
  /// @param the index of the signal path
  void markSignalPathAsPower(IndexType pathIdx) {
    _db.signalPath(pathIdx).markAsPower();
  }
  /// @brief add a new pin to a signal path
  /// @param first: the index of signal path
  /// @param second: the name of cell
  /// @param third: the name of pin
  void addPinToSignalPath(IndexType pathIdx, const std::string &cellName,
                          const std::string &pinName) {
    return _db.addPinToSignalPath(pathIdx, cellName, pinName);
  }

  /// @brief mark a net as self-symmmetric
  /// @param the net index
  void markSelfSymNet(IndexType netIdx) { _db.net(netIdx).markSelfSym(); }
  /// @brief open the functionality of virtual pin assignment
  void openVirtualPinAssignment() {
    _db.parameters().openVirtualPinAssignment();
  }
  /// @brief close the functionality of virtual pin assignment
  void closeVirtualPinAssignment() {
    _db.parameters().closeVirtualPinAssignment();
  }
  /// @brief set net to be io pin
  void markAsIoNet(IndexType netIdx) { _db.net(netIdx).setIsIo(true); }
  /// @brief remove io net mark
  void revokeIoNet(IndexType netIdx) { _db.net(netIdx).setIsIo(false); }
  /// @brief mark a net as vdd
  void markAsVddNet(IndexType netIdx) {
    _db.net(netIdx).markAsVdd();
    _db.net(netIdx).setWeight(1);
  }
  /// @brief mark a net as vss
  void markAsVssNet(IndexType netIdx) {
    _db.net(netIdx).markAsVss();
    _db.net(netIdx).setWeight(1);
  }
  /// @brief get the x coordinate of io net
  LocType iopinX(IndexType netIdx) {
    return _db.net(netIdx).virtualPinLoc().x();
  }
  /// @brief get the y coordinate of io net
  LocType iopinY(IndexType netIdx) {
    return _db.net(netIdx).virtualPinLoc().y();
  }
  /// @brief get whether the io pin is on top or bottom
  /// @return true: top or bottom
  /// @return false: left or right
  bool isIopinVertical(IndexType netIdx) const {
    return _db.net(netIdx).iopinVertical();
  }
  /// @brief Add WPE vertical spacing rule. If device width > width, it need spacing
  void addWpeVerticalSpacingRule(LocType width, LocType spacing) {
    _db.tech().addWpeVerticalSpacingRule(width, spacing);
  }
  /// @brief Add WPE horizontal spacing rule. If device length > length, it need spacing
  void addWpeHorizontalSpacingRule(LocType length, LocType spacing) {
    _db.tech().addWpeHorizontalSpacingRule(length, spacing);
  }
  /*------------------------------*/
  /* Standard output interface    */
  /*------------------------------*/
  /// @brief get the x coordinate for a cell
  /// @param  the cell index
  /// @return the x coordinate
  LocType xCellLoc(IndexType cellIdx) { return _db.cell(cellIdx).xLoc(); }
  /// @brief get the y coordinate for a cell
  /// @param  the cell index
  /// @return the y coordinate
  LocType yCellLoc(IndexType cellIdx) { return _db.cell(cellIdx).yLoc(); }
  /// @brief get the index of the cell based on name
  /// @param cell name
  /// @return the cellIdx
  IndexType cellIdxName(const std::string name);
  /// @brief set cell name
  /// @param first cellIdx
  /// @param second cell name
  void setCellName(IndexType cellIdx, const std::string name) {
    _db.setCellName(cellIdx, name);
  }
  /// @brief return cell name
  /// @param first cellIdx
  /// @return cell name
  std::string cellName(IndexType cellIdx) { return _db.cell(cellIdx).name(); }
  /// @brief return number of cells
  /// @return number of cells
  IndexType numCells() { return _db.numCells(); }
  /// @brief get the index pin
  /// @param cell index
  /// @param the pin index of the cell
  /// @return the pinIdx in database
  IntType pinIdx(IndexType cellIdx, IndexType pinCellIdx) {
    if (pinCellIdx >= _db.cell(cellIdx).numPinIdx())
      return -1;
    return _db.cell(cellIdx).pinIdx(pinCellIdx);
  }
  /* Run time */
  /// @brief get the total run time
  /// @return time in us
  decltype(auto) runtimeIdeaPlaceEx() {
    return WATCH_LOOK_RECORD_TIME("IdeaPlaceEx");
  }
  /// @brief get the the run time for global placement
  /// @return time in us
  decltype(auto) runtimeGlobalPlace() {
    return WATCH_LOOK_RECORD_TIME("NlpGPlacer");
  }
  /// @brief get the the run time for calculating objectives in global placement
  /// @return time in us
  decltype(auto) runtimeGlobalPlaceCalcObj() {
    return WATCH_LOOK_RECORD_TIME("GP_calculate_obj");
  }
  /// @brief get the the run time for calculating gradients in global placement
  /// @return time in us
  decltype(auto) runtimeGlobalPlaceCalcGrad() {
    return WATCH_LOOK_RECORD_TIME("GP_calculate_gradient");
  }
  /// @brief get the the run time for the optimization kernel in global
  /// placement
  /// @return time in us
  decltype(auto) runtimeGlobalPlaceOptmKernel() {
    return WATCH_LOOK_RECORD_TIME("GP_optimizer_kernel");
  }
  /// @brief get the the run time for the optimization kernel in global
  /// placement
  /// @return time in us
  decltype(auto) runtimeGlobalPlaceOptimize() {
    return WATCH_LOOK_RECORD_TIME("GP_optimize");
  }
  /// @brief get the the run time for updating the problem (e.g. multipliers) in
  /// global placement
  /// @return time in us
  decltype(auto) runtimeGlobalPlaceUpdateProblem() {
    return WATCH_LOOK_RECORD_TIME("GP_update_problem");
  }
  /// @brief get the the run time for legalization
  /// @return time in us
  decltype(auto) runtimeLegalization() {
    return WATCH_LOOK_RECORD_TIME("legalization");
  }
  /// @brief get the the run time for detailed placement
  /// @return time in us
  decltype(auto) runtimeDetailedPlacement() {
    return WATCH_LOOK_RECORD_TIME("detailedPlacement");
  }

  /*------------------------------*/
  /* Log option                   */
  /*------------------------------*/
  void logScreenOn() { MsgPrinter::screenOn(); }
  void logScreenOff() { MsgPrinter::screenOff(); }
  LocType hpwl() const { return _db.hpwl(); }

  /*------------------------------*/
  /* Well                         */
  /*------------------------------*/
  IndexType allocateWell() {
    IndexType wellIdx = _db.allocateWell();
    return wellIdx;
  }
  void setWellName(const IndexType wellIdx, const std::string name) {
    _db.well(wellIdx).setName(name);
  }
  void setWellShape(const IndexType wellIdx,
                    const std::vector<Point<LocType>> &vPts) {
    Polygon<LocType>::RingType ring(vPts.begin(), vPts.end());
    _db.well(wellIdx).setShape(ring);
  }
  void printWellInfo(const IndexType wellIdx) const {
    _db.well(wellIdx).printInfo();
  }
  void clearWells() { _db.clearWells(); }
  void assignCellToWell() { _db.assignCellToWellAndRemoveUnusedWell(); }
  void splitWells() { _db.splitWells(); }
  IndexType numWellRects() const {
    return _db.numWellRects();
  }
  const Box<LocType> &wellRect(IndexType idx) const {
    auto rectIdx = _db.getWellRectIdx(idx);
    return _db.well(rectIdx.first).rects().at(rectIdx.second);
  }
  /// @brief Add a well contact template
  void addVddContactTemplate(LocType xLo, LocType yLo, LocType xHi, LocType yHi, IntType weight) {
    _db.tech().addVddContactTemplate(xLo, yLo, xHi, yHi, weight);
  }
  /// @brief Set the required spacing for VDD contact
  void setVddContactRequiredSpacing(IntType require) {
    _db.tech().setVddContactRequiredSpacing(require);
  }
  /* Debug */
#ifdef DEBUG_DRAW
  void drawDebug(const std::string &filename) {
    _db.drawCellBlocks(filename);
  }
#endif
  void crfSigPathStat();

protected:
  Database _db; ///< The placement engine database
  ProximityMgr _proximityMgr;
};

PROJECT_NAMESPACE_END

#endif /// IDEAPLACE_IDEAPLACEEX_H_
