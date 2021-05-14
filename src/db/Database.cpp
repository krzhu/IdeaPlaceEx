#include "Database.h"
#include <boost/geometry/algorithms/distance.hpp>

PROJECT_NAMESPACE_BEGIN

bool Database::initCells() {
  for (IndexType cellIdx = 0; cellIdx < this->numCells(); ++cellIdx) {
    if (!this->initCell(cellIdx)) {
      return false;
    }
  }
  return true;
}

bool Database::initCell(IndexType cellIdx) {
  IndexType numLayers = this->tech().numLayers();
  _cellArray.at(cellIdx).allocateLayers(numLayers);
  return true;
}


RealType Database::area() const {
    RealType xMax = LOC_TYPE_MIN;
    RealType xMin = LOC_TYPE_MAX;
    RealType yMax = LOC_TYPE_MIN;
    RealType yMin = LOC_TYPE_MAX;
    for (const auto &cell : _cellArray) {
      const auto box = cell.cellBBoxOff();
      xMax = std::max(xMax, static_cast<RealType>(box.xh()));
      xMin = std::min(xMin, static_cast<RealType>(box.xl()));
      yMax = std::max(yMax, static_cast<RealType>(box.yh()));
      yMin = std::min(yMin, static_cast<RealType>(box.yl()));
    }
    for (const auto &well : _wellArray) {
      for (const auto &pt : well.shape().outer()) {
        xMax = std::max(xMax, static_cast<RealType>(pt.x()));
        yMax = std::max(yMax, static_cast<RealType>(pt.y()));
        xMin = std::min(xMin, static_cast<RealType>(pt.x()));
        yMin = std::min(yMin, static_cast<RealType>(pt.y()));
      }
    }
    RealType area = xMax - xMin;
    area *= (yMax - yMin);
    DBG("Area %f * %f = %f \n", yMax-yMin, xMax-xMin, area);
    return area;
}

LocType Database::hpwl() const {
  LocType hpwl = 0;
  for (const auto &net : _netArray) {
    if (net.numPinIdx() <= 1) {
      continue;
    }
    LocType xMax = LOC_TYPE_MIN;
    LocType xMin = LOC_TYPE_MAX;
    LocType yMax = LOC_TYPE_MIN;
    LocType yMin = LOC_TYPE_MAX;
    for (IndexType pinIdx : net.pinIdxArray()) {
      const auto &pin = _pinArray.at(pinIdx);
      auto pinLoc = pin.midLoc();
      IndexType cellIdx = pin.cellIdx();
      const auto &cell = _cellArray.at(cellIdx);
      auto cellLoc = cell.loc();
      auto pinFinalLoc = pinLoc + cellLoc;
      xMax = std::max(xMax, pinFinalLoc.x());
      xMin = std::min(xMin, pinFinalLoc.x());
      yMax = std::max(yMax, pinFinalLoc.y());
      yMin = std::min(yMin, pinFinalLoc.y());
    }
    hpwl += ((xMax - xMin) + (yMax - yMin)) * net.weight();
  }
  return hpwl;
}

LocType Database::hpwlWithVitualPins() const {
  LocType hpwl = 0;
  for (const auto &net : _netArray) {
    if (net.numPinIdx() <= 1) {
      continue;
    }
    LocType xMax = LOC_TYPE_MIN;
    LocType xMin = LOC_TYPE_MAX;
    LocType yMax = LOC_TYPE_MIN;
    LocType yMin = LOC_TYPE_MAX;
    for (IndexType pinIdx : net.pinIdxArray()) {
      const auto &pin = _pinArray.at(pinIdx);
      auto pinLoc = pin.midLoc();
      IndexType cellIdx = pin.cellIdx();
      const auto &cell = _cellArray.at(cellIdx);
      auto cellLoc = cell.loc();
      auto pinFinalLoc = pinLoc + cellLoc;
      xMax = std::max(xMax, pinFinalLoc.x());
      xMin = std::min(xMin, pinFinalLoc.x());
      yMax = std::max(yMax, pinFinalLoc.y());
      yMin = std::min(yMin, pinFinalLoc.y());
    }
    if (net.isIo()) {
      const auto &virPin = net.virtualPinLoc();
      xMax = std::max(xMax, virPin.x());
      xMin = std::min(xMin, virPin.x());
      yMax = std::max(yMax, virPin.y());
      yMin = std::min(yMin, virPin.y());
    }
    hpwl += ((xMax - xMin) + (yMax - yMin)) * net.weight();
  }
  return hpwl;
}

LocType Database::findSymAxis() {
#ifndef MULTI_SYM_GROUP
  LocType axis = LOC_TYPE_MIN;
#else
  LocType axis = LOC_TYPE_MIN;
  Assert(false);
#endif
  for (const auto &symGrp : _symGroups) {
    if (symGrp.numSymPairs() == 0) {
      break;
    }
    LocType symAxis;
    for (IndexType symPairIdx = 0; symPairIdx < symGrp.numSymPairs();
         ++symPairIdx) {
      const auto &symPair = symGrp.symPair(symPairIdx);
      symAxis = (cell(symPair.firstCell()).xCenter() +
                 cell(symPair.secondCell()).xCenter()) /
                2;
      /*
          ERR("Found sym. axis %d, cell %s at x= %d, cell %s at x= %d \n",
         symAxis, cell(symPair.firstCell()).name().c_str(),
         cell(symPair.firstCell()).xCenter(),
                 cell(symPair.secondCell()).name().c_str(),
         cell(symPair.secondCell()).xCenter()); const auto &cell1 =
         cell(symPair.firstCell()); const auto &cell2 =
         cell(symPair.secondCell()); DBG("cell1 xlo %d xcenter %d xhi %d, xlen
         %d \n", cell1.xLo(), cell1.xCenter(), cell1.xHi(),
         cell1.cellBBox().xLen()); DBG("cell2 xlo %d xcenter %d xhi %d, xlen %d
         \n", cell2.xLo(), cell2.xCenter(), cell2.xHi(),
         cell2.cellBBox().xLen());
          */
      break;
    }
    if (axis != LOC_TYPE_MIN) {
      if (axis != symAxis) {
        ERR("Check Symmetry: different sym axises! \n");
        return false;
      }
    } else {
      axis = symAxis;
    }
  }
  return axis;
}

bool Database::checkSym() {
#ifndef MULTI_SYM_GROUP
  LocType axis = LOC_TYPE_MIN;
#else
  LocType axis = LOC_TYPE_MIN;
  Assert(false);
#endif
  for (const auto &symGrp : _symGroups) {
    if (symGrp.numSymPairs() == 0) {
      break;
    }
    LocType symAxis;
    for (IndexType symPairIdx = 0; symPairIdx < symGrp.numSymPairs();
         ++symPairIdx) {
      const auto &symPair = symGrp.symPair(symPairIdx);
      symAxis = (cell(symPair.firstCell()).xCenter() +
                 cell(symPair.secondCell()).xCenter()) /
                2;
      /*
          ERR("Found sym. axis %d, cell %s at x= %d, cell %s at x= %d \n",
         symAxis, cell(symPair.firstCell()).name().c_str(),
         cell(symPair.firstCell()).xCenter(),
                 cell(symPair.secondCell()).name().c_str(),
         cell(symPair.secondCell()).xCenter()); const auto &cell1 =
         cell(symPair.firstCell()); const auto &cell2 =
         cell(symPair.secondCell()); DBG("cell1 xlo %d xcenter %d xhi %d, xlen
         %d \n", cell1.xLo(), cell1.xCenter(), cell1.xHi(),
         cell1.cellBBox().xLen()); DBG("cell2 xlo %d xcenter %d xhi %d, xlen %d
         \n", cell2.xLo(), cell2.xCenter(), cell2.xHi(),
         cell2.cellBBox().xLen());
          */
      break;
    }
    if (axis != LOC_TYPE_MIN) {
      if (axis != symAxis) {
        ERR("Check Symmetry: different sym axises! \n");
        return false;
      }
    } else {
      axis = symAxis;
    }
    for (IndexType symPairIdx = 0; symPairIdx < symGrp.numSymPairs();
         ++symPairIdx) {
      const auto &symPair = symGrp.symPair(symPairIdx);
      if (symAxis != (cell(symPair.firstCell()).xCenter() +
                      cell(symPair.secondCell()).xCenter()) /
                         2) {
        ERR("Check Symmetry: symPair failed. axis %d, cell %d %s at x= %d, cell "
            "%d %s at x= %d \n",
            symAxis, symPair.firstCell(), cell(symPair.firstCell()).name().c_str(),
            cell(symPair.firstCell()).xCenter(),
            symPair.secondCell(), cell(symPair.secondCell()).name().c_str(),
            cell(symPair.secondCell()).xCenter());
        const auto &cell1 = cell(symPair.firstCell());
        const auto &cell2 = cell(symPair.secondCell());
        DBG("cell1 xlo %d xcenter %d xhi %d, xlen %d \n", cell1.xLo(),
            cell1.xCenter(), cell1.xHi(), cell1.cellBBox().xLen());
        DBG("cell2 xlo %d xcenter %d xhi %d, xlen %d \n", cell2.xLo(),
            cell2.xCenter(), cell2.xHi(), cell2.cellBBox().xLen());
        return false;
      }
      if ((cell(symPair.firstCell()).yLoc() !=
           cell(symPair.secondCell()).yLoc())) {
        ERR("Check Symmetry: symPair failed. different y, cell %d %d \n",
            symAxis, cell(symPair.firstCell()).yLoc(),
            cell(symPair.secondCell()).yLoc());
        return false;
      }
    }
    for (IndexType ssIdx = 0; ssIdx < symGrp.numSelfSyms(); ++ssIdx) {
      if (cell(symGrp.selfSym(ssIdx)).xCenter() != symAxis) {
        ERR("Check Symmetry: selfSym failed. axis %d, cell idx %d at %d \n",
            symAxis, symGrp.selfSym(ssIdx),
            cell(symGrp.selfSym(ssIdx)).xCenter());
        return false;
      }
    }
  }
  return true;
}

void Database::assignCellToWellAndRemoveUnusedWell() {
  for (IndexType cellIdx = 0; cellIdx <  _cellArray.size(); ++cellIdx) {
    if (not cell(cellIdx).needWell()) { continue; }
    LocType shortestDist = LOC_TYPE_MAX;
    IndexType targetWellIdx = INDEX_TYPE_MAX;
    const Point<LocType> cellCenter = cell(cellIdx).cellBBoxOff().center();
    for (IndexType wellIdx = 0; wellIdx < _wellArray.size(); ++wellIdx) {
      auto &well = _wellArray.at(wellIdx);
      LocType dist = boost::geometry::distance(cellCenter, well.shape());
      if (dist < shortestDist) {
        targetWellIdx = wellIdx;
        shortestDist = dist;
      }
    }
    if (targetWellIdx != INDEX_TYPE_MAX) {
      _wellArray.at(targetWellIdx).addCellIdx(cellIdx);
    }
  }
  _wellArray.erase(
      std::remove_if(_wellArray.begin(), _wellArray.end(),
        [](const Well &well) { return not well.assignedWithAnyCell();}),
      _wellArray.end()
      );
  // Record well index for cells
  for (IndexType wellIdx = 0; wellIdx < _wellArray.size(); ++wellIdx) {
    const auto &well = _wellArray.at(wellIdx);
    for (IndexType cellIdx  : well.sCellIds()) {
      _cellArray.at(cellIdx).setWellIdx(wellIdx);
    }
  }
}

void Database::offsetLayout() {
  LocType xMin = LOC_TYPE_MAX;
  LocType yMin = LOC_TYPE_MAX;
  // Cell
  for (const auto &cell : _cellArray) {
    xMin = std::min(xMin, cell.xLo());
    yMin = std::min(yMin, cell.yLo());
  }
  // Well
  for (const auto &well : _wellArray) {
    for (const auto &pt : well.shape().outer()) {
      xMin = std::min(xMin, pt.x());
      yMin = std::min(yMin, pt.y());
    }
  }
  // Pin
  for (const auto &net : _netArray) {
    if (not net.isValidVirtualPin()) {
      return;
    }
    xMin = std::min(xMin, net.virtualPinLoc().x());
    yMin = std::min(yMin, net.virtualPinLoc().y());
  }
  LocType xOffset = std::max(-xMin, 0);
  LocType yOffset = std::max(-yMin, 0);
  if (xOffset == 0 and yOffset == 0) {
    return;
  }
  // Cell
  for (auto &cell : _cellArray) {
    cell.setXLoc(cell.xLoc() + xOffset);
    cell.setYLoc(cell.yLoc() + yOffset);
  }
  // Well
  for (auto &well : _wellArray) {
    for (auto &pt : well.shape().outer()) {
      pt= pt +  Point<LocType>(xOffset, yOffset);
    }
  }
  // Pin
  for (auto &net : _netArray) {
    if (not net.isValidVirtualPin()) {
      return;
    }
    net.offsetVirtualPin(Point<LocType>(xOffset, yOffset));
  }

}

/*------------------------------*/
/* Debug functions              */
/*------------------------------*/
#ifdef DEBUG_DRAW
PROJECT_NAMESPACE_END
#include "util/Polygon2Rect.h"
#include "writer/gdsii/WriteGds.h"
PROJECT_NAMESPACE_BEGIN
void Database::drawCellBlocks(const std::string &filename) {
  auto wg = std::make_shared<WriteGds>(filename);
  if (!wg->initWriter()) {
    return;
  }
  if (!wg->createLib("TOP", 2000, 1e-6 / 2000)) // Hardcoded numbers
  {
    return;
  }
  if (!wg->writeCellBgn("DEBUG")) {
    return;
  }
  // Write all the cells
  for (IndexType cellIdx = 0; cellIdx < this->numCells(); ++cellIdx) {
    const auto &cell = this->cell(cellIdx);
    Box<LocType> cellBox = cell.cellBBox();
    cellBox.enlargeBy(0);
    auto cellLoc = Point<LocType>(cell.xLoc(), cell.yLoc());
    cellBox.offsetBy(cellLoc);
    wg->writeRectangle(cellBox, cellIdx, 0);
#if 1
    // Also write pins
    for (IndexType pinIdxInCell = 0; pinIdxInCell < cell.numPinIdx();
         ++pinIdxInCell) {
      const auto &pin = this->pin(cell.pinIdx(pinIdxInCell));
      Box<LocType> pinBox = pin.shape();
      pinBox.offsetBy(cellLoc);
      wg->writeRectangle(pinBox, 100 + cellIdx, 0);
    }
#endif
  }
  // net virtual pins
  for (IndexType netIdx = 0; netIdx < numNets(); ++netIdx) {
    const auto &net = _netArray.at(netIdx);
    if (net.isValidVirtualPin()) {
      Box<LocType> box = Box<LocType>(net.virtualPinLoc(), net.virtualPinLoc());
      box.enlargeBy(100);
      wg->writeRectangle(box, 500 + netIdx, 0);
    }
  }
  // Well
  for (IndexType wellIdx = 0; wellIdx < _wellArray.size(); ++wellIdx) {
    const auto &well = _wellArray.at(wellIdx);
    std::vector<Box<int>> boxes; // Splited polygon
    if (not klib::convertPolygon2Rects(well.shape().outer(), boxes)) {
      ERR("NlpGPlacer:: cannot split well polygon! \n");
      Assert(false);
    }
    for (const auto &box : boxes) {
      wg->writeRectangle(box, wellIdx + 200, 0);
    }
  }
  // END
  wg->writeCellEnd();
  wg->endLib();
  DBG("Ideaplace: saved debugging gds to %s\n", filename.c_str());
}

#endif // DEBUG_DRAW

void Database::splitSignalPathsBySymPairs() {
  for (IndexType pathIdx = 0; pathIdx < this->vSignalPaths().size();
       ++pathIdx) {
    auto &path = this->vSignalPaths().at(pathIdx);
    // since the automatic generated current flow could include both sides of
    // differential pair, we need to split the path is seen such scenario
    std::set<IndexType> symCells; // The other side of cells of the cells have
                                  // been seen in this path
    auto judgeCellIdx = [&](IndexType cellIdx) {
      if (symCells.find(cellIdx) != symCells.end()) {
        symCells.clear();
        return false;
      }
      if (this->cell(cellIdx).hasSymPair()) {
        symCells.insert(this->cell(cellIdx).symNetIdx());
      }
      return true;
    };
    bool needBreak = false;
    IndexType breakIdx = 0;
    for (IndexType idx = 0; idx < path.vPinIdxArray().size(); ++idx) {
      const IndexType pinIdx = path.vPinIdxArray().at(idx);
      const auto &pin = this->pin(pinIdx);
      const IndexType cellIdx = pin.cellIdx();
      if (not judgeCellIdx(cellIdx)) {
        needBreak = true;
        breakIdx = idx;
        break;
      }
    }
    if (needBreak) {
      // Split the parts after breakIdx into a new path
      auto newPathIdx = this->allocateSignalPath();
      this->signalPath(newPathIdx).copySettings(this->signalPath(pathIdx));
      // copy the pins to the new path and erase from the original
      for (IndexType idx = breakIdx;
           idx < this->signalPath(pathIdx).vPinIdxArray().size(); ++idx) {
        auto tempPath = this->signalPath(pathIdx);
        this->signalPath(newPathIdx)
            .addPinIdx(this->signalPath(pathIdx).vPinIdxArray().at(idx));
      }
      this->signalPath(pathIdx).vPinIdxArray().erase(
          this->signalPath(pathIdx).vPinIdxArray().begin() + breakIdx,
          this->signalPath(pathIdx).vPinIdxArray().end());
    }
  }
}
PROJECT_NAMESPACE_END
