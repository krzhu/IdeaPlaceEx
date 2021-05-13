#include "CGLegalizer.h"
#include "constraintGraphGeneration.h"
#include "pinassign/VirtualPinAssigner.h"
#include "place/legalize/lp_legalize.hpp"
#include "place/legalize/legalize_well.hpp"

PROJECT_NAMESPACE_BEGIN


enum class Relation {
  RIGHT, LEFT, TOP, BOTTOM, // overlap or parallel run,
  TOP_RIGHT, TOP_LEFT, BOTTOM_RIGHT, BOTTOM_LEFT // Completely disjoined
};

Relation determineBoxRelation(const Box<LocType> &box1, const Box<LocType> &box2, bool isSymPair=false, bool isSelfSym=false) {
  Assert(not (isSymPair and isSelfSym));
  if (isSymPair) {
    if (box1.xLo() < box2.xLo()) {
      return Relation::RIGHT;
    }
    else {
      return Relation::LEFT;
    }
  }
  if (isSelfSym) {
    if (box1.yLo() < box2.yLo()) {
      return Relation::TOP;
    }
    else {
      return Relation::BOTTOM;
    }
  }
  if (box1.xHi() <= box2.xLo()) {
    if (box1.yHi() <= box2.yLo()) {
      return Relation::TOP_RIGHT;
    }
    else if (box1.yLo() >= box2.yHi()) {
      return Relation::BOTTOM_RIGHT;
    }
    else {
      // Parallel run on y
      return Relation::RIGHT;
    }
  }
  else if (box1.xLo() >= box2.xHi()) {
    if (box1.yHi() <= box2.yLo()) {
      return Relation::TOP_LEFT;
    }
    else if (box1.yLo() >= box2.yHi()) {
      return Relation::BOTTOM_LEFT;
    }
    else {
      // Parallel run on y
      return Relation::LEFT;
    }
  }
  else {
    // parallel run on x
    if (box1.yHi() <= box2.yLo()) {
      return Relation::TOP;
    }
    else if (box1.yLo() >= box2.yHi()) {
      return Relation::BOTTOM;
    }
    // overlap
    bool left, bottom;
    LocType overlapX, overlapY;
    if (box1.xLo() <= box2.xLo()) {
      left = false;
      overlapX = box1.xHi() - box2.xLo();
    }
    else {
      left = true;
      overlapX = box2.xHi() - box1.xLo();
    }
    if (box1.yLo() <= box2.yLo()) {
      bottom = false;
      overlapY = box1.yHi() - box2.yLo();
    }
    else {
      bottom = true;
      overlapY = box2.yHi() - box1.yLo();
    }
    if (overlapX < overlapY) {
      if (left) { return Relation::LEFT; }
      else { return Relation::RIGHT; }
    }
    else {
      if (bottom) { return Relation::BOTTOM; }
      else { return Relation::TOP; }
    }
  }
}

void CGLegalizer::legalizeWells(bool addContact) {
  WellLegalizer wellLegalizer(_db);
  if (addContact) {
    wellLegalizer.legalizeAndAddContact();
  }
  else {
    wellLegalizer.legalize();
  }
}

void CGLegalizer::generateIndividualWells(){
  WellLegalizer wellLegalizer(_db);
  wellLegalizer.generateIndividualWells();
}

bool CGLegalizer::legalize() {

  VirtualPinAssigner pinAssigner(_db);

  // this->generateConstraints();
  _db.offsetLayout();


  areaDrivenCompaction();

  if (_db.parameters().ifUsePinAssignment()) {
    pinAssigner.solveFromDB();
  }
  if (!wirelengthDrivenCompaction()) {
    INF("CG Legalizer: detailed placement fine tunning failed. Directly output "
        "legalization output. \n");
    return true;
  }
  return true;

  INF("CG Legalizer: legalization finished\n");
  return true;
}

void CGLegalizer::findCellBoundary() {
  LocType xMin = LOC_TYPE_MAX;
  LocType xMax = LOC_TYPE_MIN;
  LocType yMin = LOC_TYPE_MAX;
  LocType yMax = LOC_TYPE_MIN;
  for (IndexType cellIdx = 0; cellIdx < _db.numCells(); ++cellIdx) {
    const auto &cell = _db.cell(cellIdx);
    xMin = std::min(xMin, cell.xLo());
    xMax = std::max(xMax, cell.xHi());
    yMin = std::min(yMin, cell.yLo());
    yMax = std::max(yMax, cell.yHi());
  }
    for (const auto &well : _db.vWells()) {
      for (const auto &pt : well.shape().outer()) {
        xMax = std::max(xMax, pt.x());
        yMax = std::max(yMax, pt.y());
        xMin = std::min(xMin, pt.x());
        yMin = std::min(yMin, pt.y());
      }
    }
  _wStar = std::max(0.0, static_cast<RealType>(xMax - xMin)) + 10;
  _hStar = std::max(0.0, static_cast<RealType>(yMax - yMin)) + 10;
}

void CGLegalizer::prepare() {
  _db.offsetLayout();
}

void CGLegalizer::generateHorConstraints() {

  _db.splitWells(1);

  _hConstraints.clear();
  _vConstraints.clear();
  // Init the irredundant constraint edges
  
  auto getCellOrWellBBoc = [&](IndexType cellIdx) {
    if (cellIdx < _db.numCells()) {
      return _db.cell(cellIdx).cellBBoxOff();
    }
    auto wellIdx = _db.getWellRectIdx(cellIdx - _db.numCells());
    return _db.well(wellIdx.first).rects().at(wellIdx.second);
  };

  auto exemptVerticalMoveFunc = [&](IndexType cellIdx1, IndexType cellIdx2) {
    IndexType largestIdx = _db.numCells() + _db.numWellRects();
    if (cellIdx1 >= largestIdx or cellIdx2 >= largestIdx) {
      return false;
    }
    const auto box1 = getCellOrWellBBoc(cellIdx1);
    const auto box2 = getCellOrWellBBoc(cellIdx2);
    const auto relation = determineBoxRelation(box1, box2);
    if (relation == Relation::TOP or relation == Relation::BOTTOM) {
      return true;
    }
    return false;
  };

  SweeplineConstraintGraphGenerator sweepline(_db, _hConstraints,
                                              _vConstraints);
  sweepline.setExemptFunc(exemptVerticalMoveFunc);
  if (_wellAware) {
    sweepline.openConsiderWell();
  }
  else {
    sweepline.closeConsiderWell();
  }
  sweepline.solve();
}
void CGLegalizer::generateVerConstraints() {
  _db.splitWells(2);
  _hConstraints.clear();
  _vConstraints.clear();
  // Init the irredundant constraint edges

  SweeplineConstraintGraphGenerator sweepline(_db, _hConstraints,
                                              _vConstraints);
  if (_wellAware) {
    sweepline.openConsiderWell();
  }
  else {
    sweepline.closeConsiderWell();
  }
  sweepline.solve();
}


BoolType CGLegalizer::areaDrivenCompaction() {
  auto area = _db.area();
  auto lastArea = area;
  IndexType iter = 0;
  bool legalizedHor = false;
  bool legalizedVer = false;
  do {
    lastArea = area;
    // Horizontal
    if (iter % 2 == 0 ) {
      this->generateHorConstraints();
      INF("CG legalizer: area-driven legalize horizontal LP...\n");
      auto horSolver = lp_legalize::LpLegalizeArea<lp_legalize::LEGALIZE_HORIZONTAL_DIRECTION, lp_legalize::DO_NOT_RELAX_SYM_CONSTR>(_db, _hConstraints);
      bool horpass = horSolver.solve();
      if (horpass) {
        legalizedHor = true;
        horSolver.exportSolution();
      } else {
        if (legalizedHor and legalizedVer) {
          INF("CG legalizer: LP failed, but valid solution exists. Exit \n");
          return true;
        }
        INF("CG legalizer: LP failed! No valid solution \n");
        return false;
      }
      _db.drawCellBlocks("./debug/hor.gds");
    }
    else {
      // Vertical
      this->generateVerConstraints();
      INF("CG legalizer: area-driven legalize vertical LP...\n");
      auto verSolver = lp_legalize::LpLegalizeArea<lp_legalize::LEGALIZE_VERTICAL_DIRECTION, lp_legalize::DO_NOT_RELAX_SYM_CONSTR>(_db, _vConstraints);
      bool verpass = verSolver.solve();
      if (verpass) {
        legalizedVer = true;
        verSolver.exportSolution();
      } else {
        if (legalizedHor and legalizedVer) {
          INF("CG legalizer: LP failed, but valid solution exists. Exit \n");
          return true;
        }
        INF("CG legalizer: LP failed! No valid solution \n");
        return false;
      }
    }
    ++iter;
    area = _db.area();
  } while (lastArea != area or iter <= 2);
  return true;
}

BoolType CGLegalizer::wirelengthDrivenCompaction() {
  LocType hpwl = _db.hpwl();
  LocType lastHpwl = LOC_TYPE_MAX;
  IndexType iter = 0;
  // Find the fixed boundary for this step
  findCellBoundary();
  do {
    lastHpwl = hpwl;
    // Horizontal
    if (iter % 2 == 0) {
      this->generateHorConstraints();
      INF("CG legalizer: wirelength-driven detailed placement horizontal LP...\n");
      auto horSolver = lp_legalize::LpLegalizeWirelength<lp_legalize::LEGALIZE_HORIZONTAL_DIRECTION, lp_legalize::DO_NOT_RELAX_SYM_CONSTR>(_db, _hConstraints);
#ifdef DEBUG_LEGALIZE
      DBG("wstar for width %f \n", _wStar);
#endif
      horSolver.setBoundary(_wStar);
      bool horpass = horSolver.solve();
      if (horpass) {
        horSolver.exportSolution();
      } else {
        Assert(false);
        return false;
      }
    }
    else {
      // Vertical
      this->generateVerConstraints();
      INF("CG legalizer:  wirelength-driven detailed placement vertical LP...\n");
      auto verSolver = lp_legalize::LpLegalizeWirelength<lp_legalize::LEGALIZE_VERTICAL_DIRECTION, lp_legalize::DO_NOT_RELAX_SYM_CONSTR>(_db, _vConstraints);
      verSolver.setBoundary(_hStar);
      bool verpass = verSolver.solve();
      if (verpass) {
        verSolver.exportSolution();
      } else {
        Assert(false);
        return false;
      }
    }
    ++iter;
    hpwl = _db.hpwl();
  } while (lastHpwl != hpwl or iter <= 2);

#ifdef DEBUG_LEGALIZE
#ifdef DEBUG_DRAW
  _db.drawCellBlocks("./debug/after_dp.gds");
#endif
#endif
  return true;
}

BoolType CGLegalizer::preserveRelationCompaction(LocType extraSpacing) {
  generateRedundantConstraintGraph();
  INF("CG legalizer: Prerserve relational coodinate horizontal LP...\n");
  auto horSolver = lp_legalize::LpLegalizeWirelength<lp_legalize::LEGALIZE_HORIZONTAL_DIRECTION, lp_legalize::DO_NOT_RELAX_SYM_CONSTR>(_db, _hConstrPreserve);
  horSolver.setExtraSpacing(extraSpacing);
  if (horSolver.solve()) {
    horSolver.exportSolution();
  }
  else {
    return false;
  }
  INF("CG legalizer: Prerserve relational coodinate vetical LP...\n");
  auto verSolver = lp_legalize::LpLegalizeWirelength<lp_legalize::LEGALIZE_VERTICAL_DIRECTION, lp_legalize::DO_NOT_RELAX_SYM_CONSTR>(_db, _vConstrPreserve);
  verSolver.setExtraSpacing(extraSpacing);
  if (verSolver.solve()) {
    verSolver.exportSolution();
  }
  else {
    return false;
  }

  return true;
}
void CGLegalizer::generateRedundantConstraintGraph() {
  if (not _hConstrPreserve.edges().empty() or not _vConstrPreserve.edges().empty()) {
    return;
  }
  for (IndexType ci1 = 0; ci1 < _db.numCells(); ++ci1) {
    const Box<LocType> box1 = _db.cell(ci1).cellBBoxOff();
    for (IndexType ci2 = ci1 + 1; ci2 < _db.numCells(); ++ci2) {
      bool isSymPair = _db.cell(ci1).hasSymPair() and _db.cell(ci1).symNetIdx() == ci2;
      bool bothIsSelfSym = _db.cell(ci1).isSelfSym() and _db.cell(ci2).isSelfSym();
      const Box<LocType> box2 = _db.cell(ci2).cellBBoxOff();
      auto relation = determineBoxRelation(box1, box2, isSymPair, bothIsSelfSym);
      if (relation == Relation::LEFT) {
        _hConstrPreserve.addConstraintEdge(ci2, ci1);
      }
      else if (relation == Relation::RIGHT) {
        _hConstrPreserve.addConstraintEdge(ci1, ci2);
      }
      else if (relation == Relation::BOTTOM) {
        _vConstrPreserve.addConstraintEdge(ci2, ci1);
      }
      else if (relation == Relation::TOP) {
        _vConstrPreserve.addConstraintEdge(ci1, ci2);
      }
      else if (relation == Relation::BOTTOM_LEFT) {
        _hConstrPreserve.addConstraintEdge(ci2, ci1);
        _vConstrPreserve.addConstraintEdge(ci2, ci1);
      }
      else if (relation == Relation::BOTTOM_RIGHT) {
        _hConstrPreserve.addConstraintEdge(ci1, ci2);
        _vConstrPreserve.addConstraintEdge(ci2, ci1);
      }
      else if (relation == Relation::TOP_LEFT) {
        _hConstrPreserve.addConstraintEdge(ci2, ci1);
        _vConstrPreserve.addConstraintEdge(ci1, ci2);
      }
      else {
        Assert(relation == Relation::TOP_RIGHT);
        _vConstrPreserve.addConstraintEdge(ci1, ci2);
        _hConstrPreserve.addConstraintEdge(ci1, ci2);
      }
    }
  }
}


PROJECT_NAMESPACE_END
