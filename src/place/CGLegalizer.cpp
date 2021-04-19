#include "CGLegalizer.h"
#include "constraintGraphGeneration.h"
#include "pinassign/VirtualPinAssigner.h"
#include "place/legalize/lp_legalize.hpp"

PROJECT_NAMESPACE_BEGIN

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
  _wStar = std::max(0.0, static_cast<RealType>(xMax - xMin)) + 10;
  _hStar = std::max(0.0, static_cast<RealType>(yMax - yMin)) + 10;
}

void CGLegalizer::prepare() {
  _db.offsetLayout();
}

void CGLegalizer::generateHorConstraints() {

  _hConstraints.clear();
  _vConstraints.clear();
  // Init the irredundant constraint edges

  auto exemptSelfSymsFunc = [&](IndexType cellIdx1, IndexType cellIdx2) {
    return false;
    if (cellIdx1 >= _db.numCells()) {
      return false;
    }
    if (cellIdx2 >= _db.numCells()) {
      return false;
    }
    if (_db.cell(cellIdx1).isSelfSym() and _db.cell(cellIdx2).isSelfSym()) {
      return true;
    }
    return false;
  };

  SweeplineConstraintGraphGenerator sweepline(_db, _hConstraints,
                                              _vConstraints);
  sweepline.setExemptFunc(exemptSelfSymsFunc);
  if (_wellAware) {
    sweepline.openConsiderWell();
  }
  else {
    sweepline.closeConsiderWell();
  }
  sweepline.solve();
}
void CGLegalizer::generateVerConstraints() {
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
  this->generateHorConstraints();
  INF("CG legalizer: area-driven legalize horizontal LP...\n");
  auto horSolver = lp_legalize::LpLegalizeArea<lp_legalize::LEGALIZE_HORIZONTAL_DIRECTION, lp_legalize::DO_NOT_RELAX_SYM_CONSTR>(_db, _hConstraints);
  if (horSolver.solve()) {
    horSolver.exportSolution();
  } else {
    Assert(false);
    return false;
  }
  this->generateVerConstraints();
  INF("CG legalizer: area-driven legalize vertical LP...\n");
  auto verSolver = lp_legalize::LpLegalizeArea<lp_legalize::LEGALIZE_VERTICAL_DIRECTION, lp_legalize::DO_NOT_RELAX_SYM_CONSTR>(_db, _vConstraints);
  if (verSolver.solve()) {
    verSolver.exportSolution();
  } else {
    Assert(false);
    return false;
  }
  return true;
}

BoolType CGLegalizer::wirelengthDrivenCompaction() {
  // Find the fixed boundary for this step
  findCellBoundary();
  // Horizontal
  //this->generateHorConstraints();
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
    return false;
  }

  // Vertical
  //this->generateVerConstraints();
  INF("CG legalizer:  wirelength-driven detailed placement vertical LP...\n");
  auto verSolver = lp_legalize::LpLegalizeWirelength<lp_legalize::LEGALIZE_VERTICAL_DIRECTION, lp_legalize::DO_NOT_RELAX_SYM_CONSTR>(_db, _vConstraints);
  verSolver.setBoundary(_hStar);
  bool verpass = verSolver.solve();
  if (verpass) {
    verSolver.exportSolution();
  } else {
    return false;
  }

#ifdef DEBUG_LEGALIZE
#ifdef DEBUG_DRAW
  _db.drawCellBlocks("./debug/after_dp.gds");
#endif
#endif
  return true;
}

BoolType CGLegalizer::preserveRelationCompaction() {
  generateRedundantConstraintGraph();
  INF("CG legalizer: Prerserve relational coodinate horizontal LP...\n");
  auto horSolver = lp_legalize::LpLegalizeWirelength<lp_legalize::LEGALIZE_HORIZONTAL_DIRECTION, lp_legalize::DO_NOT_RELAX_SYM_CONSTR>(_db, _hConstraints);
  if (horSolver.solve()) {
    horSolver.exportSolution();
  }
  else {
    return false;
  }
  INF("CG legalizer: Prerserve relational coodinate vetical LP...\n");
  auto verSolver = lp_legalize::LpLegalizeWirelength<lp_legalize::LEGALIZE_VERTICAL_DIRECTION, lp_legalize::DO_NOT_RELAX_SYM_CONSTR>(_db, _vConstraints);
  if (verSolver.solve()) {
    verSolver.exportSolution();
  }
  else {
    return false;
  }

  return true;
}

enum class Relation {
  RIGHT, LEFT, TOP, BOTTOM, // overlap or parallel run,
  TOP_RIGHT, TOP_LEFT, BOTTOM_RIGHT, BOTTOM_LEFT // Completely disjoined
};

Relation determineBoxRelation(const Box<LocType> &box1, const Box<LocType> &box2) {
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

void CGLegalizer::generateRedundantConstraintGraph() {
  _hConstraints.clear();
  _vConstraints.clear();
  for (IndexType ci1 = 0; ci1 < _db.numCells(); ++ci1) {
    const Box<LocType> box1 = _db.cell(ci1).cellBBoxOff();
    for (IndexType ci2 = ci1 + 1; ci2 < _db.numCells(); ++ci2) {
      const Box<LocType> box2 = _db.cell(ci2).cellBBoxOff();
      auto relation = determineBoxRelation(box1, box2);
      if (relation == Relation::LEFT) {
        _hConstraints.addConstraintEdge(ci2, ci1);
      }
      else if (relation == Relation::RIGHT) {
        _hConstraints.addConstraintEdge(ci1, ci2);
      }
      else if (relation == Relation::BOTTOM) {
        _vConstraints.addConstraintEdge(ci2, ci1);
      }
      else if (relation == Relation::TOP) {
        _vConstraints.addConstraintEdge(ci1, ci2);
      }
      else if (relation == Relation::BOTTOM_LEFT) {
        _hConstraints.addConstraintEdge(ci2, ci1);
        _vConstraints.addConstraintEdge(ci2, ci1);
      }
      else if (relation == Relation::BOTTOM_RIGHT) {
        _hConstraints.addConstraintEdge(ci1, ci2);
        _vConstraints.addConstraintEdge(ci2, ci1);
      }
      else if (relation == Relation::TOP_LEFT) {
        _hConstraints.addConstraintEdge(ci2, ci1);
        _vConstraints.addConstraintEdge(ci1, ci2);
      }
      else {
        Assert(relation == Relation::TOP_RIGHT);
        _vConstraints.addConstraintEdge(ci1, ci2);
        _hConstraints.addConstraintEdge(ci1, ci2);
      }
    }
  }
}


PROJECT_NAMESPACE_END
