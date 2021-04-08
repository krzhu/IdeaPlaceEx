#include "CGLegalizer.h"
#include "constraintGraphGeneration.h"
#include "pinassign/VirtualPinAssigner.h"

PROJECT_NAMESPACE_BEGIN

bool CGLegalizer::legalize() {

  VirtualPinAssigner pinAssigner(_db);

  // this->generateConstraints();

  auto legalizationStopWath = WATCH_CREATE_NEW("legalization");
  legalizationStopWath->start();

  this->generateVerConstraints();
  _hStar = lpLegalization(false);
  this->generateHorConstraints();
  _wStar = lpLegalization(true);
  legalizationStopWath->stop();

  LocType xMin = LOC_TYPE_MAX;
  LocType xMax = LOC_TYPE_MIN;
  LocType yMin = LOC_TYPE_MAX;
  LocType yMax = LOC_TYPE_MIN;
  for (IndexType cellIdx = 0; cellIdx < _db.numCells(); ++cellIdx) {
    auto cellBox = _db.cell(cellIdx).cellBBoxOff();
    xMin = std::min(xMin, cellBox.xLo());
    xMax = std::max(xMax, cellBox.xHi());
    yMin = std::min(yMin, cellBox.yLo());
    yMax = std::max(yMax, cellBox.yHi());
  }
  _wStar = std::max(0.0, static_cast<RealType>(xMax - xMin)) + 10;
  _hStar = std::max(0.0, static_cast<RealType>(yMax - yMin)) + 10;
  // this->generateConstraints();

  auto dpStopWatch = WATCH_CREATE_NEW("detailedPlacement");
  dpStopWatch->start();
  if (_db.parameters().ifUsePinAssignment()) {
    pinAssigner.solveFromDB();
  }
  if (!lpDetailedPlacement()) {
    INF("CG Legalizer: detailed placement fine tunning failed. Directly output "
        "legalization output. \n");
    dpStopWatch->stop();
    return true;
  }
  if (!lpDetailedPlacement()) {
    INF("CG Legalizer: detailed placement fine tunning failed. Directly output "
        "legalization output. \n");
    dpStopWatch->stop();
    return true;
  }
  dpStopWatch->stop();
  return true;

  INF("CG Legalizer: legalization finished\n");
  return true;
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
  sweepline.solve();
}
void CGLegalizer::generateVerConstraints() {
  _hConstraints.clear();
  _vConstraints.clear();
  // Init the irredundant constraint edges

  SweeplineConstraintGraphGenerator sweepline(_db, _hConstraints,
                                              _vConstraints);
  sweepline.solve();
}


RealType CGLegalizer::lpLegalization(bool isHor) {
  RealType obj;
  if (isHor) {
    INF("CG legalizer: legalize horizontal LP...\n");
    auto solver = LpLegalizeSolver(_db, _hConstraints, isHor);
    bool optimal = solver.solve();
    if (optimal) {
      solver.exportSolution();
    } else {
      return -1;
    }
    obj = solver.evaluateObj();
    //_db.drawCellBlocks("./debug/after_legalization_hor.gds");
  } else {
    INF("CG legalizer: legalize vertical LP...\n");
    auto solver = LpLegalizeSolver(_db, _vConstraints, isHor);
    bool optimal = solver.solve();
    if (optimal) {
      solver.exportSolution();
    } else {
      return -1;
    }
    obj = solver.evaluateObj();
  }
#ifdef DEBUG_LEGALIZE
#ifdef DEBUG_DRAW
  _db.drawCellBlocks("./debug/after_legalization.gds");
#endif
#endif

  return obj;
}

bool CGLegalizer::lpDetailedPlacement() {
  // Horizontal
  this->generateHorConstraints();
  INF("CG legalizer: detailed placement horizontal LP...\n");
  auto horSolver = LpLegalizeSolver(_db, _hConstraints, true, 1, 0);
#ifdef DEBUG_LEGALIZE
  DBG("wstar for width %f \n", _wStar);
#endif
  horSolver.setWStar(_wStar);
  bool horpass = horSolver.solve();
  if (horpass) {
    horSolver.exportSolution();
  } else {
    return false;
  }

  // Vertical
  this->generateVerConstraints();
  INF("CG legalizer: detailed placement vertical LP...\n");
  auto verSolver = LpLegalizeSolver(_db, _vConstraints, false, 1, 0);
  verSolver.setWStar(_hStar);
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

PROJECT_NAMESPACE_END
