#include "CGLegalizer.h"
#include "constraintGraphGeneration.h"
#include "pinassign/VirtualPinAssigner.h"
#include "place/legalize/lp_legalize.hpp"

PROJECT_NAMESPACE_BEGIN

bool CGLegalizer::legalize() {

  VirtualPinAssigner pinAssigner(_db);

  // this->generateConstraints();
  _db.offsetLayout();

  auto legalizationStopWath = WATCH_CREATE_NEW("legalization");
  legalizationStopWath->start();

  this->generateVerConstraints();
  lpLegalization(false);
  this->generateHorConstraints();
  lpLegalization(true);
  legalizationStopWath->stop();

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


void CGLegalizer::lpLegalization(bool isHor) {
  if (isHor) {
    INF("CG legalizer: legalize horizontal LP...\n");
    auto solver = lp_legalize::LpLegalizeArea<lp_legalize::LEGALIZE_HORIZONTAL_DIRECTION, lp_legalize::DO_NOT_RELAX_SYM_CONSTR>(_db, _hConstraints);
    bool optimal = solver.solve();
    if (optimal) {
      solver.exportSolution();
    } else {
      return;
    }
  } else {
    INF("CG legalizer: legalize vertical LP...\n");
    auto solver = lp_legalize::LpLegalizeArea<lp_legalize::LEGALIZE_VERTICAL_DIRECTION, lp_legalize::DO_NOT_RELAX_SYM_CONSTR>(_db, _vConstraints);
    bool optimal = solver.solve();
    if (optimal) {
      solver.exportSolution();
    } else {
      return;
    }
  }
#ifdef DEBUG_LEGALIZE
#ifdef DEBUG_DRAW
  _db.drawCellBlocks("./debug/after_legalization.gds");
#endif
#endif

}

bool CGLegalizer::lpDetailedPlacement() {
  // Horizontal
  this->generateHorConstraints();
  INF("CG legalizer: detailed placement horizontal LP...\n");
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
  this->generateVerConstraints();
  INF("CG legalizer: detailed placement vertical LP...\n");
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


PROJECT_NAMESPACE_END
