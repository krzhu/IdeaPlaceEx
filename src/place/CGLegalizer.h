/**
 * @file CGLegalizer.h
 * @brief The "legalization" solver with constraint graph + LP
 * @author Keren Zhu
 * @date 11/25/2019
 */

#ifndef IDEAPLACE_CG_LEGALIZER_H_
#define IDEAPLACE_CG_LEGALIZER_H_

#include "place/legalize/constraint_graph.hpp"
#include "db/Database.h"

PROJECT_NAMESPACE_BEGIN


class CGLegalizer {

public:
  /// @brief Constructor
  /// @param The database of IdeaPlaceEx
  explicit CGLegalizer(Database &db) : _db(db) {}
  /// @brief legalize the design
  bool legalize();

  /// @brief prepare the legalization
  void prepare();
  /// @brief Area-driven compaction
  /// @return true, successful
  BoolType areaDrivenCompaction();
  /// @brief LP-based detailed placement. For optimizing wire length
  BoolType wirelengthDrivenCompaction();
  /// @brief compaction and preserve current coordinate relation
  /// @param Extra spacing to add to the legalization
  BoolType preserveRelationCompaction(LocType extraSpacing = -1); 
  /// @brief legalzie with-in well 
  void legalizeWells();
  /// @brief Generate individual wells and legalize
  void generateIndividualWells();
  /// @brief close the well-aware functionality
  void closeWellAware() {
    _wellAware = false;
  }
  /// @brief open the well-aware legalization
  void openWellAware() {
    _wellAware = true;
  }
private:
  /// @brief Generate the constraints (not optimal in number of constraints).
  /// Based on sweeping algorithm
  void generateHorConstraints();
  void generateVerConstraints();

  void findCellBoundary();
  /// @brief generate constraint graph to preserve the current relation
  void generateRedundantConstraintGraph();

private:
  Database &_db;             ///< The database of IdeaPlaceEx
  Constraints _hConstraints; ///< The horizontal constraint edges
  Constraints _vConstraints; ///< The vertical constraint edges
  RealType _wStar; ///< The width from the objective function of the first LP
  RealType _hStar; ///< The width from the objective function of the first LP
  BoolType _wellAware = false;
};

PROJECT_NAMESPACE_END

#endif // IDEAPLACE_CG_LEGALIZER_H_
