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
private:
  class BoxEdge {
  public:
    explicit BoxEdge(LocType coord_, IndexType cellIdx_, bool isTop_)
        : coord(coord_), cellIdx(cellIdx_), isTop(isTop_) {}
    bool operator<(const BoxEdge &rhs) const {
      if (this->coord == rhs.coord) {
        if (this->isTop == rhs.isTop) {
          return this->cellIdx < rhs.cellIdx;
        } else {
          if (!this->isTop) {
            return false;
          } else {
            return true;
          }
        }
      } else {
        return this->coord < rhs.coord;
      }
    }
    std::string toStr() const {
      std::stringstream ss;
      ss << "LocType " << coord << " cellIdx " << cellIdx << " isTop " << isTop;
      return ss.str();
    }

  public:
    LocType coord;     ///< Coordinate of the edge
    IndexType cellIdx; ///< The index of the cell
    bool isTop;        ///< True: top/right edge. False: bottom/left edge
  };

public:
  /// @brief Constructor
  /// @param The database of IdeaPlaceEx
  explicit CGLegalizer(Database &db) : _db(db) {}
  /// @brief legalize the design
  bool legalize();

private:
  /// @brief Generate the constraints (not optimal in number of constraints).
  /// Based on sweeping algorithm
  void generateHorConstraints();
  void generateVerConstraints();
  /// @brief linear programming-based legalization
  /// @param if solving horizontal or vertical
  /// @return the resulting objective function. if negative, failed
  void lpLegalization(bool isHor);
  /// @brief LP-based detailed placement. For optimizing wire length
  bool lpDetailedPlacement();

private:
  Database &_db;             ///< The database of IdeaPlaceEx
  Constraints _hConstraints; ///< The horizontal constraint edges
  Constraints _vConstraints; ///< The vertical constraint edges
  RealType _wStar; ///< The width from the objective function of the first LP
  RealType _hStar; ///< The width from the objective function of the first LP
};

PROJECT_NAMESPACE_END

#endif // IDEAPLACE_CG_LEGALIZER_H_
