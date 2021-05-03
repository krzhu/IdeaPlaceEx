/**
 * @file legalize_well.h
 * @brief Legalize the wells and in-well cells
 * @author Keren Zhu
 * @date 04/18/2021
 */

#pragma once

#include "db/Database.h"

PROJECT_NAMESPACE_BEGIN

class WellLegalizer {
public:
  explicit WellLegalizer(Database &db) : _db(db) {}
  BoolType legalize();
private:
  /// @brief add cell-edge spacing
  void legalizeCellEdgeSpacing();
  /// @brief legalize the NW spacing via patches for min step
  void legalizeMinStep();
private:
    Database &_db;
};

PROJECT_NAMESPACE_END
