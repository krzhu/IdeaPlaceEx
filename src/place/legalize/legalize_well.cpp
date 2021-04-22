#include "legalize_well.hpp"

PROJECT_NAMESPACE_BEGIN

BoolType WellLegalizer::legalize() {
  legalizeCellEdgeSpacing();
  return true;
}

void WellLegalizer::legalizeCellEdgeSpacing() {
  const LocType cellToWellEdgeSpacing = _db.parameters().cellToWellEdgeSpacing();
  for (IndexType wellIdx = 0; wellIdx < _db.vWells().size(); ++wellIdx) {
    auto &well = _db.well(wellIdx);
    // Here we generate a rectangle for each cell that satisfying the spacing
    // and then we merge them together into a new polygon
    // boost::geometry has a union_ function, but it is not very convenient
    // Here use the boost::polygon implementation
    typedef boost::polygon::property_merge_90<LocType, IntType> PropertyMergeType; //use int as property_type -> we don't care basically
    typedef boost::polygon::polygon_90_set_data<LocType> PolygonSetType; // Potentially we can use our own Polygon set implementation
    // But since we are using boost::geometry outside, it might be safe to explicitly converting boost::polygon polygon back
    typedef std::map<std::set<IntType>, PolygonSetType> PropertyMergeResultType;
    PropertyMergeType pm;
    PropertyMergeResultType result;
    pm.insert(well.shape(), 0);
    for (IndexType cellIdx : well.sCellIds()) {
      auto cellBox = _db.cell(cellIdx).cellBBoxOff();
      cellBox.expand(cellToWellEdgeSpacing);
      pm.insert(cellBox, 0);
    }
    pm.merge(result);
    std::vector<Polygon<LocType>> polyVec;
    (*result.begin()).second.get_polygons(polyVec);
    Assert(polyVec.size() == 1);
    well.setShape(polyVec[0].outer());
  }
}

PROJECT_NAMESPACE_END
