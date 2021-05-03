#include "legalize_well.hpp"

PROJECT_NAMESPACE_BEGIN

enum class JogOrient_t {
  NE = 0,
  SE = 1,
  SW = 2,
  NW = 3,
  INVALID = 4
};

bool clockwise(const Point<LocType>& p0, const Point<LocType>& p1, const Point<LocType>& p2, JogOrient_t& orient) {
  const LocType path1_deltaX = p1.x() - p0.x();
  const LocType path1_deltaY = p1.y() - p0.y();
  const LocType path2_deltaX = p2.x() - p1.x();
  const LocType path2_deltaY = p2.y() - p1.y();
  if (path1_deltaX > 0 and path1_deltaY == 0 and path2_deltaX == 0 and path2_deltaY < 0) {
    orient = JogOrient_t::NE;
    return true;
  }
  if (path1_deltaX == 0 and path1_deltaY < 0 and path2_deltaX < 0 and path2_deltaY == 0) {
    orient = JogOrient_t::SE;
    return true;
  }
  if (path1_deltaX < 0 and path1_deltaY == 0 and path2_deltaX == 0 and path2_deltaY > 0) {
    orient = JogOrient_t::SW;
    return true;
  }
  if (path1_deltaX == 0 and path1_deltaY > 0 and path2_deltaX > 0 and path2_deltaY == 0) {
    orient = JogOrient_t::NW;
    return true;
  }
  orient = JogOrient_t::INVALID;
  return false;
}

bool counterClockwise(const Point<LocType>& p0, const Point<LocType>& p1, const Point<LocType>& p2, JogOrient_t& orient) {
  const LocType path1_deltaX = p1.x() - p0.x();
  const LocType path1_deltaY = p1.y() - p0.y();
  const LocType path2_deltaX = p2.x() - p1.x();
  const LocType path2_deltaY = p2.y() - p1.y();
  if (path1_deltaX == 0 and path1_deltaY < 0 and path2_deltaX > 0 and path2_deltaY == 0) {
    orient = JogOrient_t::NE;
    return true;
  }
  if (path1_deltaX > 0 and path1_deltaY == 0 and path2_deltaX == 0 and path2_deltaY > 0) {
    orient = JogOrient_t::NW;
    return true;
  }
  if (path1_deltaX == 0 and path1_deltaY > 0 and path2_deltaX < 0 and path2_deltaY == 0) {
    orient = JogOrient_t::SW;
    return true;
  }
  if (path1_deltaX < 0 and path1_deltaY == 0 and path2_deltaX == 0 and path2_deltaY < 0) {
    orient = JogOrient_t::SE;
    return true;
  }
  orient = JogOrient_t::INVALID;
  return false;
}

/// @brief patch the polygon to eliminate concave jogs. 
/// @return if detect jog in the process
bool patchConcaveJogs(LocType minStep, Polygon<LocType>& poly) {
  bool hasJog = false;
  const auto &ring = poly.outer();
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
  pm.insert(poly, 0);
  for (IndexType j = 1; j < ring.size(); ++j) {
    const auto& pt0 = ring[j - 1];
    const auto& pt1 = ring[j];
    const auto& pt2 = j + 1 == ring.size() ? ring[0] : ring[j + 1];
    if (::klib::manhattanDistance(pt0, pt1) < minStep
        and ::klib::manhattanDistance(pt1, pt2) < minStep) {
      auto box = Box<LocType>(
                     std::min({pt0.x(), pt1.x(), pt2.x()}),
                     std::min({pt0.y(), pt1.y(), pt2.y()}),
                     std::max({pt0.x(), pt1.x(), pt2.x()}),
                     std::max({pt0.y(), pt1.y(), pt2.y()}));
      JogOrient_t orient;
      if (counterClockwise(pt0, pt1, pt2, orient)) { // concave
        hasJog = true;
        pm.insert(box, 0);
        assert(orient != JogOrient_t::INVALID);
      }
    }
  }
  pm.merge(result);
  std::vector<Polygon<LocType>> polyVec;
  (*result.begin()).second.get_polygons(polyVec);
  Assert(polyVec.size() == 1);
  poly.setOuter(polyVec[0].outer());
  return hasJog;
}

/// @brief patch the polygon to eliminate convex jogs
/// @return if detect jog in the process
bool patchConvexJogs(LocType minStep, Polygon<LocType>& poly) {
  bool hasJog = false;
  const auto &ring = poly.outer();
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
  pm.insert(poly, 0);
  for (IndexType j = 1; j < ring.size(); ++j) {
    const auto& pt0 = ring[j - 1];
    const auto& pt1 = ring[j];
    const auto& pt2 = j + 1 == ring.size() ? ring[0] : ring[j + 1];
    if (::klib::manhattanDistance(pt0, pt1) < minStep
        and ::klib::manhattanDistance(pt1, pt2) < minStep) {
      auto box = Box<LocType>(
                     std::min({pt0.x(), pt1.x(), pt2.x()}),
                     std::min({pt0.y(), pt1.y(), pt2.y()}),
                     std::max({pt0.x(), pt1.x(), pt2.x()}),
                     std::max({pt0.y(), pt1.y(), pt2.y()}));
      JogOrient_t orient;
      if (clockwise(pt0, pt1, pt2, orient)) { // convex
        hasJog = true;
        //_cir.addMaskWire(box, i);
        switch(orient) {
          case JogOrient_t::NE: box.setYH(box.yh() + minStep); break;
          case JogOrient_t::SE: box.setYL(box.yl() - minStep); break;
          case JogOrient_t::SW: box.setYL(box.yl() - minStep); break;  
          case JogOrient_t::NW: box.setYH(box.yh() + minStep); break;
          default: assert(false);
        }
        pm.insert(box, 0);
      }
    }
  }
  pm.merge(result);
  std::vector<Polygon<LocType>> polyVec;
  (*result.begin()).second.get_polygons(polyVec);
  Assert(polyVec.size() == 1);
  poly.setOuter(polyVec[0].outer());
  return hasJog;
}


/// @brief patch the polygon to eliminate u-shape spacing error
/// @return if detect jog in the process
bool patchSpacingJogs(LocType minStep, Polygon<LocType>& poly) {
  bool hasJog = false;
  const auto &ring = poly.outer();
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
  pm.insert(poly, 0);
  for (IndexType j = 1; j < ring.size(); ++j) {
    const auto& pt0 = ring[j - 1];
    const auto& pt1 = ring[j];
    const auto& pt2 = j + 1 == ring.size() ? ring[0] : ring[j + 1];
    const auto& pt3 = j + 2 == ring.size() ? ring[1] : ring[j + 2];
    if (::klib::manhattanDistance(pt1, pt2) < minStep 
        or (::klib::manhattanDistance(pt0, pt1) < minStep
          and ::klib::manhattanDistance(pt2, pt3) < minStep)) {
      auto box = Box<LocType>(
                     std::min({pt0.x(), pt1.x(), pt2.x(), pt3.x()}),
                     std::min({pt0.y(), pt1.y(), pt2.y(), pt3.y()}),
                     std::max({pt0.x(), pt1.x(), pt2.x(), pt3.x()}),
                     std::max({pt0.y(), pt1.y(), pt2.y(), pt3.y()}));
      JogOrient_t orientA, orientB;
      if (counterClockwise(pt0, pt1, pt2, orientA) and counterClockwise(pt1, pt2, pt3, orientB) ) { // concave
        hasJog = true;
        pm.insert(box, 0);
        Assert(orientA != JogOrient_t::INVALID);
        Assert(orientB != JogOrient_t::INVALID);
      }
    }
  }
  pm.merge(result);
  std::vector<Polygon<LocType>> polyVec;
  (*result.begin()).second.get_polygons(polyVec);
  Assert(polyVec.size() == 1);
  poly.setOuter(polyVec[0].outer());
  return hasJog;
}


BoolType WellLegalizer::legalize() {
  legalizeCellEdgeSpacing();
  legalizeMinStep();
  return true;
}

void WellLegalizer::legalizeMinStep() {
  if (not _db.tech().isNwellLayerSet()) {
    return;
  }
  const LocType minStep = _db.tech().spacingRule(_db.tech().nwellLayerIdx());
  for (IndexType wellIdx = 0; wellIdx < _db.vWells().size(); ++wellIdx) {
    auto &well = _db.well(wellIdx);
    while (1) {
      bool hasJog = false;
      if (patchConcaveJogs(minStep, well.shape())) {
        hasJog = true;
        continue;
      }

      if (patchConvexJogs(minStep, well.shape())) {
        hasJog = true;
        continue;
      }
      if (patchSpacingJogs(minStep, well.shape())) {
        hasJog = true;
        continue;
      }
      if (!hasJog)
        break;
    }
  }
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
