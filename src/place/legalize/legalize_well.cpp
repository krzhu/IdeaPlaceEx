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
bool patchSpacingJogs(LocType minStep, Polygon<LocType>& poly, boost::geometry::index::rtree<Box<LocType>, boost::geometry::index::rstar<8, 2>> &rtree) {
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
    if (pt0 == pt3) {
      continue;
    }
    auto box = Box<LocType>(
                   std::min({pt0.x(), pt1.x(), pt2.x(), pt3.x()}),
                   std::min({pt0.y(), pt1.y(), pt2.y(), pt3.y()}),
                   std::max({pt0.x(), pt1.x(), pt2.x(), pt3.x()}),
                   std::max({pt0.y(), pt1.y(), pt2.y(), pt3.y()}));
    std::vector<Box<LocType>> queryResults;
    rtree.query(boost::geometry::index::intersects(box), std::back_inserter(queryResults));
    if (::klib::manhattanDistance(pt1, pt2) < minStep 
        or (::klib::manhattanDistance(pt0, pt1) < minStep
          and ::klib::manhattanDistance(pt2, pt3) < minStep)
        or not queryResults.empty()) {
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

BoolType WellLegalizer::legalizeAndAddContact() {
  legalizeCellEdgeSpacing();
  legalizeContact();
  legalizeMinStep();
  return true;
}

void WellLegalizer::legalizeMinStep() {
  if (not _db.tech().isNwellLayerSet()) {
    return;
  }
  using rtree_type  = boost::geometry::index::rtree<Box<LocType>, boost::geometry::index::rstar<8, 2>>;
  const LocType minStep = _db.tech().spacingRule(_db.tech().nwellLayerIdx());
  for (IndexType wellIdx = 0; wellIdx < _db.vWells().size(); ++wellIdx) {
    auto &well = _db.well(wellIdx);
    std::vector<Box<LocType>> outWellDevicesBoxes;
    for (IndexType cellIdx = 0; cellIdx < _db.numCells(); ++cellIdx) {
      if (well.sCellIds().find(cellIdx) == well.sCellIds().end()) {
        outWellDevicesBoxes.emplace_back(_db.cell(cellIdx).cellBBoxOff());
      }
    }
    rtree_type outDeviceRtree(outWellDevicesBoxes); // Store the out of well devices shapes
    while (1) {
      bool hasJog = false;
      if (patchConcaveJogs(minStep, well.shape())) {
        hasJog = true;
        continue;
      }

      if (patchConvexJogs(std::max(minStep, _db.parameters().smallestWellPolyStep()), well.shape())) { 
        hasJog = true;
        continue;
      }
      if (patchSpacingJogs(std::max(minStep, _db.parameters().smallestWellPolyStep()), well.shape(), outDeviceRtree)) {
        hasJog = true;
        continue;
      }
      if (!hasJog)
        break;
    }
  }
}

void WellLegalizer::generateIndividualWells() {
  _db.clearWells();
  for (IndexType cellIdx = 0; cellIdx < _db.numCells(); ++cellIdx) {
    auto &cell = _db.cell(cellIdx);
    if (not cell.needWell()) {
      continue;
    }
    auto cellBox = cell.cellBBoxOff();
    auto wpeSpacing = _db.wpeSpacing(cellIdx);
    cellBox.expandX(wpeSpacing.first);
    cellBox.expandY(wpeSpacing.second);
    auto wellIdx = _db.allocateWell();
    auto &well = _db.well(wellIdx);
    well.setShape(Polygon<LocType>(cellBox).outer());
    cell.setWellIdx(wellIdx);
    well.addCellIdx(cellIdx);
  }
  legalizeContact();
  legalizeMinStep();
}

void WellLegalizer::legalizeCellEdgeSpacing() {
  for (IndexType wellIdx = 0; wellIdx < _db.vWells().size(); ++wellIdx) {
    auto &well = _db.well(wellIdx);
    // Here we generate a rectangle for each cell that satisfying the spacing
    // and then we merge them together into a new polygon
    // boost::geometry has a union_ function, but it is not very convenient
    // Here use the boost::polygon implementation
    typedef boost::polygon::property_merge_90<LocType, IntType> PropertyMergeType; //use int as property_type -> we don't care basically
    typedef boost::polygon::polygon_90_set_data<LocType> PolygonSetType; // Potentially we can use our own Polygon set implementation
    // But since we are using boost::geometry outside, it might be safe to explicitly converting boost::polygon polygon back
    typedef std::map<std::vector<IntType>, PolygonSetType> PropertyMergeResultType;
    PropertyMergeType pm;
    PropertyMergeResultType result;
    boost::geometry::correct(well.shape());
    pm.insert(well.shape(), 0);
    for (IndexType cellIdx : well.sCellIds()) {
      auto cellBox = _db.cell(cellIdx).cellBBoxOff();
      auto wpeSpacing = _db.wpeSpacing(cellIdx);
      cellBox.expandX(wpeSpacing.first);
      cellBox.expandY(wpeSpacing.second);
      pm.insert(cellBox, 0);
    }
    pm.merge(result);
    std::vector<Polygon<LocType>> polyVec;
    (*result.begin()).second.get_polygons(polyVec);
    //Assert(polyVec.size() == 1);
    if (polyVec.size() != 1) {
      ERR("WellLegalizer::Well away from cell. WellIdx %d \n", wellIdx);
      ERR("Result size %d  vec %d\n", result.size(), polyVec.size());
      well.printInfo();
      Box<LocType> bbox(LOC_TYPE_MAX, LOC_TYPE_MAX, LOC_TYPE_MIN, LOC_TYPE_MIN);
      for (IndexType cellIdx : well.sCellIds()) {
        auto cellBox = _db.cell(cellIdx).cellBBoxOff();
        auto wpeSpacing = _db.wpeSpacing(cellIdx);
        cellBox.expandX(wpeSpacing.first);
        cellBox.expandY(wpeSpacing.second);
        bbox.unionBox(cellBox);
      }
      for (const auto &pt : well.shape().outer()) {
        bbox.coverPoint(pt);
      }
      Ring<LocType> ring({
          Point<LocType>(bbox.xLo(), bbox.yLo()),
          Point<LocType>(bbox.xLo(), bbox.yHi()),
          Point<LocType>(bbox.xHi(), bbox.yHi()),
          Point<LocType>(bbox.xHi(), bbox.yLo()),
          Point<LocType>(bbox.xLo(), bbox.yLo())
          });
      well.setShape(ring);
      continue;
    }
    well.setShape(polyVec[0].outer());
  }
}

namespace _legalize_well_vdd_contact_details {


  static constexpr IndexType MAX_CONTACT_CANDIDATE = 2000;

  void topologicalSort(const std::vector<std::vector<IndexType>> &outEdges, std::vector<IndexType> &results) {
    std::vector<char> visited(outEdges.size(), false);
    std::function<void(IndexType)> util = [&](IndexType i) {
      visited.at(i) = true;
      for (IndexType idx : outEdges.at(i)) {
        if (not visited.at(idx)) {
          util(idx);
        }
      }
      results.emplace_back(i);
    };
    util(outEdges.size() - 2);  // We know this is source, don't need to iterate through all nodes
  }

  void selectVddContactCandidates(const std::vector<std::tuple<Box<LocType>, IndexType, LocType, IntType, IndexType>> &candidates, 
      std::pair<IntType, IntType> costThreshold, 
      LocType contactSpacing,
      std::vector<IndexType> &results) {
    Assert(candidates.size() > 0);

    // Figure out the valid candidates with in the cost range
    IndexType startIdx = 0;
    IndexType numCandidates = 0;
    for (IndexType idx = 0; idx < candidates.size(); ++idx) {
      const auto cost = std::get<2>(candidates.at(idx));
      if (cost < costThreshold.first) {
        ++startIdx;
      }
      if (cost > costThreshold.second) {
        numCandidates = idx - startIdx;
        Assert(idx> 0);
        break;
      }
    }
    if (numCandidates == 0) {
      numCandidates = candidates.size() - startIdx; //Upper threshold is more than max(candidates cost)
    }
    // Construct the graph
    std::vector<std::vector<IndexType>> outEdges;
    outEdges.resize(numCandidates + 2);
    const IndexType sourceIdx = numCandidates; // S
    const IndexType targetIdx = numCandidates + 1; // T
    for (IndexType idx = 0; idx < numCandidates; ++idx) {
      outEdges.at(sourceIdx).emplace_back(idx);
      outEdges.at(idx).emplace_back(targetIdx);
    }
    // Find the independent shape and give them a directed edge
    for (IndexType i = 0; i < numCandidates; ++i) {
      const auto &boxI = std::get<0>(candidates.at(i + startIdx));
      for (IndexType j = i+1; j < numCandidates; ++j) {
        const auto &boxJ = std::get<0>(candidates.at(j + startIdx));
        if (boxI.xLo() >= boxJ.xHi() + contactSpacing and boxI.yLo() >= boxJ.yHi() + contactSpacing) {
          // i is right to j
          outEdges.at(j).emplace_back(i);
        }
        else if (boxJ.xLo() >= boxI.xHi() + contactSpacing and boxJ.yLo() >= boxI.yHi() + contactSpacing) {
          // j is right to i
          outEdges.at(i).emplace_back(j);
        }
      }
    }
    // Topological sort
    std::vector<IndexType> sorted; // In reverse order
    topologicalSort(outEdges, sorted);
    Assert(sorted.front() == targetIdx);

    // Find longest path
    std::vector<IntType> dist(numCandidates + 2, -1);
    dist.at(sourceIdx) = 0;
    std::vector<IndexType> prevNode(numCandidates + 2, INDEX_TYPE_MAX);
    for (IndexType i = 0; i < sorted.size() - 1; ++i) { // sorted is in reverse order. Skip target
      IndexType idx = sorted.at(sorted.size() - 1 - i);
      IntType currentWeight = 0;
      if (i != 0) { // i=0 -> source node -> give it 0
        currentWeight = std::get<3>(candidates.at(idx + startIdx));
      }
      Assert(dist.at(idx) >= 0);
      for (IndexType outIdx : outEdges.at(idx)) {
        if (dist.at(outIdx) < dist.at(idx) + currentWeight) {
          dist.at(outIdx) = dist.at(idx) + currentWeight;
          prevNode.at(outIdx) = idx;
        }
      }
    }
    // Extract the result
    IndexType nodeIdx = targetIdx;
    while(prevNode.at(nodeIdx) != sourceIdx) {
      if (nodeIdx != targetIdx) {
        results.emplace_back(nodeIdx + startIdx);
      }
      nodeIdx = prevNode.at(nodeIdx);
    }
    results.emplace_back(nodeIdx + startIdx);
    Assert(nodeIdx != sourceIdx and nodeIdx != targetIdx);
  }

  /// @brief Try add as many contacts as possible, but make sure add at least one
  void legalizeContactByAddingAsManyContactsAsPossible(Database &db, IndexType wellIdx) {
    using rtree_type  = boost::geometry::index::rtree<Box<LocType>, boost::geometry::index::rstar<8, 2>>;
    Well &well = db.well(wellIdx);
    const LocType contactSpacing = db.tech().vddContactSpacing();
    const Polygon<LocType> poly = well.shape();
    std::vector<Box<LocType>> deviceBoxs;
    for (IndexType cellIdx : well.sCellIds()) {
      deviceBoxs.emplace_back(db.cell(cellIdx).cellBBoxOff());
    }
    rtree_type deviceRtree(deviceBoxs); // Store the shape of device
    std::vector<Box<LocType>> outWellDevicesBoxes;
    for (IndexType cellIdx = 0; cellIdx < db.numCells(); ++cellIdx) {
      if (well.sCellIds().find(cellIdx) == well.sCellIds().end()) {
        outWellDevicesBoxes.emplace_back(db.cell(cellIdx).cellBBoxOff());
      }
    }
    rtree_type outDeviceRtree(outWellDevicesBoxes); // Store the out of well devices shapes
    // Generate the contact candidates
    const auto currentWellBBox = well.boundingBox();
    const LocType genStep = std::max(db.parameters().gridStep(), 1);
    const LocType genXLo = genStep * static_cast<LocType>((currentWellBBox.xLo() - contactSpacing)/ genStep);
    const LocType genYLo = genStep * static_cast<LocType>((currentWellBBox.yLo() - contactSpacing)/ genStep);
    const LocType genXHi = genStep * static_cast<LocType>((currentWellBBox.xHi() + contactSpacing)/ genStep + 1);
    const LocType genYHi = genStep * static_cast<LocType>((currentWellBBox.yHi() + contactSpacing)/ genStep + 1);
    std::vector<std::tuple<Box<LocType>,IndexType, LocType, IntType, IndexType>> candidates; // 0=shape, 1= idx, 2= cost, 3=weight 4=templateIdx
    IndexType candidateIdx = 0;
    const LocType genStepX = std::max(genStep, (genXHi - genXLo) / 200);
    const LocType genStepY = std::max(genStep, (genYHi - genYLo) / 200);
    for (IndexType contactTemplateIdx = 0; contactTemplateIdx < db.tech().numVddContactTemplates(); ++contactTemplateIdx) {
      const Box<LocType> contactTemplate = db.tech().vddContactTemplate(contactTemplateIdx);
      for (IntType xIdx = 0; xIdx <= (genXHi - genXLo) / genStepX; ++xIdx) {
        for (IntType yIdx = 0; yIdx <= (genYHi - genYLo) / genStepX; ++yIdx) {
          Box<LocType> contactShape = contactTemplate.offsetBox(Point<LocType>(xIdx * genStepX  + genXLo, yIdx * genStepY + genYLo));
          auto queryBox = contactShape;
          queryBox.expand(contactSpacing);
          std::vector<Box<LocType>> queryResults;
          deviceRtree.query(boost::geometry::index::intersects(queryBox), std::back_inserter(queryResults));
          if (not queryResults.empty()) {
            LocType dif = std::max(queryResults.front().yHi() - queryBox.yLo(), 0);
            yIdx += dif / genStepY;
            continue;
          }
          std::vector<Polygon<LocType>> intersects;
          boost::geometry::intersection(contactShape, poly, intersects); // Get the intersection between well and contact

          LocType intersectArea = 0;
          for (const auto & inter : intersects) {
            intersectArea += boost::geometry::area(inter);
          }
          std::vector<Box<LocType>> outWellIntersects;
          outDeviceRtree.query(boost::geometry::index::intersects(queryBox), std::back_inserter(outWellIntersects)); // Get the intersection between other devices and contact
          LocType outWellDeviceIntersectArea = 0;
          for (const auto &inter : outWellIntersects) {
            outWellDeviceIntersectArea += boost::geometry::area(inter);
          }
          candidates.emplace_back(std::make_tuple(contactShape, candidateIdx, contactShape.area() - intersectArea + outWellDeviceIntersectArea, db.tech().vddContactWeight(contactTemplateIdx), contactTemplateIdx));
          ++candidateIdx;
        }
      }
    }
    AssertMsg(candidates.size() > 0, "Ok, then we have to have a loop to handle no valid candidate situation \n");
    std::stable_sort(candidates.begin(), candidates.end(), [&](const auto &i, const auto &j){ return std::get<2>(i) < std::get<2>(j);});
    for (IndexType candIdx = 0; candIdx < candidates.size(); ++candIdx) {
      std::get<1>(candidates.at(candIdx)) = candIdx;
    }
    IndexType numZeroCost = 0;
    std::pair<LocType, LocType> candidateThresholds;
    for (IndexType candIdx = 0; candIdx < candidates.size(); ++candIdx) {
      if (std::get<2>(candidates.at(candIdx)) == 0) {
        numZeroCost += 1;
      }
    }
    if (numZeroCost == 0) {
      WRN("Legalize well:: Add non zero cost contact \n");
      if (candidates.size() > MAX_CONTACT_CANDIDATE) {
        candidateThresholds = std::make_pair(0, std::get<2>(candidates.at(MAX_CONTACT_CANDIDATE)));
      }
      else {
        candidateThresholds = std::make_pair(0, LOC_TYPE_MAX);
      }
    }
    else {
      candidateThresholds = std::make_pair(0, 0);
    }
    std::vector<IndexType> selected;
    selectVddContactCandidates(candidates, candidateThresholds, contactSpacing, selected);
    Assert(selected.size() > 0);
    typedef boost::polygon::property_merge_90<LocType, IntType> PropertyMergeType; //use int as property_type -> we don't care basically
    typedef boost::polygon::polygon_90_set_data<LocType> PolygonSetType; // Potentially we can use our own Polygon set implementation
    // But since we are using boost::geometry outside, it might be safe to explicitly converting boost::polygon polygon back
    typedef std::map<std::set<IntType>, PolygonSetType> PropertyMergeResultType;
    PropertyMergeType pm;
    PropertyMergeResultType result;
    boost::geometry::correct(well.shape());
    pm.insert(well.shape(), 0);
    Box<LocType> bigBBox(LOC_TYPE_MAX, LOC_TYPE_MAX, LOC_TYPE_MIN, LOC_TYPE_MIN);
    for (IndexType idx : selected) {
      Box<LocType> box = std::get<0>(candidates.at(idx));
      const IndexType templIdx = std::get<4>(candidates.at(idx));
      const Box<LocType> &templBox = db.tech().vddContactTemplate(templIdx);
      const Point<LocType> pos = box.ll() -  templBox.ll();
      well.addVddContact(templIdx, pos);
      box.expand(contactSpacing);
      bigBBox.unionBox(box);
      pm.insert(box, 0);
    }
    pm.merge(result);
    std::vector<Polygon<LocType>> polyVec;
    (*result.begin()).second.get_polygons(polyVec);
    if (polyVec.size() != 1) {
      for (const auto &pt : well.shape().outer()) {
        bigBBox.coverPoint(pt);
      }
      Ring<LocType> ring({
          Point<LocType>(bigBBox.xLo(), bigBBox.yLo()),
          Point<LocType>(bigBBox.xLo(), bigBBox.yHi()),
          Point<LocType>(bigBBox.xHi(), bigBBox.yHi()),
          Point<LocType>(bigBBox.xHi(), bigBBox.yLo()),
          Point<LocType>(bigBBox.xLo(), bigBBox.yLo())
          });
      well.setShape(ring);
    }
    else {
      well.setShape(polyVec[0].outer());
    }
  }

} // namespace _legalize_well_vdd_contact_details

void WellLegalizer::legalizeContact() {
  if (_db.tech().numVddContactTemplates() == 0) {
    return;
  }
  for (IndexType wellIdx = 0; wellIdx < _db.vWells().size(); ++wellIdx) {
    _legalize_well_vdd_contact_details::legalizeContactByAddingAsManyContactsAsPossible(_db, wellIdx);
  }
}

PROJECT_NAMESPACE_END
