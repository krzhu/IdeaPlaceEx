#include "CGLegalizer.h"

PROJECT_NAMESPACE_BEGIN

/// @brief Find which direction is the least displacement direction to make the two boxes disjoint
/// @return 1 if moving box2 left
/// @return 2 if moving box2 right
/// @return 3 if moving box2 down
/// @return 4 if moving box2 up
IntType minOverlappingDirection(const Box<LocType> &box1, const Box<LocType> &box2)
{
    //Assert(box1.overlap(box2));
    LocType dispL = box2.xHi() - box1.xLo();
    LocType dispR = box1.xHi() - box2.xLo();
    LocType dispB = box2.yHi() - box1.yLo();
    LocType dispT = box1.yHi() - box2.yLo();
    LocType minDisp = std::min(std::min(dispL, dispR), std::min(dispB, dispT));
    if (minDisp == dispL) { return 1; }
    if (minDisp == dispR) { return 2; }
    if (minDisp == dispB) { return 3; }
    if (minDisp == dispT) { return 4; }
    return 0;
}


void CGLegalizer::generateConstraints()
{
    // Init the irredundant constraint edges
    
    // Sort the edges of cells
    std::vector<BoxEdge> hBoxEdges, vBoxEdges;
    for (IndexType cellIdx = 0; cellIdx < _db.numCells(); ++cellIdx)
    {
        const auto &cell = _db.cell(cellIdx);
        hBoxEdges.emplace_back(BoxEdge(cell.yLoc(), cellIdx, false));
        hBoxEdges.emplace_back(BoxEdge(cell.yHi(), cellIdx, true));
        vBoxEdges.emplace_back(BoxEdge(cell.xLoc(), cellIdx, false));
        vBoxEdges.emplace_back(BoxEdge(cell.xHi(), cellIdx, true));
    }
    std::sort(hBoxEdges.begin(), hBoxEdges.end());
    std::sort(vBoxEdges.begin(), vBoxEdges.end());

    // Add edges
    // Also calculate overlapAny. overlapAny[cellIdx] = true if the cell is overlapping with any other
    std::vector<bool> overlapAny(_db.numCells(), false);
    for (IndexType idx1 = 0; idx1 < _db.numCells(); ++idx1)
    {
        const auto &cell1 = _db.cell(idx1);
        auto cellBox1 = _db.cell(idx1).cellBBoxOff();
        for (IndexType idx2 = idx1 + 1; idx2 < _db.numCells(); ++idx2)
        {
            const auto &cell2 = _db.cell(idx2);
            auto cellBox2 = _db.cell(idx2).cellBBoxOff();
            // overlapAny part
            if (cellBox1.overlap(cellBox2))
            {
                AT(overlapAny, idx1) = true;
                AT(overlapAny, idx2) = true;
            }
            // Add edge part
            // First detect whether the pair are sym pair or are both self-symmetric
            bool symPairFlag = false;
            bool selfSymFlag = false;
            if (_db.cell(idx1).hasSymPair() && _db.cell(idx2).hasSymPair())
            {
                if (cell1.symPairCellIdx() == idx2 && cell2.symPairCellIdx() == idx1)
                {
                    symPairFlag = true;
                }
            }
            if (cell1.hasSelfSym() && cell2.hasSelfSym())
            {
                selfSymFlag = true;
            }
            // If sym pair add them in the following way
            if (symPairFlag)
            {
                // Only add horizontal edges
                // Assume all horizontally symmetric
                if (cell1.xLoc() <  cell2.xLoc())
                {
                    _hConstraints.addConstraintEdge(idx1, idx2, - cell1.cellBBox().xLen());
                }
                else
                {
                    _hConstraints.addConstraintEdge(idx2, idx1, - cell2.cellBBox().xLen());
                }
                // Do not proceed into other types of constraint edges
                continue;
            }
            // Both are self symmetric cells
            if (selfSymFlag)
            {
                // Only add vertical edges
                // Add edges based on the displacement distance of moving one over another
                LocType dispB = cell2.yHi() - cell1.yLoc();
                LocType dispT = cell1.yHi() - cell2.yLoc();
                if (dispB < dispT)
                {
                    _vConstraints.addConstraintEdge(idx1, idx2, - cell1.cellBBox().yLen());
                }
                else
                {
                    _vConstraints.addConstraintEdge(idx2, idx1, - cell2.cellBBox().yLen());
                }
            }
            IntType dispDirection = minOverlappingDirection(cellBox1, cellBox2);
            if (dispDirection == 1)
            {
                _hConstraints.addConstraintEdge(idx2, idx1, -cellBox2.xLen());
            }
            else if (dispDirection == 2)
            {
                _hConstraints.addConstraintEdge(idx1, idx2, -cellBox1.xLen());
            }
            else if (dispDirection == 3)
            {
                _vConstraints.addConstraintEdge(idx2, idx1, -cellBox2.yLen());
            }
            else if (dispDirection == 4)
            {
                _vConstraints.addConstraintEdge(idx1, idx2, -cellBox1.yLen());
            }
        }
    }
    
    // Execute the saved BoxEdge. 
    std::vector<IndexType> vOrders, hOrders;
    // "the right hand side candidates". The direction relation between vertices
    // The two additional vertices are virtural source and target
    std::vector<IntType> hCand(_db.numCells() + 2, -1), vCand(_db.numCells() + 2, -1);

    for (const auto &ev : hBoxEdges)
    {
        if (ev.isTop)
        {
            initIrredundantEdgesDelete(true, hOrders, hCand, ev.cellIdx, overlapAny);
        }
        else
        {
            initIrredundantEdgesInsert(true, hOrders, hCand, ev.cellIdx);
        }
    }
    for (const auto &ev : vBoxEdges)
    {
        if (ev.isTop)
        {
            initIrredundantEdgesDelete(false, vOrders, vCand, ev.cellIdx, overlapAny);
        }
        else
        {
            initIrredundantEdgesInsert(false, vOrders, vCand, ev.cellIdx);
        }
    }
#ifdef DEBUG_LEGALIZE
    DBG("Print the generated constraints... \n\n");
    for (IndexType i = 0; i < _hConstraints.edges().size(); ++i)
    {
        const auto &edge = _hConstraints.edges()[i];
        DBG("horizontal i %d %s \n", i, edge.toStr().c_str());
    }
    for (IndexType i = 0; i < _vConstraints.edges().size(); ++i)
    {
        const auto &edge = _vConstraints.edges()[i];
        DBG("vertical i %d %s \n", i, edge.toStr().c_str());
    }
#endif
    // Construct the boost constraint graphs
    constructConstraintGraphs();
    // Minimize the DAG
    dagTransitiveReduction(_hCG);
    dagTransitiveReduction(_vCG);

}

void CGLegalizer::dagTransitiveReduction(ConstraintGraph &cg)
{
    // DFS-based transitive reduction N(N+M)
    IndexType numNodes = boost::num_vertices(cg.boostGraph());
    // Build edge matrix
    Vector2D<IntType> edgeMat(numNodes, numNodes), reachable(numNodes, numNodes);
    auto edges = boost::edges(cg.boostGraph());
    for (auto it = edges.first; it != edges.second; ++it)
    {
        IntType sourceNode = boost::source(*it, cg.boostGraph());
        IntType targetNode = boost::target(*it, cg.boostGraph());
        edgeMat.at(sourceNode, targetNode) = true;
    }
    // DFS remove edge
    std::vector<bool> visited(numNodes, false);
    ConstraintGraph::IndexMap idxMap = boost::get(boost::vertex_index, cg.boostGraph());
    dfsRemoveTransitiveEdge(cg, edgeMat, _db.numCells(), visited, reachable, idxMap);
}

bool CGLegalizer::dfsRemoveTransitiveEdge(ConstraintGraph &cg, Vector2D<IntType> &edgeMat, IndexType node, 
        std::vector<bool> &visited, Vector2D<IntType> &reachable, ConstraintGraph::IndexMap &idxMap)
{
    bool hasTransitiveEdge = false;
    visited[node] = true;
    auto neighbors = boost::adjacent_vertices(boost::vertex(node, cg.boostGraph()), cg.boostGraph());
    for (; neighbors.first != neighbors.second; ++neighbors.first)
    {
        IndexType neighborNode = idxMap[*neighbors.first];
        if (!visited[neighborNode])
        {
            if (dfsRemoveTransitiveEdge(cg, edgeMat, neighborNode, visited, reachable, idxMap))
            {
                hasTransitiveEdge = true;
            }
        }
        IndexType numNodes = boost::num_vertices(cg.boostGraph());
        for (IndexType i = 0; i < numNodes; ++i)
        {
            if (reachable.at(node, i) == 0 && reachable.at(neighborNode, i) == 1)
            {
                // If has edge, remove it
                if (edgeMat.at(node, i) == 1)
                {
                    edgeMat.at(node, i) = 0;
                    boost::remove_edge(boost::vertex(node, cg.boostGraph()), 
                            boost::vertex(i, cg.boostGraph()), 
                            cg.boostGraph());
                    hasTransitiveEdge = true;
                }
                reachable.at(node, i) = 1;
            }
        }
        reachable.at(node, neighborNode) = 1;
    }
    return hasTransitiveEdge;
}

void CGLegalizer::constructConstraintGraphs()
{
    IndexType numHNodes = _db.numCells() + 2;
    IndexType numVNodes = _db.numCells() + 2;
    _hCG.allocateVertices(numHNodes);
    for (const auto &edge : _hConstraints.edges())
    {
        _hCG.addEdge(edge.source(), edge.target(), edge.weight());
    }
    _vCG.allocateVertices(numVNodes);
    for (const auto &edge : _vConstraints.edges())
    {
        _vCG.addEdge(edge.source(), edge.target(), edge.weight());
    }
}

bool horizontalOverlapping(const Box<LocType> &lhs, const Box<LocType> &rhs)
{
    if (lhs.xHi() <= rhs.xLo() || rhs.xHi() <= lhs.xLo())
    {
        return false;
    }
    return true;
}


bool verticalOverlapping(const Box<LocType> &lhs, const Box<LocType> &rhs)
{
    if (lhs.yHi() <= rhs.yLo() || rhs.yHi() <= lhs.yLo())
    {
        return false;
    }
    return true;
}

void CGLegalizer::initIrredundantEdgesInsert(bool isHor , std::vector<IndexType> &orders, 
        std::vector<IntType> &cand, IndexType cellIdx)
{
    orders.emplace_back(cellIdx);
    if (isHor)
    {
        auto compLeft = [&] (const IndexType lhs, const IndexType rhs)
        {
            return _db.cell(lhs).xLoc() < _db.cell(rhs).xLoc();
        };
        std::sort(orders.begin(), orders.end(), compLeft);
    }
    else
    {
        auto compBottom = [&] (const IndexType lhs, const IndexType rhs)
        {
            return _db.cell(lhs).yLoc() < _db.cell(rhs).yLoc();
        };
        std::sort(orders.begin(), orders.end(), compBottom);
    }
    auto it = std::find(orders.begin(), orders.end(), cellIdx);
    // First check the outgoing node
    if (it != orders.end() - 1)
    {
        // the cell has right neighbor
        const auto & cellBox1 = _db.cell(*it).cellBBoxOff();
        auto it2 = it + 1; 
        while (it2 != orders.end())
        {
            const auto &cellBox2 = _db.cell(*it2).cellBBoxOff();
            /// FIXME: I am not sure about whether the following checking should be separated?
            if (!horizontalOverlapping(cellBox1, cellBox2) || !verticalOverlapping(cellBox1, cellBox2))
            {
                break;
            }
            ++it2;
        }
        if (it2 == orders.end())
        {
            AT(cand, cellIdx) = cand.size() - 1; // Point to virtual target
        }
        else
        {
            IndexType cellIdx2 = *it2;
            AT(cand, cellIdx) = cellIdx2;
        }
    } 
    else // it != orders.end() - 1
    {
        AT(cand, cellIdx) = cand.size() - 1; // Point to virtual target
    }
    // Then check the ingoing node
    if (it != orders.begin())
    {
        // Has left neighbor
        const auto & cellBox1 = _db.cell(*it).cellBBoxOff();
        auto it2 = it - 1;
        while (it2 + 1 != orders.begin())
        {
            const auto &cellBox2 = _db.cell(*it2).cellBBoxOff();
            /// FIXME: I am not sure about whether the following checking should be separated?
            if (!horizontalOverlapping(cellBox1, cellBox2) || !verticalOverlapping(cellBox1, cellBox2))
            {
                break;
            }
            --it2;
        }
        if (it2 + 1 == orders.begin())
        {
            AT(cand, cand.size() - 1) = cellIdx; // Point virtual source to this cell
        }
        else
        {
            IndexType cellIdx2 = *it2;
            AT(cand, cellIdx2) = cellIdx; // Point it2 to it
        }
    }
    else // it != orders.begin()
    {
        AT(cand, cand.size() - 1) = cellIdx; // Point virtual source to this cell
    }
}

void CGLegalizer::initIrredundantEdgesDelete(bool isHor, std::vector<IndexType> &orders, std::vector<IntType> &cand, 
        IndexType cellIdx, std::vector<bool> &overlapAny)
{
    IntType targetIdx = static_cast<IntType>(cand.size()) - 1;
    IntType sourceIdx = static_cast<IntType>(cand.size()) - 2;
    auto it = std::find(orders.begin(), orders.end(), cellIdx);
    // First working on the right hand side
    if (it != orders.end() - 1)
    {
        // Has right neighbor
        const auto & cellBox1 = _db.cell(*it).cellBBoxOff();
        auto it2 = it + 1;
        while (it2 != orders.end())
        {
            const auto &cellBox2 = _db.cell(*it2).cellBBoxOff();
            while (!horizontalOverlapping(cellBox1, cellBox2) || !verticalOverlapping(cellBox1, cellBox2))
            {
                break;
            }
            ++it2;
        }
        if (it2 == orders.end())
        {
            // Point to the virtual target
            if (AT(cand, cellIdx) == targetIdx)
            {
                if (!isHor)
                {
                    IntType weight = - cellBox1.yLen();
                    _vConstraints.addConstraintEdge(cellIdx, targetIdx, weight);
                }
                else
                {
                    IntType weight = - cellBox1.xLen();
                    _hConstraints.addConstraintEdge(cellIdx, targetIdx, weight);
                }
            }
        }
        else // it2 == orders.end()
        {
            IndexType cellIdx2 = *it2;
            auto cellBox2 = _db.cell(*it2).cellBBoxOff();
            if (AT(cand, cellIdx) == static_cast<IntType>(cellIdx2) 
                    || AT(cand, cellIdx) == targetIdx
                    || (!horizontalOverlapping(cellBox1, cellBox2) && !verticalOverlapping(cellBox1, cellBox2))
                    )
            {
                if (isHor)
                {
                    IntType weight = - cellBox1.xLen();
                    _hConstraints.addConstraintEdge(cellIdx, cellIdx2, weight);
                }
                else
                {
                    IntType weight = - cellBox1.yLen();
                    _vConstraints.addConstraintEdge(cellIdx, cellIdx2, weight);
                }
            }
            else // AT(cand, cellIdx) ...
            {
                while (it2 != orders.end()
                        && AT(overlapAny, cellIdx2)
                        && ! (horizontalOverlapping(cellBox1, cellBox2) && verticalOverlapping(cellBox1, cellBox2)))
                {
                    if (isHor)
                    {
                        IntType weight = -cellBox1.xLen();
                        _hConstraints.addConstraintEdge(cellIdx, *it2, weight);
                    }
                    else
                    {
                        IntType weight = -cellBox2.yLen();
                        _vConstraints.addConstraintEdge(cellIdx, *it2, weight);
                    }
                    ++it2;
                    if (it2 == orders.end())
                    {
                        break;
                    }
                    cellBox2 = _db.cell(*it2).cellBBoxOff();
                }
            }
        }
    }
    else // it != orders.end() - 1
    {
        const auto & cellBox1 = _db.cell(cellIdx).cellBBoxOff();
        // Point to the virtual target
        if (AT(cand, cellIdx) == targetIdx)
        {
            if (isHor)
            {
                IntType weight = -cellBox1.xLen();
                _hConstraints.addConstraintEdge(cellIdx, targetIdx, weight);
            }
            else
            {
                IntType weight = -cellBox1.yLen();
                _vConstraints.addConstraintEdge(cellIdx, targetIdx, weight);
            }
        }
    }

    // Then check the left
    if (it != orders.begin())
    {
        const auto &cellBox1 = _db.cell(cellIdx).cellBBoxOff();
        // Has left neighbor
        auto it2 = it - 1;
        while (it2 + 1 != orders.begin())
        {
            const auto &cellBox2 = _db.cell(*it2).cellBBoxOff();
            if (!horizontalOverlapping(cellBox1, cellBox2) || !verticalOverlapping(cellBox1, cellBox2))
            {
                break;
            }
            --it2;
        }
        if (it2 + 1 == orders.begin())
        {
            // Virtual source point to it
            if (AT(cand, sourceIdx) == static_cast<IntType>(cellIdx) 
                    || AT(overlapAny, cellIdx))
            {
                // Add both h and v constraints
                // Weight = 0
                if (isHor)
                {
                    _hConstraints.addConstraintEdge(sourceIdx, cellIdx, 0);
                }
                else
                {
                    _vConstraints.addConstraintEdge(sourceIdx, cellIdx, 0);
                }
            }
        }
        else // it2 + 1 == orders.begin()
        {
            auto cellBox2 = _db.cell(*it2).cellBBoxOff();
            IndexType cellIdx2 = *it2;
            if (isHor)
            {
                IntType weight = -cellBox2.xLen();
                _hConstraints.addConstraintEdge(cellIdx2, cellIdx, weight);
            }
            else
            {
                IntType weight = -cellBox2.yLen();
                _vConstraints.addConstraintEdge(cellIdx2, cellIdx, weight);
            }
            // Ensure no overlap
            if (it2 != orders.begin())
            {
                --it2; // So that the --it2 in the while loop can be in the end
            }
            while (it2 != orders.begin())
            {
                auto cellBox2 = _db.cell(*it2).cellBBoxOff();
                if (!AT(overlapAny, *it2))
                {
                    --it2;
                    continue;
                }
                if (horizontalOverlapping(cellBox1, cellBox2) && verticalOverlapping(cellBox1, cellBox2))
                {
                    --it2;
                    continue;
                }
                if (isHor)
                {
                    IntType weight = -cellBox2.xLen();
                    _hConstraints.addConstraintEdge(cellIdx2, cellIdx, weight);
                }
                else
                {
                    IntType weight = -cellBox2.yLen();
                    _vConstraints.addConstraintEdge(cellIdx2, cellIdx, weight);
                }
                --it2;
            }
        }
    } // it != orders.begin()
    else
    {
        // Virtual source point to it
        if (AT(cand, sourceIdx) == static_cast<IntType>(cellIdx) || AT(overlapAny, cellIdx))
        {
            // Weight = 0
            if (isHor)
            {
                _hConstraints.addConstraintEdge(sourceIdx, cellIdx, 0); 
            }
            else
            {
                _vConstraints.addConstraintEdge(sourceIdx, cellIdx, 0);
            }
        }
    }
    orders.erase(it);
}


PROJECT_NAMESPACE_END
