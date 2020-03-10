#include "VirtualPinAssigner.h"
#include <lemon/list_graph.h>
#include <lemon/network_simplex.h>

PROJECT_NAMESPACE_BEGIN

bool VirtualPinAssigner::solveFromDB()
{
    reconfigureVirtualPinLocationFromDB();
    return pinAssignmentFromDB();
}

void VirtualPinAssigner::reconfigureVirtualPinLocationFromDB()
{
    Box<LocType> boundary(LOC_TYPE_MAX, LOC_TYPE_MAX, LOC_TYPE_MIN, LOC_TYPE_MIN);
    for (IndexType cellIdx = 0; cellIdx < _db.numCells(); ++cellIdx)
    {
        const auto &cell = _db.cell(cellIdx);
        boundary.unionBox(cell.cellBBoxOff());
    }
    reconfigureVirtualPinLocations(boundary);
}


bool VirtualPinAssigner::pinAssignmentFromDB()
{
    auto cellLocQueryFunc = [&] (IndexType cellIdx)
    {
        const auto &cell = _db.cell(cellIdx);
        return cell.loc();
    };
    return pinAssignment(cellLocQueryFunc);
}

IntType gcd(IntType a, IntType b)
{
    for (;;)
    {
        if (a == 0) return b;
        b %= a;
        if (b == 0) return a;
        a %= b;
    }
}

IntType lcm(IntType a, IntType b)
{
    IntType temp = gcd(a, b);
    return temp ? (a / temp * b) : 0;
}

void VirtualPinAssigner::reconfigureVirtualPinLocations(const Box<LocType> &cellsBBox)
{
    _virtualPinInterval = _db.parameters().virtualPinInterval();
    _virtualBoundaryExtension = _db.parameters().virtualBoundaryExtension();
    _boundary = cellsBBox;
    _boundary.enlargeBy(_virtualBoundaryExtension);
    LocType pinInterval = _virtualPinInterval;
    // Align to grid
    if (_db.parameters().hasGridStep())
    {
        LocType gridStep = _db.parameters().gridStep();
        LocType center = _boundary.center().x();
        LocType targetCenter = (center / gridStep) * gridStep + gridStep / 2;
        LocType targetWidth = std::max(_boundary.xHi() - targetCenter, targetCenter - _boundary.xLo());
        targetWidth = targetWidth + gridStep - (targetWidth % gridStep);
        _boundary.setXLo(targetCenter - targetWidth);
        _boundary.setYLo(_boundary.yLo() - (_boundary.yLo() % gridStep));
        _boundary.setXHi(targetCenter + targetWidth);
        _boundary.setYHi(_boundary.yHi() + gridStep -  (_boundary.yHi() % gridStep));
        pinInterval = lcm(pinInterval, gridStep);
    }
    // generate the virtual pin locations
    _virtualPins.clear();
    for (LocType x = _boundary.xLo() + pinInterval;  x < _boundary.center().x() - pinInterval / 2 ; x += pinInterval)
    {
        continue;
        LocType rightX = 2 * _boundary.center().x() - x;
        if (rightX <= x)
        {
            ERR("Pin assignment: unexpected pin location x %d center %d right x %d pinInterval %d \n", x, _boundary.center().x(), rightX, pinInterval);
            continue;
        }
        // left bottom
        _virtualPins.emplace_back(XY<LocType>(x, _boundary.yLo()));
        _virtualPins.back().setDirection(Direction2DType::SOUTH);

        // right bottom
        _virtualPins.emplace_back(XY<LocType>(rightX, _boundary.yLo()));
        _virtualPins.back().setDirection(Direction2DType::SOUTH);

        // Add to map
        _leftToRightMap[_virtualPins.size() - 2] = _virtualPins.size() - 1;

        // left top
        _virtualPins.emplace_back(XY<LocType>(x, _boundary.yHi()));
        _virtualPins.back().setDirection(Direction2DType::NORTH);

        // right top
        _virtualPins.emplace_back(XY<LocType>(rightX, _boundary.yHi()));
        _virtualPins.back().setDirection(Direction2DType::NORTH);

        // Add to map
        _leftToRightMap[_virtualPins.size() - 2] = _virtualPins.size() - 1;
    }
    for (LocType y = _boundary.yLo() + pinInterval;  y < _boundary.yHi() - pinInterval; y += pinInterval)
    {
        _virtualPins.emplace_back(XY<LocType>(_boundary.xLo(), y));
        _virtualPins.back().setDirection(Direction2DType::WEST);
        _virtualPins.emplace_back(XY<LocType>(_boundary.xHi(), y));
        _virtualPins.back().setDirection(Direction2DType::EAST);

        // Add to map
        _leftToRightMap[_virtualPins.size() - 2] = _virtualPins.size() - 1;
    }
    std::set<VirtualPin> pinSet;
    for (const auto & vp : _virtualPins)
    {
        auto findIter = pinSet.find(vp);
        Assert(findIter == pinSet.end());
        pinSet.insert(vp);
    }
}

bool VirtualPinAssigner::_networkSimplexPinAssignment(std::function<bool(IndexType)> useNetFunc,
        std::function<bool(IndexType)> usePinFunc,
        std::function<LocType(IndexType, IndexType)> netToPinCostFunc,
        std::function<void(IndexType, IndexType)> setNetToVirtualPinFunc)
{
    // Prepare the nets and pins available to the problem
    std::vector<IndexType> nets;
    for (IndexType netIdx = 0; netIdx < _db.numNets(); ++netIdx)
    {
        if (useNetFunc(netIdx))
        {
            nets.emplace_back(netIdx);
        }
    }
    std::vector<IndexType> iopins; 
    for (IndexType iopinIdx = 0; iopinIdx < _virtualPins.size(); ++iopinIdx)
    {
        if (usePinFunc(iopinIdx))
        {
            iopins.emplace_back(iopinIdx);
        }
    }
    // If # of io pins < # nets, the problem is infeasible
    IndexType numNets = nets.size();
    IndexType numSites = iopins.size();
    if (numNets > numSites)
    {
        return false;
    }

    // Build the network
    // Use lemon min cost max flow with network simplex algorithm
    lemon::ListDigraph graph;
    std::vector<lemon::ListDigraph::Node> lNodes, rNodes;
    std::vector<lemon::ListDigraph::Arc> lArcs, rArcs, mArcs;
    std::vector<std::pair<IndexType, IndexType>> mArcPairs;
    // Reserve the needed space (optional)
    lNodes.reserve(numNets); rNodes.reserve(numSites);
    lArcs.reserve(numNets); rArcs.reserve(numSites); mArcs.reserve(numNets * numSites);
    mArcPairs.reserve(numNets * numSites);

    lemon::ListDigraph::ArcMap<IntType> capLo(graph); // capacity low
    lemon::ListDigraph::ArcMap<IntType> capHi(graph); // capacity high
    lemon::ListDigraph::ArcMap<IntType> costMap(graph); // Cost map

    // Adding the nodes and arcs
    // Source and target nodes
    auto s = graph.addNode();
    auto t = graph.addNode();

    // Add arcs source->left
    for (IndexType i = 0; i < numNets; ++i)
    {
        lNodes.emplace_back(graph.addNode());
        lArcs.emplace_back(graph.addArc(s, lNodes.back()));
        capLo[lArcs.back()] = 0;
        capHi[lArcs.back()] = 1;
        costMap[lArcs.back()] = 0;
    }
    // Add arcs right->target
    for (IndexType i = 0; i < numSites; ++i)
    {
        rNodes.emplace_back(graph.addNode());
        rArcs.emplace_back(graph.addArc(rNodes.back(), t));
        capLo[rArcs.back()] = 0;
        capHi[rArcs.back()] = 1;
        costMap[rArcs.back()] = 0;
    }
    // Mid arcs
    for (IndexType l = 0; l < numNets; ++l)
    {
        for (IndexType r = 0; r < numSites; ++r)
        {
            mArcs.emplace_back(graph.addArc(lNodes[l], rNodes[r]));
            mArcPairs.emplace_back(l, r);
            capLo[mArcs.back()] = 0;
            capHi[mArcs.back()] = 1;
            costMap[mArcs.back()] = netToPinCostFunc(nets[l], iopins[r]);
        }
    }

    // Solve min cost max flow sing network simplex algorithm
    lemon::NetworkSimplex<lemon::ListDigraph, IntType> networkSimplex(graph);
    networkSimplex.stSupply(s, t, numNets);
    networkSimplex.lowerMap(capLo).upperMap(capHi).costMap(costMap);
    networkSimplex.run();
    
    // Collect the solution and export into the database
    for (IndexType i = 0; i < mArcs.size(); ++i)
    {
        if (networkSimplex.flow(mArcs[i]))
        {
            const auto &pair = mArcPairs[i];
            setNetToVirtualPinFunc(nets[pair.first], iopins[pair.second]);
        }
    }
    return true;
}

bool VirtualPinAssigner::pinAssignment(std::function<XY<LocType>(IndexType)> cellLocQueryFunc)
{

    // Calculate the current HPWLs without virtual pin
    std::vector<Box<LocType>> curNetBBox;
    curNetBBox.resize(_db.numNets());
    for (IndexType netIdx = 0; netIdx < _db.numNets(); ++netIdx)
    {
        IndexType pinIdx = _db.net(netIdx).pinIdx(0);
        XY<LocType> pinOff = _db.pin(pinIdx).midLoc();
        XY<LocType> cellLoc = cellLocQueryFunc(_db.pin(pinIdx).cellIdx());
        curNetBBox.at(netIdx) = Box<LocType>(cellLoc + pinOff, cellLoc + pinOff); // Init into an area=0 box
        // The rest pins
        for (IndexType idx = 1; idx < _db.net(netIdx).numPinIdx(); ++idx)
        {
            pinIdx = _db.net(netIdx).pinIdx(idx);
            pinOff = _db.pin(pinIdx).midLoc();
            cellLoc = cellLocQueryFunc(_db.pin(pinIdx).cellIdx());
            curNetBBox.at(netIdx).join(cellLoc + pinOff);
        }
    }

    auto calculateIncreasedHpwl = [&](IndexType netIdx, IndexType pinIdx)
    {
        const auto &virtualPinLoc = _virtualPins.at(pinIdx);
        auto difX = std::max(virtualPinLoc.x() - curNetBBox.at(netIdx).xHi(), curNetBBox.at(netIdx).xLo() - virtualPinLoc.x());
        difX = std::max(difX, 0);
        auto difY = std::max(virtualPinLoc.y() - curNetBBox.at(netIdx).yHi(), curNetBBox.at(netIdx).yLo() - virtualPinLoc.y());
        difY = std::max(difY, 0);
        return difX + difY;
    };

    // Calculate the added HPWL if adding the virtual pin
    auto directNetToPinCostFunc = [&](IndexType netIdx, IndexType virtualPinIdx)
    {
        return calculateIncreasedHpwl(netIdx, virtualPinIdx);
    };

    auto useASymNet = [&](IndexType netIdx)
    {
        if (_db.net(netIdx).isSelfSym())
        {
            return false;
        }
        if (_db.net(netIdx).hasSymNet())
        {
            return false;
        }
        return true;
    };

    auto useFreePin = [&](IndexType pinIdx)
    {
        if (_virtualPins.at(pinIdx).assigned())
        {
            return false;
        }
        return true;
    };

    auto directAssignNetToPinFunc = [&](IndexType netIdx, IndexType virtualPinIdx)
    {
        _virtualPins[virtualPinIdx].assign(netIdx);
        _db.net(netIdx).setVirtualPin(_virtualPins[virtualPinIdx]);
    };

    // THe symmetric parts
    auto symPairNetToPinCostFunc = [&](IndexType netIdx, IndexType leftPinIdx)
    {
        auto rightPinIdx = _leftToRightMap.at(leftPinIdx);
        if (_db.net(netIdx).isSelfSym())
        {
            return calculateIncreasedHpwl(netIdx, leftPinIdx) + calculateIncreasedHpwl(netIdx, rightPinIdx);
        }
        else if (_db.net(netIdx).hasSymNet())
        {
            auto otherNetIdx = _db.net(netIdx).symNetIdx();
            auto netCost0 = calculateIncreasedHpwl(netIdx, leftPinIdx) + calculateIncreasedHpwl(otherNetIdx, rightPinIdx);
            auto netCost1 = calculateIncreasedHpwl(otherNetIdx, leftPinIdx) + calculateIncreasedHpwl(netIdx, rightPinIdx);
            return std::min(netCost0, netCost1);
        }
        Assert(false);
        return 0;
    };

    auto useSymNet = [&](IndexType netIdx)
    {
        if (_db.net(netIdx).isSelfSym())
        {
            return true;
        }
        if (_db.net(netIdx).hasSymNet())
        {
            if (netIdx < _db.net(netIdx).symNetIdx())
            {
                return true;
            }
        }
        return false;
    };
    
    auto useLeftPin = [&](IndexType pinIdx)
    {
        if (_virtualPins.at(pinIdx).assigned())
        {
            return false;
        }
        if (_leftToRightMap.find(pinIdx) ==  _leftToRightMap.end())
        {
            return false;
        }
        return true;
    };

    auto symPairAssignNetToPinFunc = [&](IndexType netIdx, IndexType leftPinIdx)
    {
        auto rightPinIdx = _leftToRightMap.at(leftPinIdx);
        if (_db.net(netIdx).isSelfSym())
        {
            directAssignNetToPinFunc(netIdx, leftPinIdx);
            // TODO: also add right pin
        }
        else if (_db.net(netIdx).hasSymNet())
        {
            auto otherNetIdx = _db.net(netIdx).symNetIdx();
            auto netCost0 = calculateIncreasedHpwl(netIdx, leftPinIdx) + calculateIncreasedHpwl(otherNetIdx, rightPinIdx);
            auto netCost1 = calculateIncreasedHpwl(otherNetIdx, leftPinIdx) + calculateIncreasedHpwl(netIdx, rightPinIdx);
            if (netCost0 <= netCost1)
            {
                // net -> left. other net -> right
                directAssignNetToPinFunc(netIdx, leftPinIdx);
                directAssignNetToPinFunc(otherNetIdx, rightPinIdx);
            }
            else
            {
                // net -> right. other net -> left
                directAssignNetToPinFunc(netIdx, rightPinIdx);
                directAssignNetToPinFunc(otherNetIdx, leftPinIdx);
            }
        }
        else
        {
            Assert(false);
        }

    };

    if (!_networkSimplexPinAssignment(useSymNet, useLeftPin, symPairNetToPinCostFunc, symPairAssignNetToPinFunc))
    {
        return false;
    }

    return _networkSimplexPinAssignment(useASymNet, useFreePin, directNetToPinCostFunc, directAssignNetToPinFunc);
}

PROJECT_NAMESPACE_END
