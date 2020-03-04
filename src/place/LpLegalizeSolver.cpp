#include "CGLegalizer.h"

PROJECT_NAMESPACE_BEGIN

RealType LpLegalizeSolver::evaluateObj()
{
    return _ilpModel.evaluateObjective();
}

void LpLegalizeSolver::exportSolution()
{
    for (IndexType cellIdx = 0; cellIdx < _db.numCells(); ++cellIdx)
    {
        auto var = _ilpModel.variableSolution(_locs.at(cellIdx));
        // convert to cell original location
        if (_isHor)
        {
            _db.cell(cellIdx).setXLo(static_cast<LocType>(var) + _db.parameters().layoutOffset());
        }
        else
        {
            _db.cell(cellIdx).setYLo(static_cast<LocType>(var) + _db.parameters().layoutOffset());
        }
    }
}

bool LpLegalizeSolver::solve()
{
    // Add variables
    addIlpVars();
    // add constraints
    addIlpConstraints();
    // Configure the objective function
    configureObjFunc();
    // Solve the LP problem
    return solveLp();
}

bool LpLegalizeSolver::solveLp()
{
    _ilpModel.setObjective(_obj);
    _ilpModel.setOptimizeType(limbo::solvers::MIN);
    _params.setVerbose(2); // 2: SEVERE
    SolverType _solver(&_ilpModel);
    _optimStatus = _solver(&_params);
    if (_optimStatus == limbo::solvers::UNBOUNDED)
    {
        ERR("LP legalization solver: LP unbounded \n");
        return false;
    }
    else if (_optimStatus == limbo::solvers::OPTIMAL)
    {
        INF("LP legalization solver: LP optimal \n");
        return true;
    }
    else if (_optimStatus == limbo::solvers::INFEASIBLE)
    {
        ERR("LP legalization solver: LP infeasible \n");
        return false;
    }
    else
    {
        ERR("LP legalization solver: Unknown LP status %d \n", _optimStatus);
        return false;
    }
}

void LpLegalizeSolver::addWirelengthObj()
{
    if (_optHpwl == 1)
    {
        bool hasAtLeastOneNet = false;
        for (IndexType netIdx = 0; netIdx < _db.numNets(); ++netIdx)
        {
            if (_db.net(netIdx).numPinIdx() == 0)
            {
                continue;
            }
            if (_db.net(netIdx).numPinIdx() == 1 && !_db.net(netIdx).isValidVirtualPin())
            {
                continue;
            }
            hasAtLeastOneNet = true;
            auto weight = _db.net(netIdx).weight();
            _obj += weight * (_wlR[netIdx] - _wlL[netIdx]);
        }
        if (!hasAtLeastOneNet)
        {
            ERR("LP Legalizer:: No valid net \n");
            Assert(false);
        }
    }
}

void LpLegalizeSolver::addAreaObj()
{
    if (_optArea == 1)
    {
        _obj += _dim;
    }
}

void LpLegalizeSolver::configureObjFunc()
{
    // Wirelength
    this->addWirelengthObj();
    // Area
    this->addAreaObj();
}

IndexType LpLegalizeSolver::numVars() const
{
    auto numLocVars = _db.numCells();
    auto numHpwlVars = _db.numNets() * 2 * _optHpwl;
    auto numBoundaryVars = 1 * _optArea;
    IndexType numSymVars;
    if (_isMultipleSymGrp)
    {
        numSymVars = _db.numSymGroups();
    }
    else
    {
        numSymVars = 1;
    }
    return numLocVars + numHpwlVars + numBoundaryVars + numSymVars;
}

void LpLegalizeSolver::addLocVars()
{
    // NOTE: the _locs variables here are general location variables
    _locs.resize(_db.numCells());
    for (IndexType i = 0; i < _db.numCells(); ++i)
    {
        _locs.at(i) = _ilpModel.addVariable(0, std::numeric_limits<RealType>::max(),
                                                limbo::solvers::CONTINUOUS, 
                                                "loc_" + std::to_string(i));
    }
}

void LpLegalizeSolver::addWirelengthVars()
{
    // add HPWL variables
    if (_optHpwl == 1)
    {
        _wlL.resize(_db.numNets());
        _wlR.resize(_db.numNets());
        for (IndexType i =0; i < _db.numNets(); ++i)
        {
            _wlL.at(i) = _ilpModel.addVariable(0, std::numeric_limits<RealType>::max(),
                                                limbo::solvers::CONTINUOUS,
                                                "wll_" + std::to_string(i));
            _wlR.at(i) = _ilpModel.addVariable(0, std::numeric_limits<RealType>::max(),
                                                limbo::solvers::CONTINUOUS,
                                                "wlr_" + std::to_string(i));
        }
    }
}

void LpLegalizeSolver::addAreaVars()
{
    if (_optArea == 1)
    {
        _dim = _ilpModel.addVariable(0, std::numeric_limits<RealType>::max(),
                                    limbo::solvers::CONTINUOUS,
                                    "dim");
    }
}

void LpLegalizeSolver::addSymVars()
{
    if (_isMultipleSymGrp)
    {
        // Symmetric group axis variables
        _symLocs.resize(_db.numSymGroups());
        for (IndexType i = 0; i < _db.numSymGroups(); ++i)
        {
            _symLocs.at(i) = _ilpModel.addVariable(0,
                    std::numeric_limits<RealType>::max(),
                    limbo::solvers::CONTINUOUS,
                    "symLoc_"+std::to_string(i));
        }
    }
    else
    {
        _symLocs.resize(1);
        _symLocs[0] = _ilpModel.addVariable(0,
                    std::numeric_limits<RealType>::max(),
                    limbo::solvers::CONTINUOUS,
                    "symLoc_"+std::to_string(0));
    }
}

void LpLegalizeSolver::addIlpVars()
{
    // Calculate how many variables are there in the ILP model
    auto numVars = this->numVars();
    // Reserve the variables in the model
    _ilpModel.reserveVariables(numVars);

    // add coordinate variables
    this->addLocVars();
    // Add wire length variables
    this->addWirelengthVars();
    // Add area variables
    this->addAreaVars();
    // Add symmetric variables
    this->addSymVars();

}

void LpLegalizeSolver::addBoundaryConstraints()
{
    for (IndexType i = 0;  i < _db.numCells(); ++i)
    {
        if (_optArea == 0)
        {
            if (_isHor)
            {
                // 0 <= x_i <= W* - w_i
#ifdef DEBUG_LEGALIZE
                DBG("Add boundary constraint: x_%d <= %f - %d \n", i, _wStar, _db.cell(i).cellBBox().xLen());
#endif
                _ilpModel.addConstraint(_locs.at(i) <= _wStar - _db.cell(i).cellBBox().xLen());
            }
            else
            {
#ifdef DEBUG_LEGALIZE
                DBG("Add boundary constraint: y_%d <= %f - %d \n", i, _wStar, _db.cell(i).cellBBox().yLen());
#endif
                _ilpModel.addConstraint(_locs.at(i) <= _wStar - _db.cell(i).cellBBox().yLen());
            }
        }
        else // if (_optArea == 0)
        {
            if (_isHor)
            {
                // 0 <= x_i <= W - w_i
                _ilpModel.addConstraint(_locs.at(i) - _dim <= - _db.cell(i).cellBBox().xLen());
            }
            else
            {
                _ilpModel.addConstraint(_locs.at(i) - _dim <= - _db.cell(i).cellBBox().yLen());
            }
        }
    }
}

void LpLegalizeSolver::addTopologyConstraints()
{
    for (auto & edge : _constrains.edges())
    {
        IndexType sourceIdx = edge.source();
        IndexType targetIdx = edge.target();
        if (sourceIdx == targetIdx)
        {
            continue;
        }
        if (sourceIdx == _db.numCells() || targetIdx == _db.numCells() + 1)
        {
            // the s, t constraints
            continue;
        }
        LocType cellDim; // width or height
        auto spacingBox = _db.cellSpacing(sourceIdx, targetIdx);
        LocType spacing;
        // Force cell1 to be lower/left to the cell2. Therefore only using spacingBox.xLo() and .yLo()
        if (_isHor)
        {
            cellDim = _db.cell(sourceIdx).cellBBox().xLen();
            spacing = spacingBox.xLo();
        }
        else
        {
            cellDim = _db.cell(sourceIdx).cellBBox().yLen();
            spacing = spacingBox.yLo();
        }
        // Add the constraint 
        // x_i + w_i + spacing <= x_j
        _ilpModel.addConstraint(_locs.at(sourceIdx) - _locs.at(targetIdx) <= - cellDim - spacing);
#ifdef DEBUG_LEGALIZE
        DBG("Add spacing constrain: from %d to %d, <= -celldim %d - spacing %d = %d \n", sourceIdx, targetIdx, cellDim, spacing, -cellDim - spacing);
#endif
    }

}

void LpLegalizeSolver::addSymmetryConstraints()
{
    if (_isHor)
    {
        // Force them to be symmetric along an axis
        for (IndexType symGrpIdx = 0; symGrpIdx < _db.numSymGroups(); ++symGrpIdx)
        {
            const auto &symGrp = _db.symGroup(symGrpIdx);
            LpModelType::variable_type *symVar;
            if (_isMultipleSymGrp)
            {
                symVar = &_symLocs[symGrpIdx];
            }
            else
            {
                symVar = &_symLocs[0];
            }

            for (IndexType symPairIdx = 0; symPairIdx < symGrp.numSymPairs(); ++symPairIdx)
            {
                const auto &symPair = symGrp.symPair(symPairIdx);
                // x1 + x2 + width =  2 * symAxis
#ifdef DEBUG_LEGALIZE
                DBG("Add sym constraint. \n symGrp %d, cell %d %d \n width %d \n",
                        symGrpIdx, symPair.firstCell(), symPair.secondCell(), _db.cell(symPair.firstCell()).cellBBox().xLen());
#endif
                _ilpModel.addConstraint(_locs.at(symPair.firstCell()) 
                        + _locs.at(symPair.secondCell()) 
                        - 2 * (*symVar) 
                        == 
                        -_db.cell(symPair.firstCell()).cellBBox().xLen()); // Two cells are equal in width <- assumption
                AssertMsg(_db.cell(symPair.firstCell()).cellBBox().xLen() == _db.cell(symPair.secondCell()).cellBBox().xLen(), "cell %s and cell %s \n", _db.cell(symPair.firstCell()).name().c_str(),  _db.cell(symPair.secondCell()).name().c_str());
            }
            for (IndexType selfSymIdx = 0; selfSymIdx < symGrp.numSelfSyms(); ++selfSymIdx)
            {
                IndexType ssCellIdx = symGrp.selfSym(selfSymIdx);
                // x1 + width + x2 = 2 * symAxis
                _ilpModel.addConstraint(2 * _locs.at(ssCellIdx) - 2 * (*symVar)
                        == - _db.cell(ssCellIdx).cellBBox().xLen());
            }
        }
    }
    else
    {
        // Force they have the same y coordinate
        for (IndexType symGroupIdx = 0;  symGroupIdx < _db.numSymGroups(); ++symGroupIdx)
        {
            const auto & symGroup = _db.symGroup(symGroupIdx);
            for (IndexType symPairIdx = 0; symPairIdx < symGroup.numSymPairs(); ++symPairIdx)
            {
                const auto &symPair = symGroup.symPair(symPairIdx);
                // y_i = y_j
                _ilpModel.addConstraint(_locs.at(symPair.firstCell()) - _locs.at(symPair.secondCell()) == 0.0);
            }
        }
    }
}

void LpLegalizeSolver::addHpwlConstraints()
{
    if (_optHpwl)
    {
        for (IndexType netIdx = 0; netIdx < _db.numNets(); ++netIdx)
        {
            const auto &net = _db.net(netIdx);
            for (IndexType pinIdxInNet =0; pinIdxInNet < net.numPinIdx(); ++pinIdxInNet)
            {
                IndexType pinIdx = net.pinIdx(pinIdxInNet);
                const auto &pin = _db.pin(pinIdx);
                const auto &cell = _db.cell(pin.cellIdx());
                auto midLoc = pin.midLoc();
                XY<RealType> cellLoLoc = XY<RealType>(cell.cellBBox().xLo(), cell.cellBBox().yLo());
                midLoc -= cellLoLoc;
                if (_isHor)
                {
                    RealType loc = static_cast<RealType>(midLoc.x());
                    // wl_l <= _loc + pin_offset for all pins in the net
                    _ilpModel.addConstraint(  _wlL.at(netIdx)
                            - _locs.at(pin.cellIdx())
                            <=  loc);
                    // wl_r >= _loc + pin_offset for all pins in the net
                    _ilpModel.addConstraint(  _wlR.at(netIdx)
                            - _locs.at(pin.cellIdx()) 
                            >=  loc);
                }
                else
                {
                    RealType loc = static_cast<RealType>(midLoc.y());
                    // wl_l <= _loc + pin_offset for all pins in the net
                    _ilpModel.addConstraint( _wlL.at(netIdx)
                            - _locs.at(pin.cellIdx())
                            <=  loc);
                    // wl_r >= _loc + pin_offset for all pins in the net
                    _ilpModel.addConstraint(  _wlR.at(netIdx)
                            - _locs.at(pin.cellIdx()) 
                            >=  loc);
                }
            }
            // Wirelength with virtual pin
            if (net.isValidVirtualPin())
            {
                if (_isHor)
                {
                    RealType loc = static_cast<RealType>(net.virtualPinLoc().x());
                    // wl_l <= _loc + pin_offset for all pins in the net
                    _ilpModel.addConstraint(  _wlL.at(netIdx)
                            <=  std::max(loc, 0.0));
                    // wl_r >= _loc + pin_offset for all pins in the net
                    _ilpModel.addConstraint(  _wlR.at(netIdx)
                            >=  loc);
                }
                else
                {
                    RealType loc = static_cast<RealType>(net.virtualPinLoc().y());
                    // wl_l <= _loc + pin_offset for all pins in the net
                    _ilpModel.addConstraint( _wlL.at(netIdx)
                            <=  std::max(loc, 0.0));
                    // wl_r >= _loc + pin_offset for all pins in the net
                    _ilpModel.addConstraint( _wlR.at(netIdx)
                            >=  loc);
                }
            }
        }
    }
}

void LpLegalizeSolver::addIlpConstraints()
{
    // Add boundary constraint
    addBoundaryConstraints();
    // Add topology constraints
    addTopologyConstraints();
    // Add symmetric constraints
    addSymmetryConstraints();
    // Add HPWL constraints
    addHpwlConstraints();
}

PROJECT_NAMESPACE_END
