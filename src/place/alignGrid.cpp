#include "alignGrid.h"

PROJECT_NAMESPACE_BEGIN

template<typename T>
T floorDif(T n, T stepSize)
{
    throw std::invalid_argument("Unsupport argument types. Need to be integer");
}
template<>
IntType floorDif(IntType n, IntType stepSize)
{
    return n % stepSize;
}


template<typename T>
T ceilDif(T n, T stepSize)
{
    throw std::invalid_argument("Unsupport argument types. Need to be integer");
}
template<>
IntType ceilDif(IntType n, IntType stepSize)
{
    return stepSize - n % stepSize;
}

void GridAligner::align(LocType stepSize)
{
    _stepSize = stepSize;
#ifdef MULTI_SYM_GROUP
    Assert(0);
#else
    naiveAlign();
#endif
    adjustOffset(XY<LocType>(0, 0));
}

void GridAligner::naiveAlign()
{
    // Assume one symmetric axis
    LocType symAxis = 0;
    for (IndexType symGrpIdx = 0; symGrpIdx < _db.numSymGroups(); ++symGrpIdx)
    {
        const auto &symGrp = _db.symGroup(symAxis);
        for (IndexType symPairIdx = 0; symPairIdx < symGrp.numSymPairs(); ++symPairIdx)
        {
            const auto &symPair = symGrp.symPair(symPairIdx);
            symAxis = (_db.cell(symPair.firstCell()).xLo() + _db.cell(symPair.secondCell()).xHi()) / 2;
            goto theEnd;
        }
        for (IndexType ssIdx = 0; ssIdx < symGrp.numSelfSyms(); ++ssIdx)
        {
            const auto &ssCell = _db.cell(symGrp.selfSym(ssIdx));
            symAxis = (ssCell.xLo() + ssCell.xHi()) / 2;
            goto theEnd;
        }
    }
theEnd:
    for (IndexType cellIdx = 0; cellIdx < _db.numCells(); ++cellIdx)
    {
        auto &cell = _db.cell(cellIdx);
        if (cell.xLo() <= symAxis)
        {
            cell.setXLoc(cell.xLoc() - floorDif(cell.xLo(), _stepSize));
        }
        else
        {
            cell.setXLoc(cell.xLoc() - floorDif(cell.xLo(), _stepSize));
        }
        cell.setYLoc(cell.yLoc() -  floorDif(cell.yLo(), _stepSize));
    }
}

void GridAligner::adjustOffset(const XY<LocType> &offset)
{
    LocType xLo = LOC_TYPE_MAX;
    LocType yLo = LOC_TYPE_MAX;
    for (IndexType cellIdx = 0; cellIdx < _db.numCells(); ++cellIdx)
    {
        const auto & cell = _db.cell(cellIdx);
        xLo = std::min(xLo, cell.xLo());
        yLo = std::min(yLo, cell.yLo());
    }
    // Compute the dif
    xLo = offset.x() - xLo;
    yLo = offset.y() - yLo;
    for (IndexType cellIdx = 0; cellIdx < _db.numCells(); ++cellIdx)
    {
        auto & cell = _db.cell(cellIdx);
        cell.setXLoc(cell.xLoc() + xLo);
        cell.setYLoc(cell.yLoc() + yLo);
    }
}

PROJECT_NAMESPACE_END
