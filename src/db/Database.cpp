#include "Database.h"

PROJECT_NAMESPACE_BEGIN

bool Database::initCells()
{
    for (IndexType cellIdx = 0; cellIdx < this->numCells(); ++cellIdx)
    {
        if (!this->initCell(cellIdx))
        {
            return false;
        }
    }
    return true;
}


bool Database::initCell(IndexType cellIdx)
{
    IndexType numLayers = this->tech().numLayers();
    AT(_cellArray, cellIdx).allocateLayers(numLayers);
    return true;
}

LocType Database::hpwl() const
{
    LocType hpwl = 0;
    for (const auto & net : _netArray)
    {
        if (net.numPinIdx() <= 1)
        {
            continue;
        }
        LocType xMax = LOC_TYPE_MIN;
        LocType xMin = LOC_TYPE_MAX;
        LocType yMax = LOC_TYPE_MIN;
        LocType yMin = LOC_TYPE_MAX;
        for (IndexType pinIdx : net.pinIdxArray())
        {
            const auto &pin = AT(_pinArray, pinIdx);
            auto pinLoc = pin.midLoc();
            IndexType cellIdx = pin.cellIdx();
            const auto &cell = AT(_cellArray, cellIdx);
            auto cellLoc = cell.loc();
            auto pinFinalLoc = pinLoc + cellLoc;
            xMax = std::max(xMax, pinFinalLoc.x());
            xMin = std::min(xMin, pinFinalLoc.x());
            yMax = std::max(yMax, pinFinalLoc.y());
            yMin = std::min(yMin, pinFinalLoc.y());
        }
        hpwl = (xMax - xMin) + (yMax - yMin);
    }
    return hpwl;
}

/*------------------------------*/ 
/* Debug functions              */
/*------------------------------*/ 
#ifdef DEBUG_DRAW
PROJECT_NAMESPACE_END
#include "writer/gdsii/WriteGds.h"
PROJECT_NAMESPACE_BEGIN
void Database::drawCellBlocks(const std::string &filename)
{
    auto wg = std::make_shared<WriteGds>(filename);
    if (!wg->initWriter())
    {
        return;
    }
    if (!wg->createLib("TOP", 2000, 1e-6/2000)) // Hardcoded numbers
    {
        return;
    }
    if (!wg->writeCellBgn("DEBUG"))
    {
        return;
    }
    // Write all the cells
    for (IndexType cellIdx = 0; cellIdx < this->numCells(); ++cellIdx)
    {
        const auto &cell = this->cell(cellIdx);
        Box<LocType> cellBox = cell.cellBBox();
        cellBox.enlargeBy(0);
        XY<LocType> cellLoc = XY<LocType>(cell.xLoc(), cell.yLoc());
        cellBox.offsetBy(cellLoc);
        wg->writeRectangle(cellBox, cellIdx, 0);
        // Also wrute pins
        for (IndexType pinIdxInCell = 0; pinIdxInCell < cell.numPinIdx(); ++pinIdxInCell)
        {
            const auto &pin = this->pin(cell.pinIdx(pinIdxInCell));
            Box<LocType> pinBox = pin.shape();
            pinBox.offsetBy(cellLoc);
            wg->writeRectangle(pinBox, 100 + cellIdx, 0);
        }
    }
    // END
    wg->writeCellEnd();
    wg->endLib();
    DBG("Database::%s: debug cell block saved in %s \n", __FUNCTION__, filename.c_str());
}

#endif //DEBUG_DRAW
PROJECT_NAMESPACE_END
