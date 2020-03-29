#include "signalPathMgr.h"

PROJECT_NAMESPACE_BEGIN

SigPathMgr::SigPathMgr(Database &db)
    : _db(db)
{
    this->decomposeSignalPaths();
}

void SigPathMgr::decomposeSignalPaths()
{
    for (const auto &path : _db.vSignalPaths())
    {
        for (IndexType idx = 0; idx < path.vPinIdxArray().size() - 3; ++idx)
        {
            const IndexType pinIdx1 = path.vPinIdxArray().at(idx);   const auto &pin1 = _db.pin(pinIdx1); const IndexType cellIdx1 = pin1.cellIdx();
            const IndexType pinIdx2 = path.vPinIdxArray().at(idx+1); const auto &pin2 = _db.pin(pinIdx2); const IndexType cellIdx2 = pin2.cellIdx();
            const IndexType pinIdx3 = path.vPinIdxArray().at(idx+2); const auto &pin3 = _db.pin(pinIdx3); const IndexType cellIdx3 = pin3.cellIdx();
            const IndexType pinIdx4 = path.vPinIdxArray().at(idx+3); const auto &pin4 = _db.pin(pinIdx4); const IndexType cellIdx4 = pin4.cellIdx();

            bool pin23Same = cellIdx2 == cellIdx3;
            bool pin12Diff = cellIdx1 != cellIdx2;
            bool pin34Diff = cellIdx3 != cellIdx4;
            bool pin14Diff = cellIdx1 != cellIdx4;
            if (pin23Same and pin12Diff and pin34Diff and pin14Diff)
            {
                // valid segment
                _segs.emplace_back(SigPathSeg(pinIdx1, pinIdx2, pinIdx3, pinIdx4));
                DBG("SigPathMgr: add cell %d pin %d -> cell %d pin %d -> cell %d pin %d -> cell %d -< pin %d \n", cellIdx1, pinIdx1, cellIdx2, pinIdx2, cellIdx3, pinIdx3, cellIdx4, pinIdx4);
            }
        }
    }
}

PROJECT_NAMESPACE_END
