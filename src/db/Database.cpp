#include "Database.h"

PROJECT_NAMESPACE_BEGIN

bool Database::initCells()
{
    IndexType numLayers = this->tech().numLayers();
    for (auto &cell : _cellArray)
    {
        cell.allocateLayers(numLayers);
    }
    return true;
}

PROJECT_NAMESPACE_END
