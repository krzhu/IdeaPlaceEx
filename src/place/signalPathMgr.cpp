#include "signalPathMgr.h"

PROJECT_NAMESPACE_BEGIN

SigPathMgr::SigPathMgr(Database &db)
    : _db(db)
{
    this->decomposeSignalPaths();
}

void SigPathMgr::decomposeSignalPaths()
{
}

PROJECT_NAMESPACE_END
