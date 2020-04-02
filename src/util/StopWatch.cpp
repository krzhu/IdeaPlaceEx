#include "StopWatch.hpp"

namespace klib
{
    std::vector<std::uint64_t> StopWatchMgr::_us = std::vector<std::uint64_t>(1, 0);
    std::unordered_map<std::string, std::uint32_t> StopWatchMgr::_nameToIdxMap;
    StopWatch StopWatchMgr::_watch = StopWatch(0); 

    std::shared_ptr<StopWatch> StopWatchMgr::createNewStopWatch(const std::string &name) 
    {
        auto idx = _us.size();
        _us.emplace_back(std::numeric_limits<uint64_t>::max());
        _nameToIdxMap[name] = idx;
        return std::make_shared<StopWatch>(StopWatch(idx));
    }
    void StopWatchMgr::quickStart()
    {
        _watch.clear();
        _watch.start();
    }
    std::uint64_t StopWatchMgr::quickEnd()
    {
        _watch.stop();
        return _watch.record();
    }
}
