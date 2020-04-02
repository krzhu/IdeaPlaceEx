#include "NlpGPlacer.h"
#include "place/signalPathMgr.h"


PROJECT_NAMESPACE_BEGIN

using namespace nt;

template<typename nlp_settings>
IntType NlpGPlacerBase<nlp_settings>::solve()
{
    this->initProblem();
    this->initPlace();
    this->initOperators();
    this->constructTasks();
    this->optimize();
    return 0;
}

template<typename nlp_settings>
void NlpGPlacerBase<nlp_settings>::optimize()
{
    _wrapObjAllTask.run();
    DBG("obj: %f %f %f %f %f %f \n", _obj, _objHpwl, _objOvl, _objOob, _objAsym, _objCos);
}

template<typename nlp_settings>
void NlpGPlacerBase<nlp_settings>::initOptimizationKernelMembers()
{
    _stopCondition = stop_condition_trait::construct(*this);
}

template<typename nlp_settings>
void NlpGPlacerBase<nlp_settings>::initProblem()
{
    initHyperParams();
    initBoundaryParams();
    initVariables();
}

template<typename nlp_settings>
void NlpGPlacerBase<nlp_settings>::initHyperParams()
{
    _alpha = NLP_WN_CONJ_ALPHA;
}

template<typename nlp_settings>
void NlpGPlacerBase<nlp_settings>::initBoundaryParams()
{
    auto maxWhiteSpace = _db.parameters().maxWhiteSpace();
    // Total cell area
    RealType totalCellArea = static_cast<RealType>(_db.calculateTotalCellArea());
    _scale = sqrt(100 / (totalCellArea));
    _totalCellArea = 100;

    // Placement Boundary
    if (_db.parameters().isBoundaryConstraintSet())
    {
        // If the constraint is set in the database, follow it.
        const auto &bb = _db.parameters().boundaryConstraint();
        _boundary.setXLo(static_cast<RealType>(bb.xLo()) * _scale);
        _boundary.setYLo(static_cast<RealType>(bb.yLo()) * _scale);
        _boundary.setXHi(static_cast<RealType>(bb.xHi()) * _scale);
        _boundary.setYHi(static_cast<RealType>(bb.yHi()) * _scale);
    }
    else
    {
        // If the constraint is not set, calculate a rough boundry with 1 aspect ratio
        RealType aspectRatio = 0.85;
        const RealType tolerentArea = _totalCellArea * (1 + maxWhiteSpace);
        const RealType xLen = std::sqrt(tolerentArea * aspectRatio);
        const RealType yLen= tolerentArea / xLen;
        _boundary.set(-xLen  / 2 , - yLen / 2 , xLen / 2 , yLen / 2 );
        INF("NlpWnconj::%s: automatical set boundary to be %s \n", __FUNCTION__, _boundary.toStr().c_str());
    }

    _totalCellArea = 0;
    for (IndexType cellIdx = 0; cellIdx < _db.numCells(); ++cellIdx)
    {
        auto bbox = _db.cell(cellIdx).cellBBox();
        _totalCellArea +=  bbox.xLen() * _scale * bbox.yLen() * _scale;
    }

    // Default sym axis is at the middle
    _defaultSymAxis = (_boundary.xLo() + _boundary.xHi()) / 2;
}

template<typename nlp_settings>
void NlpGPlacerBase<nlp_settings>::initVariables()
{
    // The number of nlp problem variables
    _numCells = _db.numCells();
    IntType size = _db.numCells() * 2 + _db.numSymGroups();
    _pl.resize(size);
    _plx = std::make_shared<EigenMap>(EigenMap(_pl.data(), _numCells));
    _ply = std::make_shared<EigenMap>(EigenMap(_pl.data() + _numCells, _db.numCells()));
    _sym = std::make_shared<EigenMap>(EigenMap(_pl.data() + 2* _numCells, _db.numSymGroups()));
#ifndef MULTI_SYM_GROUP
    (*_sym)(0) = _defaultSymAxis; // Set the default symmtry axisx`
#endif
}

template<typename nlp_settings>
void NlpGPlacerBase<nlp_settings>::initPlace()
{
    auto initPlace = init_place_trait::construct(*this);
    init_place_trait::initPlace(initPlace, *this);
}

template<typename nlp_settings>
void NlpGPlacerBase<nlp_settings>::initOperators()
{
    auto getAlphaFunc = [&]()
    {
        return _alpha;
    };
    auto getLambdaFuncOvr = [&]()
    {
        return 1.0;
    };
    auto getLambdaFuncBoundary = [&]()
    {
        return 1.0;
    };
    auto getLambdaFuncHpwl = [&]()
    {
        return 1.0;
    };
    auto getLambdaFuncAsym = [&]()
    {
        return 1.0;
    };
    auto getLambdaFuncCosine = [&]()
    {
        return 1.0;
    };
    auto getVarFunc = [&] (IndexType cellIdx, Orient2DType orient)
    {
        return _pl(plIdx(cellIdx, orient));
    };

    auto calculatePinOffset = [&](IndexType pinIdx)
    {
        const auto &pin = _db.pin(pinIdx);
        IndexType cellIdx = pin.cellIdx();
        const auto &cell = _db.cell(cellIdx);
        // Get the cell location from the input arguments
        XY<RealType> midLoc = XY<RealType>(pin.midLoc().x(), pin.midLoc().y()) * _scale;
        XY<RealType> cellLoLoc = XY<RealType>(cell.cellBBox().xLo(), cell.cellBBox().yLo()) * _scale;
        return midLoc - cellLoLoc;
    };
    // Hpwl
    for (IndexType netIdx = 0; netIdx < _db.numNets(); ++netIdx)
    {
        const auto &net = _db.net(netIdx);
        _hpwlOps.emplace_back(nlp_hpwl_type(getAlphaFunc, getLambdaFuncHpwl));
        auto &op = _hpwlOps.back();
        op.setWeight(net.weight());
        for (IndexType idx = 0; idx < net.numPinIdx(); ++idx)
        {
            // Get the pin location referenced to the cell
            IndexType pinIdx = net.pinIdx(idx);
            auto pinLoc = calculatePinOffset(pinIdx);
            op.addVar(_db.pin(pinIdx).cellIdx(), pinLoc.x(), pinLoc.y());
        }
        op.setGetVarFunc(getVarFunc);
    }
    // Pair-wise cell overlapping
    for (IndexType cellIdxI = 0; cellIdxI < _db.numCells(); ++cellIdxI)
    {
        const auto cellBBoxI = _db.cell(cellIdxI).cellBBox();
        for (IndexType cellIdxJ = cellIdxI + 1; cellIdxJ < _db.numCells(); ++cellIdxJ)
        {
            const auto cellBBoxJ = _db.cell(cellIdxJ).cellBBox();
            _ovlOps.emplace_back(nlp_ovl_type(
                        cellIdxI,
                        cellBBoxI.xLen() * _scale,
                        cellBBoxI.yLen() * _scale,
                        cellIdxJ,
                        cellBBoxJ.xLen() * _scale,
                        cellBBoxJ.yLen() * _scale,
                        getAlphaFunc,
                        getLambdaFuncOvr
                        ));
            _ovlOps.back().setGetVarFunc(getVarFunc);
        }
    }
    // Out of boundary
    for (IndexType cellIdx = 0; cellIdx < _db.numCells(); ++cellIdx)
    {
        const auto &cellBBox = _db.cell(cellIdx).cellBBox();
        _oobOps.emplace_back(nlp_oob_type(
                    cellIdx,
                    cellBBox.xLen() * _scale,
                    cellBBox.yLen() * _scale,
                    &_boundary,
                    getAlphaFunc,
                    getLambdaFuncBoundary
                    ));
        _oobOps.back().setGetVarFunc(getVarFunc);
    }
    // Asym
    for (IndexType symGrpIdx = 0; symGrpIdx < _db.numSymGroups(); ++symGrpIdx)
    {
        const auto &symGrp = _db.symGroup(symGrpIdx);
        _asymOps.emplace_back(nlp_asym_type(symGrpIdx, getLambdaFuncAsym));
        for (const auto &symPair : symGrp.vSymPairs())
        {
            IndexType cellIdxI = symPair.firstCell();
            IndexType cellIdxJ = symPair.secondCell();
            RealType widthI = _db.cell(cellIdxI).cellBBox().xLen() * _scale;
            _asymOps.back().addSymPair(cellIdxI, cellIdxJ, widthI);
        }
        for (const auto &ssCellIdx : symGrp.vSelfSyms())
        {
            RealType width = _db.cell(ssCellIdx).cellBBox().xLen() * _scale;
            _asymOps.back().addSelfSym(ssCellIdx, width);
        }
        _asymOps.back().setGetVarFunc(getVarFunc);
    }
    // Signal path
    SigPathMgr pathMgr(_db);
    for (const auto &seg : pathMgr.vSegList())
    {
        IndexType sPinIdx = seg.beginPinFirstSeg();
        IndexType midPinIdxA = seg.endPinFirstSeg();
        IndexType midPinIdxB = seg.beginPinSecondSeg();
        IndexType tPinIdx = seg.endPinSecondSeg();

        const auto &sPin = _db.pin(sPinIdx);
        IndexType sCellIdx = sPin.cellIdx();
        const auto &mPinA = _db.pin(midPinIdxA);
        IndexType mCellIdx = mPinA.cellIdx();
        const auto &tPin = _db.pin(tPinIdx);
        IndexType tCellIdx = tPin.cellIdx();

        auto sOffset = calculatePinOffset(sPinIdx);
        auto midOffsetA = calculatePinOffset(midPinIdxA);
        auto midOffsetB = calculatePinOffset(midPinIdxB);
        auto tOffset = calculatePinOffset(tPinIdx);
        _cosOps.emplace_back(sCellIdx, sOffset,
                mCellIdx, midOffsetA, midOffsetB,
                tCellIdx, tOffset,
                getLambdaFuncCosine);
        _cosOps.back().setGetVarFunc(getVarFunc);
    }
}


template<typename nlp_settings>
void NlpGPlacerBase<nlp_settings>::writeOut()
{
    // find the min value
    RealType minX =1e10; 
    RealType minY = 1e10;
    for (IndexType cellIdx = 0; cellIdx < _db.numCells(); ++cellIdx)
    {
        if ((*_plx)(cellIdx) < minX)
        {
            minX = (*_plx)(cellIdx);
        }
        if ((*_ply)(cellIdx) < minY)
        {
            minY = (*_ply)(cellIdx);
        }
    }
    // Dump the cell locations to database
    for (IndexType cellIdx = 0; cellIdx < _db.numCells(); ++cellIdx)
    {
        auto & cell = _db.cell(cellIdx);
        LocType xLo = ::klib::autoRound<LocType>(((*_plx)(cellIdx) - minX) / _scale + _db.parameters().layoutOffset());
        LocType yLo = ::klib::autoRound<LocType>(((*_ply)(cellIdx) - minY) / _scale + _db.parameters().layoutOffset());
        _db.cell(cellIdx).setXLoc(xLo - cell.cellBBox().xLo());
        _db.cell(cellIdx).setYLoc(yLo - cell.cellBBox().yLo());
    }
}

template<typename nlp_settings>
void NlpGPlacerBase<nlp_settings>::constructTasks()
{
    constructObjTasks();
    constructStopConditionTask();
}

template<typename nlp_settings>
void NlpGPlacerBase<nlp_settings>::constructObjTasks()
{
    constructObjectiveCalculationTasks();
    constructSumObjTasks();
#ifdef DEBUG_SINGLE_THREAD_GP
    constructWrapObjTask();
#endif
}

template<typename nlp_settings>
void NlpGPlacerBase<nlp_settings>::constructObjectiveCalculationTasks()
{
    using EvaObjTask = EvaObjTask<nlp_numerical_type>;
    for (const auto &hpwl : _hpwlOps)
    {
        auto eva = [&]() { return diff::placement_differentiable_traits<nlp_hpwl_type>::evaluate(hpwl);};
        _evaHpwlTasks.emplace_back(Task<EvaObjTask>(EvaObjTask(eva)));
    }
    for (const auto &ovl : _ovlOps)
    {
        auto eva = [&]() { return diff::placement_differentiable_traits<nlp_ovl_type>::evaluate(ovl);};
        _evaOvlTasks.emplace_back(Task<EvaObjTask>(EvaObjTask(eva)));
    }
    for (const auto &oob : _oobOps)
    {
        auto eva = [&]() { return diff::placement_differentiable_traits<nlp_oob_type>::evaluate(oob);};
        _evaOobTasks.emplace_back(Task<EvaObjTask>(EvaObjTask(eva)));
    }
    for (const auto &asym : _asymOps)
    {
        auto eva = [&]() { return diff::placement_differentiable_traits<nlp_asym_type>::evaluate(asym);};
        _evaAsymTasks.emplace_back(Task<EvaObjTask>(EvaObjTask(eva)));
    }
    for (const auto &cos : _cosOps)
    {
        auto eva = [&]() { return diff::placement_differentiable_traits<nlp_cos_type>::evaluate(cos);};
        _evaCosTasks.emplace_back(Task<EvaObjTask>(EvaObjTask(eva)));
    }
}

template<typename nlp_settings>
void NlpGPlacerBase<nlp_settings>::constructSumObjTasks()
{
    auto hpwl = [&]() 
    {
        _objHpwl = 0.0;
        for (const auto &eva : _evaHpwlTasks)
        {
            _objHpwl += eva.taskData().obj();
        }
    };
    _sumObjHpwlTask = Task<FuncTask>(FuncTask(hpwl));
    auto ovl = [&]() 
    {
        _objOvl = 0.0;
        for (const auto &eva : _evaOvlTasks)
        {
            _objOvl += eva.taskData().obj();
        }
    };
    _sumObjOvlTask = Task<FuncTask>(FuncTask(ovl));
    auto oob = [&]()
    {
        _objOob = 0.0;
        for (const auto &eva : _evaOobTasks)
        {
            _objOob += eva.taskData().obj();
        }
    };
    _sumObjOobTask = Task<FuncTask>(FuncTask(oob));
    auto asym = [&]()
    {
        _objAsym = 0.0;
        for (const auto &eva : _evaAsymTasks)
        {
            _objAsym += eva.taskData().obj();
        }
    };
    _sumObjAsymTask = Task<FuncTask>(FuncTask(asym));
    auto cos = [&]()
    {
        _objCos = 0.0;
        for (const auto &eva : _evaCosTasks)
        {
            _objCos += eva.taskData().obj();
        }
    };
    _sumObjCosTask = Task<FuncTask>(FuncTask(cos));
    auto all = [&]()
    {
        _obj = 0.0;
        _obj += _objHpwl;
        _obj += _objOvl;
        _obj += _objOob;
        _obj += _objAsym;
        _obj += _objCos;
    };
    _sumObjAllTask = Task<FuncTask>(FuncTask(all));
}

#ifdef DEBUG_SINGLE_THREAD_GP
template<typename nlp_settings>
void NlpGPlacerBase<nlp_settings>::constructWrapObjTask()
{
    auto hpwl = [&]()
    {
        for (auto &eva : _evaHpwlTasks)
        {
            eva.run();
        }
        _sumObjHpwlTask.run();
    };
    _wrapObjHpwlTask = Task<FuncTask>(FuncTask(hpwl));
    auto ovl = [&]()
    {
        for (auto &eva : _evaOvlTasks)
        {
            eva.run();
        }
        _sumObjOvlTask.run();
    };
    _wrapObjOvlTask = Task<FuncTask>(FuncTask(ovl));
    auto oob = [&]()
    {
        for (auto &eva : _evaOobTasks)
        {
            eva.run();
        }
        _sumObjOobTask.run();
    };
    _wrapObjOobTask = Task<FuncTask>(FuncTask(oob));
    auto asym = [&]()
    {
        for (auto &eva : _evaAsymTasks)
        {
            eva.run();
        }
        _sumObjAsymTask.run();
    };
    _wrapObjAsymTask = Task<FuncTask>(FuncTask(asym));
    auto cos = [&]()
    {
        for (auto &eva : _evaCosTasks)
        {
            eva.run();
        }
        _sumObjCosTask.run();
    };
    _wrapObjCosTask = Task<FuncTask>(FuncTask(cos));
    auto all = [&]()
    {
        _wrapObjHpwlTask.run();
        _wrapObjOvlTask.run();
        _wrapObjOobTask.run();
        _wrapObjAsymTask.run();
        _wrapObjCosTask.run();
        _sumObjAllTask.run();
    };
    _wrapObjAllTask = Task<FuncTask>(FuncTask(all));
}
#endif //DEBUG_SINGLE_THREAD_GP

template<typename nlp_settings>
void NlpGPlacerBase<nlp_settings>::constructOptimizationKernelTasks()
{
    constructStopConditionTask();
}

template<typename nlp_settings>
void NlpGPlacerBase<nlp_settings>::constructStopConditionTask()
{
    auto stopCondition = [&]()
    {
        return stop_condition_trait::stopPlaceCondition(_stopCondition, *this);
    };
    _checkStopConditionTask = Task<ConditionTask>(ConditionTask(stopCondition));
}


/* FirstOrder */

template<typename nlp_settings>
void NlpGPlacerFirstOrder<nlp_settings>::optimize()
{
    tf::Executor exe; 
    this->_wrapObjAllTask.regTask(this->_taskflow);
    this->_wrapCalcGradTask.regTask(this->_taskflow);
    exe.run(this->_taskflow).wait();
    std::cout<<this->_pl <<std::endl;
}

template<typename nlp_settings>
void NlpGPlacerFirstOrder<nlp_settings>::initProblem()
{
    this->initHyperParams();
    this->initBoundaryParams();
    this->initVariables();
    this->initFirstOrderGrad();
}

template<typename nlp_settings>
void NlpGPlacerFirstOrder<nlp_settings>::initFirstOrderGrad()
{
    const IntType size = this->_db.numCells() * 2 + this->_db.numSymGroups();
    _grad.resize(size);
    _gradHpwl.resize(size);
    _gradOvl.resize(size);
    _gradOob.resize(size);
    _gradAsym.resize(size);
    _gradCos.resize(size);
}

template<typename nlp_settings>
void NlpGPlacerFirstOrder<nlp_settings>::constructTasks()
{
    this->constructObjTasks();
    this->constructStopConditionTask();
    constructFirstOrderTasks();
}

template<typename nlp_settings>
void NlpGPlacerFirstOrder<nlp_settings>::constructFirstOrderTasks()
{
    constructCalcPartialsTasks();
    constructUpdatePartialsTasks();
    constructClearGradTasks();
    constructSumGradTask();
#ifdef DEBUG_SINGLE_THREAD_GP
    constructWrapCalcGradTask();
#endif
}

template<typename nlp_settings>
void NlpGPlacerFirstOrder<nlp_settings>::constructCalcPartialsTasks()
{
    using Hpwl = CalculateOperatorPartialTask<nlp_hpwl_type, EigenVector>;
    using Ovl = CalculateOperatorPartialTask<nlp_ovl_type, EigenVector>;
    using Oob = CalculateOperatorPartialTask<nlp_oob_type, EigenVector>;
    using Asym = CalculateOperatorPartialTask<nlp_asym_type, EigenVector>;
    using Cos = CalculateOperatorPartialTask<nlp_cos_type, EigenVector>;
    for (auto &hpwlOp : this->_hpwlOps)
    {
        _calcHpwlPartialTasks.emplace_back(Task<Hpwl>(Hpwl(&hpwlOp)));
    }
    for (auto &ovlOp : this->_ovlOps)
    {
        _calcOvlPartialTasks.emplace_back(Task<Ovl>(Ovl(&ovlOp)));
    }
    for (auto &oobOp : this->_oobOps)
    {
        _calcOobPartialTasks.emplace_back(Task<Oob>(Oob(&oobOp)));
    }
    for (auto &asymOp : this->_asymOps)
    {
        _calcAsymPartialTasks.emplace_back(Task<Asym>(Asym(&asymOp)));
    }
    for (auto &cosOp : this->_cosOps)
    {
        _calcCosPartialTasks.emplace_back(Task<Cos>(Cos(&cosOp)));
    }
}

template<typename nlp_settings>
void NlpGPlacerFirstOrder<nlp_settings>::constructUpdatePartialsTasks()
{
    using Hpwl = UpdateGradientFromPartialTask<nlp_hpwl_type, EigenVector>;
    using Ovl = UpdateGradientFromPartialTask<nlp_ovl_type, EigenVector>;
    using Oob = UpdateGradientFromPartialTask<nlp_oob_type, EigenVector>;
    using Asym = UpdateGradientFromPartialTask<nlp_asym_type, EigenVector>;
    using Cos = UpdateGradientFromPartialTask<nlp_cos_type, EigenVector>;
    auto getIdxFunc = [&](IndexType cellIdx, Orient2DType orient) { return this->plIdx(cellIdx, orient); }; // wrapper the convert cell idx to pl idx
    for (auto &hpwl : _calcHpwlPartialTasks)
    {
        _updateHpwlPartialTasks.emplace_back(Task<Hpwl>(Hpwl(hpwl.taskDataPtr(), &_gradHpwl, getIdxFunc)));
    }
    for (auto &ovl : _calcOvlPartialTasks)
    {
        _updateOvlPartialTasks.emplace_back(Task<Ovl>(Ovl(ovl.taskDataPtr(), &_gradOvl, getIdxFunc)));
    }
    for (auto &oob : _calcOobPartialTasks)
    {
        _updateOobPartialTasks.emplace_back(Task<Oob>(Oob(oob.taskDataPtr(), &_gradOob, getIdxFunc)));
    }
    for (auto &asym : _calcAsymPartialTasks)
    {
        _updateAsymPartialTasks.emplace_back(Task<Asym>(Asym(asym.taskDataPtr(), &_gradAsym, getIdxFunc)));
    }
    for (auto &cos : _calcCosPartialTasks)
    {
        _updateCosPartialTasks.emplace_back(Task<Cos>(Cos(cos.taskDataPtr(), &_gradCos, getIdxFunc)));
    }
}

template<typename nlp_settings>
void NlpGPlacerFirstOrder<nlp_settings>::constructClearGradTasks()
{
    _clearGradTask = Task<FuncTask>(FuncTask([&]() { _grad.setZero(); }));
    _clearHpwlGradTask = Task<FuncTask>(FuncTask([&]() { _gradHpwl.setZero(); }));
    _clearOvlGradTask = Task<FuncTask>(FuncTask([&]() { _gradOvl.setZero(); }));
    _clearOobGradTask = Task<FuncTask>(FuncTask([&]() { _gradOob.setZero(); }));
    _clearAsymGradTask = Task<FuncTask>(FuncTask([&]() { _gradAsym.setZero(); }));
    _clearCosGradTask = Task<FuncTask>(FuncTask([&]() { _gradCos.setZero(); }));
}

template<typename nlp_settings>
void NlpGPlacerFirstOrder<nlp_settings>::constructSumGradTask()
{
    _sumGradTask = Task<FuncTask>(FuncTask([&](){ _grad = _gradHpwl + _gradOvl + _gradOob + _gradAsym + _gradCos; }));
}

#ifdef DEBUG_SINGLE_THREAD_GP
template<typename nlp_settings>
void NlpGPlacerFirstOrder<nlp_settings>::constructWrapCalcGradTask()
{
    auto calcGrad = [&]()
    {
        _clearGradTask.run();
        _clearHpwlGradTask.run();
        _clearOvlGradTask.run();
        _clearOobGradTask.run();
        _clearAsymGradTask.run();
        _clearCosGradTask.run();
        for (auto & calc : _calcHpwlPartialTasks) { calc.run(); }
        for (auto & update : _updateHpwlPartialTasks) { update.run(); }
        for (auto & calc : _calcOvlPartialTasks) { calc.run(); }
        for (auto & update : _updateOvlPartialTasks) { update.run(); }
        for (auto & calc : _calcOobPartialTasks) { calc.run(); }
        for (auto & update : _updateOobPartialTasks) { update.run(); }
        for (auto & calc : _calcAsymPartialTasks) { calc.run(); }
        for (auto & update : _updateAsymPartialTasks) { update.run(); }
        for (auto & calc : _calcCosPartialTasks) { calc.run(); }
        for (auto & update : _updateCosPartialTasks) { update.run(); }
        _sumGradTask.run();
    };
    _wrapCalcGradTask = Task<FuncTask>(FuncTask(calcGrad));
}
#endif


// declare the default settings for linking
template class NlpGPlacerBase<nlp::nlp_default_settings>;
template class NlpGPlacerFirstOrder<nlp::nlp_default_settings>;

PROJECT_NAMESPACE_END
