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
    initOptimizationKernelMembers();
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
    _numVariables = size;
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
void NlpGPlacerBase<nlp_settings>::alignToSym()
{
    for (const auto &symGrp : _db.vSymGrpArray())
    {
        for (const auto &symPair : symGrp.vSymPairs())
        {
            IndexType cell1 = symPair.firstCell();
            IndexType cell2 = symPair.secondCell();
            _pl(plIdx(cell2, Orient2DType::HORIZONTAL)) = 2 * _defaultSymAxis - _pl(plIdx(cell1, Orient2DType::HORIZONTAL)) + _db.cell(cell1).cellBBox().xLen() * _scale;
            _pl(plIdx(cell2, Orient2DType::VERTICAL)) = _pl(plIdx(cell1, Orient2DType::VERTICAL));
        }
        for (const auto cellIdx : symGrp.vSelfSyms())
        {
            _pl(plIdx(cellIdx, Orient2DType::HORIZONTAL)) = _defaultSymAxis + _db.cell(cellIdx).cellBBox().xLen() * _scale ;
        }
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
    constructWrapObjTask();
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
        #pragma omp parallel for schedule(static)
        for (IndexType idx = 0; idx < _evaHpwlTasks.size(); ++idx)
        {
            _evaHpwlTasks[idx].run();
        }
        _sumObjHpwlTask.run();
    };
    _wrapObjHpwlTask = Task<FuncTask>(FuncTask(hpwl));
    auto ovl = [&]()
    {
        #pragma omp parallel for schedule(static)
        for (IndexType idx = 0; idx < _evaOvlTasks.size(); ++idx)
        {
            _evaOvlTasks[idx].run();
        }
        _sumObjOvlTask.run();
    };
    _wrapObjOvlTask = Task<FuncTask>(FuncTask(ovl));
    auto oob = [&]()
    {
        #pragma omp parallel for schedule(static)
        for (IndexType idx = 0; idx < _evaOobTasks.size(); ++idx)
        {
            _evaOobTasks[idx].run();
        }
        _sumObjOobTask.run();
    };
    _wrapObjOobTask = Task<FuncTask>(FuncTask(oob));
    auto asym = [&]()
    {
        #pragma omp parallel for schedule(static)
        for (IndexType idx = 0; idx < _evaAsymTasks.size(); ++idx)
        {
            _evaAsymTasks[idx].run();
        }
        _sumObjAsymTask.run();
    };
    _wrapObjAsymTask = Task<FuncTask>(FuncTask(asym));
    auto cos = [&]()
    {
        #pragma omp parallel for schedule(static)
        for (IndexType idx = 0; idx < _evaCosTasks.size(); ++idx)
        {
            _evaCosTasks[idx].run();
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
        return stop_condition_trait::stopPlaceCondition(*this, _stopCondition);
    };
    _checkStopConditionTask = Task<ConditionTask>(ConditionTask(stopCondition));
}

#ifdef IDEAPLACE_TASKFLOR_FOR_GRAD_OBJ_

template<typename nlp_settings>
void NlpGPlacerBase<nlp_settings>::regEvaHpwlTaskflow(tf::Taskflow &tfFlow)
{
    _sumObjHpwlTask.regTask(tfFlow);
    for (IndexType idx = 0; idx < _hpwlOps.size(); ++idx)
    {
        _evaHpwlTasks[idx].regTask(tfFlow);
        _evaHpwlTasks[idx].precede(_sumObjHpwlTask);
    }
}

template<typename nlp_settings>
void NlpGPlacerBase<nlp_settings>::regEvaOvlTaskflow(tf::Taskflow &tfFlow)
{
    _sumObjOvlTask.regTask(tfFlow);
    for (IndexType idx = 0; idx < _ovlOps.size(); ++idx)
    {
        _evaOvlTasks[idx].regTask(tfFlow);
        _evaOvlTasks[idx].precede(_sumObjOvlTask);
    }
}

template<typename nlp_settings>
void NlpGPlacerBase<nlp_settings>::regEvaOobTaskflow(tf::Taskflow &tfFlow)
{
    _sumObjOobTask.regTask(tfFlow);
    for (IndexType idx = 0; idx < _oobOps.size(); ++idx)
    {
        _evaOobTasks[idx].regTask(tfFlow);
        _evaOobTasks[idx].precede(_sumObjOobTask);
    }
}

template<typename nlp_settings>
void NlpGPlacerBase<nlp_settings>::regEvaAsymTaskflow(tf::Taskflow &tfFlow)
{
    _sumObjAsymTask.regTask(tfFlow);
    for (IndexType idx = 0; idx < _asymOps.size(); ++idx)
    {
        _evaAsymTasks[idx].regTask(tfFlow);
        _evaAsymTasks[idx].precede(_sumObjAsymTask);
    }
}

template<typename nlp_settings>
void NlpGPlacerBase<nlp_settings>::regEvaCosTaskflow(tf::Taskflow &tfFlow)
{
    _sumObjCosTask.regTask(tfFlow);
    for (IndexType idx = 0; idx < _cosOps.size(); ++idx)
    {
        _evaCosTasks[idx].regTask(tfFlow);
        _evaCosTasks[idx].precede(_sumObjCosTask);
    }
}

template<typename nlp_settings>
void NlpGPlacerBase<nlp_settings>::regEvaAllObjTaskflow(tf::Taskflow &tfFlow)
{
    _sumObjAllTask.regTask(tfFlow);
    this->regEvaHpwlTaskflow(tfFlow);
    this->regEvaOvlTaskflow(tfFlow);
    this->regEvaOobTaskflow(tfFlow);
    this->regEvaAsymTaskflow(tfFlow);
    this->regEvaCosTaskflow(tfFlow);
    _sumObjHpwlTask.precede(_sumObjAllTask);
    _sumObjOvlTask.precede(_sumObjAllTask);
    _sumObjOobTask.precede(_sumObjAllTask);
    _sumObjAsymTask.precede(_sumObjAllTask);
    _sumObjCosTask.precede(_sumObjAllTask);
}

#endif //IDEAPLACE_TASKFLOR_FOR_GRAD_OBJ_

/* FirstOrder */

template<typename nlp_settings>
void NlpGPlacerFirstOrder<nlp_settings>::optimize()
{
    WATCH_QUICK_START();
    // setting up the multipliers
    this->_wrapObjAllTask.run();
    _wrapCalcGradTask.run();
    optm_type optm;
    mult_type multiplier = mult_trait::construct(*this);
    mult_trait::init(*this, multiplier);
    IntType iter = 0;
    do
    {
        std::string debugGdsFilename  = "./debug/";
        debugGdsFilename += "gp_iter_" + std::to_string(iter)+".gds";
        base_type::drawCurrentLayout(debugGdsFilename);
        DBG("iter %d \n", iter);
        optm_trait::optimize(*this, optm);
        mult_trait::update(*this, multiplier);
        DBG("obj %f hpwl %f ovl %f oob %f asym %f cos %f \n", this->_obj, this->_objHpwl, this->_objOvl, this->_objOob, this->_objAsym, this->_objCos);
        ++iter;
    } while (not base_type::stop_condition_trait::stopPlaceCondition(*this, this->_stopCondition));
    auto end = WATCH_QUICK_END();
    //std::cout<<"grad"<<"\n"<< _grad <<std::endl;
    std::cout<<"time "<< end / 1000 <<" ms" <<std::endl;
    this->writeOut();
}

template<typename nlp_settings>
void NlpGPlacerFirstOrder<nlp_settings>::initProblem()
{
    this->initHyperParams();
    this->initBoundaryParams();
    this->initVariables();
    this->initOptimizationKernelMembers();
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
    constructWrapCalcGradTask();
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
    _sumHpwlGradTask = Task<FuncTask>(FuncTask([&](){ for (auto &upd : _updateHpwlPartialTasks){ upd.run(); }}));
    _sumOvlGradTask = Task<FuncTask>(FuncTask([&](){ for (auto &upd : _updateOvlPartialTasks){ upd.run(); }}));
    _sumOobGradTask = Task<FuncTask>(FuncTask([&](){ for (auto &upd : _updateOobPartialTasks){ upd.run(); }}));
    _sumAsymGradTask = Task<FuncTask>(FuncTask([&](){ for (auto &upd : _updateAsymPartialTasks){ upd.run(); }}));
    _sumCosGradTask = Task<FuncTask>(FuncTask([&](){ for (auto &upd : _updateCosPartialTasks){ upd.run(); }}));
    _sumGradTask = Task<FuncTask>(FuncTask([&](){ _grad = _gradHpwl + _gradOvl + _gradOob + _gradAsym + _gradCos; }));
}

template<typename nlp_settings>
void NlpGPlacerFirstOrder<nlp_settings>::constructWrapCalcGradTask()
{
    auto calcGradLambda = [&]()
    {
        _clearGradTask.run();
        _clearHpwlGradTask.run();
        _clearOvlGradTask.run();
        _clearOobGradTask.run();
        _clearAsymGradTask.run();
        _clearCosGradTask.run();
        #pragma omp parallel for schedule(static)
        for (IndexType i = 0; i < _calcHpwlPartialTasks.size(); ++i ) { _calcHpwlPartialTasks[i].run(); }
        for (IndexType i = 0; i < _updateHpwlPartialTasks.size(); ++i ) { _updateHpwlPartialTasks[i].run(); }
        #pragma omp parallel for schedule(static)
        for (IndexType i = 0; i < _calcOvlPartialTasks.size(); ++i ) { _calcOvlPartialTasks[i].run(); }
        for (IndexType i = 0; i < _updateOvlPartialTasks.size(); ++i ) { _updateOvlPartialTasks[i].run(); }
        #pragma omp parallel for schedule(static)
        for (IndexType i = 0; i < _calcOobPartialTasks.size(); ++i ) { _calcOobPartialTasks[i].run(); }
        for (IndexType i = 0; i < _updateOobPartialTasks.size(); ++i ) { _updateOobPartialTasks[i].run(); }
        #pragma omp parallel for schedule(static)
        for (IndexType i = 0; i < _calcAsymPartialTasks.size(); ++i ) { _calcAsymPartialTasks[i].run(); }
        for (IndexType i = 0; i < _updateAsymPartialTasks.size(); ++i ) { _updateAsymPartialTasks[i].run(); }
        #pragma omp parallel for schedule(static)
        for (IndexType i = 0; i < _calcCosPartialTasks.size(); ++i ) { _calcCosPartialTasks[i].run(); }
        for (IndexType i = 0; i < _updateCosPartialTasks.size(); ++i ) { _updateCosPartialTasks[i].run(); }
        _sumGradTask.run();
    };
    _wrapCalcGradTask = Task<FuncTask>(FuncTask(calcGradLambda));
}

#ifdef IDEAPLACE_TASKFLOR_FOR_GRAD_OBJ_

template<typename nlp_settings>
void NlpGPlacerFirstOrder<nlp_settings>::regCalcHpwlGradTaskFlow(tf::Taskflow &tfFlow)
{
    IntType opSize = this->_hpwlOps.size();
    _clearHpwlGradTask.regTask(tfFlow);
    _sumHpwlGradTask.regTask(tfFlow);
    auto taskIntEndCalc = tfFlow.parallel_for(0, opSize, 1, [&] (const int idx)
            {
                _calcHpwlPartialTasks[idx].run();
            });
    _clearHpwlGradTask.precede(taskIntEndCalc.first);
    //taskIntEndCalc.second.precede(_sumHpwlGradTask.tfTask());
    taskIntEndCalc.second.precede(_endGradCalcTask.tfTask());
    _endGradCalcTask.precede(_sumHpwlGradTask);
}


template<typename nlp_settings>
void NlpGPlacerFirstOrder<nlp_settings>::regCalcOvlGradTaskFlow(tf::Taskflow &tfFlow)
{
    IntType opSize = this->_ovlOps.size();
    _clearOvlGradTask.regTask(tfFlow);
    _sumOvlGradTask.regTask(tfFlow);
    auto taskIntEndCalc = tfFlow.parallel_for(0, opSize, 1, [&] (const int idx)
            {
                _calcOvlPartialTasks[idx].run();
            });
    _clearOvlGradTask.precede(taskIntEndCalc.first);
    //taskIntEndCalc.second.precede(_sumOvlGradTask.tfTask());
    taskIntEndCalc.second.precede(_endGradCalcTask.tfTask());
    _endGradCalcTask.precede(_sumOvlGradTask);
}


template<typename nlp_settings>
void NlpGPlacerFirstOrder<nlp_settings>::regCalcOobGradTaskFlow(tf::Taskflow &tfFlow)
{
    IntType opSize = this->_oobOps.size();
    _clearOobGradTask.regTask(tfFlow);
    _sumOobGradTask.regTask(tfFlow);
    auto taskIntEndCalc = tfFlow.parallel_for(0, opSize, 1, [&] (const int idx)
            {
                _calcOobPartialTasks[idx].run();
            });
    _clearOobGradTask.precede(taskIntEndCalc.first);
    //taskIntEndCalc.second.precede(_sumOobGradTask.tfTask());
    taskIntEndCalc.second.precede(_endGradCalcTask.tfTask());
    _endGradCalcTask.precede(_sumOobGradTask);
}

template<typename nlp_settings>
void NlpGPlacerFirstOrder<nlp_settings>::regCalcAsymGradTaskFlow(tf::Taskflow &tfFlow)
{
    IntType opSize = this->_asymOps.size();
    _clearAsymGradTask.regTask(tfFlow);
    _sumAsymGradTask.regTask(tfFlow);
    auto taskIntEndCalc = tfFlow.parallel_for(0, opSize, 1, [&] (const int idx)
            {
                _calcAsymPartialTasks[idx].run();
            });
    _clearAsymGradTask.precede(taskIntEndCalc.first);
    //taskIntEndCalc.second.precede(_sumAsymGradTask.tfTask());
    taskIntEndCalc.second.precede(_endGradCalcTask.tfTask());
    _endGradCalcTask.precede(_sumAsymGradTask);
}

template<typename nlp_settings>
void NlpGPlacerFirstOrder<nlp_settings>::regCalcCosGradTaskFlow(tf::Taskflow &tfFlow)
{
    IntType opSize = this->_cosOps.size();
    _clearCosGradTask.regTask(tfFlow);
    _sumCosGradTask.regTask(tfFlow);
    auto taskIntEndCalc = tfFlow.parallel_for(0, opSize, 1, [&] (const int idx)
            {
                _calcCosPartialTasks[idx].run();
            });
    _clearCosGradTask.precede(taskIntEndCalc.first);
    //taskIntEndCalc.second.precede(_sumCosGradTask.tfTask());
    taskIntEndCalc.second.precede(_endGradCalcTask.tfTask());
    _endGradCalcTask.precede(_sumCosGradTask);
}

template<typename nlp_settings>
void NlpGPlacerFirstOrder<nlp_settings>::regCalcGradForLoop(tf::Taskflow &tfFlow)
{
    //_sumHpwlGradTask.regTask(tfFlow);
    //_sumOvlGradTask.regTask(tfFlow);
    //_sumOobGradTask.regTask(tfFlow);
    //_sumAsymGradTask.regTask(tfFlow);
    //_sumCosGradTask.regTask(tfFlow);
    _beginGradCalcTask.regTask(tfFlow);
    _endGradCalcTask.regTask(tfFlow);
    for (auto & op : _calcHpwlPartialTasks)
    {
        op.regTask(tfFlow);
        _beginGradCalcTask.precede(op);
        op.precede(_endGradCalcTask);
    }
    for (auto & op : _calcOvlPartialTasks)
    {
        op.regTask(tfFlow);
        _beginGradCalcTask.precede(op);
        op.precede(_endGradCalcTask);
    }
    for (auto & op : _calcOobPartialTasks)
    {
        op.regTask(tfFlow);
        _beginGradCalcTask.precede(op);
        op.precede(_endGradCalcTask);
    }
    for (auto & op : _calcAsymPartialTasks)
    {
        op.regTask(tfFlow);
        _beginGradCalcTask.precede(op);
        op.precede(_endGradCalcTask);
    }
    for (auto & op : _calcCosPartialTasks)
    {
        op.regTask(tfFlow);
        _beginGradCalcTask.precede(op);
        op.precede(_endGradCalcTask);
    }
}

template<typename nlp_settings>
void NlpGPlacerFirstOrder<nlp_settings>::regCalcAllGradTaskFlow(tf::Taskflow &tfFlow)
{
    //_sumGradTask.regTask(tfFlow);
    //_endGradCalcTask.regTask(tfFlow);
    //this->regCalcHpwlGradTaskFlow(tfFlow);
    //this->regCalcOvlGradTaskFlow(tfFlow);
    //this->regCalcOobGradTaskFlow(tfFlow);
    //this->regCalcAsymGradTaskFlow(tfFlow);
    //this->regCalcCosGradTaskFlow(tfFlow);
    //_sumHpwlGradTask.precede(_sumGradTask);
    //_sumOvlGradTask.precede(_sumGradTask);
    //_sumOobGradTask.precede(_sumGradTask);
    //_sumAsymGradTask.precede(_sumGradTask);
    //_sumCosGradTask.precede(_sumGradTask);
    
    this->regCalcGradForLoop(tfFlow);
}
#endif // IDEAPLACE_TASKFLOR_FOR_GRAD_OBJ_


#ifdef DEBUG_GR
#ifdef DEBUG_DRAW
PROJECT_NAMESPACE_END
#include "writer/gdsii/WriteGds.h"
PROJECT_NAMESPACE_BEGIN
template<typename nlp_settings>
void NlpGPlacerBase<nlp_settings>::drawCurrentLayout(const std::string &filename)
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
    for (IndexType cellIdx = 0; cellIdx < _db.numCells(); ++cellIdx)
    {
        const auto &cell = _db.cell(cellIdx);
        Box<LocType> cellBox = cell.cellBBox();
        LocType xLo = static_cast<LocType>(_pl(plIdx(cellIdx, Orient2DType::HORIZONTAL)) / _scale);
        LocType yLo = static_cast<LocType>(_pl(plIdx(cellIdx, Orient2DType::VERTICAL)) / _scale);
        Box<LocType> cellShape = Box<LocType>(xLo, yLo, xLo + cellBox.xLen(), yLo + cellBox.yLen());
        wg->writeRectangle(cellShape, cellIdx, 0);
    }
    LocType boundaryXLo = static_cast<LocType>(_boundary.xLo() / _scale);
    LocType boundaryXHi = static_cast<LocType>(_boundary.xHi() / _scale);
    LocType boundaryYLo = static_cast<LocType>(_boundary.yLo() / _scale);
    LocType boundaryYHi = static_cast<LocType>(_boundary.yHi() / _scale);
    Box<LocType> scaleBoundary(boundaryXLo, boundaryYLo, boundaryXHi, boundaryYHi);
    wg->writeRectangle(scaleBoundary, 666, 0);
    // END
    wg->writeCellEnd();
    wg->endLib();
    DBG("Database::%s: debug cell block saved in %s \n", __FUNCTION__, filename.c_str());
}
#endif
#endif
PROJECT_NAMESPACE_END


PROJECT_NAMESPACE_BEGIN

// declare the default settings for linking
template class NlpGPlacerBase<nlp::nlp_default_settings>;
template class NlpGPlacerFirstOrder<nlp::nlp_default_settings>;

PROJECT_NAMESPACE_END
