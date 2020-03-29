/**
 * @file NlpWnconj.h
 * @brief The global placement solver with non-linear optimization + wnlib conjugated gradient
 * @author Keren Zhu
 * @date 10/12/2019
 */

#ifndef IDEAPLACE_NLP_WNCONJ_H_
#define IDEAPLACE_NLP_WNCONJ_H_

#include "db/Database.h"
#include "wnconj.h"
#include <math.h>
#include <numeric>
#include <functional>
#include "pinassign/VirtualPinAssigner.h"
#include "different.h"

PROJECT_NAMESPACE_BEGIN

/// @class IDEAPLACE::NlpWnconj
/// @brief the non-linear optimization for global placement with wnlib conjugated gradient solver
class NlpWnconj
{
    typedef RealType nlp_coordinate_type;
    typedef RealType nlp_numerical_type;
    typedef diff::LseHpwlDifferentiable<nlp_numerical_type, nlp_coordinate_type> nlp_hpwl_type;
    typedef diff::CellPairOverlapPenaltyDifferentiable<nlp_numerical_type, nlp_coordinate_type> nlp_ovl_type;
    typedef diff::CellOutOfBoundaryPenaltyDifferentiable<nlp_numerical_type, nlp_coordinate_type> nlp_oob_type;
    typedef diff::AsymmetryDifferentiable<nlp_numerical_type, nlp_coordinate_type> nlp_asym_type;
    typedef diff::CosineDatapathDifferentiable<nlp_numerical_type, nlp_coordinate_type> nlp_cos_type;
    struct OpIdxType
    {
        OpIdxType(IndexType idx_, diff::OpEnumType type_) : idx(idx_), type(type_) {}
        IndexType idx;
        diff::OpEnumType type;
    };
    public:
        /// @brief default constructor
        explicit NlpWnconj(Database &db) : _db(db) , _assigner(db) {}
        /// @brief run the NLP placer
        /// @return if successful
        bool solve();
        /// @brief set tough model
        void setToughMode(bool set) { _toughModel = set; }
    private:
        /*------------------------------*/ 
        /* Flow                         */
        /*------------------------------*/ 
        /// @brief initializing the variables
        /// @return if successful
        bool initVars();
        /// @brief nonlinear optimization kernel
        /// @return if successful
        bool nlpKernel();
        /// @brief write out the solution to the database
        /// @return if successful
        bool writeOut();
        /// @brief clean up the pointers
        /// @return if successful
        bool cleanup();
        /*------------------------------*/ 
        /* Supporting functions         */
        /*------------------------------*/ 
        /// @brief calculate the total overlap area
        /// @return the total overlap area
        RealType totalOvlArea(RealType *values);
        /// @brief get total asymmetric distance
        /// @return the total asymmetric distance
        RealType totalAsymDist() { return 0; }
        /// @brief evalute the current solution. Specifically, calculating _curOvlRatio, _curOOBRatio, _curAsymDist
        void evaluteSolution();
    public:
        /*------------------------------*/ 
        /* For WNLIB interface          */
        /*------------------------------*/ 
        /// @brief the objective function
        /// @param the pointer to the current variables
        /// @return the evaluated objective function
        RealType objFunc(RealType *values);
        /// @brief the gradient function
        /// @param first: the pointer to the gradients
        /// @param second: the pointer to the values
        void gradFunc(RealType *grad, RealType *values);
    private:
        /*------------------------------*/ 
        /* Update penalty multiplier    */
        /*------------------------------*/ 
        /// @brief calculate the step size
        RealType stepSize();
        /// @brief update penalty terms
        /// @return true if the global routing should be terminated
        bool updateMultipliers();
        bool updateMultipliers2();

        void alignSym();

        /* Pin assignment */
        void assignPin()
        {
            if (!_db.parameters().ifUsePinAssignment()) { return; }
            auto cellLocQueryFunc = [&] (IndexType cellIdx)
            {
                RealType x = _solutionVect[cellIdx * 2];
                RealType y = _solutionVect[cellIdx * 2 + 1];
                LocType xLoc = ::klib::autoRound<LocType>(x / _scale);
                LocType yLoc = ::klib::autoRound<LocType>(y / _scale);
                return XY<LocType>(xLoc, yLoc);
            };
            RealType xLo = REAL_TYPE_MAX; RealType xHi = REAL_TYPE_MIN; RealType yLo = REAL_TYPE_MAX; RealType yHi = REAL_TYPE_MIN;
            for (IndexType cellIdx = 0; cellIdx < _db.numCells(); ++cellIdx)
            {
                if (_solutionVect[cellIdx *2] < xLo)
                {
                    xLo = _solutionVect[cellIdx * 2];
                }
                if (_solutionVect[cellIdx *2 + 1] < yLo)
                {
                    yLo = _solutionVect[cellIdx *2 + 1];
                }
                if (_solutionVect[cellIdx *2] > xHi)
                {
                    xHi = _solutionVect[cellIdx * 2];
                }
                if (_solutionVect[cellIdx *2 + 1] < yHi)
                {
                    yHi = _solutionVect[cellIdx *2 + 1];
                }
            }
            xLo = std::min(xLo, _boundary.xLo());
            xHi = std::max(xHi, _boundary.xHi());
            yLo = std::min(yLo, _boundary.xLo());
            yHi = std::max(yHi, _boundary.yHi());

            LocType xLoBox = ::klib::autoRound<LocType>(xLo / _scale);
            LocType yLoBox = ::klib::autoRound<LocType>(yLo / _scale);
            LocType xHiBox = ::klib::autoRound<LocType>(xHi / _scale);
            LocType yHiBox = ::klib::autoRound<LocType>(yHi / _scale);
            _assigner.reconfigureVirtualPinLocations(Box<LocType>(xLoBox, yLoBox, xHiBox, yHiBox));
            if ( _assigner.pinAssignment(cellLocQueryFunc))
            {
                // Update the hpwl operators
                for (IndexType netIdx = 0; netIdx < _db.numNets(); ++netIdx)
                {
                    const auto &net = _db.net(netIdx);
                    if (net.isValidVirtualPin())
                    {
                        XY<RealType> virtualPinLoc = XY<RealType>(net.virtualPinLoc().x(), net.virtualPinLoc().y());
                        virtualPinLoc *= _scale;
                        _hpwlOps[netIdx].setVirtualPin(virtualPinLoc.x(), virtualPinLoc.y());
                    }
                }
            }
        }

#ifdef DEBUG_GR
#ifdef DEBUG_DRAW
        void drawCurrentLayout(const std::string &filename, RealType * values);
#endif
#endif

        /* Init operators */
        void initOperators();

        
    private:
        Database &_db; ///< The placement engine database
        int _code = -1; ///< The wnlib status of the completed search. WN_SUCCESS WN_SUBOPTIMAL WN_UNBOUNDED
        RealType _valMin = 0; ///< The wnlib objective function value at the final solution
        RealType *_solutionVect = nullptr; ///< Loaded as beginning points and return containing the final solution
        int _len = -1; ///< The number of variables in "_solutionVect"
        RealType _totalCellArea = 0; ///< The total cell area
        RealType _curOvlRatio = 1; ///< The current overlapping ratio
        RealType _curOOBRatio = 1; ///< The current out of boundry ratio
        RealType _curAsymDist = 1; ///< The current asymmetric distance
        RealType _lambda1; ///< The coefficient for x and y overlap penalty 
        RealType _lambda2; ///< The cofficient for out of boundry penalty
        RealType _lambda3; ///< The coefficient for wire length penalty
        RealType _lambda4; ///< The coefficient for asymmetry penalty
        RealType _maxWhiteSpace;
        RealType _alpha; ///< Used in objective function
        Box<RealType> _boundary; ///< The boundary constraint for the placement
        RealType _overlapThreshold = NLP_WN_CONJ_OVERLAP_THRESHOLD; ///< Threshold for whether increase penalty for overlapping penalty
        RealType _oobThreshold = NLP_WN_CONJ_OOB_THRESHOLD; ///< The threshold for wehther increasing the penalty for out of boundry
        RealType _asymThreshold = NLP_WN_CONJ_ASYM_THRESHOLD; ///< The threshold for whether increasing the penalty for asymmetry
        RealType _scale = 0.01;
        RealType _fOverlap = 0;
        RealType _fOOB = 0;
        RealType _fHpwl = 0;
        RealType _fAsym = 0;
        RealType _epsilon = NLP_WN_CONJ_EPISLON;///< For updating penalty multipliers
        IndexType _iter = 0; ///< Current iteration
        RealType _tao = 0.5; ///< The exponential decay factor for step size
        IndexType _maxIter = NLP_WN_CONJ_DEFAULT_MAX_ITER; ///< The maximum iterations
        IndexType _innerIter = 0; ///< The iteration in the inner non-linear optimization problem
        RealType _defaultSymAxis = 0;
        bool _toughModel = false; ///< Whether try hard to be feasible

        VirtualPinAssigner _assigner;
        /* NLP differentiable operators */
        std::vector<OpIdxType> _ops; ///< The operator index and enum type struct. Record the types and indices
        std::vector<nlp_hpwl_type> _hpwlOps; ///< The HPWL cost 
        std::vector<nlp_ovl_type> _ovlOps; ///< The cell pair overlapping penalty operators
        std::vector<nlp_oob_type> _oobOps; ///< The cell out of boundary penalty operators 
        std::vector<nlp_asym_type> _asymOps; ///< The asymmetric penalty operators
        std::vector<nlp_cos_type> _cosOps;

};

inline RealType NlpWnconj::objFunc(RealType *values)
{
    auto getVarFunc = [&] (IndexType cellIdx, Orient2DType orient)
    {
        if (orient == Orient2DType::HORIZONTAL)
        {
            return values[2 * cellIdx];
        }
        else if (orient == Orient2DType::VERTICAL)
        {
            return values[2 * cellIdx + 1];
        }
        else
        {
#ifdef MULTI_SYM_GROUP
            return values[2 * _db.numCells() + cellIdx];
#else
            return _defaultSymAxis;
#endif
        }
    };

    // Initial the objective to be 0 and add the non-zero to it
    nlp_numerical_type result = 0;
    _fOverlap = 0;
    _fOOB = 0;
    _fHpwl = 0;
    _fAsym = 0;
    std::vector<nlp_numerical_type> ovl, oob, hpwl, asym, cos;
    ovl.resize(_ovlOps.size());
    oob.resize(_oobOps.size());
    hpwl.resize(_hpwlOps.size());
    asym.resize(_asymOps.size());
    cos.resize(_cosOps.size());


    #pragma omp parallel for schedule(static)
    for (IndexType idx = 0; idx < _ops.size(); ++idx)
    {
        const auto &opIdx = _ops[idx];
        const auto &type = opIdx.type;
        if (type == diff::OpEnumType::ovl)
        {
            ovl[opIdx.idx] = diff::placement_differentiable_traits<nlp_ovl_type>::evaluate(_ovlOps[opIdx.idx], getVarFunc);
        }
        else if (type == diff::OpEnumType::oob)
        {
            oob[opIdx.idx] = diff::placement_differentiable_traits<nlp_oob_type>::evaluate(_oobOps[opIdx.idx], getVarFunc);
        }
        else if (type == diff::OpEnumType::hpwl)
        {
            hpwl[opIdx.idx] = diff::placement_differentiable_traits<nlp_hpwl_type>::evaluate(_hpwlOps[opIdx.idx], getVarFunc);
        }
        else if (type == diff::OpEnumType::asym)
        {
            asym[opIdx.idx] = diff::placement_differentiable_traits<nlp_asym_type>::evaluate(_asymOps[opIdx.idx], getVarFunc);
        }
        else
        {
            cos[opIdx.idx] = diff::placement_differentiable_traits<nlp_cos_type>::evaluate(_cosOps[opIdx.idx], getVarFunc);
        }
    }

    _fOverlap = std::accumulate(ovl.begin(), ovl.end(), 0.0);
    _fOOB = std::accumulate(oob.begin(), oob.end(), 0.0);
    _fHpwl = std::accumulate(hpwl.begin(), hpwl.end(), 0.0);
    _fAsym = std::accumulate(asym.begin(), asym.end(), 0.0);

    result = _fOverlap + _fOOB + _fHpwl + _fAsym;


    return result;
}

inline void NlpWnconj::gradFunc(RealType *grad, RealType *values)
{
    auto getVarFunc = [&] (IndexType cellIdx, Orient2DType orient)
    {
        if (orient == Orient2DType::HORIZONTAL)
        {
            return values[2 * cellIdx];
        }
        else if (orient == Orient2DType::VERTICAL)
        {
            return values[2 * cellIdx + 1];
        }
        else
        {
#ifdef MULTI_SYM_GROUP
            return values[2 * _db.numCells() + cellIdx];
#else
            return _defaultSymAxis;
#endif
        }
    };
    
    struct Partial
    {
        Partial() { value = 0.0; idx = 0; }
        Partial(nlp_numerical_type val, IndexType i) : value(val), idx(i) {}
        nlp_numerical_type value;
        IndexType idx;
    };


    struct OpGrad
    {
        OpGrad(IndexType numCells) { _ds.resize(0); _numCells = numCells; }
        void addPartial(nlp_numerical_type value, IndexType cellIdx, Orient2DType orient)
        {
            if (orient == Orient2DType::HORIZONTAL)
            {
                _ds.emplace_back(Partial(value, 2 * cellIdx));
            }
            else if (orient == Orient2DType::VERTICAL)
            {
                _ds.emplace_back(Partial(value, 2 * cellIdx + 1));
            }
#ifdef MULTI_SYM_GROUP
            else
            {
                _ds.emplace_back(Partial(value, 2 * _numCells + cellIdx));
            }
#endif
        }
        std::vector<Partial> _ds;
        IndexType _numCells = 0;
    };

    std::vector<OpGrad> _grads(_ops.size(), _db.numCells());

    _innerIter++;

    #pragma omp parallel for schedule(static)
    for (IndexType idx = 0; idx < _ops.size(); ++idx)
    {
        auto f = 
            std::bind(&OpGrad::addPartial, &(_grads[idx]), 
                    std::placeholders::_1, std::placeholders::_2, std::placeholders::_3);
        const auto &opIdx = _ops[idx];
        const auto &type = opIdx.type;
        if (type == diff::OpEnumType::ovl)
        {
            diff::placement_differentiable_traits<nlp_ovl_type>::accumlateGradient(_ovlOps[opIdx.idx], getVarFunc, f);
        }
        else if (type == diff::OpEnumType::oob)
        {
            diff::placement_differentiable_traits<nlp_oob_type>::accumlateGradient(_oobOps[opIdx.idx], getVarFunc, f);
        }
        else if (type == diff::OpEnumType::hpwl)
        {
            diff::placement_differentiable_traits<nlp_hpwl_type>::accumlateGradient(_hpwlOps[opIdx.idx], getVarFunc, f);
        }
        else if (type == diff::OpEnumType::asym)
        {
            diff::placement_differentiable_traits<nlp_asym_type>::accumlateGradient(_asymOps[opIdx.idx], getVarFunc, f);
        }
        else
        {
            diff::placement_differentiable_traits<nlp_cos_type>::accumlateGradient(_cosOps[opIdx.idx], getVarFunc, f);
        }
    }

    for (const auto & part : _grads)
    {
        for (const auto &partial : part._ds)
        {
            grad[partial.idx] += partial.value;
        }
    }

#ifdef DEBUG_GR
#ifdef DEBUG_DRAW
    //std::stringstream ss;
    //ss<< "./debug/gr_iter_" << _iter  << "_"<< _innerIter << ".gds";
    //drawCurrentLayout(ss.str(), values);
#endif
#endif
    
}

inline RealType NlpWnconj::totalOvlArea(RealType *values)
{
    // naive implementation
    // TODO: better implementation
    RealType area = 0;
    for (IndexType cellIdxI = 0; cellIdxI < _db.numCells(); ++cellIdxI)
    {
        const auto &bboxI = _db.cell(cellIdxI).cellBBox();
        for (IndexType cellIdxJ = cellIdxI + 1; cellIdxJ < _db.numCells(); ++cellIdxJ)
        {
            const auto &bboxJ = _db.cell(cellIdxJ).cellBBox();
            // Values arrangement x0, y0, x1, y1...
            RealType xLoI = values[2 * cellIdxI]; RealType xHiI = xLoI + bboxI.xLen() * _scale;
            RealType xLoJ = values[2 * cellIdxJ]; RealType xHiJ = xLoJ + bboxJ.xLen() * _scale;
            RealType yLoI = values[2 * cellIdxI + 1]; RealType yHiI = yLoI + bboxI.yLen() * _scale;
            RealType yLoJ = values[2 * cellIdxJ + 1]; RealType yHiJ = yLoJ + bboxJ.yLen() * _scale;
            // max (min(xHiI - xLoJ, xHiJ - xLoI), 0)
            RealType var1X = xHiI - xLoJ;
            RealType var2X = xHiJ - xLoI;
            RealType overlapX = std::min(var1X, var2X);
            overlapX = std::max(overlapX, 0.0);
            // max (min(yHiI - yLoJ, yHiJ - yLoI), 0)
            RealType var1Y = yHiI - yLoJ;
            RealType var2Y = yHiJ - yLoI;
            RealType overlapY = std::min(var1Y, var2Y);
            overlapY = std::max(overlapY, 0.0);
            if (overlapX * overlapY >= 0.000001)
            area +=  overlapX * overlapY;
        }
    }
    return area;
}

/// @brief evalute the current solution. Specifically, calculating _curOvlRatio, _curOOBRatio, _curAsymDist
inline void NlpWnconj::evaluteSolution()
{
    RealType totalOvlArea = this->totalOvlArea(_solutionVect);
    _curOvlRatio = totalOvlArea / (_totalCellArea );//* _scale * _scale);
    // out of boundary
    RealType totalOOBarea = 0.0;
    for (IndexType cellIdx = 0; cellIdx < _db.numCells(); ++cellIdx)
    {
        const auto &bbox = _db.cell(cellIdx).cellBBox();
        RealType xLo = _solutionVect[2 * cellIdx]; RealType xHi = xLo + bbox.xLen() * _scale;
        RealType yLo = _solutionVect[2 * cellIdx + 1]; RealType yHi = yLo + bbox.yLen() * _scale;
        /*
        RealType varX = 0;
        varX += std::max(-xLo + _boundary.xLo(), 0.0);
        varX += std::max(xHi - _boundary.xHi(), 0.0);
        RealType varY = 0;
        varY += std::max(_boundary.yLo() - yLo, 0.0);
        varY += std::max(yHi - _boundary.yHi(), 0.0);
        totalOOBarea +=  varX * varY;
        */
        RealType varX = 0;
        varX += std::max(-xLo + _boundary.xLo(), 0.0);
        varX += std::max(xHi - _boundary.xHi(), 0.0);
        RealType varY = 0;
        varY += std::max(_boundary.yLo() - yLo, 0.0);
        varY += std::max(yHi - _boundary.yHi(), 0.0);
        totalOOBarea +=  (varX - _boundary.xLen()) * (varY - _boundary.yLen());
        //DBG("total oob area %f \n", totalOOBarea);
    }
    _curOOBRatio = totalOOBarea / ( _totalCellArea);
    // Asymmetry
    _curAsymDist = 0;
    for (IndexType symGrpIdx = 0; symGrpIdx < _db.numSymGroups(); ++symGrpIdx)
    {
        const auto &symGrp = _db.symGroup(symGrpIdx);
#ifdef MULTI_SYM_GROUP
        RealType symAxis = _solutionVect[2 * _db.numCells() + symGrpIdx];
#else
        RealType symAxis = _defaultSymAxis;
#endif
        for (IndexType symPairIdx = 0; symPairIdx < symGrp.numSymPairs(); ++symPairIdx)
        {
            const auto &symPair = symGrp.symPair(symPairIdx);
            IndexType cell1 = symPair.firstCell();
            IndexType cell2 = symPair.secondCell();
            RealType cellWidth = _db.cell(cell1).cellBBox().xLen() * _scale;
            _curAsymDist += std::abs(_solutionVect[2 * cell1 + 1] - _solutionVect[2 * cell2 + 1]);
            _curAsymDist += std::abs(_solutionVect[2 * cell1] + _solutionVect[2 * cell2] + cellWidth
                    - 2 * symAxis);
        }
        for (IndexType selfSymIdx = 0; selfSymIdx < symGrp.numSelfSyms(); ++selfSymIdx)
        {
            IndexType cellIdx = symGrp.selfSym(selfSymIdx);
            RealType cellWidth = _db.cell(cellIdx).cellBBox().xLen() * _scale;
            _curAsymDist += std::abs(_solutionVect[2 * cellIdx] +  cellWidth / 2 - symAxis);
        }
    }
}
PROJECT_NAMESPACE_END

#endif
