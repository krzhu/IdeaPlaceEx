#include "NlpWnconj.h"
#include <cstdlib>
#include <functional>

PROJECT_NAMESPACE_BEGIN

bool NlpWnconj::solve()
{
    // Init the variables used
    if (!this->initVars())
    {
        ERR("NLP wnlib conj::%s initialize variables failed! \n", __FUNCTION__);
        Assert(false);
        return false;
    }
    // Run the NLP kernel
    if (!this->nlpKernel())
    {
        ERR("NLP wnlib conj::%s Optimization kernel failed! \n", __FUNCTION__);
        Assert(false);
        return false;
    }
    // Dump the results into the database
    if (!this->writeOut())
    {
        ERR("NLP wnlib conj::%s Outputing results failed! \n", __FUNCTION__);
        Assert(false);
        return false;
    }

    // Clean up the variables
    if (!this->cleanup())
    {
        ERR("NLP wnlib conj::%s cleanning up the variables failed! \n", __FUNCTION__);
        Assert(false);
        return false;
    }
    INF("Nlp:: global placement finished \n");
    return true;
}

bool NlpWnconj::writeOut()
{
    // find the min value
    RealType minX =1e10; 
    RealType minY = 1e10;
    for (IndexType cellIdx = 0; cellIdx < _db.numCells(); ++cellIdx)
    {
        if (_solutionVect[cellIdx *2] < minX)
        {
            minX = _solutionVect[cellIdx * 2];
        }
        if (_solutionVect[cellIdx *2 + 1] < minY)
        {
            minY = _solutionVect[cellIdx *2 + 1];
        }
    }
    // Dump the cell locations to database
    for (IndexType cellIdx = 0; cellIdx < _db.numCells(); ++cellIdx)
    {
        auto & cell = _db.cell(cellIdx);
        LocType xLo = ::klib::autoRound<LocType>((_solutionVect[cellIdx * 2] - minX) / _scale);
        LocType yLo = ::klib::autoRound<LocType>((_solutionVect[cellIdx * 2 + 1] - minY) / _scale);
        _db.cell(cellIdx).setXLoc(xLo + cell.cellBBox().xLo());
        _db.cell(cellIdx).setYLoc(yLo + cell.cellBBox().yLo());
    }
    return true;
}

RealType NlpWnconj::stepSize()
{
    RealType eps = _epsilon * ( exp( _iter * _tao));
    RealType obj = this->objFunc(_solutionVect); // also calculate fOOB etc.
    RealType violate = _fOverlap + _fOOB + _fAsym;
    return eps * ( obj - 0.0) / violate; // 0.0 is better to be replaced by lower bound baseline
}

bool NlpWnconj::updateMultipliers()
{
    auto mu = stepSize();
#ifdef DEBUG_GR
    DBG("\n\niter %d mu %f lambda1 %f lambda2 %f lambda4 %f\n", _iter, mu,  _lambda1, _lambda2, _lambda4);
#endif
    RealType violate = _fOverlap  + _fOOB + _fAsym + _fHpwl ;
    _lambda1 = _lambda1  + mu * _fOverlap  / violate;
    _lambda2 = _lambda2 + mu * _fOOB  / violate;
    _lambda4 = _lambda4 + mu *_fAsym  / violate;
    if (_lambda1 + _lambda2 + _lambda4 > NLP_WN_MAX_PENALTY)
    {
        _lambda1 *= NLP_WN_REDUCE_PENALTY;
        _lambda2 *= NLP_WN_REDUCE_PENALTY;
        _lambda4 *= NLP_WN_REDUCE_PENALTY;
    }
#ifdef DEBUG_GR
    DBG("\n\niter %d after mu %f lambda1 %f lambda2 %f lambda4 %f\n", _iter, mu,  _lambda1, _lambda2, _lambda4);
#endif
    if (_toughModel)
    {
        _lambdaMaxOverlap += mu ;
    }
    return false;
}

bool NlpWnconj::initVars()
{
    // The number of nlp problem variables
    _len = 2 * static_cast<int>(_db.numCells()) + static_cast<int>(_db.numSymGroups());
    // Allocate the soluction vec
    _solutionVect = (double*)malloc(sizeof(double) * _len);

    // The penalty coefficients
    _lambda1 = LAMBDA_1Init;
    _lambda2 = LAMBDA_2Init;
    _lambda3 = LAMBDA_3Init;
    _lambda4 = LAMBDA_4Init;
    _lambdaMaxOverlap = 0;

    // max white space
    _maxWhiteSpace = NLP_WN_CONJ_DEFAULT_MAX_WHITE_SPACE;

    if (_toughModel)
    {
        INF("NLP global placement: trying hard mode \n");
        _lambda1 = 1;
        _lambda2 = 0;
        _lambda3 = 0;
        _lambda4 = 4;
        _lambdaMaxOverlap = LAMBDA_MAX_OVERLAP_Init;
        _maxWhiteSpace = 50;
        _maxIter = 32;
    }

    // Other static variables
    _alpha = NLP_WN_CONJ_ALPHA;

    // Total cell area
    RealType totalCellArea = static_cast<double>(_db.calculateTotalCellArea());
    _scale = sqrt(100 / (totalCellArea));
    //_totalCellArea = static_cast<double>(_db.calculateTotalCellArea()) * _scale * _scale;
    _totalCellArea = 100;

    // Placement Boundary
    if (_db.parameters().isBoundaryConstraintSet())
    {
        // If the constraint is set in the database, follow it.
        const auto &bb = _db.parameters().boundaryConstraint();
        _boundary.setXLo(static_cast<double>(bb.xLo()) * _scale);
        _boundary.setYLo(static_cast<double>(bb.yLo()) * _scale);
        _boundary.setXHi(static_cast<double>(bb.xHi()) * _scale);
        _boundary.setYHi(static_cast<double>(bb.yHi()) * _scale);
    }
    else
    {
        // If the constraint is not set, calculate a rough boundry with 1 aspect ratio
        double aspectRatio = 1;
        double xLo = 0; double yLo = 0; 
        double tolerentArea = _totalCellArea * (1 + _maxWhiteSpace);
        double xHi = std::sqrt(tolerentArea * aspectRatio);
        double yHi = tolerentArea / xHi;
        _boundary.set(xLo , yLo , xHi , yHi );
        INF("NlpWnconj::%s: automatical set boundary to be %s \n", __FUNCTION__, _boundary.toStr().c_str());
        Box<LocType> bb = Box<LocType>(static_cast<LocType>(_boundary.xLo()),
                        static_cast<LocType>(_boundary.yLo()),
                        static_cast<LocType>(_boundary.xHi()),
                        static_cast<LocType>(_boundary.yHi())
                        );
        _db.parameters().setBoundaryConstraint(bb);
        INF("NlpWnconj::initVars: add boundary constraints as calculated \n");
    }

    // Default sym axis is at the middle
    _defaultSymAxis = (_boundary.xLo() + _boundary.xHi()) / 2;


    /*
    // initialize the initial value of the solution
    // Biying mark a FIXME here
    for (IntType idx = 0; idx < _len; ++idx)
    {
        _solutionVect[idx] = static_cast<double>(std::rand() %20);
    }

    // Ensure the initial coordinate are different for symmetric pairs and shape coordinate
    // TODO: symmetric group
    */
    
    // My alternative intialization, just give the variables fixed linear value
    for (IntType idx = 0; idx < _len; ++idx)
    {
        double value = (static_cast<double>(idx) * _boundary.xHi() + _boundary.xLo()) / static_cast<double>(_len);
        _solutionVect[idx] = value;
        //_solutionVect[idx] = _defaultSymAxis;
    }
#if 1
    if (1)
    {
        for (IndexType symGrpIdx = 0; symGrpIdx < _db.numSymGroups(); ++symGrpIdx)
        {
            const auto &symGrp = _db.symGroup(symGrpIdx);
            for (IndexType symPairIdx = 0; symPairIdx < symGrp.numSymPairs(); ++ symPairIdx)
            {
                const auto &symPair = symGrp.symPair(symPairIdx);
                IndexType cell1 = symPair.firstCell();
                IndexType cell2 = symPair.secondCell();
                RealType xLo1 = _solutionVect[2 * cell1];
                RealType yLo1 = _solutionVect[2 * cell1 + 1];
                RealType xHi1 = xLo1 + _db.cell(cell1).cellBBox().xLen() * _scale;
                RealType xLo2 = _defaultSymAxis * 2 -xHi1;
                _solutionVect[2 * cell2] = xLo2;
                _solutionVect[2 * cell2 + 1] = yLo1;
            }
            for (IndexType selfSymIdx = 0; selfSymIdx < symGrp.numSelfSyms(); ++selfSymIdx)
            {
                const auto selfSym = symGrp.selfSym(selfSymIdx);
                _solutionVect[2 * selfSym] =  _defaultSymAxis / 2;
            }
        }
    }
#endif

    // Calculate tao
    // target = exp( tao * max_iter)
    _tao = log(NLP_WN_CONJ_EXP_DECAY_TARGET) / _maxIter;
    return true;
}

// The following is a trick to avoid using static member functions
// The wnlib using functional pointer to obj function and grad function, which needs to be static if naively implemented
extern NlpWnconj *nlpPtr;
double objFuncWrapper(double *vec)
{
    return nlpPtr->objFunc(vec);
}

void gradFuncWrapper(double *vec1, double *vec2)
{
    return nlpPtr->gradFunc(vec1, vec2);
}

bool NlpWnconj::nlpKernel()
{
    // Iteratively solve NLP untial total overlap is less than the threshold
    _iter = 0; // Current number of iterations


    double *initial_coord_x0s;
    initial_coord_x0s = (double*)malloc(sizeof(double) * _len);
    for (int i = 0; i < _len; ++i)
    {
        initial_coord_x0s[i] = 1.0;
    }

    // Iteratively solving NLP
    while (_iter < _maxIter)
    {
        wn_conj_gradient_method(&_code, &_valMin, _solutionVect, _len, objFuncWrapper, gradFuncWrapper, 1000);
        if (_toughModel)
        {
            for (IndexType symGrpIdx = 0; symGrpIdx < _db.numSymGroups(); ++symGrpIdx)
            {
                const auto &symGrp = _db.symGroup(symGrpIdx);
                for (IndexType symPairIdx = 0; symPairIdx < symGrp.numSymPairs(); ++ symPairIdx)
                {
                    const auto &symPair = symGrp.symPair(symPairIdx);
                    IndexType cell1 = symPair.firstCell();
                    IndexType cell2 = symPair.secondCell();
                    RealType xLo1 = _solutionVect[2 * cell1];
                    RealType yLo1 = _solutionVect[2 * cell1 + 1];
                    RealType xHi1 = xLo1 + _db.cell(cell1).cellBBox().xLen() * _scale;
                    RealType xLo2 = _defaultSymAxis * 2 -xHi1;
                    _solutionVect[2 * cell2] = xLo2;
                    _solutionVect[2 * cell2 + 1] = yLo1;
                }
                for (IndexType selfSymIdx = 0; selfSymIdx < symGrp.numSelfSyms(); ++selfSymIdx)
                {
                    const auto selfSym = symGrp.selfSym(selfSymIdx);
                    _solutionVect[2 * selfSym] =  _defaultSymAxis / 2;
                }
            }
        }
#ifdef DEBUG_GR
        //wn_conj_direction_method(&_code, &_valMin, _solutionVect, initial_coord_x0s, _len, objFuncWrapper, 1000);
        double objective = objFunc(_solutionVect);
        DBG("NlpWnconj::%s: objective %f at iter %d \n", __FUNCTION__, objective, _iter);
        for (IndexType cellIdx = 0; cellIdx < _db.numCells(); ++cellIdx)
        {
            DBG("cell %d x %f y %f \n", cellIdx, _solutionVect[cellIdx * 2], _solutionVect[cellIdx *2 + 1]);
        }
#endif
        if (updateMultipliers())
        {
            INF("Global placement terminates \n");
            break;
        }
#if 0
        if (_curOvlRatio > _overlapThreshold)
        {
            DBG("Overlap check failed \n");
            _lambda1 = _lambda1 * 2; // Double
            breakFlag = false;
        }
        if (_curOOBRatio > _oobThreshold)
        {
            DBG("OOB check failed \n");
            _lambda2 = _lambda2 * 2; // Double
            breakFlag = false;
        }
        if (0)
        //if (_curAsymDist > _asymThreshold)
        {
            DBG("asym check failed \n");
            _lambda4 = _lambda4 * 2; // Double
            breakFlag = false;
        }
        // Break if all criterials are met
        if (_code == 0 && breakFlag)
        {
            break;
        }
#endif
        _iter++;
    }
    return true;
}

bool NlpWnconj::cleanup()
{
    delete(_solutionVect);
    _solutionVect = nullptr;
    return true;
}

PROJECT_NAMESPACE_END
