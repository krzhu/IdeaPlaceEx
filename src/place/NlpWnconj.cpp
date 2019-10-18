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

    // Clean up the variables
    if (!this->cleanup())
    {
        ERR("NLP wnlib conj::%s cleanning up the variables failed! \n", __FUNCTION__);
        Assert(false);
        return false;
    }
    return true;
}

bool NlpWnconj::initVars()
{
    // The number of nlp problem variables
    _len = 2 * static_cast<int>(_db.numCells()) + static_cast<int>(_db.numSymGroups());
    // Allocate the soluction vec
    _solutionVect = (double*)malloc(sizeof(double) * _len);

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
        double value = static_cast<double>(idx) * 20.0 / static_cast<double>(_len);
        _solutionVect[idx] = value;
    }

    // The penalty coefficients
    _lambda1 = LAMBDA_1Init;
    _lambda2 = LAMBDA_2Init;
    _lambda3 = LAMBDA_3Init;
    _lambda4 = LAMBDA_4Init;

    // Other static variables
    _alpha = NLP_WN_CONJ_ALPHA;

    // Total cell area
    _totalCellArea = static_cast<double>(_db.calculateTotalCellArea());

    // Placement Boundary
    if (_db.parameters().isBoundaryConstraintSet())
    {
        // If the constraint is set in the database, follow it.
        const auto &bb = _db.parameters().boundaryConstraint();
        _boundary.setXLo(static_cast<double>(bb.xLo()));
        _boundary.setYLo(static_cast<double>(bb.yLo()));
        _boundary.setXHi(static_cast<double>(bb.xHi()));
        _boundary.setYHi(static_cast<double>(bb.yHi()));
    }
    else
    {
        // If the constraint is not set, calculate a rough boundry with 1 aspect ratio
        double aspectRatio = 1;
        double maxWhiteSpace = 0.2;
        double xLo = 0; double yLo = 0; 
        double tolerentArea = _totalCellArea * (1 + maxWhiteSpace);
        double xHi = std::sqrt(tolerentArea * aspectRatio);
        double yHi = tolerentArea / xHi;
        _boundary.set(xLo, yLo, xHi, yHi);
    }

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
    double test [] = {1,2 };
    objFuncWrapper(test);
    return true;
    // Iteratively solve NLP untial total overlap is less than the threshold
    size_t numIter = 0; // Current number of iterations
    size_t maxIter = 10; // The maximum number of iterations

    // Iteratively solving NLP
    while (numIter < maxIter)
    {
        wn_conj_gradient_method(&_code, &_valMin, _solutionVect, _len, objFuncWrapper, gradFuncWrapper, 1000);
        bool breakFlag = true; // Whether all criterials are met
        // Evaluate the current solution, calculate _curOvlRatio, _curOOBratio and _curAsymDist
        this->evaluteSolution();
        // Increase penalty for overlapping if the overlapping ratio is larger than the threshold
        if (_curOvlRatio > _overlapThreshold)
        {
            _lambda1 = _lambda1 * 2; // Double
            breakFlag = false;
        }
        if (_curOOBRatio > _oobThreshold)
        {
            _lambda2 = _lambda2 * 2; // Double
            breakFlag = false;
        }
        if (_curAsymDist > _asymThreshold)
        {
            _lambda4 = _lambda4 * 2; // Double
            breakFlag = false;
        }
        // Break if all criterials are met
        if (_code == 0 && breakFlag)
        {
            break;
        }
    }
    return true;
}

bool NlpWnconj::cleanup()
{
    free(_solutionVect);
    _solutionVect = nullptr;
    return true;
}

PROJECT_NAMESPACE_END
