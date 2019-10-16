#include "NlpWnconj.h"
#include <cstdlib>

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

    return true;
}

bool NlpWnconj::nlpKernel()
{
    // Iteratively solve NLP untial total overlap is less than the threshold
    double curOvlRatio = std::numeric_limits<double>::max(); // Current overlap ratio
    size_t numIter = 0; // Current number of iterations
    size_t maxIter = 10; // The maximum number of iterations

    // Iteratively solving NLP
    while (numIter < maxIter)
    {
        wn_conj_gradient_method(&_code, &_valMin, _solutionVect, _len, this->objFunc, this->gradFunc, 1000);
        double ovlArea = this->totalOvlArea();
        curOvlRatio = ovlArea; // TODO
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
