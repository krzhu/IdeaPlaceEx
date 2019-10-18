/**
 * @file NlpWnconj.h
 * @brief The "global placement" solver with non-linear optimization + wnlib conjugated gradient
 * @author Keren Zhu
 * @date 10/12/2019
 */

#ifndef IDEAPLACE_NLP_WNCONJ_H_
#define IDEAPLACE_NLP_WNCONJ_H_

#include "db/Database.h"
#include "wnconj.h"
#include <math.h>

PROJECT_NAMESPACE_BEGIN

namespace NlpWnconjDetails
{
    /// @brief calculate the log sum exp, used to smooth min max function
    /// @return the calculated log sum exp
    inline double logSumExp(double var1, double var2, double alpha)
    {
        return alpha * log(exp(var1 / alpha) + exp(var2 / alpha));
    }

    /// @brief Some smooth function from Biying I cannot understand
    /// @return the smoothed result
    inline double logSumExp0(double var, double alpha)
    {
        // modify it to handle overflow issue
        // avg > 709.8 => overflow; avg < -709.8 => underflow
        if (var / alpha < 709.8)
        {
            return alpha * log(exp(var / alpha) + 1);
        }
        else
        {
            // -var + 32
            return alpha * log(exp(32 / alpha) + exp((-var + 32)/alpha)) + var - 32;
        }
    }
}

/// @class IDEAPLACE::NlpWnconj
/// @brief the non-linear optimization for global placement with wnlib conjugated gradient solver
class NlpWnconj
{
    public:
        /// @brief default constructor
        explicit NlpWnconj(Database &db) : _db(db) {}
        /// @brief run the NLP placer
        /// @return if successful
        bool solve();
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
        /// @brief clean up the pointers
        /// @return if successful
        bool cleanup();
        /*------------------------------*/ 
        /* Supporting functions         */
        /*------------------------------*/ 
        /// @brief calculate the total overlap area
        /// @return the total overlap area
        double totalOvlArea(double *values) { return 0; }
        /// @brief get total asymmetric distance
        /// @return the total asymmetric distance
        double totalAsymDist() {}
        /// @brief evalute the current solution. Specifically, calculating _curOvlRatio, _curOOBRatio, _curAsymDist
        void evaluteSolution() {}
        /*------------------------------*/ 
        /* Math                         */
        /*------------------------------*/ 
        /// @brief calculate the smoothed HPWL of each net
        /// @param first: the array of the current values
        /// @param second: the net index
        /// @param third: alpha parameter in the smooth function
        double logSumExpHPWL(double *values, IndexType netIdx, double alpha)
        {
            double xmax = 0, xmin = 0, ymax = 0, ymin = 0;
            const auto &net = _db.net(netIdx);
            for (IndexType idx = 0; idx < net.numPinIdx(); ++idx)
            {
                // Get the pin location referenced to the cell
                IndexType pinIdx = net.pinIdx(idx);
                const auto &pin = _db.pin(pinIdx);
                IndexType cellIdx = pin.cellIdx();
                // Get the cell location from the input arguments
                XY<LocType> cellLoc = XY<LocType>(values[cellIdx * 2], values[cellIdx * 2 + 1]);
                XY<LocType> pinLoc = cellLoc + pin.midLoc();
                xmax += exp(pinLoc.x() / alpha);
                xmin += exp(-pinLoc.x() / alpha);
                ymax += exp(pinLoc.y() / alpha);
                ymin += exp(-pinLoc.y() / alpha);
            }
            xmax = log(xmax);
            xmin = log(xmin);
            ymax = log(ymax);
            ymin = log(ymin);
            return alpha * (xmax + xmin + ymax + ymin);
        }
    public:
        /*------------------------------*/ 
        /* For WNLIB interface          */
        /*------------------------------*/ 
        /// @brief the objective function
        /// @param the pointer to the current variables
        /// @return the evaluated objective function
        double objFunc(double *values);
        /// @brief the gradient function
        /// @param first: the pointer to the gradients
        /// @param second: the pointer to the values
        void gradFunc(double *grad, double *values);
        
    private:
        Database &_db; ///< The placement engine database
        int _code = -1; ///< The wnlib status of the completed search. WN_SUCCESS WN_SUBOPTIMAL WN_UNBOUNDED
        double _valMin = 0; ///< The wnlib objective function value at the final solution
        double *_solutionVect = nullptr; ///< Loaded as beginning points and return containing the final solution
        int _len = -1; ///< The number of variables in "_solutionVect"
        double _totalCellArea = 0; ///< The total cell area
        double _curOvlRatio = 1; ///< The current overlapping ratio
        double _curOOBRatio = 1; ///< The current out of boundry ratio
        double _curAsymDist = 1; ///< The current asymmetric distance
        double _lambda1; ///< The coefficient for x and y overlap penalty 
        double _lambda2; ///< The cofficient for out of boundry penalty
        double _lambda3; ///< The coefficient for wire length penalty
        double _lambda4; ///< The coefficient for asymmetry penalty
        double _alpha; ///< Used in objective function
        Box<double> _boundary; ///< The boundary constraint for the placement
        double _overlapThreshold = NLP_WN_CONJ_OVERLAP_THRESHOLD; ///< Threshold for whether increase penalty for overlapping penalty
        double _oobThreshold = NLP_WN_CONJ_OOB_THRESHOLD; ///< The threshold for wehther increasing the penalty for out of boundry
        double _asymThreshold = NLP_WN_CONJ_ASYM_THRESHOLD; ///< The threshold for whether increasing the penalty for asymmetry
};

inline double NlpWnconj::objFunc(double *values)
{
    // Initial the objective to be 0 and add the non-zero to it
    double result = 0;
    // log-sum-exp model for overlap penalty
    IndexType numCells = _db.numCells();
    // Calculate the cell-wise overlapping area penalty
    // TODO: improve the codes below from O(n^2) to O(nlogn)
    for (IndexType cellIdxI = 0; cellIdxI < numCells; ++cellIdxI)
    {
        const auto &bboxI = _db.cell(cellIdxI).cellBBox();
        for (IndexType cellIdxJ = cellIdxI + 1; cellIdxJ < numCells; ++cellIdxJ)
        {
            const auto &bboxJ = _db.cell(cellIdxJ).cellBBox();
            // Values arrangement x0, y0, x1, y1...
            double xLoI = values[2 * cellIdxI]; double xHiI = xLoI + bboxI.xLen();
            double xLoJ = values[2 * cellIdxJ]; double xHiJ = xLoJ + bboxJ.xLen();
            double yLoI = values[2 * cellIdxI + 1]; double yHiI = yLoI + bboxI.yLen();
            double yLoJ = values[2 * cellIdxJ + 1]; double yHiJ = yLoJ + bboxJ.yLen();
            // max (min(xHiI - xLoJ, xHiJ - xLoI), 0)
            double var1X = xHiI - xLoJ;
            double var2X = xHiJ - xLoI;
            double overlapX = NlpWnconjDetails::logSumExp(var1X, var2X, -_alpha);
            overlapX = NlpWnconjDetails::logSumExp0(overlapX, _alpha);
            // max (min(yHiI - yLoJ, yHiJ - yLoI), 0)
            double var1Y = yHiI - yLoJ;
            double var2Y = yHiJ - yLoI;
            double overlapY = NlpWnconjDetails::logSumExp(var1Y, var2Y, -_alpha);
            overlapY = NlpWnconjDetails::logSumExp0(overlapY, _alpha);
            // Add to the objective
            result += _lambda1 * overlapX * overlapY;
        }
    }

    // Out of boundary penalty
    for (IndexType cellIdx = 0; cellIdx < numCells; ++cellIdx)
    {
        const auto &bbox = _db.cell(cellIdx).cellBBox();
        double xLo = values[2 * cellIdx]; double xHi = xLo + bbox.xLen();
        double yLo = values[2 * cellIdx + 1]; double yHi = yLo + bbox.yLen();
        // max (-xLo +boundary_xLo, 0), max(xHi - boundary_xHi, 0)
        double obXLo = NlpWnconjDetails::logSumExp0(_boundary.xLo() - xLo, _alpha);
        double obXHi = NlpWnconjDetails::logSumExp0(xHi - _boundary.xHi(), _alpha);
        // also y
        double obYLo = NlpWnconjDetails::logSumExp0(_boundary.yLo() - yLo, _alpha);
        double obYHi = NlpWnconjDetails::logSumExp0(yHi - _boundary.yLo(), _alpha);
        // Add to the objective
        result += _lambda2 * (obXLo + obXHi + obYLo + obYHi);
    }

    // Wire length penalty
    for (IndexType netIdx = 0; netIdx < _db.numNets(); ++netIdx)
    {
        double smoothHPWL = this->logSumExpHPWL(values, netIdx, _alpha);
        result += _lambda3 * _db.net(netIdx).weight() * smoothHPWL;
    }
    
    // ASYMMETRY
    // TODO
    
    return result;
}

inline void NlpWnconj::gradFunc(double *grad, double *values)
{
}

PROJECT_NAMESPACE_END

#endif
