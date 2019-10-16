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

PROJECT_NAMESPACE_BEGIN

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
        /// @brief the objective function
        /// @param the pointer to the current variables
        /// @return the evaluated objective function
        static double objFunc(double *values) { return 0; }
        /// @brief the gradient function
        /// @param first: the pointer to the gradients
        /// @param second: the pointer to the values
        static void gradFunc(double *grad, double *values) {}
        /// @brief calculate the total overlap area
        /// @return the total overlap area
        double totalOvlArea() { return 0; }
    private:
        Database &_db; ///< The placement engine database
        int _code = -1; ///< The wnlib status of the completed search. WN_SUCCESS WN_SUBOPTIMAL WN_UNBOUNDED
        double _valMin = 0; ///< The wnlib objective function value at the final solution
        double *_solutionVect = nullptr; ///< Loaded as beginning points and return containing the final solution
        int _len = -1; ///< The number of variables in "_solutionVect"

};

PROJECT_NAMESPACE_END

#endif
