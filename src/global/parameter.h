/**
 * @file parameter.h
 * @brief Define some hyperparameters
 * @author Keren Zhu
 * @date 09/30/2019
 */

#ifndef IDEAPLACE_PARAMETER_H_
#define IDEAPLACE_PARAMETER_H_

#include "type.h"

PROJECT_NAMESPACE_BEGIN

/* NLP wn conj */
constexpr double LAMBDA_1Init = 1; ///< The initial value for x/y overlapping penalty
constexpr double LAMBDA_2Init = 1; ///< The initial value for out of boundry penalty
constexpr double LAMBDA_3Init = 16; ///< The initial value for wirelength penalty
constexpr double LAMBDA_4Init = 1; ///< The initial value for asymmetric penalty
constexpr double NLP_WN_CONJ_OVERLAP_THRESHOLD = 0.05; ///< The threshold for whether increase the penalty for overlapping
constexpr double NLP_WN_CONJ_OOB_THRESHOLD = 0.05; ///< The threshold for wehther increasing the penalty for out of boundry
constexpr double NLP_WN_CONJ_ASYM_THRESHOLD = 0.5; ///< The threshodl for whether increase the penalty for asymmetry
constexpr double NLP_WN_CONJ_ALPHA = 1; ///< "alpha" should be a very small value. Used in objective function
PROJECT_NAMESPACE_END

#endif ///IDEAPLACE_PARAMETER_H_
