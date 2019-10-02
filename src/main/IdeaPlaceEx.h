/**
 * @file IdeaPlaceEx.h
 * @brief wrapper of Everything
 * @author Keren Zhu
 * @date 10/02/2019
 */

#ifndef IDEAPLACE_IDEAPLACEEX_H_
#define IDEAPLACE_IDEAPLACEEX_H_

#include "db/Database.h"

PROJECT_NAMESPACE_BEGIN

/// @class IDEAPLACE::IdeaPlaceEx
/// @brief the main wrapper for the placement engine
class IdeaPlaceEx
{
    protected:
        Database _db; ///< The placement engine database 
};

PROJECT_NAMESPACE_END

#endif ///IDEAPLACE_IDEAPLACEEX_H_
