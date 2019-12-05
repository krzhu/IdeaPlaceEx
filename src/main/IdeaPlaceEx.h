/**
 * @file IdeaPlaceEx.h
 * @brief wrapper of Everything
 * @author Keren Zhu
 * @date 10/02/2019
 */

#ifndef IDEAPLACE_IDEAPLACEEX_H_
#define IDEAPLACE_IDEAPLACEEX_H_

#include "db/Database.h"
/* Solver */
#include "place/CGLegalizer.h"
#include "place/NlpWnconj.h" //< This stupid package must be included here

PROJECT_NAMESPACE_BEGIN


/// @class IDEAPLACE::IdeaPlaceEx
/// @brief the main wrapper for the placement engine
class IdeaPlaceEx
{
    public:
        /// @brief default constructor
        explicit IdeaPlaceEx() = default;
        /// @brief the file-based input
        /// @param the system arguments
        /// @return if the parsing is successful
        bool parseFileBased(int argc, char** argv);
        /// @brief run the placement algorithm
        /// @return whether the placement is successful
        bool solve();
        /// @brief the file-based output
        /// @param the system arguments
        /// @return if the writing is successful
        bool outputFileBased(int argc, char** argv);
    protected:
        Database _db; ///< The placement engine database 
};

PROJECT_NAMESPACE_END

#endif ///IDEAPLACE_IDEAPLACEEX_H_
