/**
 * @file Constraints.h
 * @brief The placement database data structure
 * @author Keren Zhu
 * @date 10/02/2019
 */

#ifndef IDEAPLACE_DATABASE_CONSTRAINTS_H_
#define IDEAPLACE_DATABASE_CONSTRAINTS_H_

#include "global/global.h"

PROJECT_NAMESPACE_BEGIN

/// @brief the data structures for symmetric pairs
class SymPair
{
    public:
        /// @brief default constructor
        explicit SymPair() = default;
        /// @brief constructor
        /// @param one cell index
        /// @param another cell index
        explicit SymPair(IndexType cellIdx1, IndexType cellIdx2) { this->setCellIdx(cellIdx1, cellIdx2); }
        /// @brief set the cell indices. The order stored is always left is lower in index
        /// @param one cell index
        /// @param another cell index
        void setCellIdx(IndexType cellIdx1, IndexType cellIdx2)
        {
            if (cellIdx1 < cellIdx2)
            {
                _firstCell = cellIdx1;
                _secondCell = cellIdx2;
            }
            else
            {
                _secondCell = cellIdx1;
                _firstCell = cellIdx2;
            }
        }
        /// @brief get the first cell index
        /// @return the first cell index
        IndexType firstCell() const { Assert(_firstCell != INDEX_TYPE_MAX);  return _firstCell; }
        /// @brief get the second cell index
        /// @return the second cell index
        IndexType secondCell() const { Assert(_secondCell != INDEX_TYPE_MAX);  return _secondCell; }
    private:
        IndexType _firstCell = INDEX_TYPE_MAX; ///< The first cell
        IndexType _secondCell = INDEX_TYPE_MAX; ///< The second cell
};

PROJECT_NAMESPACE_END

#endif //IDEAPLACE_DATABASE_CONSTRAINTS_H_
