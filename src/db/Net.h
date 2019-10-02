/**
 * @file Net.h
 * @brief The placement net data structure
 * @author Keren Zhu
 * @date 09/30/2019
 */

#ifndef IDEAPLACE_NET_H_
#define IDEAPLACE_NET_H_

#include "global/global.h"

PROJECT_NAMESPACE_BEGIN

/// @class IDEAPLACE::Net
/// @brief the net class
class Net
{
    public:
        /// @brief default constructor
        explicit Net() = default;
        /*------------------------------*/ 
        /* Getters                      */
        /*------------------------------*/ 
        /// @brief get the weight of the net
        /// @return the weight of this net
        IntType weight() const { return _weight; }
        /*------------------------------*/ 
        /* Setters                      */
        /*------------------------------*/ 
        /// @brief set the weight of the net
        /// @param the weight of the net
        void setWeight(IntType weight) { _weight = weight; }
        /*------------------------------*/ 
        /* Vector operations            */
        /*------------------------------*/ 
        /// @brief get the number of pins
        /// @return the number of pins this cell has
        IndexType numPinIdx() const { return _pinIdxArray.size(); }
        /// @brief add one pin to the cell
        /// @param the index of pin that want to be added
        void addPin(IndexType pinIdx) { _pinIdxArray.emplace_back(pinIdx); }
        /// @brief get the database pin index
        /// @param the index of pin array
        /// @return the databse pin index
        IndexType pinIdx(IndexType idx) const { return AT(_pinIdxArray, idx); }
    private:
        std::vector<IndexType> _pinIdxArray; ///< The index to the pins belonging to the net
        IntType _weight = 1; ///< The weight of this net
};

PROJECT_NAMESPACE_END

#endif /// IDEAPLACE_NET_H_
