/**
 * @file Database.h
 * @brief The placement database data structure
 * @author Keren Zhu
 * @date 10/02/2019
 */

#ifndef IDEAPLACE_DATABASE_H_
#define IDEAPLACE_DATABASE_H_

#include "Cell.h"
#include "Net.h"
#include "Pin.h"
#include "Tech.h"

PROJECT_NAMESPACE_BEGIN

/// @class IDEAPLACE::Database
/// @brief the database class of the placement engine
class Database
{
    public:
        /// @brief default database
        explicit Database() = default;
        /*------------------------------*/ 
        /* Initialization               */
        /*------------------------------*/ 
        /// @brief initializing cells. Allocate the correct number of layers
        /// @return if successful
        bool initCells();
        /*------------------------------*/ 
        /* Getters                      */
        /*------------------------------*/ 
        /// @brief get the technology-wrapper of the placement
        /// @return the technology
        const Tech & tech() const { return _tech; }
        /// @brief get the technology-wrapper of the placement
        /// @return the technology
        Tech & tech() { return _tech; }
        /*------------------------------*/ 
        /* Vector operations            */
        /*------------------------------*/ 
        /// @brief get the number of symmetric groups
        /// @return the number of the symmetric groups
        IndexType numSymGroups() const { return 0; }
        /// @brief get the number of cells
        /// @return the number of cells
        IndexType numCells() const { return _cellArray.size(); }
        /// @brief get one cell from the array
        /// @param the index of the cell
        /// @return the cell of the index
        const Cell & cell(IndexType cellIdx) const { return AT(_cellArray, cellIdx); }
        /// @brief get one cell from the array
        /// @param the index of the cell
        /// @return the cell of the index
        Cell & cell(IndexType cellIdx) { return AT(_cellArray, cellIdx); }
        /// @brief allocate a new cell in the array
        /// @return the index of the new cell
        IndexType allocateCell() { _cellArray.emplace_back(Cell()); return _cellArray.size() - 1; }
        /// @brief get the number of nets
        /// @param the index of the nets
        /// @return the net of the index
        IndexType numNets() const { return _netArray.size(); }
        /// @brief get one net from the array
        /// @param the index of the net
        /// @return the net of the index
        const Net & net(IndexType netIdx) const { return AT(_netArray, netIdx); }
        /// @brief get one net from the array
        /// @param the index of the net
        /// @return the net of the index
        Net & net(IndexType netIdx) { return AT(_netArray, netIdx); }
        /// @brief allocate a new net in the array
        /// @return index of the new net
        IndexType allocateNet() { _netArray.emplace_back(Net()); return _netArray.size() -1; }
        /// @brief get the number of pins
        /// @return the number of pins
        IndexType numPins() const { return _pinArray.size(); }
        /// @brief get one pin from the array
        /// @param the index of the pin
        /// @return the pin of the index
        const Pin & pin(IndexType pinIdx) const { return AT(_pinArray, pinIdx); }
        /// @brief get one pin from the array
        /// @param the index of the pin
        /// @return the pin of the index
        Pin & pin(IndexType pinIdx) { return AT(_pinArray, pinIdx); }
        /// @brief allocate a new pin in the array
        /// @return the index of the pin
        IndexType allocatePin() { _pinArray.emplace_back(Pin()); return _pinArray.size() - 1; }
        /*------------------------------*/ 
        /* Supporting functions         */
        /*------------------------------*/ 
        /// @brief calcuate the total cell area
        /// @return the total cell area
        LocType calculateTotalCellArea() const;

    private:
        std::vector<Cell> _cellArray; ///< The cells of the placement problem
        std::vector<Net> _netArray; ///< The nets of the placement problem
        std::vector<Pin> _pinArray; ///< The pins of the placement problem
        Tech _tech; ///< The tech information
};

inline LocType Database::calculateTotalCellArea() const
{
}

PROJECT_NAMESPACE_END

#endif /// IDEAPLACE_DATABASE_H_
