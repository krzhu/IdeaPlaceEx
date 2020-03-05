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


class VirtualPin
{
    public:
        VirtualPin() = default;
        VirtualPin(const XY<LocType> &loc) : _loc(loc) {}
        const XY<LocType> & loc() const { return _loc; }
        XY<LocType> & loc() { return _loc; }
        LocType x() const { return _loc.x(); }
        LocType y() const { return _loc.y(); }
        IndexType cellIdx() const { return _cellIdx; }
        bool assigned() const { return _cellIdx !=  INDEX_TYPE_MAX; }
        void free() { _cellIdx = INDEX_TYPE_MAX; }
        void assign(IndexType cellIdx) { _cellIdx = cellIdx; }
        void setDirection(Direction2DType dir) { _dir = dir; }
        Direction2DType direction() const { return _dir; }

        bool operator<(const VirtualPin &pin) const 
        {
            return _loc < pin._loc;
        }
    private:
        XY<LocType> _loc;
        IndexType _cellIdx = INDEX_TYPE_MAX;
        Direction2DType _dir; ///< The location on the placement boundary
};
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
        /// @brief get the name for the net
        /// @return the name of this net
        const std::string & name() const { return _name; }
        /// @brief get whether this pin is an IO pin
        bool isIo() const { return _isIo; }
        /// @brief get the virtual pin location
        /// @return the location for the virtual pin
        const XY<LocType> &virtualPinLoc() const { return _virtualPin.loc(); }
        /// @brief get whether need to consider the virtual pin: If not IO net, or if no vitual pin assigned
        bool isValidVirtualPin() const { return _isIo && _virtualPin.assigned(); }
        /// @brief get whether this net is a dummy net
        /// @return whether this net is a dummy net
        bool isDummyNet() const { return _isDummy; }
        /// @brief get whether the io pin is on top or bottom
        /// @return true: top or bottom
        /// @return false: left or right
        bool iopinVertical() const { return _virtualPin.direction() == Direction2DType::NORTH or _virtualPin.direction() == Direction2DType::SOUTH; }
        /*------------------------------*/ 
        /* Setters                      */
        /*------------------------------*/ 
        /// @brief set the weight of the net
        /// @param the weight of the net
        void setWeight(IntType weight) { _weight = weight; }
        /// @brief set the name of the net
        /// @param the name of the net
        void setName(const std::string &name) { _name = name; }
        /// @brief set whether this net is an IO net
        void setIsIo(bool isIo) { _isIo =isIo; }
        /// @brief set the virtual pin location of this net
        void setVirtualPin(const VirtualPin &virtualPinLocation) { _virtualPin = virtualPinLocation; }
        /// @brief invalidate the virtual pin
        void invalidateVirtualPin() { _virtualPin.free(); }
        /// @brief mark this net as a dummy net
        void markAsDummyNet() { _isDummy = true; }
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
        IndexType pinIdx(IndexType idx) const { return _pinIdxArray.at(idx); }
        const std::vector<IndexType> & pinIdxArray() const { return _pinIdxArray; }
    private:
        std::string _name = ""; ///< The name for the net
        std::vector<IndexType> _pinIdxArray; ///< The index to the pins belonging to the net
        IntType _weight = 1; ///< The weight of this net
        bool _isIo = false; ///< Whether thisnet is an IO net 
        VirtualPin _virtualPin; ///< The virtual pin location
        bool _isDummy = false; ///< Whether this is a dummy net
};

PROJECT_NAMESPACE_END

#endif /// IDEAPLACE_NET_H_
