/**
 * @file Net.h
 * @brief The placement net data structure
 * @author Keren Zhu
 * @date 09/30/2019
 */

#ifndef IDEAPLACE_NET_H_
#define IDEAPLACE_NET_H_

#include "global/global.h"
#include "util/box.hpp"

PROJECT_NAMESPACE_BEGIN


class VirtualPin
{
    public:
        VirtualPin() = default;
        VirtualPin(const Point<LocType> &loc) : _loc(loc) {}
        const Point<LocType> & loc() const { return _loc; }
        Point<LocType> & loc() { return _loc; }
        LocType x() const { return _loc.x(); }
        LocType y() const { return _loc.y(); }
        IndexType netIdx() const { return _netIdx; }
        bool assigned() const { return _netIdx !=  INDEX_TYPE_MAX; }
        void free() { _netIdx = INDEX_TYPE_MAX; }
        void assign(IndexType netIdx) { _netIdx = netIdx; }
        void setDirection(Direction2DType dir) { _dir = dir; }
        Direction2DType direction() const { return _dir; }

        bool operator<(const VirtualPin &pin) const 
        {
            return _loc < pin._loc;
        }
        std::string toStr() const
        {
            std::stringstream ss;
            ss << _loc.toStr() <<" net: " << _netIdx;
            return ss.str();
        }
    private:
        Point<LocType> _loc;
        IndexType _netIdx = INDEX_TYPE_MAX;
        Direction2DType _dir; ///< The location on the placement boundary
};

enum class IoPinAssignment
{
    UNDEFINED = 0,
    LEFT = 1,
    RIGHT = 2
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
        const Point<LocType> &virtualPinLoc() const { return _virtualPin.loc(); }
        /// @brief get whether need to consider the virtual pin: If not IO net, or if no vitual pin assigned
        bool isValidVirtualPin() const { return (_isIo or _isVdd or _isVss) && _virtualPin.assigned(); }
        /// @brief get whether this net is a dummy net
        /// @return whether this net is a dummy net
        bool isDummyNet() const { return _isDummy; }
        /// @brief get whether the io pin is on top or bottom
        /// @return true: top or bottom
        /// @return false: left or right
        bool iopinVertical() const { return _virtualPin.direction() == Direction2DType::NORTH or _virtualPin.direction() == Direction2DType::SOUTH; }
        /// @brief the net index that is sym pair of this one
        /// @return the net index that is sym pair of this one
        IndexType symNetIdx() const { return _symNetIdx; }
        /// @brief get whether this net is self-symmetric
        /// @return whether this net is self-symmetric
        bool isSelfSym() const { return _isSelfSym; }
        /// @brief get whether this net has symmetric net
        bool hasSymNet() const { return _symNetIdx != INDEX_TYPE_MAX; }
        /// @brief get whether this net is vdd
        bool isVdd() const { return _isVdd; }
        /// @brief get whether this net is vss
        bool isVss() const { return _isVss; }
        /// @brief whether this should be on the left of the symnet pair
        bool isLeftSym() const { return _isLeftSym; }
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
        /// @brief set the symmetric pair net index of this one
        /// @param the net index that this one should be symmetric to
        void setSymNet(IndexType symNet, bool isLeftSym) { _symNetIdx = symNet; _isLeftSym = isLeftSym; }
        /// @brief mark this net as self symmetric
        void markSelfSym() { _isSelfSym = true; }
        /// @brief revoke the mark of this net is self symmetric
        void revokeSelfSym() { _isSelfSym = false; }
        /// @brief mark this net as vdd
        void markAsVdd() { _isVdd = true; _isVss = false; _isIo = false; }
        /// @brief mark this net as vss
        void markAsVss() { _isVss = true; _isVdd = false; _isIo = false; }
        /// @brief assign to left
        void fpIoPinAssignLeft() { _fpAssignment = IoPinAssignment::LEFT; }
        /// @brief assign to right
        void fpIoPinAssignRight() { _fpAssignment = IoPinAssignment::RIGHT; }
        /// @brief clear io pin assignment status from floorplan
        void clearFpIoPinAssign() { _fpAssignment = IoPinAssignment::UNDEFINED; }
        /// @brief get whether this net is assigned to particular location in the floorplan process
        /// @return true: this is either assigned to left or right. false: no
        bool isAssignedInFpIoPin() const { return _fpAssignment != IoPinAssignment::UNDEFINED; }
        /// @brief get whether this net is assigned to left
        /// @return true: it is assigned to left. false: no
        bool isLeftAssignedInFpIoPin() const { return _fpAssignment == IoPinAssignment::LEFT; }
        /// @brief get whether this net is assigned to right
        /// @return true: it is assigned to right. false: no
        bool isRightAssignedInFpIoPin() const { return _fpAssignment == IoPinAssignment::RIGHT; }
        /// @brief get whether this net has external bounding box set
        /// @return true: is set. false: not
        bool isExternalBBoxSet() const { return _externalBBox.valid(); }
        /// @brief get the external net bounding box
        /// @return the external net bounding box
        const Box<LocType> & externalBBox() const { return _externalBBox; }
        /// @brief set the external bounding box
        /// @param a box
        void setExternalBBox(const Box<LocType> &bbox) { _externalBBox = bbox; }
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
        IndexType _symNetIdx = INDEX_TYPE_MAX; ///< The symmetric pair of this net. If INDEX_TYPE_MAX, it does not have a sym pair
        bool _isLeftSym = true; ///< Whethet this should be on the left in sym net pair
        bool _isSelfSym = false; ///< Whether this net is self-symmetric
        bool _isVdd = false; ///< whether this net is vdd
        bool _isVss = false; ///< Whether this net is vss
        IoPinAssignment _fpAssignment = IoPinAssignment::UNDEFINED; ///< The assignment from the floorplan
        Box<LocType> _externalBBox = Box<LocType>(1,1, -1, -1); ///< The bounding box for the external nets
};

PROJECT_NAMESPACE_END

#endif /// IDEAPLACE_NET_H_
