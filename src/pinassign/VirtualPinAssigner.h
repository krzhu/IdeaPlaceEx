/**
 * @file VirtualPinAssigner.h
 * @brief Kernel for assigning virtual pins on the placement boundary
 * @author Keren Zhu
 * @date 02/21/2020
 */

#ifndef IDEAPLACE_VIRTUAL_PIN_ASSIGNMER_H_
#define IDEAPLACE_VIRTUAL_PIN_ASSIGNMER_H_

#include "db/Database.h"

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
    private:
        XY<LocType> _loc;
        IndexType _cellIdx = INDEX_TYPE_MAX;
};

/// @class IDEAPLACE::VirtualPinAssigner
/// @brief The kernel for assigning virtual pins
class VirtualPinAssigner
{
    public:
        explicit VirtualPinAssigner(Database &db) : _db(db) 
        {
            _virtualPinInterval = db.parameters().virtualPinInterval();
            _virtualBoundaryExtension = db.parameters().virtualBoundaryExtension();
        }
        /* Kernal interface */
        /// @brief cnfigure the virtual boundary based on databse
        void reconfigureVirtualPinLocationFromDB();
        /// @brief reconfigure the virtual boundary and pin locations
        void reconfigureVirtualPinLocations(const Box<LocType> &cellsBBox);
        /// @brief solve pin assignment from information from DB
        bool pinAssignmentFromDB();
        /// @brief short cut to solve the problem from databse information
        bool solveFromDB();
        /// @brief solve the pin assignment problem. Export the solution to the database
        /// @param a function for query cell location
        bool pinAssignment(std::function<XY<LocType>(IndexType)> cellLocQueryFunc);
        /* Parameter settings */
        /// @brief set the extension distance of placement boundary to cell boundary
        void setVirtualBoundaryExtension(LocType ex) { _virtualBoundaryExtension = ex; }
        /// @brief set the interval between pins
        void setVirtualPinInterval(LocType in) { _virtualPinInterval = in; }
    private:
        bool symPinAssign(std::function<XY<LocType>(IndexType)> cellLocQueryFunc,
                std::function<LocType(IndexType, IndexType)> edgeCostFunc);
    private:
        Database &_db; ///< The placement database
        Box<LocType> _boundary; ///< The virtual boundary of the placement
        std::vector<VirtualPin> _virtualPins; ///< The locations for virtual pins
        LocType _virtualBoundaryExtension = -1; ///< The extension to placement cell bounding box
        LocType _virtualPinInterval = -1; ///< The interval between virtual pins
        std::vector<IndexType> _leftVirtualPins;
        std::vector<IndexType> _rightVirtualPins;
};

PROJECT_NAMESPACE_END

#endif //IDEAPLACE_VIRTUAL_PIN_ASSIGNMER_H_
