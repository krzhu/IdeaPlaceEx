/**
 * @file Parameters.h
 * @brief The placement parameters
 * @author Keren Zhu
 * @date 10/16/2019
 */

#ifndef IDEAPLACE_PARAMETERS_H_
#define IDEAPLACE_PARAMETERS_H_

#include "global/global.h"
#include "util/Box.h"

PROJECT_NAMESPACE_BEGIN

/// @class IDEAPLACE::Parameters
/// @brief The placement eingie parameters
class Parameters
{
    public:
        /// @brief default constrcutor
        explicit Parameters() = default;
        /*------------------------------*/ 
        /* Set the parameters           */
        /*------------------------------*/ 
        /// @brief set the boundry constraints
        /// @param the placement boundry
        void setBoundaryConstraint(const Box<LocType> &boundaryConstraint) { _boundaryConstraint = boundaryConstraint; }
        /// @brief open the functionality of virtual pin assignment
        void openVirtualPinAssignment() { _ifUsePinAssignment = true; }
        /// @brief close the functionality of virtual pin assignment 
        void closeVirtualPinAssignment() { _ifUsePinAssignment = false; }
        /// @brief set the number of threads
        void setNumThreads(IndexType numThreads) { _numThreads = numThreads; }
        /// @brief set the grid step constraint
        void setGridStep(LocType gridStep) { _gridStep = gridStep; } 
        /// @brief set the virutal boundary extension. The extension of boundary to io pin with respect to the cell placement
        void setVirtualBoundaryExtension(LocType virtualBoundaryExtension) { _virtualBoundaryExtension = virtualBoundaryExtension; _layoutOffset = 2 * virtualBoundaryExtension; }
        /// @brief set the pin interval 
        void setVirtualPinInterval(LocType virtualPinInterval) { _virtualPinInterval  = virtualPinInterval; }
        /*------------------------------*/ 
        /* Query the parameters         */
        /*------------------------------*/ 
        /// @brief whether the boundry constraint is set
        /// @return whether the boundry constraint is set
        bool isBoundaryConstraintSet() const { return _boundaryConstraint.valid(); }
        /// @brief get the boundry constraint
        /// @return the boundry constraint
        const Box<LocType> & boundaryConstraint() const { return _boundaryConstraint; }
        /// @brief get whether to use the virtual pin assignment functionality
        bool ifUsePinAssignment() const { return _ifUsePinAssignment; }
        /// @brief get the number of thread
        IndexType numThreads() const { return _numThreads; }
        /// @brief get the grid step
        LocType gridStep() const { return _gridStep; }
        /// @brief get wether there is grid step constraint
        bool hasGridStep() const { return _gridStep > 0; }
        /// @brief get the extension of virtual boundary of cell placement boundary
        LocType virtualBoundaryExtension() const { return _virtualBoundaryExtension; }
        /// @brief get the interval of virtual io pins
        LocType virtualPinInterval() const { return _virtualPinInterval; }
        /// @brief get the layout offset
        LocType layoutOffset() const { return _layoutOffset; }
    private:
        Box<LocType> _boundaryConstraint = Box<LocType>(LOC_TYPE_MAX, LOC_TYPE_MAX, LOC_TYPE_MIN, LOC_TYPE_MIN);
        bool _ifUsePinAssignment = true; ///< If do pin assignment
        IndexType _numThreads = 10;
        LocType _gridStep = -1;
        LocType _virtualBoundaryExtension = 200; ///< The extension of current virtual boundary to the bounding box of placement
        LocType _virtualPinInterval = 400; ///< The interval between each virtual pin
        LocType _layoutOffset = 1000; ///< The default offset for the placement

};

PROJECT_NAMESPACE_END

#endif ///IDEAPLACE_PARAMETERS_H_
