/**
 * @file IdeaPlaceEx.cpp
 * @brief The pybind11 interface for the core placement engine api
 * @author Keren Zhu
 * @date 12/16/2019
 */

#include "main/IdeaPlaceEx.h"
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;

void initIdeaPlaceExAPI(py::module &m) {
  py::class_<PROJECT_NAMESPACE::IdeaPlaceEx>(m, "IdeaPlaceEx")
      .def(py::init<>())
      .def("solve", &PROJECT_NAMESPACE::IdeaPlaceEx::solve, "Solve the problem")
      .def("endPlace", &PROJECT_NAMESPACE::IdeaPlaceEx::endPlace, "End and wrapping up the placement, should not be used together with solve()")
      .def("initGlobalPlacer",
           &PROJECT_NAMESPACE::IdeaPlaceEx::initGlobalPlacer,
           "Initialize a global placer")
      .def("initLegalizer",
          &PROJECT_NAMESPACE::IdeaPlaceEx::initLegalizer,
          "Initialize a legalizer")
      .def("alignToGrid", &PROJECT_NAMESPACE::IdeaPlaceEx::alignToGrid,
           "Align the placement to grid")
      .def("numThreads", &PROJECT_NAMESPACE::IdeaPlaceEx::setNumThreads,
           "Set number of threads")
      .def("setGridStep", &PROJECT_NAMESPACE::IdeaPlaceEx::setGridStep, "Set the grid step" )
      .def("readTechSimpleFile",
           &PROJECT_NAMESPACE::IdeaPlaceEx::readTechSimpleFile,
           "Internal usage: Read in the techsimple file")
      .def("readPinFile", &PROJECT_NAMESPACE::IdeaPlaceEx::readPinFile,
           "Internal usage: Read in the .pin file")
      .def("readConnectionFile",
           &PROJECT_NAMESPACE::IdeaPlaceEx::readConnectionFile,
           "Internal usage: Read in the .connection file")
      .def("readNetWgtFile", &PROJECT_NAMESPACE::IdeaPlaceEx::readNetWgtFile,
           "Internal usage: Read in the .netwgt file")
      .def("readSymFile", &PROJECT_NAMESPACE::IdeaPlaceEx::readSymFile,
           "Internal usage: Read in the .sym file")
      .def("readGdsLayout", &PROJECT_NAMESPACE::IdeaPlaceEx::readGdsLayout,
           "Internal Usage: Read in a gds for cell, the name of the GDS cell "
           "need to match the cell name")
      .def("readSymNetFile", &PROJECT_NAMESPACE::IdeaPlaceEx::readSymNetFile,
           "Internal Usage: Read in a .symnet file")
      .def("readSigpathFile", &PROJECT_NAMESPACE::IdeaPlaceEx::readSigpathFile,
           "Internal Usage: Read in a .sigpath file")
      .def("addGdsLayer", &PROJECT_NAMESPACE::IdeaPlaceEx::addGdsLayer,
           py::arg("cellIdx") = PROJECT_NAMESPACE::INDEX_TYPE_MAX,
           "Add a gds to a cell")
      .def("finishAddingGdsLayer",
           &PROJECT_NAMESPACE::IdeaPlaceEx::finishAddingGdsLayer,
           "Finish the gds layer adding, trigger a init function")
      .def("setIntraLayerSpacing",
          &PROJECT_NAMESPACE::IdeaPlaceEx::setIntraLayerSpacing)
      .def("setInterLayerSpacing",
          &PROJECT_NAMESPACE::IdeaPlaceEx::setInterLayerSpacing)
      .def("setNwellLayerIdx",
          &PROJECT_NAMESPACE::IdeaPlaceEx::setNwellLayerIdx)
      .def("allocateCell", &PROJECT_NAMESPACE::IdeaPlaceEx::allocateCell,
           "Allocate a new cell, return the index of the cell")
      .def("setCellName", &PROJECT_NAMESPACE::IdeaPlaceEx::setCellName,
           "Set the name of a cell")
      .def("setCellFlip", &PROJECT_NAMESPACE::IdeaPlaceEx::setCellFlip,
           "Flip a cell")
      .def("allocatePin", &PROJECT_NAMESPACE::IdeaPlaceEx::allocatePin,
           "Allocate a new pin. Return the index of the pin")
      .def("addPinShape", &PROJECT_NAMESPACE::IdeaPlaceEx::addPinShape,
           "Add a shape to a pin")
      .def("setPinName", &PROJECT_NAMESPACE::IdeaPlaceEx::setPinName,
           "Set the name of a pin")
      .def("allocateNet", &PROJECT_NAMESPACE::IdeaPlaceEx::allocateNet,
           "Allocate a new net. Return the index of the net")
      .def("setNetName", &PROJECT_NAMESPACE::IdeaPlaceEx::setNetName,
           "Set the name of a net")
      .def("setNetExternalBBox",
           &PROJECT_NAMESPACE::IdeaPlaceEx::setNetExternalBBox,
           "Set the external net bouding box of the circuit")
      .def("addPinToNet", &PROJECT_NAMESPACE::IdeaPlaceEx::addPinToNet,
           "Add a pin to a net")
      .def("setNetWgt", &PROJECT_NAMESPACE::IdeaPlaceEx::setNetWgt,
           "Set the weight of a net")
      .def("fpIoPinAssignLeft",
           &PROJECT_NAMESPACE::IdeaPlaceEx::fpIoPinAssignLeft,
           "set a net should be assigned to the left for IO pin")
      .def("fpIoPinAssignRight",
           &PROJECT_NAMESPACE::IdeaPlaceEx::fpIoPinAssignRight,
           "set a net should be assigned to the right for IO pin")
      .def("clearFpIoPinAssign",
           &PROJECT_NAMESPACE::IdeaPlaceEx::clearFpIoPinAssign,
           "Clear a net's IO Pin assignment")
      .def("allocateSymGrp", &PROJECT_NAMESPACE::IdeaPlaceEx::allocateSymGrp,
           "Allocate a new symmetric group")
      .def("addSymPair", &PROJECT_NAMESPACE::IdeaPlaceEx::addSymPair,
           "Add a new symmetric pair to a symmetric group")
      .def("addSelfSym", &PROJECT_NAMESPACE::IdeaPlaceEx::addSelfSym,
           "Add a self-symmetric constraint to a symmetric group")
      .def("addCellShape", &PROJECT_NAMESPACE::IdeaPlaceEx::addCellShape,
           "Add a shape to a cell")
      .def("numCells", &PROJECT_NAMESPACE::IdeaPlaceEx::numCells,
           "Get the number of cells")
      .def("cellIdxName", &PROJECT_NAMESPACE::IdeaPlaceEx::cellIdxName,
           "Get the index based on cell name")
      .def("cellName", &PROJECT_NAMESPACE::IdeaPlaceEx::cellName,
           "Get the cell name")
      .def("pinIdx", &PROJECT_NAMESPACE::IdeaPlaceEx::pinIdx,
           "Get the index based on pin name")
      .def("openVirtualPinAssignment",
           &PROJECT_NAMESPACE::IdeaPlaceEx::openVirtualPinAssignment,
           "Open the virtual pin assignment functionality")
      .def("closeVirtualPinAssignment",
           &PROJECT_NAMESPACE::IdeaPlaceEx::closeVirtualPinAssignment,
           "Close the virtual pin assignment functionality")
      .def("openFastMode", &PROJECT_NAMESPACE::IdeaPlaceEx::openFastMode,
           "Open the placement fast mode. Notice that the constraints might be "
           "waived")
      .def("closeFastMode", &PROJECT_NAMESPACE::IdeaPlaceEx::closeFastMode,
           "Close the placement fast mode.")
      .def("setBoundaryConstraint",
           &PROJECT_NAMESPACE::IdeaPlaceEx::setBoundaryConstraint,
           "Set the boundary constraints. params: xLo, yLo, xHi, yHi")
      .def("setIoPinBoundaryExtension",
           &PROJECT_NAMESPACE::IdeaPlaceEx::setIoPinBoundaryExtension,
           "Set the extension of io pin locations to the boundary of cell "
           "placements")
      .def("setIoPinInterval",
           &PROJECT_NAMESPACE::IdeaPlaceEx::setIoPinInterval,
           "Set the minimum interval of io pins")
      .def("markIoNet", &PROJECT_NAMESPACE::IdeaPlaceEx::markAsIoNet,
           "Mark a net as IO net")
      .def("revokeIoNet", &PROJECT_NAMESPACE::IdeaPlaceEx::revokeIoNet,
           "Revoke IO net flag on a net")
      .def("markAsVddNet", &PROJECT_NAMESPACE::IdeaPlaceEx::markAsVddNet,
           "Mark a net as VDD")
      .def("markAsVssNet", &PROJECT_NAMESPACE::IdeaPlaceEx::markAsVssNet,
           "Mark a net as VSS")
      .def("addSymNetPair", &PROJECT_NAMESPACE::IdeaPlaceEx::addSymNetPair,
           "Add a pair of symmetric nets")
      .def("markSelfSym", &PROJECT_NAMESPACE::IdeaPlaceEx::markSelfSymNet,
           "Mark a net as self symmetric")
      .def("allocateProximityGroup",
           &PROJECT_NAMESPACE::IdeaPlaceEx::allocateProximityGroup,
           "Allocate a new proximity group")
      .def("addCellToProximityGroup",
           &PROJECT_NAMESPACE::IdeaPlaceEx::addCellToProximityGroup,
           "Add a cell to a proximity group")
      .def("allocateSignalPath",
           &PROJECT_NAMESPACE::IdeaPlaceEx::allocateSignalPath,
           "Allocate a new signal path")
      .def("markSignalPathAsPower",
           &PROJECT_NAMESPACE::IdeaPlaceEx::markSignalPathAsPower,
           "Mark a signal path as power")
      .def("addPinToSignalPath",
           &PROJECT_NAMESPACE::IdeaPlaceEx::addPinToSignalPath,
           "add a pin to a signal path")
      .def("iopinX", &PROJECT_NAMESPACE::IdeaPlaceEx::iopinX,
           "Get the x coordinate of a net")
      .def("iopinY", &PROJECT_NAMESPACE::IdeaPlaceEx::iopinY,
           "Get the y coordinate of a net")
      .def("isIoPinVertical", &PROJECT_NAMESPACE::IdeaPlaceEx::isIopinVertical,
           "true if io pins are on top or bottom")
      .def("xCellLoc", &PROJECT_NAMESPACE::IdeaPlaceEx::xCellLoc,
           "Get x coordinate of a cell location")
      .def("yCellLoc", &PROJECT_NAMESPACE::IdeaPlaceEx::yCellLoc,
           "Get y coordinate of a cell location")
      .def("runtimeIdeaPlaceEx",
           &PROJECT_NAMESPACE::IdeaPlaceEx::runtimeIdeaPlaceEx,
           "Get the runtime for the Ideaplace")
      .def("runtimeGlobalPlace",
           &PROJECT_NAMESPACE::IdeaPlaceEx::runtimeGlobalPlace,
           "Get the time used for global placement")
      .def("runtimeGlobalPlaceCalcObj",
           &PROJECT_NAMESPACE::IdeaPlaceEx::runtimeGlobalPlaceCalcObj,
           "Get the time used for calculating the objectives in global "
           "placement")
      .def(
          "runtimeGlobalPlaceCalcGrad",
          &PROJECT_NAMESPACE::IdeaPlaceEx::runtimeGlobalPlaceCalcGrad,
          "Get the time used for calculating the gradients in global placement")
      .def("runtimeGlobalPlaceOptmKernel",
           &PROJECT_NAMESPACE::IdeaPlaceEx::runtimeGlobalPlaceOptmKernel,
           "Get the time used for optimizer kernel in the global placement")
      .def("runtimeGlobalPlaceOptimize",
           &PROJECT_NAMESPACE::IdeaPlaceEx::runtimeGlobalPlaceOptimize,
           "Get the time used for global placement optimize routine")
      .def("runtimeGlobalPlaceUpdateProblem",
           &PROJECT_NAMESPACE::IdeaPlaceEx::runtimeGlobalPlaceUpdateProblem,
           "Get the time used for updating the problem in gloobal placement")
      .def("runtimeLegalization",
           &PROJECT_NAMESPACE::IdeaPlaceEx::runtimeLegalization,
           "Get the time used for legalization")
      .def("runtimeDetailedPlacement",
           &PROJECT_NAMESPACE::IdeaPlaceEx::runtimeDetailedPlacement,
           "Get the time used for detailed placement")
      .def("logScreenOn", &PROJECT_NAMESPACE::IdeaPlaceEx::logScreenOn)
      .def("logScreenOff", &PROJECT_NAMESPACE::IdeaPlaceEx::logScreenOff)
      .def("hpwl", &PROJECT_NAMESPACE::IdeaPlaceEx::hpwl, "Get HPWL")
      // Well
      .def("allocateWell", &PROJECT_NAMESPACE::IdeaPlaceEx::allocateWell,
           "Allocate a new well")
      .def("setWellName", &PROJECT_NAMESPACE::IdeaPlaceEx::setWellName,
           "Set the name of a well")
      .def("setWellShape", &PROJECT_NAMESPACE::IdeaPlaceEx::setWellShape,
           "Set the shape (polygon) of a well")
      .def("printWellInfo", &PROJECT_NAMESPACE::IdeaPlaceEx::printWellInfo,
           "Print the info of a well")
      .def("clearWells", &PROJECT_NAMESPACE::IdeaPlaceEx::clearWells, 
          "Clear the current wells")
      .def("setCellWellType", &PROJECT_NAMESPACE::IdeaPlaceEx::setCellWellType, "Set the well type of the cell")
      .def("setNumWellTypes", &PROJECT_NAMESPACE::IdeaPlaceEx::setNumWellTypes, "Set the number of well types to be considered")
      .def("setCellFingerChannelWidth", &PROJECT_NAMESPACE::IdeaPlaceEx::setCellFingerChannelWidth, "Set the finger width for the device")
      .def("setCellFingerChannelLength", &PROJECT_NAMESPACE::IdeaPlaceEx::setCellFingerChannelLength, "Set the finger length for the device")
      .def("assignCellToWell", &PROJECT_NAMESPACE::IdeaPlaceEx::assignCellToWell, "Assign cells to wells")
      .def("splitWells", &PROJECT_NAMESPACE::IdeaPlaceEx::splitWells, "Split the wells into rectangles")
      .def("numWells", &PROJECT_NAMESPACE::IdeaPlaceEx::numWells, "Get the number of well ")
      .def("numWellRects", &PROJECT_NAMESPACE::IdeaPlaceEx::numWellRects, "Get the number of well rects")
      .def("wellRect", &PROJECT_NAMESPACE::IdeaPlaceEx::wellRect, "Get one rect from the well")
      .def("setVddContactRequiredSpacing", &PROJECT_NAMESPACE::IdeaPlaceEx::setVddContactRequiredSpacing, "Set the required spacing for VDD contacts")
      .def("addVddContactTemplate", &PROJECT_NAMESPACE::IdeaPlaceEx::addVddContactTemplate, "Add a new VDD contact template")
      .def("addWpeHorizontalSpacingRule", &PROJECT_NAMESPACE::IdeaPlaceEx::addWpeHorizontalSpacingRule, "Add a WPE horizontal spacing rule: a pair of length and spacing")
      .def("addWpeVerticalSpacingRule", &PROJECT_NAMESPACE::IdeaPlaceEx::addWpeVerticalSpacingRule, "Add a WPE vertical spacing rule: a pair of width and spacing")
      .def("numContacts", &PROJECT_NAMESPACE::IdeaPlaceEx::numContacts)
      .def("contactTemplate", &PROJECT_NAMESPACE::IdeaPlaceEx::contactTemplate)
      .def("contactPosition", &PROJECT_NAMESPACE::IdeaPlaceEx::contactPosition)
      .def("individualWell", &PROJECT_NAMESPACE::IdeaPlaceEx::individualWell, "Only for comparsion purpose")
#ifdef DEBUG_DRAW
      .def("debugDraw", &PROJECT_NAMESPACE::IdeaPlaceEx::drawDebug, "Draw the debug layout")
#endif // ifdef DEBUG_DRAW
      ;
}

void initPointAPI(py::module &m) {
  py::class_<PROJECT_NAMESPACE::Point<int>>(m, "Point")
      .def(py::init<int, int>())
      .def("set", &PROJECT_NAMESPACE::Point<int>::set)
      .def("setX", &PROJECT_NAMESPACE::Point<int>::setX)
      .def("setY", &PROJECT_NAMESPACE::Point<int>::setY)
      .def("setXY", &PROJECT_NAMESPACE::Point<int>::setXY)
      .def("shift", &PROJECT_NAMESPACE::Point<int>::shift)
      .def("shiftX", &PROJECT_NAMESPACE::Point<int>::shiftX)
      .def("shiftY", &PROJECT_NAMESPACE::Point<int>::shiftY)
      .def("shiftXY", &PROJECT_NAMESPACE::Point<int>::shiftXY)
      .def("x", &PROJECT_NAMESPACE::Point<int>::x)
      .def("y", &PROJECT_NAMESPACE::Point<int>::y)
      .def("rotate90", &PROJECT_NAMESPACE::Point<int>::rotate90)
      .def("rotate180", &PROJECT_NAMESPACE::Point<int>::rotate180)
      .def("flipX", &PROJECT_NAMESPACE::Point<int>::flipX)
      .def("flipY", &PROJECT_NAMESPACE::Point<int>::flipY);
}


void initBoxAPI(py::module &m) {
  using LocType = PROJECT_NAMESPACE::LocType;
  using BOX = PROJECT_NAMESPACE::Box<LocType>;
  py::class_<BOX>(m, "Box")
      .def(py::init<LocType, LocType, LocType, LocType>())
      .def("xLo", &BOX::xLo)
      .def("xHi", &BOX::xHi)
      .def("yLo", &BOX::yLo)
      .def("yHi", &BOX::yHi)
      ;
}
