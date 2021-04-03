/**
 * @file gplacerApi.cpp
 * @brief The pybind11 interface for the core placement engine api
 * @author Keren Zhu
 * @date 12/16/2019
 */

#include "main/IdeaPlaceEx.h"
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;

void initGPlacerApi(py::module &m) {
  using GP = PROJECT_NAMESPACE::IdeaPlaceEx::GPlacer;
  py::class_<GP>(m, "GPlacer")
      .def(py::init<PROJECT_NAMESPACE::Database &>())
      .def("writeLocs", &GP::writeLocs, "Write the cell locations in database")
      .def("readLocs", &GP::readLocs, "Read the cell locations in database")
      .def("reinitWellOperators", &GP::reinitWellOperators, "Re-initialize the well-related operators and tasks")
      .def("prepareWellAwarePlace", &GP::prepareWellAwarePlace, "Prepare the well-aware placement problem")
      .def("stepOptmIter", &GP::stepOptmIter, "Run an iteration of optimization problem")
      .def("overlapAreaRatio", &GP::overlapAreaRatio, "Calculate the current overlapping area ratio over total cell area")
      ;
}
