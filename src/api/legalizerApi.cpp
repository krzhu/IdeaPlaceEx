/**
 * @file legalizer.cpp
 * @brief The pybind11 interface for the legalizer engine api
 * @author Keren Zhu
 * @date 04/12/2021
 */

#include "place/CGLegalizer.h"
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;

void initLegalizerApi(py::module &m) {
  using LEGAL = PROJECT_NAMESPACE::CGLegalizer;
  py::class_<LEGAL>(m, "Legalizer")
      .def(py::init<PROJECT_NAMESPACE::Database &>())
      .def("prepare", &LEGAL::prepare, "Prepare the legalizer")
      .def("preserveRelationCompaction", &LEGAL::preserveRelationCompaction, "Compaction that preserve the current relations")
      .def("areaDrivenCompaction", &LEGAL::areaDrivenCompaction, "Area-driven legalization")
      .def("wirelengthDrivenCompaction", &LEGAL::wirelengthDrivenCompaction, "Wirelength driven legalization")
      .def("openWellAware", &LEGAL::openWellAware, "Open well-aware legalization")
      .def("closeWellAware", &LEGAL::closeWellAware, "Close well-aware legalization")
      ;
}
