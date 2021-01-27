/**
 * @file api.cpp
 * @brief The top level cpp for initialize the pybind module
 * @author Keren Zhu
 * @date 12/16/2019
 */

#include <pybind11/pybind11.h>
#include <pybind11/stl_bind.h>
#include "global/global.h"

namespace py = pybind11;

void initIdeaPlaceExAPI(py::module &);
void initPointAPI(py::module&);

PYBIND11_MAKE_OPAQUE(std::vector<PROJECT_NAMESPACE::IndexType>);
PYBIND11_MAKE_OPAQUE(std::vector<std::pair<PROJECT_NAMESPACE::LocType, PROJECT_NAMESPACE::LocType>>);

PYBIND11_MODULE(IdeaPlaceExPy, m)
{
  initIdeaPlaceExAPI(m);
  initPointAPI(m);
}
