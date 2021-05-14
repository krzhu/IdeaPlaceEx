#include "Cell.h"
#include "util/Polygon2Rect.h"

PROJECT_NAMESPACE_BEGIN

void Well::splitIntoRects(IntType mode) {
  limbo::geometry::slicing_orientation_2d slicing_orient = limbo::geometry::HOR_VER_SLICING;
  if (mode == 1) {
    slicing_orient = limbo::geometry::HORIZONTAL_SLICING;
  }
  if (mode == 2) {
    slicing_orient = limbo::geometry::VERTICAL_SLICING;
  }
  _splitedRects.clear();
  if (not klib::convertPolygon2Rects(shape().outer(), _splitedRects, slicing_orient)) {
    ERR("Well:: cannot split well polygon! \n");
    _splitedRects.emplace_back(this->boundingBox());
  }
}

PROJECT_NAMESPACE_END
