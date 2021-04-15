#include "Cell.h"
#include "util/Polygon2Rect.h"

PROJECT_NAMESPACE_BEGIN

void Well::splitIntoRects() {
  _splitedRects.clear();
  if (not klib::convertPolygon2Rects(shape().outer(), _splitedRects)) {
    ERR("Well:: cannot split well polygon! \n");
    Assert(false);
  }
}

PROJECT_NAMESPACE_END
