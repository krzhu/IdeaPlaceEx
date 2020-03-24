/**
 * @file init_different.h
 * @brief The methods to init differentiation operators
 * @author Keren Zhu
 * @date 03/24/2020
 */

#ifndef IDEAPLACE_INIT_DIFFERENT_H_
#define IDEAPLACE_INIT_DIFFERENT_H_

#include "db/Database.h"
#include "different.h"

PROJECT_NAMESPACE_BEGIN

namespace diff
{

    enum class OpEnumType
    {
        hpwl, ovl, oob, asym
    };
} // namespace diff

PROJECT_NAMESPACE_END

#endif //IDEAPLACE_INIT_DIFFERENT_H_
