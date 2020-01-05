
/*
    pbrt source code is Copyright(c) 1998-2016
                        Matt Pharr, Greg Humphreys, and Wenzel Jakob.

    This file is part of pbrt.

    Redistribution and use in source and binary forms, with or without
    modification, are permitted provided that the following conditions are
    met:

    - Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer.

    - Redistributions in binary form must reproduce the above copyright
      notice, this list of conditions and the following disclaimer in the
      documentation and/or other materials provided with the distribution.

    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
    IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
    TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
    PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
    HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
    SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
    LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
    DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
    THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
    (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
    OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

 */

#if defined(_MSC_VER)
#define NOMINMAX
#pragma once
#endif

#ifndef PBRT_CORE_FILTER_H
#define PBRT_CORE_FILTER_H

// core/filter.h*
#include "pbrt.h"
#include "geometry.h"

namespace pbrt {

// Filter Declarations
class Filter {
  public:
    // Filter Interface
    virtual ~Filter();
    Filter(const Vector2f &radius)
        : radius(radius), invRadius(Vector2f(1 / radius.x, 1 / radius.y)) {}

    // 参考 P473 式 7.12, 返回采样点在重构时对应的权重
    // 过滤器以(0, 0)为原点, 系统保证输入的 p 在半径距离 radius 内, 所以不会对 p 做检查
    // PBRT 使用的过滤器都是径向对称的, 也就是 x 和 y 可以分开算
    virtual Float Evaluate(const Point2f &p) const = 0;

    // Filter Public Data
    const Vector2f radius, invRadius;
};

}  // namespace pbrt

#endif  // PBRT_CORE_FILTER_H
