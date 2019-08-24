
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

#ifndef PBRT_CORE_QUATERNION_H
#define PBRT_CORE_QUATERNION_H

/*

## Section2.9.1 Quaternions

### 四元数的性质

### 四元数和旋转矩阵之间的转换

There is a useful relationship between unit quaternions (quaternions whose components satisfy $\left.x^{2}+y^{2}+z^{2}+w^{2}=1\right)$ and the space of rotations in $\mathbf{R}^{3} :$ specifically, a rotation of angle 2$\theta$ about a unit axis $\hat{\mathbf{v}}$ can be mapped to a unit quaternion $(\hat{\mathbf{v}} \sin \theta, \cos \theta),$ in which case the following quaternion product is equivalent to applying the rotation to a point $\mathrm{p}$ expressed in homogeneous coordinate form: $$ \mathrm{p}^{\prime}=\mathrm{qpq}^{-1} $$
对于绕轴 $\hat{\mathbf{v}}$ 旋转 2$\theta$ 角度的变换, 可以用单位四元数 $(\hat{\mathbf{v}} \sin \theta, \cos \theta),$ 来表示
Furthermore, the product of several rotation quaternions produces another quaternion that is equivalent to applying the rotations in sequence.
同时四元数的旋转也是可结合的

## Section2.9.2 Quaternion Interpolation

### slerp 的计算

*/

// core/quaternion.h*
#include "pbrt.h"
#include "stringprint.h"
#include "geometry.h"

namespace pbrt {

// Quaternion Declarations
struct Quaternion {
    // Quaternion Public Methods
    // 初始化为单位四元数
    Quaternion() : v(0, 0, 0), w(1) {}
    Quaternion &operator+=(const Quaternion &q) {
        v += q.v;
        w += q.w;
        return *this;
    }
    friend Quaternion operator+(const Quaternion &q1, const Quaternion &q2) {
        Quaternion ret = q1;
        return ret += q2;
    }
    Quaternion &operator-=(const Quaternion &q) {
        v -= q.v;
        w -= q.w;
        return *this;
    }
    Quaternion operator-() const {
        Quaternion ret;
        ret.v = -v;
        ret.w = -w;
        return ret;
    }
    friend Quaternion operator-(const Quaternion &q1, const Quaternion &q2) {
        Quaternion ret = q1;
        return ret -= q2;
    }
    Quaternion &operator*=(Float f) {
        v *= f;
        w *= f;
        return *this;
    }
    Quaternion operator*(Float f) const {
        Quaternion ret = *this;
        ret.v *= f;
        ret.w *= f;
        return ret;
    }
    Quaternion &operator/=(Float f) {
        v /= f;
        w /= f;
        return *this;
    }
    Quaternion operator/(Float f) const {
        Quaternion ret = *this;
        ret.v /= f;
        ret.w /= f;
        return ret;
    }

    // 从四元数生成旋转矩阵
    Transform ToTransform() const;
    // 从旋转矩阵生成四元数
    Quaternion(const Transform &t);

    friend std::ostream &operator<<(std::ostream &os, const Quaternion &q) {
        os << StringPrintf("[ %f, %f, %f, %f ]", q.v.x, q.v.y, q.v.z,
                           q.w);
        return os;
    }

    // Quaternion Public Data 
    // $\mathbf{q}=\left(\mathbf{q}_{x}, \mathbf{q}_{y}, \mathbf{q}_{z}, \mathbf{q}_{w}\right)$ or as $\mathbf{q}=\left(\mathbf{q}_{x y z}, \mathbf{q}_{w}\right)$
    Vector3f v;
    Float w;
};

// 用四元数进行球面线性插值(spherical linear interpolation)
// Slerp 具有转矩最小, 角速度恒定的优点
Quaternion Slerp(Float t, const Quaternion &q1, const Quaternion &q2);

// Quaternion Inline Functions
inline Quaternion operator*(Float f, const Quaternion &q) { return q * f; }

// 这里做的不是四元数的乘法, 而是内积
// $\mathbf{q} \cdot \mathbf{q}^{\prime}=x x^{\prime}+y y^{\prime}+z z^{\prime}+w w^{\prime}$
inline Float Dot(const Quaternion &q1, const Quaternion &q2) {
    return Dot(q1.v, q2.v) + q1.w * q2.w;
}

inline Quaternion Normalize(const Quaternion &q) {
    return q / std::sqrt(Dot(q, q));
}

}  // namespace pbrt

#endif  // PBRT_CORE_QUATERNION_H
