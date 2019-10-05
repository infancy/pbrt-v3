
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

#ifndef PBRT_CORE_GEOMETRY_H
#define PBRT_CORE_GEOMETRY_H

/*

from bottom to top:

pbrt.h(每个头文件都包含了 pbrt.h, 不妨忽略掉它)
 \    \ 
  \    stringprint.h
   \      /   \
    \    /     \
  geometry.h quaternion.h
      \      /
       \    /
     transform.h
        ...

*/

/*

## Section2.1 Coordinate Systems

$$ \mathbf{v}=s_{1} \mathbf{v}_{1}+\cdots+s_{n} \mathbf{v}_{n} $$
$$ \mathrm{p}=\mathrm{p}_{0}+s_{1} \mathbf{v}_{1}+\cdots+s_{n} \mathbf{v}_{n} $$

PBRT 使用左手系

*/

/*

Overview:

- isNaN

- class
    - Vector2
    - Vector3
    - Point2
    - Point3
    - Normal3
    - Bounds2
    - Bounds3
    - Ray
    - RayDifferential

- function
    - Vector3
    - Vector2
    - Point3
    - Point2
    - Normal3
    - Bounds3
    - Bounds2
    - others

*/

// core/geometry.h*
#include "pbrt.h"
#include "stringprint.h"
#include <iterator>

namespace pbrt {

template <typename T>
inline bool isNaN(const T x) {
    return std::isnan(x);
}
template <>
inline bool isNaN(const int x) {
    return false;
}

// Vector Declarations
template <typename T>
class Vector2 {
  public:
    // Vector2 Public Methods
    Vector2() { x = y = 0; }
    Vector2(T xx, T yy) : x(xx), y(yy) { DCHECK(!HasNaNs()); }
    bool HasNaNs() const { return isNaN(x) || isNaN(y); }
    explicit Vector2(const Point2<T> &p);
    explicit Vector2(const Point3<T> &p);
#ifndef NDEBUG
    // The default versions of these are fine for release builds; for debug
    // we define them so that we can add the Assert checks.
    // 在 release 中用编译器生成的, debug 中为了断言测试用手写版本的
    Vector2(const Vector2<T> &v) {
        DCHECK(!v.HasNaNs());
        x = v.x;
        y = v.y;
    }
    Vector2<T> &operator=(const Vector2<T> &v) {
        DCHECK(!v.HasNaNs());
        x = v.x;
        y = v.y;
        return *this;
    }
#endif  // !NDEBUG

    Vector2<T> operator+(const Vector2<T> &v) const {
        DCHECK(!v.HasNaNs());
        return Vector2(x + v.x, y + v.y);
    }

    Vector2<T> &operator+=(const Vector2<T> &v) {
        DCHECK(!v.HasNaNs());
        x += v.x;
        y += v.y;
        return *this;
    }
    Vector2<T> operator-(const Vector2<T> &v) const {
        DCHECK(!v.HasNaNs());
        return Vector2(x - v.x, y - v.y);
    }

    Vector2<T> &operator-=(const Vector2<T> &v) {
        DCHECK(!v.HasNaNs());
        x -= v.x;
        y -= v.y;
        return *this;
    }
    // 就这么直接比较好吗
    bool operator==(const Vector2<T> &v) const { return x == v.x && y == v.y; }
    bool operator!=(const Vector2<T> &v) const { return x != v.x || y != v.y; }
    template <typename U>
    Vector2<T> operator*(U f) const {
        return Vector2<T>(f * x, f * y);
    }

    template <typename U>
    Vector2<T> &operator*=(U f) {
        DCHECK(!isNaN(f));
        x *= f;
        y *= f;
        return *this;
    }
    template <typename U>
    Vector2<T> operator/(U f) const {
        CHECK_NE(f, 0);
        Float inv = (Float)1 / f;
        return Vector2<T>(x * inv, y * inv);
    }

    template <typename U>
    Vector2<T> &operator/=(U f) {
        CHECK_NE(f, 0);
        Float inv = (Float)1 / f;
        x *= inv;
        y *= inv;
        return *this;
    }
    Vector2<T> operator-() const { return Vector2<T>(-x, -y); }
    T operator[](int i) const {
        DCHECK(i >= 0 && i <= 1);
        if (i == 0) return x;
        return y;
    }

    T &operator[](int i) {
        DCHECK(i >= 0 && i <= 1);
        if (i == 0) return x;
        return y;
    }
    Float LengthSquared() const { return x * x + y * y; }
    Float Length() const { return std::sqrt(LengthSquared()); }

    // Vector2 Public Data
    T x, y;
};

template <typename T>
inline std::ostream &operator<<(std::ostream &os, const Vector2<T> &v) {
    os << "[ " << v.x << ", " << v.y << " ]";
    return os;
}

template <>
inline std::ostream &operator<<(std::ostream &os, const Vector2<Float> &v) {
    os << StringPrintf("[ %f, %f ]", v.x, v.y);
    return os;
}

template <typename T>
class Vector3 {
  public:
    // Vector3 Public Methods
    T operator[](int i) const {
        DCHECK(i >= 0 && i <= 2);
        if (i == 0) return x;
        if (i == 1) return y;
        return z;
    }
    T &operator[](int i) {
        DCHECK(i >= 0 && i <= 2);
        if (i == 0) return x;
        if (i == 1) return y;
        return z;
    }
    Vector3() { x = y = z = 0; }
    Vector3(T x, T y, T z) : x(x), y(y), z(z) { DCHECK(!HasNaNs()); }
    bool HasNaNs() const { return isNaN(x) || isNaN(y) || isNaN(z); }
    explicit Vector3(const Point3<T> &p);
#ifndef NDEBUG
    // The default versions of these are fine for release builds; for debug
    // we define them so that we can add the Assert checks.
    Vector3(const Vector3<T> &v) {
        DCHECK(!v.HasNaNs());
        x = v.x;
        y = v.y;
        z = v.z;
    }

    Vector3<T> &operator=(const Vector3<T> &v) {
        DCHECK(!v.HasNaNs());
        x = v.x;
        y = v.y;
        z = v.z;
        return *this;
    }
#endif  // !NDEBUG
    Vector3<T> operator+(const Vector3<T> &v) const {
        DCHECK(!v.HasNaNs());
        return Vector3(x + v.x, y + v.y, z + v.z);
    }
    Vector3<T> &operator+=(const Vector3<T> &v) {
        DCHECK(!v.HasNaNs());
        x += v.x;
        y += v.y;
        z += v.z;
        return *this;
    }
    Vector3<T> operator-(const Vector3<T> &v) const {
        DCHECK(!v.HasNaNs());
        return Vector3(x - v.x, y - v.y, z - v.z);
    }
    Vector3<T> &operator-=(const Vector3<T> &v) {
        DCHECK(!v.HasNaNs());
        x -= v.x;
        y -= v.y;
        z -= v.z;
        return *this;
    }
    bool operator==(const Vector3<T> &v) const {
        return x == v.x && y == v.y && z == v.z;
    }
    bool operator!=(const Vector3<T> &v) const {
        return x != v.x || y != v.y || z != v.z;
    }
    template <typename U>
    Vector3<T> operator*(U s) const {
        return Vector3<T>(s * x, s * y, s * z);
    }
    template <typename U>
    Vector3<T> &operator*=(U s) {
        DCHECK(!isNaN(s));
        x *= s;
        y *= s;
        z *= s;
        return *this;
    }
    template <typename U>
    Vector3<T> operator/(U f) const {
        CHECK_NE(f, 0);
        Float inv = (Float)1 / f;
        return Vector3<T>(x * inv, y * inv, z * inv);
    }

    template <typename U>
    Vector3<T> &operator/=(U f) {
        CHECK_NE(f, 0);
        Float inv = (Float)1 / f;
        x *= inv;
        y *= inv;
        z *= inv;
        return *this;
    }
    Vector3<T> operator-() const { return Vector3<T>(-x, -y, -z); }
    Float LengthSquared() const { return x * x + y * y + z * z; }
    // 用于归一化, 计算两点距离等
    Float Length() const { return std::sqrt(LengthSquared()); }
    explicit Vector3(const Normal3<T> &n);

    // Vector3 Public Data
    T x, y, z;
};

template <typename T>
inline std::ostream &operator<<(std::ostream &os, const Vector3<T> &v) {
    os << "[ " << v.x << ", " << v.y << ", " << v.z << " ]";
    return os;
}

template <>
inline std::ostream &operator<<(std::ostream &os, const Vector3<Float> &v) {
    os << StringPrintf("[ %f, %f, %f ]", v.x, v.y, v.z);
    return os;
}

typedef Vector2<Float> Vector2f;
typedef Vector2<int> Vector2i;
typedef Vector3<Float> Vector3f;
typedef Vector3<int> Vector3i;

// Point Declarations
template <typename T>
class Point2 {
  public:
    // Point2 Public Methods
    // 允许从 Point3 到 Point2 的显式转换
    explicit Point2(const Point3<T> &p) : x(p.x), y(p.y) { DCHECK(!HasNaNs()); }
    Point2() { x = y = 0; }
    Point2(T xx, T yy) : x(xx), y(yy) { DCHECK(!HasNaNs()); }

    template <typename U>
    explicit Point2(const Point2<U> &p) {
        x = (T)p.x;
        y = (T)p.y;
        DCHECK(!HasNaNs());
    }

    template <typename U>
    explicit Point2(const Vector2<U> &p) {
        x = (T)p.x;
        y = (T)p.y;
        DCHECK(!HasNaNs());
    }

    template <typename U>
    explicit operator Vector2<U>() const {
        return Vector2<U>(x, y);
    }

#ifndef NDEBUG
    Point2(const Point2<T> &p) {
        DCHECK(!p.HasNaNs());
        x = p.x;
        y = p.y;
    }

    Point2<T> &operator=(const Point2<T> &p) {
        DCHECK(!p.HasNaNs());
        x = p.x;
        y = p.y;
        return *this;
    }
#endif  // !NDEBUG
    Point2<T> operator+(const Vector2<T> &v) const {
        DCHECK(!v.HasNaNs());
        return Point2<T>(x + v.x, y + v.y);
    }

    Point2<T> &operator+=(const Vector2<T> &v) {
        DCHECK(!v.HasNaNs());
        x += v.x;
        y += v.y;
        return *this;
    }
    Vector2<T> operator-(const Point2<T> &p) const {
        DCHECK(!p.HasNaNs());
        return Vector2<T>(x - p.x, y - p.y);
    }

    Point2<T> operator-(const Vector2<T> &v) const {
        DCHECK(!v.HasNaNs());
        return Point2<T>(x - v.x, y - v.y);
    }
    Point2<T> operator-() const { return Point2<T>(-x, -y); }
    Point2<T> &operator-=(const Vector2<T> &v) {
        DCHECK(!v.HasNaNs());
        x -= v.x;
        y -= v.y;
        return *this;
    }
    Point2<T> &operator+=(const Point2<T> &p) {
        DCHECK(!p.HasNaNs());
        x += p.x;
        y += p.y;
        return *this;
    }
    Point2<T> operator+(const Point2<T> &p) const {
        DCHECK(!p.HasNaNs());
        return Point2<T>(x + p.x, y + p.y);
    }
    template <typename U>
    Point2<T> operator*(U f) const {
        return Point2<T>(f * x, f * y);
    }
    template <typename U>
    Point2<T> &operator*=(U f) {
        x *= f;
        y *= f;
        return *this;
    }
    template <typename U>
    Point2<T> operator/(U f) const {
        CHECK_NE(f, 0);
        Float inv = (Float)1 / f;
        return Point2<T>(inv * x, inv * y);
    }
    template <typename U>
    Point2<T> &operator/=(U f) {
        CHECK_NE(f, 0);
        Float inv = (Float)1 / f;
        x *= inv;
        y *= inv;
        return *this;
    }
    T operator[](int i) const {
        DCHECK(i >= 0 && i <= 1);
        if (i == 0) return x;
        return y;
    }

    T &operator[](int i) {
        DCHECK(i >= 0 && i <= 1);
        if (i == 0) return x;
        return y;
    }
    bool operator==(const Point2<T> &p) const { return x == p.x && y == p.y; }
    bool operator!=(const Point2<T> &p) const { return x != p.x || y != p.y; }
    bool HasNaNs() const { return isNaN(x) || isNaN(y); }

    // Point2 Public Data
    T x, y;
};

template <typename T>
inline std::ostream &operator<<(std::ostream &os, const Point2<T> &v) {
    os << "[ " << v.x << ", " << v.y << " ]";
    return os;
}

template <>
inline std::ostream &operator<<(std::ostream &os, const Point2<Float> &v) {
    os << StringPrintf("[ %f, %f ]", v.x, v.y);
    return os;
}

template <typename T>
class Point3 {
  public:
    // Point3 Public Methods
    Point3() { x = y = z = 0; }
    Point3(T x, T y, T z) : x(x), y(y), z(z) { DCHECK(!HasNaNs()); }
    //允许从 Point3<U> 到 Point3<T> 的显式转化
    template <typename U>
    explicit Point3(const Point3<U> &p)
        : x((T)p.x), y((T)p.y), z((T)p.z) {
        DCHECK(!HasNaNs());
    }
    template <typename U>
    explicit operator Vector3<U>() const {
        return Vector3<U>(x, y, z);
    }
#ifndef NDEBUG
    Point3(const Point3<T> &p) {
        DCHECK(!p.HasNaNs());
        x = p.x;
        y = p.y;
        z = p.z;
    }

    Point3<T> &operator=(const Point3<T> &p) {
        DCHECK(!p.HasNaNs());
        x = p.x;
        y = p.y;
        z = p.z;
        return *this;
    }
#endif  // !NDEBUG
    // Point ± Vector, Point - Point 的操作要有趣一点
    Point3<T> operator+(const Vector3<T> &v) const {
        DCHECK(!v.HasNaNs());
        return Point3<T>(x + v.x, y + v.y, z + v.z);
    }
    Point3<T> &operator+=(const Vector3<T> &v) {
        DCHECK(!v.HasNaNs());
        x += v.x;
        y += v.y;
        z += v.z;
        return *this;
    }
    Vector3<T> operator-(const Point3<T> &p) const {
        DCHECK(!p.HasNaNs());
        return Vector3<T>(x - p.x, y - p.y, z - p.z);
    }
    Point3<T> operator-(const Vector3<T> &v) const {
        DCHECK(!v.HasNaNs());
        return Point3<T>(x - v.x, y - v.y, z - v.z);
    }
    Point3<T> &operator-=(const Vector3<T> &v) {
        DCHECK(!v.HasNaNs());
        x -= v.x;
        y -= v.y;
        z -= v.z;
        return *this;
    }
    // 虽然 Point + Point 在数学上是没有意义的, 但计算其权重(weight)之和是有意义的
    Point3<T> &operator+=(const Point3<T> &p) {
        DCHECK(!p.HasNaNs());
        x += p.x;
        y += p.y;
        z += p.z;
        return *this;
    }
    Point3<T> operator+(const Point3<T> &p) const {
        DCHECK(!p.HasNaNs());
        return Point3<T>(x + p.x, y + p.y, z + p.z);
    }
    template <typename U>
    Point3<T> operator*(U f) const {
        return Point3<T>(f * x, f * y, f * z);
    }
    template <typename U>
    Point3<T> &operator*=(U f) {
        x *= f;
        y *= f;
        z *= f;
        return *this;
    }
    template <typename U>
    Point3<T> operator/(U f) const {
        CHECK_NE(f, 0);
        Float inv = (Float)1 / f;
        return Point3<T>(inv * x, inv * y, inv * z);
    }
    template <typename U>
    Point3<T> &operator/=(U f) {
        CHECK_NE(f, 0);
        Float inv = (Float)1 / f;
        x *= inv;
        y *= inv;
        z *= inv;
        return *this;
    }
    T operator[](int i) const {
        DCHECK(i >= 0 && i <= 2);
        if (i == 0) return x;
        if (i == 1) return y;
        return z;
    }

    T &operator[](int i) {
        DCHECK(i >= 0 && i <= 2);
        if (i == 0) return x;
        if (i == 1) return y;
        return z;
    }
    bool operator==(const Point3<T> &p) const {
        return x == p.x && y == p.y && z == p.z;
    }
    bool operator!=(const Point3<T> &p) const {
        return x != p.x || y != p.y || z != p.z;
    }
    bool HasNaNs() const { return isNaN(x) || isNaN(y) || isNaN(z); }
    Point3<T> operator-() const { return Point3<T>(-x, -y, -z); }

    // Point3 Public Data
    T x, y, z;
};

template <typename T>
inline std::ostream &operator<<(std::ostream &os, const Point3<T> &v) {
    os << "[ " << v.x << ", " << v.y << ", " << v.z << " ]";
    return os;
}

template <>
inline std::ostream &operator<<(std::ostream &os, const Point3<Float> &v) {
    os << StringPrintf("[ %f, %f, %f ]", v.x, v.y, v.z);
    return os;
}

typedef Point2<Float> Point2f;
typedef Point2<int> Point2i;
typedef Point3<Float> Point3f;
typedef Point3<int> Point3i;

// Normal Declarations
// A surface normal (or just normal) is a vector that is perpendicular to a surface at a particular position. 
// It can be defined as the cross product of any two nonparallel vectors that are tangent to the surface at a point.
// because normals are defined in terms of their relationship to a particular surface, they behave differently than vectors in some situations, particularly when applying transformations.
// a normal cannot be added to a point, and one cannot take the cross product of two normals. 
// Note that, in an unfortunate turn of terminology, normals are not necessarily normalized.
template <typename T>
class Normal3 {
  public:
    // Normal3 Public Methods
    Normal3() { x = y = z = 0; }
    Normal3(T xx, T yy, T zz) : x(xx), y(yy), z(zz) { DCHECK(!HasNaNs()); }
    Normal3<T> operator-() const { return Normal3(-x, -y, -z); }
    Normal3<T> operator+(const Normal3<T> &n) const {
        DCHECK(!n.HasNaNs());
        return Normal3<T>(x + n.x, y + n.y, z + n.z);
    }

    Normal3<T> &operator+=(const Normal3<T> &n) {
        DCHECK(!n.HasNaNs());
        x += n.x;
        y += n.y;
        z += n.z;
        return *this;
    }
    Normal3<T> operator-(const Normal3<T> &n) const {
        DCHECK(!n.HasNaNs());
        return Normal3<T>(x - n.x, y - n.y, z - n.z);
    }

    Normal3<T> &operator-=(const Normal3<T> &n) {
        DCHECK(!n.HasNaNs());
        x -= n.x;
        y -= n.y;
        z -= n.z;
        return *this;
    }
    bool HasNaNs() const { return isNaN(x) || isNaN(y) || isNaN(z); }
    template <typename U>
    Normal3<T> operator*(U f) const {
        return Normal3<T>(f * x, f * y, f * z);
    }

    template <typename U>
    Normal3<T> &operator*=(U f) {
        x *= f;
        y *= f;
        z *= f;
        return *this;
    }
    template <typename U>
    Normal3<T> operator/(U f) const {
        CHECK_NE(f, 0);
        Float inv = (Float)1 / f;
        return Normal3<T>(x * inv, y * inv, z * inv);
    }

    template <typename U>
    Normal3<T> &operator/=(U f) {
        CHECK_NE(f, 0);
        Float inv = (Float)1 / f;
        x *= inv;
        y *= inv;
        z *= inv;
        return *this;
    }
    Float LengthSquared() const { return x * x + y * y + z * z; }
    Float Length() const { return std::sqrt(LengthSquared()); }

#ifndef NDEBUG
    Normal3<T>(const Normal3<T> &n) {
        DCHECK(!n.HasNaNs());
        x = n.x;
        y = n.y;
        z = n.z;
    }

    Normal3<T> &operator=(const Normal3<T> &n) {
        DCHECK(!n.HasNaNs());
        x = n.x;
        y = n.y;
        z = n.z;
        return *this;
    }
#endif  // !NDEBUG
    explicit Normal3<T>(const Vector3<T> &v) : x(v.x), y(v.y), z(v.z) {
        DCHECK(!v.HasNaNs());
    }
    bool operator==(const Normal3<T> &n) const {
        return x == n.x && y == n.y && z == n.z;
    }
    bool operator!=(const Normal3<T> &n) const {
        return x != n.x || y != n.y || z != n.z;
    }

    T operator[](int i) const {
        DCHECK(i >= 0 && i <= 2);
        if (i == 0) return x;
        if (i == 1) return y;
        return z;
    }

    T &operator[](int i) {
        DCHECK(i >= 0 && i <= 2);
        if (i == 0) return x;
        if (i == 1) return y;
        return z;
    }

    // Normal3 Public Data
    T x, y, z;
};

template <typename T>
inline std::ostream &operator<<(std::ostream &os, const Normal3<T> &v) {
    os << "[ " << v.x << ", " << v.y << ", " << v.z << " ]";
    return os;
}

template <>
inline std::ostream &operator<<(std::ostream &os, const Normal3<Float> &v) {
    os << StringPrintf("[ %f, %f, %f ]", v.x, v.y, v.z);
    return os;
}

typedef Normal3<Float> Normal3f;

// PBRT 用的包围盒是轴对齐包围盒 AABB(axis-aligned bounding boxes, section2.6), 而非有向包围盒 OBB(oriented bounding boxes)
// AABB 可以用一个顶点加三个长度, 或者两个顶点的方式实现, PBRT 用的是后者

// Bounds2 主要是用在划分胶片平面(File, section7.9)上
// Bounds Declarations
template <typename T>
class Bounds2 {
  public:
    // Bounds2 Public Methods
    Bounds2() {
        T minNum = std::numeric_limits<T>::lowest();
        T maxNum = std::numeric_limits<T>::max();
        pMin = Point2<T>(maxNum, maxNum);
        pMax = Point2<T>(minNum, minNum);
    }
    explicit Bounds2(const Point2<T> &p) : pMin(p), pMax(p) {}
    Bounds2(const Point2<T> &p1, const Point2<T> &p2) {
        pMin = Point2<T>(std::min(p1.x, p2.x), std::min(p1.y, p2.y));
        pMax = Point2<T>(std::max(p1.x, p2.x), std::max(p1.y, p2.y));
    }
    template <typename U>
    explicit operator Bounds2<U>() const {
        return Bounds2<U>((Point2<U>)pMin, (Point2<U>)pMax);
    }

    Vector2<T> Diagonal() const { return pMax - pMin; }
    T Area() const {
        Vector2<T> d = pMax - pMin;
        return (d.x * d.y);
    }
    int MaximumExtent() const {
        Vector2<T> diag = Diagonal();
        if (diag.x > diag.y)
            return 0;
        else
            return 1;
    }
    inline const Point2<T> &operator[](int i) const {
        DCHECK(i == 0 || i == 1);
        return (i == 0) ? pMin : pMax;
    }
    inline Point2<T> &operator[](int i) {
        DCHECK(i == 0 || i == 1);
        return (i == 0) ? pMin : pMax;
    }
    bool operator==(const Bounds2<T> &b) const {
        return b.pMin == pMin && b.pMax == pMax;
    }
    bool operator!=(const Bounds2<T> &b) const {
        return b.pMin != pMin || b.pMax != pMax;
    }
    Point2<T> Lerp(const Point2f &t) const {
        return Point2<T>(pbrt::Lerp(t.x, pMin.x, pMax.x),
                         pbrt::Lerp(t.y, pMin.y, pMax.y));
    }
    Vector2<T> Offset(const Point2<T> &p) const {
        Vector2<T> o = p - pMin;
        if (pMax.x > pMin.x) o.x /= pMax.x - pMin.x;
        if (pMax.y > pMin.y) o.y /= pMax.y - pMin.y;
        return o;
    }
    void BoundingSphere(Point2<T> *c, Float *rad) const {
        *c = (pMin + pMax) / 2;
        *rad = Inside(*c, *this) ? Distance(*c, pMax) : 0;
    }
    friend std::ostream &operator<<(std::ostream &os, const Bounds2<T> &b) {
        os << "[ " << b.pMin << " - " << b.pMax << " ]";
        return os;
    }

    // Bounds2 Public Data
    Point2<T> pMin, pMax;
};

// Bounds3 主要是用在空间加速结构, 即层次包围盒(bounding volume hierarchy, section4.3)上
template <typename T>
class Bounds3 {
  public:
    // Bounds3 Public Methods
    Bounds3() {
        T minNum = std::numeric_limits<T>::lowest();
        T maxNum = std::numeric_limits<T>::max();
        // 这种"错误"的设置使得之后所有的操作(如 Union() )都是正确的
        pMin = Point3<T>(maxNum, maxNum, maxNum);
        pMax = Point3<T>(minNum, minNum, minNum);
    }
    explicit Bounds3(const Point3<T> &p) : pMin(p), pMax(p) {}
    Bounds3(const Point3<T> &p1, const Point3<T> &p2)
        : pMin(std::min(p1.x, p2.x), std::min(p1.y, p2.y),
               std::min(p1.z, p2.z)),
          pMax(std::max(p1.x, p2.x), std::max(p1.y, p2.y),
               std::max(p1.z, p2.z)) {}
    const Point3<T> &operator[](int i) const;
    Point3<T> &operator[](int i);
    bool operator==(const Bounds3<T> &b) const {
        return b.pMin == pMin && b.pMax == pMax;
    }
    bool operator!=(const Bounds3<T> &b) const {
        return b.pMin != pMin || b.pMax != pMax;
    }
    // 按索引 i(0~7) 返回整个包围盒的八个角之一
    // PBRT 使用左手系
    // TODO
    Point3<T> Corner(int corner) const {
        DCHECK(corner >= 0 && corner < 8);
        return Point3<T>((*this)[(corner & 1)].x,
                         (*this)[(corner & 2) ? 1 : 0].y,
                         (*this)[(corner & 4) ? 1 : 0].z);
    }
    // 对角线
    Vector3<T> Diagonal() const { return pMax - pMin; }
    // 六个面的表面积
    T SurfaceArea() const {
        Vector3<T> d = Diagonal();
        return 2 * (d.x * d.y + d.x * d.z + d.y * d.z);
    }
    // 体积
    T Volume() const {
        Vector3<T> d = Diagonal();
        return d.x * d.y * d.z;
    }
    // 最长的一条边
    int MaximumExtent() const {
        Vector3<T> d = Diagonal();
        if (d.x > d.y && d.x > d.z)
            return 0;
        else if (d.y > d.z)
            return 1;
        else
            return 2;
    }
    // 返回从 pMin 到 pMax 的插值结果
    Point3<T> Lerp(const Point3f &t) const {
        return Point3<T>(pbrt::Lerp(t.x, pMin.x, pMax.x),
                         pbrt::Lerp(t.y, pMin.y, pMax.y),
                         pbrt::Lerp(t.z, pMin.z, pMax.z));
    }
    
    // 返回点 p 相对于 pMin, pMax 的偏移量
    // 假设点 p 在直线(pMin, pMax)上, 当 p == pMin 时, 返回(0, 0, 0); 当 p == pMax 时, 返回(1, 1, 1)
    Vector3<T> Offset(const Point3<T> &p) const {
        Vector3<T> o = p - pMin;
        if (pMax.x > pMin.x) o.x /= pMax.x - pMin.x;
        if (pMax.y > pMin.y) o.y /= pMax.y - pMin.y;
        if (pMax.z > pMin.z) o.z /= pMax.z - pMin.z;
        return o;
    }
    // 返回包含 bounds3 的球体
    void BoundingSphere(Point3<T> *center, Float *radius) const {
        *center = (pMin + pMax) / 2;
        *radius = Inside(*center, *this) ? Distance(*center, pMax) : 0;
    }
    template <typename U>
    explicit operator Bounds3<U>() const {
        return Bounds3<U>((Point3<U>)pMin, (Point3<U>)pMax);
    }
    bool IntersectP(const Ray &ray, Float *hitt0 = nullptr,
                    Float *hitt1 = nullptr) const;
    inline bool IntersectP(const Ray &ray, const Vector3f &invDir,
                           const int dirIsNeg[3]) const;
    friend std::ostream &operator<<(std::ostream &os, const Bounds3<T> &b) {
        os << "[ " << b.pMin << " - " << b.pMax << " ]";
        return os;
    }

    // Bounds3 Public Data
    Point3<T> pMin, pMax;
};

typedef Bounds2<Float> Bounds2f;
typedef Bounds2<int> Bounds2i;
typedef Bounds3<Float> Bounds3f;
typedef Bounds3<int> Bounds3i;

// 当使用如 for (Point2i pixel : tileBounds) 这样的范围 for 循环时需要相应的迭代器支持
// 在遍历胶片平面 imageTile 时会用到
class Bounds2iIterator : public std::forward_iterator_tag {
  public:
    Bounds2iIterator(const Bounds2i &b, const Point2i &pt)
        : p(pt), bounds(&b) {}
    Bounds2iIterator operator++() {
        advance();
        return *this;
    }
    Bounds2iIterator operator++(int) {
        Bounds2iIterator old = *this;
        advance();
        return old;
    }
    bool operator==(const Bounds2iIterator &bi) const {
        return p == bi.p && bounds == bi.bounds;
    }
    bool operator!=(const Bounds2iIterator &bi) const {
        return p != bi.p || bounds != bi.bounds;
    }

    Point2i operator*() const { return p; }

  private:
    void advance() {
        ++p.x;
        if (p.x == bounds->pMax.x) {
            p.x = bounds->pMin.x;
            ++p.y;
        }
    }
    Point2i p;
    const Bounds2i *bounds;
};

// Ray Declarations
// $ \mathrm{r}(t)=\mathrm{o} + t \mathrm{d} \quad 0 \leq t < \infty $
class Ray {
  public:
    // Ray Public Methods
    Ray() : tMax(Infinity), time(0.f), medium(nullptr) {}

    Ray(const Point3f &o, const Vector3f &d, Float tMax = Infinity,
        Float time = 0.f, const Medium *medium = nullptr)
        : o(o), d(d), tMax(tMax), time(time), medium(medium) {}

    // 给定 t, 将光线从位置 Origin 沿方向 Direction 移动 t 个单位到新位置
    Point3f operator()(Float t) const { return o + d * t; }
    bool HasNaNs() const { return (o.HasNaNs() || d.HasNaNs() || isNaN(tMax)); }

    friend std::ostream &operator<<(std::ostream &os, const Ray &r) {
        os << "[o=" << r.o << ", d=" << r.d << ", tMax=" << r.tMax
           << ", time=" << r.time << "]";
        return os;
    }

    // Ray Public Data
    // 起点和方向, 为了方便表示成 o & d
    Point3f o;
    // d 是单位向量吗???
    Vector3f d;
    // as parameters to ray–object intersection testing routines, which will record the offsets to the closest intersection in tMax.
    // tMax 代表了 t 的最大值, 限制了 ray 的范围
    mutable Float tMax;
    // Each ray has a time value associated with it. In scenes with animated objects, the rendering system constructs a representation of the scene at the appropriate time for each ray
    // 当场景包含运动物体时, 系统会通过 time 为这条光线构造一个合适的场景表示
    Float time;
    // Finally, each ray records the medium containing its origin. The Medium class, introduced in Section 11.3, encapsulates the (potentially spatially varying) properties of media such as a foggy atmosphere, smoke, or scattering liquids like milk or shampoo. Associating this information with rays makes it possible for other parts of the system to account correctly for the effect of rays passing from one medium to another.
    // 
    // const Medium* 意味着只访问资源, 不管理 medium 的生命周期
    const Medium *medium;
};

// 光线微分, 主要用于纹理采样时的反走样(anti-aliasing)操作
// 需要先理解纹理采样, 走样, 过滤, 反走样等概念, 否则可以先跳过 RayDifferential
// 详细解释可以参考《全局光照技术》6.5.1节和 https://blog.csdn.net/suian0424/article/details/81485563
class RayDifferential : public Ray {
  public:
    // RayDifferential Public Methods
    RayDifferential() { hasDifferentials = false; }
    RayDifferential(const Point3f &o, const Vector3f &d, Float tMax = Infinity,
                    Float time = 0.f, const Medium *medium = nullptr)
        : Ray(o, d, tMax, time, medium) {
        hasDifferentials = false;
    }
    RayDifferential(const Ray &ray) : Ray(ray) { hasDifferentials = false; }

    bool HasNaNs() const {
        return Ray::HasNaNs() ||
               (hasDifferentials &&
                (rxOrigin.HasNaNs() || ryOrigin.HasNaNs() ||
                 rxDirection.HasNaNs() || ryDirection.HasNaNs()));
    }

    // 通过估计样本距离 $\mathrm{s}$ 对辅助光线进行(更大幅度的)偏移, 以...
    void ScaleDifferentials(Float s) {
        rxOrigin = o + (rxOrigin - o) * s;
        ryOrigin = o + (ryOrigin - o) * s;
        rxDirection = d + (rxDirection - d) * s;
        ryDirection = d + (ryDirection - d) * s;
    }

    friend std::ostream &operator<<(std::ostream &os, const RayDifferential &r) {
        os << "[ " << (Ray &)r << " has differentials: " <<
            (r.hasDifferentials ? "true" : "false") << ", xo = " << r.rxOrigin <<
            ", xd = " << r.rxDirection << ", yo = " << r.ryOrigin << ", yd = " <<
            r.ryDirection;
        return os;
    }

    // RayDifferential Public Data
    bool hasDifferentials;
    // 主光线有两条辅助光线
    // These extra rays represent camera rays offset by one sample in the $x$ and $y$ direction from the main ray on the film plane.
    // 它们是通过胶片平面上主光线的起点, 向右方和上方分别偏移一个像素单位得到的
    // By determining the area that these three rays project to on an object being shaded, the Texture can estimate an area to average over for proper antialiasing.
    // 通过确定三条光线所投射到的对象上着色区域的大小, 纹理对象可以估计出用于反走样的区域大小(一般是取这个区域上贴图颜色的均值)
    Point3f rxOrigin, ryOrigin;
    Vector3f rxDirection, ryDirection;
};

// Geometry Inline Functions
// 都是先高维版本(Vector3, Point3...), 再低维版本(Vector2, Point2...)

// Vector

template <typename T>
inline Vector3<T>::Vector3(const Point3<T> &p)
    : x(p.x), y(p.y), z(p.z) {
    DCHECK(!HasNaNs());
}

template <typename T, typename U>
inline Vector3<T> operator*(U s, const Vector3<T> &v) {
    return v * s;
}
template <typename T>
Vector3<T> Abs(const Vector3<T> &v) {
    return Vector3<T>(std::abs(v.x), std::abs(v.y), std::abs(v.z));
}

// $$ \mathbf{v} \cdot \mathbf{w} = \|\mathbf{v}\| \|\mathbf{w}\| \cos \theta   \tag{2.1} $$
// It immediately follows from Equation (2.1) that if $v$ and $w$ are unit vectors, their dot product is the cosine of the angle between them.
template <typename T>
inline T Dot(const Vector3<T> &v1, const Vector3<T> &v2) {
    DCHECK(!v1.HasNaNs() && !v2.HasNaNs());
    return v1.x * v2.x + v1.y * v2.y + v1.z * v2.z;
}

template <typename T>
inline T AbsDot(const Vector3<T> &v1, const Vector3<T> &v2) {
    DCHECK(!v1.HasNaNs() && !v2.HasNaNs());
    return std::abs(Dot(v1, v2));
}

// the (result of) cross product $\|\mathbf{v} \times \mathbf{w}\|$ is a vector that is perpendicular to both of them. 
// Given orthogonal vectors $v$ and $w$, then $\|\mathbf{v} \times \mathbf{w}\|$ is defined to be a vector such that $(\mathbf{V}, \mathbf{W}, \mathbf{V} \times \mathbf{W})$ form an orthogonal coordinate system

// $$ \|\mathbf{v} \times \mathbf{w}\| = \|\mathbf{v}\| \|\mathbf{w}\| |\sin \theta| $$
// An important implication of this is that the cross product of two perpendicular unit vectors is itself a unit vector. 
// Note also that the result of the cross product is a degenerate vector if $v$ and $w$ are parallel.

template <typename T>
inline Vector3<T> Cross(const Vector3<T> &v1, const Vector3<T> &v2) {
    DCHECK(!v1.HasNaNs() && !v2.HasNaNs());
    // Using extra precision for 32-bit floating-point values here protects against error from catastrophic cancellation, a type of floating-point error that can happen when subtracting two values that are very close together.
    // 在这里使用 double 的拓展精度可以避免一种灾难性的错误, 即在减去两个非常接近的值时可能发生的浮点错误
    double v1x = v1.x, v1y = v1.y, v1z = v1.z;
    double v2x = v2.x, v2y = v2.y, v2z = v2.z;
    //  ( y1 * z2 - z1 * y2 )
    // -( x1 * z2 - z1 * x2 ) 
    //  ( x1 * y2 - y1 * x2 )
    return Vector3<T>((v1y * v2z) - (v1z * v2y), (v1z * v2x) - (v1x * v2z),
                      (v1x * v2y) - (v1y * v2x));
}

template <typename T>
inline Vector3<T> Cross(const Vector3<T> &v1, const Normal3<T> &v2) {
    DCHECK(!v1.HasNaNs() && !v2.HasNaNs());
    double v1x = v1.x, v1y = v1.y, v1z = v1.z;
    double v2x = v2.x, v2y = v2.y, v2z = v2.z;
    return Vector3<T>((v1y * v2z) - (v1z * v2y), (v1z * v2x) - (v1x * v2z),
                      (v1x * v2y) - (v1y * v2x));
}

template <typename T>
inline Vector3<T> Cross(const Normal3<T> &v1, const Vector3<T> &v2) {
    DCHECK(!v1.HasNaNs() && !v2.HasNaNs());
    double v1x = v1.x, v1y = v1.y, v1z = v1.z;
    double v2x = v2.x, v2y = v2.y, v2z = v2.z;
    return Vector3<T>((v1y * v2z) - (v1z * v2y), (v1z * v2x) - (v1x * v2z),
                      (v1x * v2y) - (v1y * v2x));
}

template <typename T>
inline Vector3<T> Normalize(const Vector3<T> &v) {
    return v / v.Length();
}

template <typename T>
T MinComponent(const Vector3<T> &v) {
    return std::min(v.x, std::min(v.y, v.z));
}

template <typename T>
T MaxComponent(const Vector3<T> &v) {
    return std::max(v.x, std::max(v.y, v.z));
}

// 返回最大维度的索引, 用 0, 1, 2 指代 x, y, z
template <typename T>
int MaxDimension(const Vector3<T> &v) {
    return (v.x > v.y) ? ((v.x > v.z) ? 0 : 2) : ((v.y > v.z) ? 1 : 2);
}

// 逐分量的 min, max 操作
template <typename T>
Vector3<T> Min(const Vector3<T> &p1, const Vector3<T> &p2) {
    return Vector3<T>(std::min(p1.x, p2.x), std::min(p1.y, p2.y),
                      std::min(p1.z, p2.z));
}

template <typename T>
Vector3<T> Max(const Vector3<T> &p1, const Vector3<T> &p2) {
    return Vector3<T>(std::max(p1.x, p2.x), std::max(p1.y, p2.y),
                      std::max(p1.z, p2.z));
}

// 重排列
// poor swizzle(e.g. v.xxy, v.zzz)
// TODO: fay::vec3().(_x, _y, _x)
template <typename T>
Vector3<T> Permute(const Vector3<T> &v, int x, int y, int z) {
    return Vector3<T>(v[x], v[y], v[z]);
}

// Because the cross product of two vectors is orthogonal to both, we can apply the cross product two times to get a set of three orthogonal vectors for the coordinate system
// Note that the two vectors generated by this technique are unique only up to a rotation about the given vector.
// 对于给定的 v1, 只有在限定相对它的旋转角度时, 生成的 v2, v3 才是唯一的
template <typename T>
inline void CoordinateSystem(const Vector3<T> &v1, Vector3<T> *v2,
                             Vector3<T> *v3) {
    if (std::abs(v1.x) > std::abs(v1.y))
        // 该实现假定 v1 已经是归一化的
        // 通过 Vector3<T>(-v1.z, 0, v1.x) 构造出垂直于 v1 的 v2, 然后对 v2 归一化
        *v2 = Vector3<T>(-v1.z, 0, v1.x) / std::sqrt(v1.x * v1.x + v1.z * v1.z);
    else
        *v2 = Vector3<T>(0, v1.z, -v1.y) / std::sqrt(v1.y * v1.y + v1.z * v1.z);

    *v3 = Cross(v1, *v2);
}

template <typename T>
Vector2<T>::Vector2(const Point2<T> &p)
    : x(p.x), y(p.y) {
    DCHECK(!HasNaNs());
}

template <typename T>
Vector2<T>::Vector2(const Point3<T> &p)
    : x(p.x), y(p.y) {
    DCHECK(!HasNaNs());
}

template <typename T, typename U>
inline Vector2<T> operator*(U f, const Vector2<T> &v) {
    return v * f;
}
template <typename T>
inline Float Dot(const Vector2<T> &v1, const Vector2<T> &v2) {
    DCHECK(!v1.HasNaNs() && !v2.HasNaNs());
    return v1.x * v2.x + v1.y * v2.y;
}

template <typename T>
inline Float AbsDot(const Vector2<T> &v1, const Vector2<T> &v2) {
    DCHECK(!v1.HasNaNs() && !v2.HasNaNs());
    return std::abs(Dot(v1, v2));
}

template <typename T>
inline Vector2<T> Normalize(const Vector2<T> &v) {
    return v / v.Length();
}
template <typename T>
Vector2<T> Abs(const Vector2<T> &v) {
    return Vector2<T>(std::abs(v.x), std::abs(v.y));
}

// Point

template <typename T>
inline Float Distance(const Point3<T> &p1, const Point3<T> &p2) {
    return (p1 - p2).Length();
}

template <typename T>
inline Float DistanceSquared(const Point3<T> &p1, const Point3<T> &p2) {
    return (p1 - p2).LengthSquared();
}

template <typename T, typename U>
inline Point3<T> operator*(U f, const Point3<T> &p) {
    DCHECK(!p.HasNaNs());
    return p * f;
}

// 计算两点的线性插值
template <typename T>
Point3<T> Lerp(Float t, const Point3<T> &p0, const Point3<T> &p1) {
    return (1 - t) * p0 + t * p1;
}

// 以下五个操作都是逐分量的
template <typename T>
Point3<T> Min(const Point3<T> &p1, const Point3<T> &p2) {
    return Point3<T>(std::min(p1.x, p2.x), std::min(p1.y, p2.y),
                     std::min(p1.z, p2.z));
}

template <typename T>
Point3<T> Max(const Point3<T> &p1, const Point3<T> &p2) {
    return Point3<T>(std::max(p1.x, p2.x), std::max(p1.y, p2.y),
                     std::max(p1.z, p2.z));
}

template <typename T>
Point3<T> Floor(const Point3<T> &p) {
    return Point3<T>(std::floor(p.x), std::floor(p.y), std::floor(p.z));
}

template <typename T>
Point3<T> Ceil(const Point3<T> &p) {
    return Point3<T>(std::ceil(p.x), std::ceil(p.y), std::ceil(p.z));
}

template <typename T>
Point3<T> Abs(const Point3<T> &p) {
    return Point3<T>(std::abs(p.x), std::abs(p.y), std::abs(p.z));
}

template <typename T>
inline Float Distance(const Point2<T> &p1, const Point2<T> &p2) {
    return (p1 - p2).Length();
}

template <typename T>
inline Float DistanceSquared(const Point2<T> &p1, const Point2<T> &p2) {
    return (p1 - p2).LengthSquared();
}

template <typename T, typename U>
inline Point2<T> operator*(U f, const Point2<T> &p) {
    DCHECK(!p.HasNaNs());
    return p * f;
}

template <typename T>
Point2<T> Floor(const Point2<T> &p) {
    return Point2<T>(std::floor(p.x), std::floor(p.y));
}

template <typename T>
Point2<T> Ceil(const Point2<T> &p) {
    return Point2<T>(std::ceil(p.x), std::ceil(p.y));
}

template <typename T>
Point2<T> Lerp(Float t, const Point2<T> &v0, const Point2<T> &v1) {
    return (1 - t) * v0 + t * v1;
}

template <typename T>
Point2<T> Min(const Point2<T> &pa, const Point2<T> &pb) {
    return Point2<T>(std::min(pa.x, pb.x), std::min(pa.y, pb.y));
}

template <typename T>
Point2<T> Max(const Point2<T> &pa, const Point2<T> &pb) {
    return Point2<T>(std::max(pa.x, pb.x), std::max(pa.y, pb.y));
}

template <typename T>
Point3<T> Permute(const Point3<T> &p, int x, int y, int z) {
    return Point3<T>(p[x], p[y], p[z]);
}

// Normal3

template <typename T, typename U>
inline Normal3<T> operator*(U f, const Normal3<T> &n) {
    return Normal3<T>(f * n.x, f * n.y, f * n.z);
}

template <typename T>
inline Normal3<T> Normalize(const Normal3<T> &n) {
    return n / n.Length();
}

template <typename T>
inline Vector3<T>::Vector3(const Normal3<T> &n)
    : x(n.x), y(n.y), z(n.z) {
    DCHECK(!n.HasNaNs());
}

template <typename T>
inline T Dot(const Normal3<T> &n1, const Vector3<T> &v2) {
    DCHECK(!n1.HasNaNs() && !v2.HasNaNs());
    return n1.x * v2.x + n1.y * v2.y + n1.z * v2.z;
}

template <typename T>
inline T Dot(const Vector3<T> &v1, const Normal3<T> &n2) {
    DCHECK(!v1.HasNaNs() && !n2.HasNaNs());
    return v1.x * n2.x + v1.y * n2.y + v1.z * n2.z;
}

template <typename T>
inline T Dot(const Normal3<T> &n1, const Normal3<T> &n2) {
    DCHECK(!n1.HasNaNs() && !n2.HasNaNs());
    return n1.x * n2.x + n1.y * n2.y + n1.z * n2.z;
}

template <typename T>
inline T AbsDot(const Normal3<T> &n1, const Vector3<T> &v2) {
    DCHECK(!n1.HasNaNs() && !v2.HasNaNs());
    return std::abs(n1.x * v2.x + n1.y * v2.y + n1.z * v2.z);
}

template <typename T>
inline T AbsDot(const Vector3<T> &v1, const Normal3<T> &n2) {
    DCHECK(!v1.HasNaNs() && !n2.HasNaNs());
    return std::abs(v1.x * n2.x + v1.y * n2.y + v1.z * n2.z);
}

template <typename T>
inline T AbsDot(const Normal3<T> &n1, const Normal3<T> &n2) {
    DCHECK(!n1.HasNaNs() && !n2.HasNaNs());
    return std::abs(n1.x * n2.x + n1.y * n2.y + n1.z * n2.z);
}

// Faceforward 执行法线和向量的翻转操作
// 注意第一个参数是要翻转的参数, 第二个是指向的目标
template <typename T>
inline Normal3<T> Faceforward(const Normal3<T> &n, const Vector3<T> &v) {
    return (Dot(n, v) < 0.f) ? -n : n;
}

template <typename T>
inline Normal3<T> Faceforward(const Normal3<T> &n, const Normal3<T> &n2) {
    return (Dot(n, n2) < 0.f) ? -n : n;
}

template <typename T>
inline Vector3<T> Faceforward(const Vector3<T> &v, const Vector3<T> &v2) {
    return (Dot(v, v2) < 0.f) ? -v : v;
}

template <typename T>
inline Vector3<T> Faceforward(const Vector3<T> &v, const Normal3<T> &n2) {
    return (Dot(v, n2) < 0.f) ? -v : v;
}

template <typename T>
Normal3<T> Abs(const Normal3<T> &v) {
    return Normal3<T>(std::abs(v.x), std::abs(v.y), std::abs(v.z));
}

template <typename T>
inline const Point3<T> &Bounds3<T>::operator[](int i) const {
    DCHECK(i == 0 || i == 1);
    return (i == 0) ? pMin : pMax;
}

template <typename T>
inline Point3<T> &Bounds3<T>::operator[](int i) {
    DCHECK(i == 0 || i == 1);
    return (i == 0) ? pMin : pMax;
}

// Bounds

// 返回两者的并集 b∪p
template <typename T>
Bounds3<T> Union(const Bounds3<T> &b, const Point3<T> &p) {
    Bounds3<T> ret;
    ret.pMin = Min(b.pMin, p);
    ret.pMax = Max(b.pMax, p);
    return ret;
}

template <typename T>
Bounds3<T> Union(const Bounds3<T> &b1, const Bounds3<T> &b2) {
    Bounds3<T> ret;
    ret.pMin = Min(b1.pMin, b2.pMin);
    ret.pMax = Max(b1.pMax, b2.pMax);
    return ret;
}

// 返回两者的交集 b ∩ p
template <typename T>
Bounds3<T> Intersect(const Bounds3<T> &b1, const Bounds3<T> &b2) {
    // Important: assign to pMin/pMax directly and don't run the Bounds2()
    // constructor, since it takes min/max of the points passed to it.  In
    // turn, that breaks returning an invalid bound for the case where we
    // intersect non-overlapping bounds (as we'd like to happen).
    Bounds3<T> ret;
    ret.pMin = Max(b1.pMin, b2.pMin);
    ret.pMax = Min(b1.pMax, b2.pMax);
    return ret;
}

// 两个包围盒是否有交叠的部分(看 Figure2.9 会比较容易理解)
template <typename T>
bool Overlaps(const Bounds3<T> &b1, const Bounds3<T> &b2) {
    bool x = (b1.pMax.x >= b2.pMin.x) && (b1.pMin.x <= b2.pMax.x);
    bool y = (b1.pMax.y >= b2.pMin.y) && (b1.pMin.y <= b2.pMax.y);
    bool z = (b1.pMax.z >= b2.pMin.z) && (b1.pMin.z <= b2.pMax.z);
    return (x && y && z);
}

// 点是否在包围盒内部或边界上
template <typename T>
bool Inside(const Point3<T> &p, const Bounds3<T> &b) {
    return (p.x >= b.pMin.x && p.x <= b.pMax.x && p.y >= b.pMin.y &&
            p.y <= b.pMax.y && p.z >= b.pMin.z && p.z <= b.pMax.z);
}

// 点是否在包围盒内部
template <typename T>
bool InsideExclusive(const Point3<T> &p, const Bounds3<T> &b) {
    return (p.x >= b.pMin.x && p.x < b.pMax.x && p.y >= b.pMin.y &&
            p.y < b.pMax.y && p.z >= b.pMin.z && p.z < b.pMax.z);
}

// 向外扩张(或向内收缩)
template <typename T, typename U>
inline Bounds3<T> Expand(const Bounds3<T> &b, U delta) {
    return Bounds3<T>(b.pMin - Vector3<T>(delta, delta, delta),
                      b.pMax + Vector3<T>(delta, delta, delta));
}

// Minimum squared distance from point to box; returns zero if point is
// inside.
// 计算点 p 到包围盒的最小距离的平方值
// 若 p 在包围盒内部, 则返回 0
template <typename T, typename U>
inline Float DistanceSquared(const Point3<T> &p, const Bounds3<U> &b) {
    // TODO: remove fay::max
    Float dx = std::max({Float(0), b.pMin.x - p.x, p.x - b.pMax.x});
    Float dy = std::max({Float(0), b.pMin.y - p.y, p.y - b.pMax.y});
    Float dz = std::max({Float(0), b.pMin.z - p.z, p.z - b.pMax.z});
    return dx * dx + dy * dy + dz * dz;
}

template <typename T, typename U>
inline Float Distance(const Point3<T> &p, const Bounds3<U> &b) {
    return std::sqrt(DistanceSquared(p, b));
}

inline Bounds2iIterator begin(const Bounds2i &b) {
    return Bounds2iIterator(b, b.pMin);
}

inline Bounds2iIterator end(const Bounds2i &b) {
    // Normally, the ending point is at the minimum x value and one past
    // the last valid y value.
    Point2i pEnd(b.pMin.x, b.pMax.y);
    // However, if the bounds are degenerate, override the end point to
    // equal the start point so that any attempt to iterate over the bounds
    // exits out immediately.
    if (b.pMin.x >= b.pMax.x || b.pMin.y >= b.pMax.y)
        pEnd = b.pMin;
    return Bounds2iIterator(b, pEnd);
}

template <typename T>
Bounds2<T> Union(const Bounds2<T> &b, const Point2<T> &p) {
    Bounds2<T> ret;
    ret.pMin = Min(b.pMin, p);
    ret.pMax = Max(b.pMax, p);
    return ret;
}

template <typename T>
Bounds2<T> Union(const Bounds2<T> &b, const Bounds2<T> &b2) {
    Bounds2<T> ret;
    ret.pMin = Min(b.pMin, b2.pMin);
    ret.pMax = Max(b.pMax, b2.pMax);
    return ret;
}

template <typename T>
Bounds2<T> Intersect(const Bounds2<T> &b1, const Bounds2<T> &b2) {
    // Important: assign to pMin/pMax directly and don't run the Bounds2()
    // constructor, since it takes min/max of the points passed to it.  In
    // turn, that breaks returning an invalid bound for the case where we
    // intersect non-overlapping bounds (as we'd like to happen).
    Bounds2<T> ret;
    ret.pMin = Max(b1.pMin, b2.pMin);
    ret.pMax = Min(b1.pMax, b2.pMax);
    return ret;
}

template <typename T>
bool Overlaps(const Bounds2<T> &ba, const Bounds2<T> &bb) {
    bool x = (ba.pMax.x >= bb.pMin.x) && (ba.pMin.x <= bb.pMax.x);
    bool y = (ba.pMax.y >= bb.pMin.y) && (ba.pMin.y <= bb.pMax.y);
    return (x && y);
}

template <typename T>
bool Inside(const Point2<T> &pt, const Bounds2<T> &b) {
    return (pt.x >= b.pMin.x && pt.x <= b.pMax.x && pt.y >= b.pMin.y &&
            pt.y <= b.pMax.y);
}

template <typename T>
bool InsideExclusive(const Point2<T> &pt, const Bounds2<T> &b) {
    return (pt.x >= b.pMin.x && pt.x < b.pMax.x && pt.y >= b.pMin.y &&
            pt.y < b.pMax.y);
}

template <typename T, typename U>
Bounds2<T> Expand(const Bounds2<T> &b, U delta) {
    return Bounds2<T>(b.pMin - Vector2<T>(delta, delta),
                      b.pMax + Vector2<T>(delta, delta));
}

/*

## Section 3.1.2 Ray–Bounds Intersections
https://zhuanlan.zhihu.com/p/42397911

把轴对齐包围盒看作是三对分别和 xyz 轴垂直的平面, 依次计算 ray 到每对平面的距离, 得到三组区间 : 
$\left[t_{\text {near}}, t_{\text {far}}\right]_{x},
 \left[t_{\text {near}}, t_{\text {far}}\right]_{y},
 \left[t_{\text {near}}, t_{\text {far}}\right]_{z}$

当三组区间存在重叠时, 即 $\max \left\{t_{near}\right\} < \min \left\{t_{far}\right\}$, 存在交点:
$t_{0}=\max \left\{t_{near}\right\}, 
 t_{1}=\min \left\{t_{far}\right\}$

若 ray 在包围盒内部, 则总是相交, 返回 hitt0 = 0, hitt1 仍为 $\min \left\{t_{far}\right\}$

对于 Num / 0 = inf/-inf 和 0 / 0 = NaN 的情况, 按正常执行路径可以处理, 无需特殊考虑

*/

// 若相交则通过 hit_t0 和 hit_t1 返回从 ray 到近, 远两个交点的距离(通过 ray 的 operator() 操作可以很容易得到对应的位置)
template <typename T>
inline bool Bounds3<T>::IntersectP(const Ray &ray, Float *hitt0,
                                   Float *hitt1) const {
    Float t0 = 0, t1 = ray.tMax;

    for (int i = 0; i < 3; ++i) {
        // Update interval for _i_th bounding box slab
		// $0 = a\left(\mathrm{o}_{x}+t \mathbf{d}_{x}\right)+b\left(\mathrm{o}_{y}+t \mathbf{d}_{y}\right)+c\left(\mathrm{o}_{z}+t \mathbf{d}_{z}\right)+d$
		// => $t = \frac{-d-((a, b, c) \cdot \mathrm{o})}{((a, b, c) \cdot \mathbf{d})}$
		// 因为每组平面都是垂直于坐标轴的, abc 三个参数中有两个为 0, 计算结果简化为 
		// $t_{i}=\frac{x_{i}-\mathrm{O}_{x}}{\mathrm{d}_{x}}$
        Float invRayDir = 1 / ray.d[i];
        Float tNear = (pMin[i] - ray.o[i]) * invRayDir;
        Float tFar = (pMax[i] - ray.o[i]) * invRayDir;

        // Update parametric interval from slab intersection $t$ values
        if (tNear > tFar) std::swap(tNear, tFar);

        // Update _tFar_ to ensure robust ray--bounds intersection
        tFar *= 1 + 2 * gamma(3);

        t0 = tNear > t0 ? tNear : t0;
        t1 = tFar < t1 ? tFar : t1;

        if (t0 > t1) return false;
    }

    if (hitt0) *hitt0 = t0;
    if (hitt1) *hitt1 = t1;

    return true;
}

/*

This intersection test is at the heart of traversing the BVHAccel acceleration structure, which is introduced in Section 4.3. 
Because so many ray–bounding box intersection tests are performed while traversing the BVH tree, 
we found that this optimized method provided approximately a 15% performance improvement in overall rendering time
compared to using the Bounds3::IntersectP() variant that didn’t take the precomputed direction reciprocals and signs.
只用在 bvh.cpp 里, Ray–Bounds Intersections 是其核心之一, 通过预先计算 invDir 和 dirIsNeg 带来了 15% 的性能提升 

*/
template <typename T>
inline bool Bounds3<T>::IntersectP(const Ray &ray, const Vector3f &invDir,
                                   const int dirIsNeg[3]) const {
    const Bounds3f &bounds = *this;
    // Check for ray intersection against $x$ and $y$ slabs
    Float tMin  = (bounds[    dirIsNeg[0]].x - ray.o.x) * invDir.x;
    Float tMax  = (bounds[1 - dirIsNeg[0]].x - ray.o.x) * invDir.x;
    Float tyMin = (bounds[    dirIsNeg[1]].y - ray.o.y) * invDir.y;
    Float tyMax = (bounds[1 - dirIsNeg[1]].y - ray.o.y) * invDir.y;

    // Update _tMax_ and _tyMax_ to ensure robust bounds intersection
    tMax  *= 1 + 2 * gamma(3);
    tyMax *= 1 + 2 * gamma(3);
    if (tMin > tyMax || tyMin > tMax) 
		return false;

    if (tyMin > tMin) tMin = tyMin;
    if (tyMax < tMax) tMax = tyMax;

    // Check for ray intersection against $z$ slab
    Float tzMin = (bounds[    dirIsNeg[2]].z - ray.o.z) * invDir.z;
    Float tzMax = (bounds[1 - dirIsNeg[2]].z - ray.o.z) * invDir.z;

    // Update _tzMax_ to ensure robust bounds intersection
    tzMax *= 1 + 2 * gamma(3);
    if (tMin > tzMax || tzMin > tMax) return false;
    if (tzMin > tMin) tMin = tzMin;
    if (tzMax < tMax) tMax = tzMax;
    return (tMin < ray.tMax) && (tMax > 0);
}

inline Point3f OffsetRayOrigin(const Point3f &p, const Vector3f &pError,
                               const Normal3f &n, const Vector3f &w) {
    Float d = Dot(Abs(n), pError);
#ifdef PBRT_FLOAT_AS_DOUBLE
    // We have tons of precision; for now bump up the offset a bunch just
    // to be extra sure that we start on the right side of the surface
    // (In case of any bugs in the epsilons code...)
    d *= 1024.;
#endif
    Vector3f offset = d * Vector3f(n);
    if (Dot(w, n) < 0) offset = -offset;
    Point3f po = p + offset;
    // Round offset point _po_ away from _p_
    for (int i = 0; i < 3; ++i) {
        if (offset[i] > 0)
            po[i] = NextFloatUp(po[i]);
        else if (offset[i] < 0)
            po[i] = NextFloatDown(po[i]);
    }
    return po;
}

inline Vector3f SphericalDirection(Float sinTheta, Float cosTheta, Float phi) {
    return Vector3f(sinTheta * std::cos(phi), sinTheta * std::sin(phi),
                    cosTheta);
}

inline Vector3f SphericalDirection(Float sinTheta, Float cosTheta, Float phi,
                                   const Vector3f &x, const Vector3f &y,
                                   const Vector3f &z) {
    return sinTheta * std::cos(phi) * x + sinTheta * std::sin(phi) * y +
           cosTheta * z;
}

inline Float SphericalTheta(const Vector3f &v) {
    return std::acos(Clamp(v.z, -1, 1));
}

inline Float SphericalPhi(const Vector3f &v) {
    Float p = std::atan2(v.y, v.x);
    return (p < 0) ? (p + 2 * Pi) : p;
}

}  // namespace pbrt

#endif  // PBRT_CORE_GEOMETRY_H
