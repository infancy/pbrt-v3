
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

#ifndef PBRT_CORE_TRANSFORM_H
#define PBRT_CORE_TRANSFORM_H

/*

## Section2.7 Transformations

In general, a transformation  is a mapping from points to points and from vectors to vectors

Linear, Continuous, One-to-one and invertibe
线性, 连续且壹壹对应,

there are two different (and incompatible!) ways that a matrix can be interpreted:
矩阵(的变换)可以被解释为两种不同的方式???:

Transformation of the frame: Given a point, the matrix could express how to compute a new point in the same frame that represents the transformation of the original point (e.g., by translating it in some direction).
空间不动, 将点从空间中的一个位置映射到另一个位置

Transformation from one frame to another: A matrix can express the coordinates of a point or vector in a new frame in terms of the coordinates in the original frame.
点的位置相对于原空间不动, 对原空间进行变换

后者意味着, 我们可以先对顶点所处的空间进行变换, 再计算顶点在变换后的空间中的位置(这样的表达有点问题, 看懂下面两个式子就行了)

$$\mathbf{M} \mathbf{p} = \mathbf{M} \mathbf{p}_{0} + s_{1} \mathbf{M} \mathbf{v}_{1} + \cdots + s_{n} \mathbf{M} \mathbf{v}_{n}$$

$$\mathbf{M} \mathbf{v} = s_{1} \mathbf{M} \mathbf{v}_{1} + \cdots + s_{n} \mathbf{M} \mathbf{v}_{n}$$



### Section2.7.1 Homogeneous Coordinates

Given a frame defined by $\left(\mathrm{p}_{\mathrm{o}}, \mathbf{V}_{\mathbf{1}}, \mathbf{V}_{\mathbf{2}}, \mathbf{v}_{\mathbf{3}}\right)$, there is ambiguity between the representation of 
a point $\left(\mathrm{p}_{x}, \mathrm{p}_{y}, \mathrm{p}_{z}\right)$ and a vector $\left(\mathbf{v}_{x}, \mathbf{v}_{y}, \mathbf{v}_{z}\right)$ with the same coordinates

使用点和向量的**齐次坐标(Homogeneous coordinate)表示**可以区分这种二义性, 也就是引入第四个分量(权重)

Converting homogeneous points into ordinary points entails dividing the first three components by the weight:
$$(x, y, z, w) \rightarrow \left(\frac{x}{w}, \frac{y}{w}, \frac{z}{w}\right)$$



In general, by characterizing how **the basis(frame) is transformed**, we know how any point or vector specified in terms of that basis is transformed. 

Because points and vectors in the current coordinate system are expressed in terms of the current coordinate system’s frame, 
applying the transformation to them directly is equivalent to 
applying the transformation to the current coordinate system’s basis and finding their coordinates in terms of the transformed basis.

PBRT 没有显式的使用齐次坐标, 而是将其隐藏在每次变换的具体实现中(在变换时, 先将点, 向量等转换成齐次坐标的表示进行变换, 再变换回普通的形式)



### Section2.7.3 Translations

变换的几个特性使得平移, 放缩和旋转各自具有一些性质: ...

#### 平移

因为在变换时使用了齐次坐标表示, 平移只会作用于点, 而不会对向量起作用

#### 放缩

#### 围绕 x, z, y 等基向量的旋转

In general, we can define an arbitrary axis from the origin in any direction and then rotate around that axis by a given angle. 
The most common rotations of this type are around the $x$, $y$, and $z$ coordinate axes. 
We will write these rotations as $\mathbf{R}_{x}(\theta), \mathbf{R}_{y}(\theta)$, and so on. The rotation around an arbitrary axis $(x, y, z)$ is denoted by $\mathbf{R}_{(x, y, z)}(\theta)$.

旋转矩阵的逆矩阵比较有意思, 等于其转置矩阵, 这样的矩阵又叫做正交矩阵(orthogonal matrix) 
its first three columns (or rows) are all normalized and orthogonal to each other. 
Fortunately, the transpose is much easier to compute than a full matrix inverse.

#### 围绕任意轴的旋转

首先将 $\mathbf_{V}$ 分解成平行于旋转轴 $a$ 的向量 $\mathbf_{V}_{parallel}$ 和垂直于 $a$ 的向量 $\mathbf_{V}_{bot} $
然后对后者进行旋转得到 $\mathbf_{V}_{bot}^{\prime}$, 
则 $\mathbf{V}^{\prime} = \mathbf_{V}_{parallel} + \mathbf_{V}_{bot}^{\prime} $

Figure2.12 的 $\mathbf_{V}$ 和 $\mathbf{V}^{\prime}$ 应该交换下, $\mathbf{v}_{\mathbf{c}}$ 也没有标识出来

To convert this to a rotation matrix, we apply this formula to the basis vectors $(1, 0, 0)$, $(0, 1, 0)$, and $(0, 0, 1)$ to get the values of the rows of the matrix.

todo: 看看 glm 的实现

#### lookat 矩阵

*/

// core/transform.h*
#include "pbrt.h"
#include "stringprint.h"
#include "geometry.h"
#include "quaternion.h"

namespace pbrt {

// Matrix4x4 的解释在 Section A.5.3 中

// Matrix4x4 Declarations
struct Matrix4x4 {
    // Matrix4x4 Public Methods
    // 默认初始化的矩阵是一个单位矩阵(identity matrix), 总是执行恒等变换(identity transformation)
    Matrix4x4() {
        m[0][0] = m[1][1] = m[2][2] = m[3][3] = 1.f;
        m[0][1] = m[0][2] = m[0][3] = m[1][0] = m[1][2] = m[1][3] = m[2][0] =
            m[2][1] = m[2][3] = m[3][0] = m[3][1] = m[3][2] = 0.f;
    }
    Matrix4x4(Float mat[4][4]);
    Matrix4x4(Float t00, Float t01, Float t02, Float t03, Float t10, Float t11,
              Float t12, Float t13, Float t20, Float t21, Float t22, Float t23,
              Float t30, Float t31, Float t32, Float t33);

    bool operator==(const Matrix4x4 &m2) const {
        for (int i = 0; i < 4; ++i)
            for (int j = 0; j < 4; ++j)
                if (m[i][j] != m2.m[i][j]) return false;
        return true;
    }
    bool operator!=(const Matrix4x4 &m2) const {
        for (int i = 0; i < 4; ++i)
            for (int j = 0; j < 4; ++j)
                if (m[i][j] != m2.m[i][j]) return true;
        return false;
    }

    friend Matrix4x4 Transpose(const Matrix4x4 &);

    void Print(FILE *f) const {
        fprintf(f, "[ ");
        for (int i = 0; i < 4; ++i) {
            fprintf(f, "  [ ");
            for (int j = 0; j < 4; ++j) {
                fprintf(f, "%f", m[i][j]);
                if (j != 3) fprintf(f, ", ");
            }
            fprintf(f, " ]\n");
        }
        fprintf(f, " ] ");
    }

    static Matrix4x4 Mul(const Matrix4x4 &m1, const Matrix4x4 &m2) {
        Matrix4x4 r;
        for (int i = 0; i < 4; ++i)
            for (int j = 0; j < 4; ++j)
                r.m[i][j] = m1.m[i][0] * m2.m[0][j] + m1.m[i][1] * m2.m[1][j] +
                            m1.m[i][2] * m2.m[2][j] + m1.m[i][3] * m2.m[3][j];
        return r;
    }
    friend Matrix4x4 Inverse(const Matrix4x4 &);

    friend std::ostream &operator<<(std::ostream &os, const Matrix4x4 &m) {
        // clang-format off
        os << StringPrintf("[ [ %f, %f, %f, %f ] "
                             "[ %f, %f, %f, %f ] "
                             "[ %f, %f, %f, %f ] "
                             "[ %f, %f, %f, %f ] ]",
                           m.m[0][0], m.m[0][1], m.m[0][2], m.m[0][3],
                           m.m[1][0], m.m[1][1], m.m[1][2], m.m[1][3],
                           m.m[2][0], m.m[2][1], m.m[2][2], m.m[2][3],
                           m.m[3][0], m.m[3][1], m.m[3][2], m.m[3][3]);
        // clang-format on
        return os;
    }

    // PBRT 的矩阵是行主序(row-major)的, 所以元素 $m_{i, j}$ 等价于 m[i][j]
    Float m[4][4];
};

// Transform Declarations
class Transform {
  public:
    // Transform Public Methods
    Transform() {}
    Transform(const Float mat[4][4]) {
        m = Matrix4x4(mat[0][0], mat[0][1], mat[0][2], mat[0][3], mat[1][0],
                      mat[1][1], mat[1][2], mat[1][3], mat[2][0], mat[2][1],
                      mat[2][2], mat[2][3], mat[3][0], mat[3][1], mat[3][2],
                      mat[3][3]);
        mInv = Inverse(m);
    }
    Transform(const Matrix4x4 &m) : m(m), mInv(Inverse(m)) {}
    Transform(const Matrix4x4 &m, const Matrix4x4 &mInv) : m(m), mInv(mInv) {}

    void Print(FILE *f) const;

    // Transform 的逆变换, 只需要简单的交换 m 和 mInv 就可以了
    friend Transform Inverse(const Transform &t) {
        return Transform(t.mInv, t.m);
    }
    friend Transform Transpose(const Transform &t) {
        return Transform(Transpose(t.m), Transpose(t.mInv));
    }

    bool operator==(const Transform &t) const {
        return t.m == m && t.mInv == mInv;
    }
    bool operator!=(const Transform &t) const {
        return t.m != m || t.mInv != mInv;
    }
    bool operator<(const Transform &t2) const {
        for (int i = 0; i < 4; ++i)
            for (int j = 0; j < 4; ++j) {
                if (m.m[i][j] < t2.m.m[i][j]) return true;
                if (m.m[i][j] > t2.m.m[i][j]) return false;
            }
        return false;
    }

    bool IsIdentity() const {
        return (m.m[0][0] == 1.f && m.m[0][1] == 0.f && m.m[0][2] == 0.f &&
                m.m[0][3] == 0.f && m.m[1][0] == 0.f && m.m[1][1] == 1.f &&
                m.m[1][2] == 0.f && m.m[1][3] == 0.f && m.m[2][0] == 0.f &&
                m.m[2][1] == 0.f && m.m[2][2] == 1.f && m.m[2][3] == 0.f &&
                m.m[3][0] == 0.f && m.m[3][1] == 0.f && m.m[3][2] == 0.f &&
                m.m[3][3] == 1.f);
    }

    const Matrix4x4 &GetMatrix() const { return m; }
    const Matrix4x4 &GetInverseMatrix() const { return mInv; }

    // It’s useful to be able to test if a transformation has a scaling term in it; an easy way to do this is to transform the three coordinate axes and see if any of their lengths are appreciably different from one.
    // 一种简单的方式, 是对三个基向量进行变换, 看向量长度是否发生了变化(在平移, 放缩和旋转中, 只有放缩会影响到长度)
    bool HasScale() const {
        Float la2 = (*this)(Vector3f(1, 0, 0)).LengthSquared();
        Float lb2 = (*this)(Vector3f(0, 1, 0)).LengthSquared();
        Float lc2 = (*this)(Vector3f(0, 0, 1)).LengthSquared();
#define NOT_ONE(x) ((x) < .999f || (x) > 1.001f)
        return (NOT_ONE(la2) || NOT_ONE(lb2) || NOT_ONE(lc2));
#undef NOT_ONE
    }

    template <typename T>
    inline Point3<T> operator()(const Point3<T> &p) const;
    template <typename T>
    inline Vector3<T> operator()(const Vector3<T> &v) const;
    template <typename T>
    inline Normal3<T> operator()(const Normal3<T> &) const;
    inline Ray operator()(const Ray &r) const;
    inline RayDifferential operator()(const RayDifferential &r) const;
    Bounds3f operator()(const Bounds3f &b) const;
    SurfaceInteraction operator()(const SurfaceInteraction &si) const;

    Transform operator*(const Transform &t2) const;

    bool SwapsHandedness() const;

    template <typename T>
    inline Point3<T> operator()(const Point3<T> &pt,
                                Vector3<T> *absError) const;
    template <typename T>
    inline Point3<T> operator()(const Point3<T> &p, const Vector3<T> &pError,
                                Vector3<T> *pTransError) const;
    template <typename T>
    inline Vector3<T> operator()(const Vector3<T> &v,
                                 Vector3<T> *vTransError) const;
    template <typename T>
    inline Vector3<T> operator()(const Vector3<T> &v, const Vector3<T> &vError,
                                 Vector3<T> *vTransError) const;
    inline Ray operator()(const Ray &r, Vector3f *oError,
                          Vector3f *dError) const;
    inline Ray operator()(const Ray &r, const Vector3f &oErrorIn,
                          const Vector3f &dErrorIn, Vector3f *oErrorOut,
                          Vector3f *dErrorOut) const;

    friend std::ostream &operator<<(std::ostream &os, const Transform &t) {
        os << "t=" << t.m << ", inv=" << t.mInv;
        return os;
    }

  private:
    // Transform Private Data

    // 为了避免重复计算, 还存储了 m 的逆矩阵 mInv
    // 这样存储一个 Transform_float 就需要 128 bytes, 比较耗费空间
    // 考虑到场景中大量 Transform 是相同的, PBRT 在 Shape 对象中就只存储到 Transform 的指针
    // 这虽然带来了限制(丧失了灵活性), 但 Transform 大多在创建场景时就确定下来了
    Matrix4x4 m, mInv;
    friend class AnimatedTransform;
    friend struct Quaternion;
};

Transform Translate(const Vector3f &delta);
Transform Scale(Float x, Float y, Float z);
Transform RotateX(Float theta);
Transform RotateY(Float theta);
Transform RotateZ(Float theta);
Transform Rotate(Float theta, const Vector3f &axis);
Transform LookAt(const Point3f &pos, const Point3f &look, const Vector3f &up);
Transform Orthographic(Float znear, Float zfar);
Transform Perspective(Float fov, Float znear, Float zfar);
bool SolveLinearSystem2x2(const Float A[2][2], const Float B[2], Float *x0,
                          Float *x1);

// Transform Inline Functions
template <typename T>
inline Point3<T> Transform::operator()(const Point3<T> &p) const {
    // 这里顶点的齐次坐标是隐式表达的
    T x = p.x, y = p.y, z = p.z;
    T xp = m.m[0][0] * x + m.m[0][1] * y + m.m[0][2] * z + m.m[0][3];
    T yp = m.m[1][0] * x + m.m[1][1] * y + m.m[1][2] * z + m.m[1][3];
    T zp = m.m[2][0] * x + m.m[2][1] * y + m.m[2][2] * z + m.m[2][3];
    T wp = m.m[3][0] * x + m.m[3][1] * y + m.m[3][2] * z + m.m[3][3];

    CHECK_NE(wp, 0);

    // 最后要除以 $w$ 以返回普通形式的坐标
    // Chapter6 的投影变换需要这种除法
    if (wp == 1)
        return Point3<T>(xp, yp, zp);
    else
        return Point3<T>(xp, yp, zp) / wp;
}

template <typename T>
inline Vector3<T> Transform::operator()(const Vector3<T> &v) const {
    // 忽略位移
    T x = v.x, y = v.y, z = v.z;
    return Vector3<T>(m.m[0][0] * x + m.m[0][1] * y + m.m[0][2] * z,
                      m.m[1][0] * x + m.m[1][1] * y + m.m[1][2] * z,
                      m.m[2][0] * x + m.m[2][1] * y + m.m[2][2] * z);
}

/*

$$
\begin{aligned} 0 
&=(\mathbf{n}^{\prime})^{T}        \mathbf{t}^{\prime} \\
&=(\mathbf{S}      \mathbf{n})^{T} \mathbf{M} \mathbf{t} \\ 
&=(\mathbf{n})^{T} \mathbf{S}^{T}  \mathbf{M} \mathbf{t} 
\end{aligned}
$$

This condition holds if $\mathbf{S}^{T} \mathbf{M}=\mathbf{I}$, the identity matrix. 
Therefore, $\mathbf{S}^{T}=\mathbf{M}^{-1}$, and so $\mathbf{S}=\left(\mathbf{M}^{-1}\right)^{T}$

需要使用变换矩阵 m 的逆矩阵 mInv 的转置矩阵 $\left(\mathbf{M}^{-1}\right)^{T}$ 乘以法线 n
*/
template <typename T>
inline Normal3<T> Transform::operator()(const Normal3<T> &n) const {
    T x = n.x, y = n.y, z = n.z;
    return Normal3<T>(mInv.m[0][0] * x + mInv.m[1][0] * y + mInv.m[2][0] * z,
                      mInv.m[0][1] * x + mInv.m[1][1] * y + mInv.m[2][1] * z,
                      mInv.m[0][2] * x + mInv.m[1][2] * y + mInv.m[2][2] * z);
}

inline Ray Transform::operator()(const Ray &r) const {
    Vector3f oError;
    Point3f o = (*this)(r.o, &oError);
    Vector3f d = (*this)(r.d);

    // todo: Section 3.9.4
    // Offset ray origin to edge of error bounds and compute _tMax_
    Float lengthSquared = d.LengthSquared();
    Float tMax = r.tMax;
    if (lengthSquared > 0) {
        Float dt = Dot(Abs(d), oError) / lengthSquared;
        o += d * dt;
        tMax -= dt;
    }
    return Ray(o, d, tMax, r.time, r.medium);
}

inline RayDifferential Transform::operator()(const RayDifferential &r) const {
    Ray tr = (*this)(Ray(r));
    RayDifferential ret(tr.o, tr.d, tr.tMax, tr.time, tr.medium);
    ret.hasDifferentials = r.hasDifferentials;
    ret.rxOrigin = (*this)(r.rxOrigin);
    ret.ryOrigin = (*this)(r.ryOrigin);
    ret.rxDirection = (*this)(r.rxDirection);
    ret.ryDirection = (*this)(r.ryDirection);
    return ret;
}

template <typename T>
inline Point3<T> Transform::operator()(const Point3<T> &p,
                                       Vector3<T> *pError) const {
    T x = p.x, y = p.y, z = p.z;
    // Compute transformed coordinates from point _pt_
    T xp = (m.m[0][0] * x + m.m[0][1] * y) + (m.m[0][2] * z + m.m[0][3]);
    T yp = (m.m[1][0] * x + m.m[1][1] * y) + (m.m[1][2] * z + m.m[1][3]);
    T zp = (m.m[2][0] * x + m.m[2][1] * y) + (m.m[2][2] * z + m.m[2][3]);
    T wp = (m.m[3][0] * x + m.m[3][1] * y) + (m.m[3][2] * z + m.m[3][3]);

    // Compute absolute error for transformed point
    T xAbsSum = (std::abs(m.m[0][0] * x) + std::abs(m.m[0][1] * y) +
                 std::abs(m.m[0][2] * z) + std::abs(m.m[0][3]));
    T yAbsSum = (std::abs(m.m[1][0] * x) + std::abs(m.m[1][1] * y) +
                 std::abs(m.m[1][2] * z) + std::abs(m.m[1][3]));
    T zAbsSum = (std::abs(m.m[2][0] * x) + std::abs(m.m[2][1] * y) +
                 std::abs(m.m[2][2] * z) + std::abs(m.m[2][3]));
    *pError = gamma(3) * Vector3<T>(xAbsSum, yAbsSum, zAbsSum);
    CHECK_NE(wp, 0);
    if (wp == 1)
        return Point3<T>(xp, yp, zp);
    else
        return Point3<T>(xp, yp, zp) / wp;
}

template <typename T>
inline Point3<T> Transform::operator()(const Point3<T> &pt,
                                       const Vector3<T> &ptError,
                                       Vector3<T> *absError) const {
    T x = pt.x, y = pt.y, z = pt.z;
    T xp = (m.m[0][0] * x + m.m[0][1] * y) + (m.m[0][2] * z + m.m[0][3]);
    T yp = (m.m[1][0] * x + m.m[1][1] * y) + (m.m[1][2] * z + m.m[1][3]);
    T zp = (m.m[2][0] * x + m.m[2][1] * y) + (m.m[2][2] * z + m.m[2][3]);
    T wp = (m.m[3][0] * x + m.m[3][1] * y) + (m.m[3][2] * z + m.m[3][3]);
    absError->x =
        (gamma(3) + (T)1) *
            (std::abs(m.m[0][0]) * ptError.x + std::abs(m.m[0][1]) * ptError.y +
             std::abs(m.m[0][2]) * ptError.z) +
        gamma(3) * (std::abs(m.m[0][0] * x) + std::abs(m.m[0][1] * y) +
                    std::abs(m.m[0][2] * z) + std::abs(m.m[0][3]));
    absError->y =
        (gamma(3) + (T)1) *
            (std::abs(m.m[1][0]) * ptError.x + std::abs(m.m[1][1]) * ptError.y +
             std::abs(m.m[1][2]) * ptError.z) +
        gamma(3) * (std::abs(m.m[1][0] * x) + std::abs(m.m[1][1] * y) +
                    std::abs(m.m[1][2] * z) + std::abs(m.m[1][3]));
    absError->z =
        (gamma(3) + (T)1) *
            (std::abs(m.m[2][0]) * ptError.x + std::abs(m.m[2][1]) * ptError.y +
             std::abs(m.m[2][2]) * ptError.z) +
        gamma(3) * (std::abs(m.m[2][0] * x) + std::abs(m.m[2][1] * y) +
                    std::abs(m.m[2][2] * z) + std::abs(m.m[2][3]));
    CHECK_NE(wp, 0);
    if (wp == 1.)
        return Point3<T>(xp, yp, zp);
    else
        return Point3<T>(xp, yp, zp) / wp;
}

template <typename T>
inline Vector3<T> Transform::operator()(const Vector3<T> &v,
                                        Vector3<T> *absError) const {
    T x = v.x, y = v.y, z = v.z;
    absError->x =
        gamma(3) * (std::abs(m.m[0][0] * v.x) + std::abs(m.m[0][1] * v.y) +
                    std::abs(m.m[0][2] * v.z));
    absError->y =
        gamma(3) * (std::abs(m.m[1][0] * v.x) + std::abs(m.m[1][1] * v.y) +
                    std::abs(m.m[1][2] * v.z));
    absError->z =
        gamma(3) * (std::abs(m.m[2][0] * v.x) + std::abs(m.m[2][1] * v.y) +
                    std::abs(m.m[2][2] * v.z));
    return Vector3<T>(m.m[0][0] * x + m.m[0][1] * y + m.m[0][2] * z,
                      m.m[1][0] * x + m.m[1][1] * y + m.m[1][2] * z,
                      m.m[2][0] * x + m.m[2][1] * y + m.m[2][2] * z);
}

template <typename T>
inline Vector3<T> Transform::operator()(const Vector3<T> &v,
                                        const Vector3<T> &vError,
                                        Vector3<T> *absError) const {
    T x = v.x, y = v.y, z = v.z;
    absError->x =
        (gamma(3) + (T)1) *
            (std::abs(m.m[0][0]) * vError.x + std::abs(m.m[0][1]) * vError.y +
             std::abs(m.m[0][2]) * vError.z) +
        gamma(3) * (std::abs(m.m[0][0] * v.x) + std::abs(m.m[0][1] * v.y) +
                    std::abs(m.m[0][2] * v.z));
    absError->y =
        (gamma(3) + (T)1) *
            (std::abs(m.m[1][0]) * vError.x + std::abs(m.m[1][1]) * vError.y +
             std::abs(m.m[1][2]) * vError.z) +
        gamma(3) * (std::abs(m.m[1][0] * v.x) + std::abs(m.m[1][1] * v.y) +
                    std::abs(m.m[1][2] * v.z));
    absError->z =
        (gamma(3) + (T)1) *
            (std::abs(m.m[2][0]) * vError.x + std::abs(m.m[2][1]) * vError.y +
             std::abs(m.m[2][2]) * vError.z) +
        gamma(3) * (std::abs(m.m[2][0] * v.x) + std::abs(m.m[2][1] * v.y) +
                    std::abs(m.m[2][2] * v.z));
    return Vector3<T>(m.m[0][0] * x + m.m[0][1] * y + m.m[0][2] * z,
                      m.m[1][0] * x + m.m[1][1] * y + m.m[1][2] * z,
                      m.m[2][0] * x + m.m[2][1] * y + m.m[2][2] * z);
}

inline Ray Transform::operator()(const Ray &r, Vector3f *oError,
                                 Vector3f *dError) const {
    Point3f o = (*this)(r.o, oError);
    Vector3f d = (*this)(r.d, dError);
    Float tMax = r.tMax;
    Float lengthSquared = d.LengthSquared();
    if (lengthSquared > 0) {
        Float dt = Dot(Abs(d), *oError) / lengthSquared;
        o += d * dt;
        //        tMax -= dt;
    }
    return Ray(o, d, tMax, r.time, r.medium);
}

inline Ray Transform::operator()(const Ray &r, const Vector3f &oErrorIn,
                                 const Vector3f &dErrorIn, Vector3f *oErrorOut,
                                 Vector3f *dErrorOut) const {
    Point3f o = (*this)(r.o, oErrorIn, oErrorOut);
    Vector3f d = (*this)(r.d, dErrorIn, dErrorOut);
    Float tMax = r.tMax;
    Float lengthSquared = d.LengthSquared();
    if (lengthSquared > 0) {
        Float dt = Dot(Abs(d), *oErrorOut) / lengthSquared;
        o += d * dt;
        //        tMax -= dt;
    }
    return Ray(o, d, tMax, r.time, r.medium);
}

/*

## Section2.9 Animating Transformations

The approach used for transformation interpolation in pbrt is based on matrix decomposition - given an arbitrary transformation matrix $\mathbf{M}$ , we decompose it into a concatenation of scale (S), rotation (R), and translation (T) transformations, $$ \mathrm{M}=\mathrm{SRT} $$ where each of those
components is independently interpolated and then the composite interpolated matrix is found by multiplying the three interpolated matrices together.

### 2.9.3 AnimatedTransform Implementation

### 2.9.4 Bounding Moving Bounding Boxes

TODO: 暂时跳过, 以后再看

*/

// AnimatedTransform Declarations
class AnimatedTransform {
  public:
    // AnimatedTransform Public Methods
    AnimatedTransform(const Transform *startTransform, Float startTime,
                      const Transform *endTransform, Float endTime);

    // decomposes the given composite transformation matrices into scaling, rotation, and translation components
    static void Decompose(const Matrix4x4 &m, Vector3f *T, Quaternion *R,
                          Matrix4x4 *S);

    // computes the interpolated transformation matrix at a given time
    void Interpolate(Float time, Transform *t) const;

    // Ray 自带 time 参数
    Ray operator()(const Ray &r) const;
    RayDifferential operator()(const RayDifferential &r) const;
    Point3f operator()(Float time, const Point3f &p) const;
    Vector3f operator()(Float time, const Vector3f &v) const;

    bool HasScale() const {
        return startTransform->HasScale() || endTransform->HasScale();
    }

    // 返回整个动画过程中, 包围盒移动所"触及"过的空间(包围盒的足迹)
    Bounds3f MotionBounds(const Bounds3f &b) const;

  private:
    // 计算包围盒一点在动画中所触及过的空间
    Bounds3f BoundPointMotion(const Point3f &p) const;

  private:
    // AnimatedTransform Private Data
    const Transform *startTransform, *endTransform;
    const Float startTime, endTime;
    const bool actuallyAnimated;

    Vector3f T[2];
    Quaternion R[2];
    Matrix4x4 S[2];

    bool hasRotation;

    // 用于 BoundPointMotion 函数的计算
    struct DerivativeTerm {
        DerivativeTerm() {}
        DerivativeTerm(Float c, Float x, Float y, Float z)
            : kc(c), kx(x), ky(y), kz(z) {}
        Float kc, kx, ky, kz;
        Float Eval(const Point3f &p) const {
            return kc + kx * p.x + ky * p.y + kz * p.z;
        }
    };
    DerivativeTerm c1[3], c2[3], c3[3], c4[3], c5[3];
};

}  // namespace pbrt

#endif  // PBRT_CORE_TRANSFORM_H
