
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

#ifndef PBRT_CORE_INTERACTION_H
#define PBRT_CORE_INTERACTION_H

// core/interaction.h*
#include "pbrt.h"
#include "geometry.h"
#include "transform.h"
#include "medium.h"
#include "material.h"

namespace pbrt {

// generating from a ray–shape intersection or from a ray passing through participating media
// interaction 记录 ray 与 surface/media 通过相交测试后生成的局部交点上的各种信息, 用于后续着色计算, 次级光线生成等

// Interaction Declarations
struct Interaction {
    // Interaction Public Methods
    // A number of Interaction constructors are available; depending on what sort of interaction is being constructed and what sort of information about it is relevant
    Interaction() : time(0) {}
    Interaction(const Point3f &p, const Normal3f &n, const Vector3f &pError,
                const Vector3f &wo, Float time,
                const MediumInterface &mediumInterface)
        : p(p),
          time(time),
          pError(pError),
          wo(Normalize(wo)),
          n(n),
          mediumInterface(mediumInterface) {}

    bool IsSurfaceInteraction() const { return n != Normal3f(); }

    // 在交点 p 上, 根据出射方向 d 或另一交点 p2/it 生成新的光线
    Ray SpawnRay(const Vector3f &d) const {
        Point3f o = OffsetRayOrigin(p, pError, n, d);
        return Ray(o, d, Infinity, time, GetMedium(d));
    }
    Ray SpawnRayTo(const Point3f &p2) const {
        Point3f origin = OffsetRayOrigin(p, pError, n, p2 - p);
        Vector3f d = p2 - p;
        return Ray(origin, d, 1 - ShadowEpsilon, time, GetMedium(d));
    }
    Ray SpawnRayTo(const Interaction &it) const {
        Point3f origin = OffsetRayOrigin(p, pError, n, it.p - p);
        Point3f target = OffsetRayOrigin(it.p, it.pError, it.n, origin - it.p);
        Vector3f d = target - origin;
        return Ray(origin, d, 1 - ShadowEpsilon, time, GetMedium(d));
    }

    Interaction(const Point3f &p, const Vector3f &wo, Float time,
                const MediumInterface &mediumInterface)
        : p(p), time(time), wo(wo), mediumInterface(mediumInterface) {}
    Interaction(const Point3f &p, Float time,
                const MediumInterface &mediumInterface)
        : p(p), time(time), mediumInterface(mediumInterface) {}

    bool IsMediumInteraction() const { return !IsSurfaceInteraction(); }

    const Medium *GetMedium(const Vector3f &w) const {
        return Dot(w, n) > 0 ? mediumInterface.outside : mediumInterface.inside;
    }
    const Medium *GetMedium() const {
        CHECK_EQ(mediumInterface.inside, mediumInterface.outside);
        return mediumInterface.inside;
    }

    // Interaction Public Data
    Point3f p;       // （世界空间中）交点的位置
    Float time;
    // For interactions where the point  was computed by ray intersection, some floating-point error is generally present in the p value. 
    // pError gives a conservative bound on this error; it’s  for points in participating media. 
    // See Section 3.9 for more on pbrt’s approach to managing floating-point error and in particular Section 3.9.4 for how this bound is computed for various shapes.
    Vector3f pError; // 累积的浮点数绝对误差

    // prev   light
    // -----  -----
    //   ^      ^
    //    \    /
    //  wo \  / wi
    //      \/
    //    ------
    //    isect
    Vector3f wo;     // 出射光线的方向(只在 ray–shape intersection 时才会生成这个参数???)
    Normal3f n;      // 交点处的表面法线
    MediumInterface mediumInterface; // 交点所处于的 medium
};

class MediumInteraction : public Interaction {
  public:
    // MediumInteraction Public Methods
    MediumInteraction() : phase(nullptr) {}
    MediumInteraction(const Point3f &p, const Vector3f &wo, Float time,
                      const Medium *medium, const PhaseFunction *phase)
        : Interaction(p, wo, time, medium), phase(phase) {}

    bool IsValid() const { return phase != nullptr; }

    // MediumInteraction Public Data
    const PhaseFunction *phase;
};

// The geometry of a particular point on a surface (often a position found by intersecting a ray against the surface) is represented by a SurfaceInteraction. 
// Having this abstraction lets most of the system work with points on surfaces without needing to consider the particular type of geometric shape the points lie on;
// SurfaceInteraction 隔离了后续 lighting, shading 等操作与具体 shape 的联系
// SurfaceInteraction Declarations
class SurfaceInteraction : public Interaction {
  public:
    // SurfaceInteraction Public Methods
    SurfaceInteraction() {}
    SurfaceInteraction(const Point3f &p, const Vector3f &pError,
                       const Point2f &uv, const Vector3f &wo,
                       const Vector3f &dpdu, const Vector3f &dpdv,
                       const Normal3f &dndu, const Normal3f &dndv, Float time,
                       const Shape *sh,
                       int faceIndex = 0);

    // 设置不同于表面原始几何结构, 用于着色的 shading 结构
    void SetShadingGeometry(const Vector3f &dpdu, const Vector3f &dpdv,
                            const Normal3f &dndu, const Normal3f &dndv,
                            bool orientationIsAuthoritative);

    // 计算交点上的散射方程, 结果会赋值给 this->bsdf, this->bssrdf 等, 用于后续的着色计算
    void ComputeScatteringFunctions(
        const RayDifferential &ray, MemoryArena &arena,
        bool allowMultipleLobes = false,
        TransportMode mode = TransportMode::Radiance);

    // compute information about the projected size of the surface area around the intersection on the image plane 
    // for use in texture antialiasing
    void ComputeDifferentials(const RayDifferential &r) const;

    // 如果交点处于一个自发光的图元上, 出射光 = 自发光 + 反射光, Le(light emit) 函数计算自发光的大小
    Spectrum Le(const Vector3f &w) const;

    // all of the shapes that pbrt supports do have at least a local parametric description - that for some range of $(u, v) $ values, points on the surface are given by some function $f$ such that $\mathrm{ p} = f(u, v) $
    // The parametric partial derivatives of the surface, $\partial \mathrm{p} / \partial u$ and $\partial \mathrm{p} / \partial v$, lie in the tangent plane but are not necessarily orthogonal.
    // The surface normal $\mathbf{n}$ is given by the cross product of $\partial \mathrm{p} / \partial u$ and $\partial \mathrm{p} / \partial v .$ 
    // The vectors $\partial \mathbf{n} / \partial u$ and $\partial \mathbf{n} / \partial v$ record the differential change in surface normal as we move $u$ and $v$ along the surface.
    // pbrt 使用的所有 shape 至少有一种表面参数化的形式

    // SurfaceInteraction Public Data
    Point2f uv;				// 表面参数化的 UV 坐标, 用于纹理采样, 生成法线等
	Vector3f dpdu, dpdv;	// 交点 p 的参数偏导数(表面的变化率), 位于点 p 的切平面上(这俩不需要正交)
	Normal3f dndu, dndv;	// 法线 n 的参数偏导数(法线的变化率)
    const Shape *shape = nullptr;

    // SurfaceInteraction stores a second instance of a surface normal and the various partial derivatives to represent possibly perturbed values of these quantities as can be generated by bump mapping or interpolated per-vertex normals with triangles.
	// 存储由凹凸纹理或三角形网格逐顶点法线插值得到的着色法线等值
    struct {
        Normal3f n;	// 着色法线
        Vector3f dpdu, dpdv;
        Normal3f dndu, dndv;
    } shading;

    const Primitive *primitive = nullptr; // 交点所处的 primitive, 用来获取交点的材质信息

    BSDF *bsdf = nullptr;
    BSSRDF *bssrdf = nullptr;

    // ???
    mutable Vector3f dpdx, dpdy;
    mutable Float dudx = 0, dvdx = 0, dudy = 0, dvdy = 0;

    // Added after book publication. Shapes can optionally provide a face
    // index with an intersection point for use in Ptex texture lookups.
    // If Ptex isn't being used, then this value is ignored.
    int faceIndex = 0;
};

}  // namespace pbrt

#endif  // PBRT_CORE_INTERACTION_H
