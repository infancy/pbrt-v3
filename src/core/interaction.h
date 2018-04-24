
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

// Interaction 用于记录 ray 与 surface/media 相交测试后所生成交点的各种信息，
// 用于后续着色计算、生成次级光线等

// Interaction Declarations
struct Interaction {
    // Interaction Public Methods
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
    Point3f p;     // （世界空间中）交点的位置
    Float time;
    Vector3f pError;  // 累积的浮点数绝对误差
    Vector3f wo;   // 入射光线的方向
    Normal3f n;    // 交点的法线
    MediumInterface mediumInterface;
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
    void SetShadingGeometry(const Vector3f &dpdu, const Vector3f &dpdv,
                            const Normal3f &dndu, const Normal3f &dndv,
                            bool orientationIsAuthoritative);
    void ComputeScatteringFunctions(
        const RayDifferential &ray, MemoryArena &arena,
        bool allowMultipleLobes = false,
        TransportMode mode = TransportMode::Radiance);
    void ComputeDifferentials(const RayDifferential &r) const;
    Spectrum Le(const Vector3f &w) const;

    // SurfaceInteraction Public Data
	Point2f uv;				// 基于表面参数化的UV坐标
	Vector3f dpdu, dpdv;	// 参数偏微分，位于切平面上，可生成法线
	Normal3f dndu, dndv;	// 表面法线变化的的偏微分
    const Shape *shape = nullptr;

	// 存储由凹凸纹理或三角形网格逐顶点法线插值得到的着色法线等值
    struct {
        Normal3f n;		// 着色法线
        Vector3f dpdu, dpdv;
        Normal3f dndu, dndv;
    } shading;

    const Primitive *primitive = nullptr;
    BSDF *bsdf = nullptr;
    BSSRDF *bssrdf = nullptr;
    mutable Vector3f dpdx, dpdy;
    mutable Float dudx = 0, dvdx = 0, dudy = 0, dvdy = 0;

    // Added after book publication. Shapes can optionally provide a face
    // index with an intersection point for use in Ptex texture lookups.
    // If Ptex isn't being used, then this value is ignored.
    int faceIndex = 0;
};

}  // namespace pbrt

#endif  // PBRT_CORE_INTERACTION_H
