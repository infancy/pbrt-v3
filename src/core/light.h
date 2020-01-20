
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

#ifndef PBRT_CORE_LIGHT_H
#define PBRT_CORE_LIGHT_H

// core/light.h*
#include "pbrt.h"
#include "memory.h"
#include "interaction.h"

namespace pbrt {

// LightFlags Declarations
enum class LightFlags : int 
{
    DeltaPosition = 1,
    DeltaDirection = 2,

    Area = 4,
    Infinite = 8
};

inline bool IsDeltaLight(int flags) 
{
    return flags & (int)LightFlags::DeltaPosition ||
           flags & (int)LightFlags::DeltaDirection;
}

// Light Declarations
class Light {
  public:
    // Light Interface
    virtual ~Light();
    Light(int flags, const Transform &LightToWorld,
          const MediumInterface &mediumInterface, int nSamples = 1);

    // 传入交点 ref，在光源上选取一点，计算选到该点的概率密度, 该点到交点的方向 wi, 入射辐射度（incident Radiance） Li 及可见性
    // 当光源是面积光源时，还需传入一个[0,1]^2范围的采样点，对面积光源上一点进行采样并记录相应**概率密度值pdf**
    virtual Spectrum Sample_Li(const Interaction &ref, const Point2f &u,
                               Vector3f *wi, Float *pdf,
                               VisibilityTester *vis) const = 0;
    // pdf of position、pdf of direction
    // 在光源上取一点 pos，在 pos 上取一个方向 dir，生成 ray
    virtual Float Pdf_Li(const Interaction &ref, const Vector3f &wi) const = 0;

    // 光源向外辐射的功率, 具体实现都只计算了近似值
    virtual Spectrum Power() const = 0;
    //在开始渲染前记录场景的一些特征，如 DistanceLight 会计算场景包围盒
    virtual void Preprocess(const Scene &scene) {}

    // 光源向 -r 方向发射的辐射度
    // 当视点发出的光线未与场景中物体相交时，计算光源朝视点方向的辐射度
    virtual Spectrum Le(const RayDifferential &r) const;

    // 从 ref->wi->light 的概率
    virtual Spectrum Sample_Le(const Point2f &u1, const Point2f &u2, Float time,
                               Ray *ray, Normal3f *nLight, Float *pdfPos,
                               Float *pdfDir) const = 0;
    // 传入 ray 和 light_normal，计算 pdfPos，pdfDir
    virtual void Pdf_Le(const Ray &ray, const Normal3f &nLight, Float *pdfPos,
                        Float *pdfDir) const = 0;

    // Light Public Data
    const int flags; // LightFlags, 在计算蒙特卡洛积分时需要区分该光源是不是 delta 分布的
    const int nSamples;  // 光源上采样点的数量, 用于区域(面积)光源
    const MediumInterface mediumInterface; // 光源所处的介质

  protected:
    // Light Protected Data
    const Transform LightToWorld, WorldToLight;
};

class VisibilityTester {
  public:
    VisibilityTester() {}
    // VisibilityTester Public Methods
    VisibilityTester(const Interaction &p0, const Interaction &p1)
        : p0(p0), p1(p1) {}

    // 某些光线传输方法需要这两个点
    const Interaction &P0() const { return p0; }
    const Interaction &P1() const { return p1; }


    bool Unoccluded(const Scene &scene) const;

    // transmittance
    // 处理当两点在介质中时的辐射度
    Spectrum Tr(const Scene &scene, Sampler &sampler) const;

  private:
    Interaction p0, p1;
};

class AreaLight : public Light {
  public:
    // AreaLight Interface
    AreaLight(const Transform &LightToWorld, const MediumInterface &medium,
              int nSamples);

    //根据光源表面一点与表面法线计算出射方向上的发射辐射度Li
    virtual Spectrum L(const Interaction &intr, const Vector3f &w) const = 0;
};

}  // namespace pbrt

#endif  // PBRT_CORE_LIGHT_H
