
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


    // P716, Figure 12.5
    // 传入(世界空间中的)交点 ref，先假设 ref 到光源间无遮挡, 
    // 计算光源传递给 ref 的入射辐射度 Li 和交点 ref 指向光源的方向 wi; 如果 Li > 0, 则继续构造可见性测试 vis
    // 光线跟踪是一类点采样算法, 当光源是面积光源时，
    // 则需要在面积光源上选取一个采样点(通过传入的随机值 u) 及相应的**概率密度值 pdf**, 通过蒙特卡洛积分来计算整个光源到达 ref 的 Li
    virtual Spectrum Sample_Li(const Interaction &ref, const Point2f &u,
                               Vector3f *wi, Float *pdf,
                               VisibilityTester *vis) const = 0;
    // pdf of position、pdf of direction
    // 在光源上取一点 pos，在 pos 上取一个方向 dir，生成 ray
    virtual Float Pdf_Li(const Interaction &ref, const Vector3f &wi) const = 0;


    // 光源向外辐射的总功率, 因为使用的范围很有限, 具体实现都只计算了近似值
    virtual Spectrum Power() const = 0;
    //在开始渲染前记录场景的一些特征，如 DistanceLight 会通过整个场景的包围盒来计算自己的功率大小
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
    // 这几个成员使用比较频繁, 所以作为常量公开了
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

    // 不考虑参与介质, 计算两点间的可见性
    bool Unoccluded(const Scene &scene) const;
    // 同时考虑阻挡物和参与介质, 计算两点间传递的辐射度
    Spectrum Tr(const Scene &scene, Sampler &sampler) const;

  private:
    Interaction p0, p1;
};

class AreaLight : public Light {
  public:
    // AreaLight Interface
    AreaLight(const Transform &LightToWorld, const MediumInterface &medium,
              int nSamples);

    // 给定区域光源表面上一点 intr, 计算出射方向 wi 上的发射辐射度 Li
    // Spectrum SurfaceInteraction::Le(const Vector3f &w) const {
    //     const AreaLight *area = primitive->GetAreaLight();
    //     return area ? area->L(*this, w) : Spectrum(0.f); }
    virtual Spectrum L(const Interaction &intr, const Vector3f &w) const = 0;
};

}  // namespace pbrt

#endif  // PBRT_CORE_LIGHT_H
