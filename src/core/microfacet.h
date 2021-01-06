
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

#ifndef PBRT_CORE_MICROFACET_H
#define PBRT_CORE_MICROFACET_H

// core/microfacet.h*
#include "pbrt.h"
#include "geometry.h"

namespace pbrt {

// 由微面元构成的表面可视为高度场, 面元分布可用统计方式来描述

// MicrofacetDistribution Declarations
class MicrofacetDistribution {
  public:
    // MicrofacetDistribution Public Methods
    virtual ~MicrofacetDistribution();

    // 微面元分布函数/法线分布函数
    // return the differential area of microfacets oriented with the given normal vector $omega$
    // 参考 P540 Figure8.16 的曲线
    // 1. 表面粗糙度越大时, 曲线越平缓
    // 2. 给定法线方向(其实是 wi 和 wo 的半法线)和表面法线夹角越大时, 朝向该方向的微面元越少, 反射光线越少
    virtual Float D(const Vector3f &wh) const = 0;

    // 几何遮挡函数的辅助函数
    virtual Float Lambda(const Vector3f &w) const = 0;

    // P541, the Smith’s masking-shadowing function G1 gives the fraction of microfacets with normal ω_h that are visible from direction ω(wo or wi)
    // 大多数情况下微面元的可见性和其朝向是独立(无关)的, 所以这里省略了其法线 ω_h
    // 在大部分方向上(除了接近垂直于表面法线的方向), 其返回值都接近 1
    Float G1(const Vector3f &w) const 
    {
        // if (Dot(w, wh) * CosTheta(w) < 0.) return 0.;
        return 1 / (1 + Lambda(w)); // P541 G_1(ω) 的简化形式
    }

    // 双向几何遮挡函数, 描述了朝向半矢量法线的微面元中, 可以**同时被入射方向和反射方向**观察的面元的比例
    // Smith’s masking-shadowing function G 假设入射遮蔽(shadowing)和出射遮挡(masking)是不相关的, 两个部分被拆分开独立计算, G = G1(wi) * G1(wo)
    // A more accurate model can be derived assuming that microfacet visibility is more likely the higher up a given point on a microfacet is.
    // This assumption leads to the model: 
    virtual Float G(const Vector3f &wo, const Vector3f &wi) const 
    {
        return 1 / (1 + Lambda(wo) + Lambda(wi));
    }

    // 根据出射方向 wo 采样法线
    virtual Vector3f Sample_wh(const Vector3f &wo, const Point2f &u) const = 0;

    Float Pdf(const Vector3f &wo, const Vector3f &wh) const;

    virtual std::string ToString() const = 0;

  protected:
    // MicrofacetDistribution Protected Methods
    // P808, 默认只采样可见区域
    MicrofacetDistribution(bool sampleVisibleArea)
        : sampleVisibleArea(sampleVisibleArea) {}

    // MicrofacetDistribution Protected Data
    const bool sampleVisibleArea;
};

inline std::ostream &operator<<(std::ostream &os,
                                const MicrofacetDistribution &md) {
    os << md.ToString();
    return os;
}

class BeckmannDistribution : public MicrofacetDistribution {
  public:
    // BeckmannDistribution Public Methods
    // P540, roughness 从 [0, 1] 更符合直觉, 方便用户调节, 这个函数负责把 roughness 从 [0, 1] 映射到 ???
    static Float RoughnessToAlpha(Float roughness) {
        roughness = std::max(roughness, (Float)1e-3);
        Float x = std::log(roughness);
        return 1.62142f + 0.819955f * x + 0.1734f * x * x +
               0.0171201f * x * x * x + 0.000640711f * x * x * x * x;
    }
    BeckmannDistribution(Float alphax, Float alphay, bool samplevis = true)
        : MicrofacetDistribution(samplevis),
          alphax(std::max(Float(0.001), alphax)),
          alphay(std::max(Float(0.001), alphay)) {}

    Float D(const Vector3f &wh) const;
    Vector3f Sample_wh(const Vector3f &wo, const Point2f &u) const;
    std::string ToString() const;

  private:
    // BeckmannDistribution Private Methods
    Float Lambda(const Vector3f &w) const;

    // BeckmannDistribution Private Data
    const Float alphax, alphay; // 控制各向异性分布的椭圆形状
};

class TrowbridgeReitzDistribution : public MicrofacetDistribution {
  public:
    // TrowbridgeReitzDistribution Public Methods
    static inline Float RoughnessToAlpha(Float roughness);
    TrowbridgeReitzDistribution(Float alphax, Float alphay,
                                bool samplevis = true)
        : MicrofacetDistribution(samplevis),
          alphax(std::max(Float(0.001), alphax)),
          alphay(std::max(Float(0.001), alphay)) {}

    Float D(const Vector3f &wh) const;
    Vector3f Sample_wh(const Vector3f &wo, const Point2f &u) const;
    std::string ToString() const;

  private:
    // TrowbridgeReitzDistribution Private Methods
    Float Lambda(const Vector3f &w) const;

    // TrowbridgeReitzDistribution Private Data
    const Float alphax, alphay;
};

// MicrofacetDistribution Inline Methods
inline Float TrowbridgeReitzDistribution::RoughnessToAlpha(Float roughness) {
    roughness = std::max(roughness, (Float)1e-3);
    Float x = std::log(roughness);
    return 1.62142f + 0.819955f * x + 0.1734f * x * x + 0.0171201f * x * x * x +
           0.000640711f * x * x * x * x;
}

}  // namespace pbrt

#endif  // PBRT_CORE_MICROFACET_H
