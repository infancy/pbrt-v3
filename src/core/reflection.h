
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

#ifndef PBRT_CORE_REFLECTION_H
#define PBRT_CORE_REFLECTION_H

// core/reflection.h*
#include "pbrt.h"
#include "geometry.h"
#include "microfacet.h"
#include "shape.h"
#include "spectrum.h"

namespace pbrt {

// Reflection Declarations

#pragma region 着色坐标系下使用的一些辅助函数

Float FrDielectric(Float cosThetaI, Float etaI, Float etaT);
Spectrum FrConductor(Float cosThetaI, const Spectrum &etaI,
                     const Spectrum &etaT, const Spectrum &k);

// P512
// 传入的 w 都是归一化后的
// BSDF Inline Functions
// P510
inline Float CosTheta(const Vector3f &w) { return w.z; }
inline Float Cos2Theta(const Vector3f &w) { return w.z * w.z; }
inline Float AbsCosTheta(const Vector3f &w) { return std::abs(w.z); }

inline Float Sin2Theta(const Vector3f &w)  { return std::max((Float)0, (Float)1 - Cos2Theta(w)); }
inline Float SinTheta(const Vector3f &w) { return std::sqrt(Sin2Theta(w)); }

inline Float TanTheta(const Vector3f &w) { return SinTheta(w) / CosTheta(w); }
inline Float Tan2Theta(const Vector3f &w) { return Sin2Theta(w) / Cos2Theta(w); }


// P511
inline Float CosPhi(const Vector3f &w) 
{
    Float sinTheta = SinTheta(w);
    return (sinTheta == 0) ? 1 : Clamp(w.x / sinTheta, -1, 1);
}

inline Float SinPhi(const Vector3f &w) 
{
    Float sinTheta = SinTheta(w);
    return (sinTheta == 0) ? 0 : Clamp(w.y / sinTheta, -1, 1);
}

inline Float Cos2Phi(const Vector3f &w) { return CosPhi(w) * CosPhi(w); }

inline Float Sin2Phi(const Vector3f &w) { return SinPhi(w) * SinPhi(w); }

inline Float CosDPhi(const Vector3f &wa, const Vector3f &wb)
{
    return Clamp(
        (wa.x * wb.x + wa.y * wb.y) / std::sqrt((wa.x * wa.x + wa.y * wa.y) *
                                                (wb.x * wb.x + wb.y * wb.y)),
        -1, 1);
}



// P526
inline Vector3f Reflect(const Vector3f &wo, const Vector3f &n) 
{
    return -wo + 2 * Dot(wo, n) * n;
}

// P531
// eta = etaI/etaT
inline bool Refract(const Vector3f &wi, const Normal3f &n, Float eta,
                    Vector3f *wt) 
{
    // Compute $\cos \theta_\roman{t}$ using Snell's law
    Float cosThetaI = Dot(n, wi);
    Float sin2ThetaI = std::max(Float(0), Float(1 - cosThetaI * cosThetaI));
    Float sin2ThetaT = eta * eta * sin2ThetaI;

    // Handle total internal reflection for transmission
    if (sin2ThetaT >= 1) return false;
    Float cosThetaT = std::sqrt(1 - sin2ThetaT);
    *wt = eta * -wi + (eta * cosThetaI - cosThetaT) * Vector3f(n);
    return true;
}

inline bool SameHemisphere(const Vector3f &w, const Vector3f &wp) { return w.z * wp.z > 0; }

inline bool SameHemisphere(const Vector3f &w, const Normal3f &wp) { return w.z * wp.z > 0; }

#pragma endregion



// BSDF Declarations
enum BxDFType 
{
    BSDF_REFLECTION   = 1 << 0,
    BSDF_TRANSMISSION = 1 << 1,

    BSDF_DIFFUSE      = 1 << 2, // 漫反射
    BSDF_GLOSSY       = 1 << 3, // 光泽反射
    BSDF_SPECULAR     = 1 << 4, // 完美镜面反射, delta 分布

    BSDF_ALL = BSDF_DIFFUSE | BSDF_GLOSSY | BSDF_SPECULAR | BSDF_REFLECTION |
               BSDF_TRANSMISSION,
};

struct FourierBSDFTable {
    // FourierBSDFTable Public Data
    Float eta;
    int mMax;
    int nChannels;
    int nMu;
    Float *mu;
    int *m;
    int *aOffset;
    Float *a;
    Float *a0;
    Float *cdf;
    Float *recip;

    // FourierBSDFTable Public Methods
    static bool Read(const std::string &filename, FourierBSDFTable *table);
    const Float *GetAk(int offsetI, int offsetO, int *mptr) const {
        *mptr = m[offsetO * nMu + offsetI];
        return a + aOffset[offsetO * nMu + offsetI];
    }
    bool GetWeightsAndOffset(Float cosTheta, int *offset,
                             Float weights[4]) const;
};



#pragma region BSDF / BxDF / ScaledBxDF

// 组合了多种 BxDF 的效果
class BSDF {
  public:
    // BSDF Public Methods
    BSDF(const SurfaceInteraction &si, Float eta = 1)
        : eta(eta),
          ng(si.n),
          ns(si.shading.n),
          ss(Normalize(si.shading.dpdu)),
          ts(Cross(ns, ss)) {}

    void Add(BxDF *b) {
        CHECK_LT(nBxDFs, MaxBxDFs);
        bxdfs[nBxDFs++] = b;
    }
    int NumComponents(BxDFType flags = BSDF_ALL) const;

    Vector3f WorldToLocal(const Vector3f &v) const 
    {
        return Vector3f(Dot(v, ss), Dot(v, ts), Dot(v, ns));
    }
    Vector3f LocalToWorld(const Vector3f &v) const 
    {
        // P574, 因为 stn 是正交矩阵，其逆矩阵即为转置矩阵
        return Vector3f(ss.x * v.x + ts.x * v.y + ns.x * v.z,
                        ss.y * v.x + ts.y * v.y + ns.y * v.z,
                        ss.z * v.x + ts.z * v.y + ns.z * v.z);
    }

    Spectrum f(const Vector3f &woW, const Vector3f &wiW,
               BxDFType flags = BSDF_ALL) const;

    Spectrum rho(int nSamples, const Point2f *samples1, const Point2f *samples2,
                 BxDFType flags = BSDF_ALL) const;
    Spectrum rho(const Vector3f &wo, int nSamples, const Point2f *samples,
                 BxDFType flags = BSDF_ALL) const;

    // 在多重重要性采样中, 对 BSDF 进行采样
    Spectrum Sample_f(const Vector3f &wo, Vector3f *wi, const Point2f &u,
                      Float *pdf, BxDFType type = BSDF_ALL,
                      BxDFType *sampledType = nullptr) const;
    Float Pdf(const Vector3f &wo, const Vector3f &wi,
              BxDFType flags = BSDF_ALL) const;

    std::string ToString() const;

    // BSDF Public Data
    const Float eta; // 表面相对折射率，只用于半透明物体

  private:
    // BSDF Private Methods
    // 析构函数是私有成员, 保证不被意外销毁(实践上这个函数不会在任何地方被调用, 而是在 MemoryArena 对象析构时直接清理掉内存)
    ~BSDF() {}

    // BSDF Private Data
    // ns,ss,ts构成了着色坐标系, 所有着色计算都是在这个坐标系下进行的
    // 在这个凹凸贴图(Bump Map)后, 着色坐标系会和原始的几何坐标系有差异
    const Normal3f ns, ng; // shading-normal/geometry-normal
    const Vector3f ss, ts; 

    // holds a set of BxDFs whose contributions are summed to give the full scattering function
    int nBxDFs = 0;
    static PBRT_CONSTEXPR int MaxBxDFs = 8;
    BxDF *bxdfs[MaxBxDFs];

    // MixMaterial 是为唯一需要直接访问 bxdfs 的对象, 这里没有另外加接口, 而是声明成友元
    friend class MixMaterial;
};

inline std::ostream &operator<<(std::ostream &os, const BSDF &bsdf) {
    os << bsdf.ToString();
    return os;
}

// BxDF Declarations
// BRDF/BTDF 的公共抽象接口
class BxDF {
  public:
    // BxDF Interface
    virtual ~BxDF() {}
    BxDF(BxDFType type) : type(type) {}

    // has_reflection, has_transmission, has_specular...
    bool MatchesFlags(BxDFType t) const { return (type & t) == type; }

    // P514, 注意返回的是 Spectrum,  每个分量代表了这个光谱上的...
    virtual Spectrum f(const Vector3f &wo, const Vector3f &wi) const = 0;

    // prev   light
    // -----  -----
    //   ^      ^
    //    \    /
    //  wo \  / wi
    //      \/
    //    ------
    //    isect
    // P806, 并非所有 BxDF 的 f(wo, wi) 函数都可用, 例如镜子, 玻璃, 水面的完美镜面分布(delta 分布, 参考 Chapter7.1 的狄拉克函数)
    // 这个时候就需要用 Sample_f, 给定出射方向 wo, 计算入射光方向 wi 并返回对应的 f(wo, wi)
    // sample 用于在光源上随机采样, 对应于概率密度 pdf, 这两者在非 delta 分布时才会用到
    virtual Spectrum Sample_f(const Vector3f &wo, Vector3f *wi,
                              const Point2f &sample, Float *pdf,
                              BxDFType *sampledType = nullptr) const;

    // P514, P815, P830, 针对某些无法通过闭式计算反射率的 BxDF，可用 rho_hd 来估算（samples[] 用于蒙特卡洛积分）
    // hemisphere_direction_reflectance
    virtual Spectrum rho(const Vector3f &wo, int nSamples,
                         const Point2f *samples) const;
    // hemisphere_hemisphere_reflectance
    virtual Spectrum rho(int nSamples, const Point2f *samples1,
                         const Point2f *samples2) const;

    // P807, 如果派生类覆盖了 Sample_f, 那么也需要覆盖 Pdf
    virtual Float Pdf(const Vector3f &wo, const Vector3f &wi) const;

    virtual std::string ToString() const = 0;

    // BxDF Public Data
	// 某些光线传输算法需要对 BRDF 和 BTDF 进行区分，所以加入 type 成员
    const BxDFType type;
};

inline std::ostream &operator<<(std::ostream &os, const BxDF &bxdf) {
    os << bxdf.ToString();
    return os;
}

// 对持有的 bxdf 进行放缩, 用在 MixMaterial 中
class ScaledBxDF : public BxDF {
  public:
    // ScaledBxDF Public Methods
    ScaledBxDF(BxDF *bxdf, const Spectrum &scale)
        : BxDF(BxDFType(bxdf->type)), bxdf(bxdf), scale(scale) {}
    Spectrum rho(const Vector3f &w, int nSamples,
                 const Point2f *samples) const {
        return scale * bxdf->rho(w, nSamples, samples);
    }
    Spectrum rho(int nSamples, const Point2f *samples1,
                 const Point2f *samples2) const {
        return scale * bxdf->rho(nSamples, samples1, samples2);
    }
    Spectrum f(const Vector3f &wo, const Vector3f &wi) const;
    Spectrum Sample_f(const Vector3f &wo, Vector3f *wi, const Point2f &sample,
                      Float *pdf, BxDFType *sampledType) const;
    Float Pdf(const Vector3f &wo, const Vector3f &wi) const;
    std::string ToString() const;

  private:
    BxDF *bxdf;
    Spectrum scale;
};

#pragma endregion


#pragma region Fresenl / Specular / FresnelSpecular

class Fresnel {
  public:
    // Fresnel Interface
    virtual ~Fresnel();
    // 给定入射光线和表面法线夹角的 cos 值, 计算 Fresenl 反射系数
    virtual Spectrum Evaluate(Float cosI) const = 0;
    virtual std::string ToString() const = 0;
};

inline std::ostream &operator<<(std::ostream &os, const Fresnel &f) {
    os << f.ToString();
    return os;
}

class FresnelConductor : public Fresnel {
  public:
    // FresnelConductor Public Methods
    Spectrum Evaluate(Float cosThetaI) const;
    FresnelConductor(const Spectrum &etaI, const Spectrum &etaT,
                     const Spectrum &k)
        : etaI(etaI), etaT(etaT), k(k) {}
    std::string ToString() const;

  private:
    Spectrum etaI, etaT, k;
};

class FresnelDielectric : public Fresnel {
  public:
    // FresnelDielectric Public Methods
    Spectrum Evaluate(Float cosThetaI) const;
    FresnelDielectric(Float etaI, Float etaT) : etaI(etaI), etaT(etaT) {}
    std::string ToString() const;

  private:
    Float etaI, etaT;
};

// 完全反射光线, 可以用在 mirror 这样的材质上
class FresnelNoOp : public Fresnel {
  public:
    Spectrum Evaluate(Float) const { return Spectrum(1.); }
    std::string ToString() const { return "[ FresnelNoOp ]"; }
};



class SpecularReflection : public BxDF {
  public:
    // SpecularReflection Public Methods
    SpecularReflection(const Spectrum &R, Fresnel *fresnel)
        : BxDF(BxDFType(BSDF_REFLECTION | BSDF_SPECULAR)),
          R(R),
          fresnel(fresnel) {}

    Spectrum f(const Vector3f &wo, const Vector3f &wi) const 
    {
        return Spectrum(0.f); // 因为 PBRT 在渲染阶段对这种 delta 分布的 BxDF 做了特殊处理(使用 Sample_f 来采样), 所以这里总是让其采样到的概率为 0
    }
    Spectrum Sample_f(const Vector3f &wo, Vector3f *wi, const Point2f &sample,
                      Float *pdf, BxDFType *sampledType) const;
    Float Pdf(const Vector3f &wo, const Vector3f &wi) const { return 0; }
    std::string ToString() const;

  private:
    // SpecularReflection Private Data
    const Spectrum R; // used to scale the reflected color
    const Fresnel *fresnel;
};

class SpecularTransmission : public BxDF {
  public:
    // SpecularTransmission Public Methods
    SpecularTransmission(const Spectrum &T, Float etaA, Float etaB,
                         TransportMode mode)
        : BxDF(BxDFType(BSDF_TRANSMISSION | BSDF_SPECULAR)),
          T(T),
          etaA(etaA),
          etaB(etaB),
          fresnel(etaA, etaB),
          mode(mode) {}

    Spectrum f(const Vector3f &wo, const Vector3f &wi) const {
        return Spectrum(0.f);
    }
    Spectrum Sample_f(const Vector3f &wo, Vector3f *wi, const Point2f &sample,
                      Float *pdf, BxDFType *sampledType) const;
    Float Pdf(const Vector3f &wo, const Vector3f &wi) const { return 0; }
    std::string ToString() const;

  private:
    // SpecularTransmission Private Data
    const Spectrum T; // transmission scale factor
    const Float etaA, etaB; // 入射介质, 折射介质的折射率
    const FresnelDielectric fresnel; // 只有两种绝缘体之间才有折射现象
    const TransportMode mode;
};


// 包含了镜面反射和镜面折射效果, 两种效果的比例由 Fresnel 方程来控制
class FresnelSpecular : public BxDF {
  public:
    // FresnelSpecular Public Methods
    FresnelSpecular(const Spectrum &R, const Spectrum &T, Float etaA,
                    Float etaB, TransportMode mode)
        : BxDF(BxDFType(BSDF_REFLECTION | BSDF_TRANSMISSION | BSDF_SPECULAR)),
          R(R),
          T(T),
          etaA(etaA),
          etaB(etaB),
          mode(mode) {}
    Spectrum f(const Vector3f &wo, const Vector3f &wi) const {
        return Spectrum(0.f);
    }
    Spectrum Sample_f(const Vector3f &wo, Vector3f *wi, const Point2f &u,
                      Float *pdf, BxDFType *sampledType) const;
    Float Pdf(const Vector3f &wo, const Vector3f &wi) const { return 0; } // 对完美镜面的 delta 分布有特殊处理, 这里干脆设为 0 了
    std::string ToString() const;

  private:
    // FresnelSpecular Private Data
    const Spectrum R, T;
    const Float etaA, etaB;
    const TransportMode mode;
};

#pragma endregion



#pragma region Lambertian

// models a perfect diffuse surface that scatters incident illumination equally in all directions.
class LambertianReflection : public BxDF {
  public:
    // LambertianReflection Public Methods
    LambertianReflection(const Spectrum &R)
        : BxDF(BxDFType(BSDF_REFLECTION | BSDF_DIFFUSE)), R(R) {}
    Spectrum f(const Vector3f &wo, const Vector3f &wi) const;

    Spectrum rho(const Vector3f &, int, const Point2f *) const { return R; }
    Spectrum rho(int, const Point2f *, const Point2f *) const { return R; }

    std::string ToString() const;

  private:
    // LambertianReflection Private Data
    const Spectrum R;
};

class LambertianTransmission : public BxDF {
  public:
    // LambertianTransmission Public Methods
    LambertianTransmission(const Spectrum &T)
        : BxDF(BxDFType(BSDF_TRANSMISSION | BSDF_DIFFUSE)), T(T) {}
    Spectrum f(const Vector3f &wo, const Vector3f &wi) const;

    Spectrum rho(const Vector3f &, int, const Point2f *) const { return T; }
    Spectrum rho(int, const Point2f *, const Point2f *) const { return T; }

    Spectrum Sample_f(const Vector3f &wo, Vector3f *wi, const Point2f &u,
                      Float *pdf, BxDFType *sampledType) const;
    Float Pdf(const Vector3f &wo, const Vector3f &wi) const;
    std::string ToString() const;

  private:
    // LambertianTransmission Private Data
    Spectrum T;
};

#pragma endregion



class OrenNayar : public BxDF {
  public:
    // OrenNayar Public Methods
    Spectrum f(const Vector3f &wo, const Vector3f &wi) const;

    // sigma, the standard deviation of the microfacet orientation angle(???), 为 0 时等同于 LambertianReflection
    OrenNayar(const Spectrum &R, Float sigma)
        : BxDF(BxDFType(BSDF_REFLECTION | BSDF_DIFFUSE)), R(R) 
    {
        // P535, 原始的分布函数只有一个参数 sigma, 这里的 A/B 是预计算的中间量
        sigma = Radians(sigma);
        Float sigma2 = sigma * sigma;
        A = 1.f - (sigma2 / (2.f * (sigma2 + 0.33f)));
        B = 0.45f * sigma2 / (sigma2 + 0.09f);
    }
    std::string ToString() const;

  private:
    // OrenNayar Private Data
    const Spectrum R;
    Float A, B;
};



// Torrance–Sparrow model(Cook-Torrance BRDF)
class MicrofacetReflection : public BxDF {
  public:
    // MicrofacetReflection Public Methods
    MicrofacetReflection(const Spectrum &R,
                         MicrofacetDistribution *distribution, Fresnel *fresnel)
        : BxDF(BxDFType(BSDF_REFLECTION | BSDF_GLOSSY)),
          R(R),
          distribution(distribution),
          fresnel(fresnel) {}

    Spectrum f(const Vector3f &wo, const Vector3f &wi) const;
    Spectrum Sample_f(const Vector3f &wo, Vector3f *wi, const Point2f &u,
                      Float *pdf, BxDFType *sampledType) const;
    Float Pdf(const Vector3f &wo, const Vector3f &wi) const;
    std::string ToString() const;

  private:
    // MicrofacetReflection Private Data
    const Spectrum R;
    const MicrofacetDistribution *distribution;
    const Fresnel *fresnel;
};

class MicrofacetTransmission : public BxDF {
  public:
    // MicrofacetTransmission Public Methods
    MicrofacetTransmission(const Spectrum &T,
                           MicrofacetDistribution *distribution, Float etaA,
                           Float etaB, TransportMode mode)
        : BxDF(BxDFType(BSDF_TRANSMISSION | BSDF_GLOSSY)),
          T(T),
          distribution(distribution),
          etaA(etaA),
          etaB(etaB),
          fresnel(etaA, etaB),
          mode(mode) {}

    Spectrum f(const Vector3f &wo, const Vector3f &wi) const;
    Spectrum Sample_f(const Vector3f &wo, Vector3f *wi, const Point2f &u,
                      Float *pdf, BxDFType *sampledType) const;
    Float Pdf(const Vector3f &wo, const Vector3f &wi) const;
    std::string ToString() const;

  private:
    // MicrofacetTransmission Private Data
    const Spectrum T;
    const MicrofacetDistribution *distribution;
    const Float etaA, etaB;
    const FresnelDielectric fresnel;
    const TransportMode mode;
};



// Fresnel Incidence Effects, models a diffuse underlying surface with a glossy specular surface above it
// 当入射角度接近表面法线时, 大多数光线被透射并漫反射; 当接近垂直于法线的方向时, 则主要发生光泽反射
class FresnelBlend : public BxDF {
  public:
    // FresnelBlend Public Methods
    FresnelBlend(const Spectrum &Rd, const Spectrum &Rs,
                 MicrofacetDistribution *distrib);

    Spectrum f(const Vector3f &wo, const Vector3f &wi) const;

    // Fresnel 的近似形式, cosTheta 是入射光线 wi 和半法线 wh 的夹角
    Spectrum SchlickFresnel(Float cosTheta) const {
        auto pow5 = [](Float v) { return (v * v) * (v * v) * v; };
        return Rs + pow5(1 - cosTheta) * (Spectrum(1.) - Rs);
    }
    Spectrum Sample_f(const Vector3f &wi, Vector3f *sampled_f, const Point2f &u,
                      Float *pdf, BxDFType *sampledType) const;
    Float Pdf(const Vector3f &wo, const Vector3f &wi) const;
    std::string ToString() const;

  private:
    // FresnelBlend Private Data
    const Spectrum Rd, Rs; //  diffuse and specular reflectance
    MicrofacetDistribution *distribution;
};



class FourierBSDF : public BxDF {
  public:
    // FourierBSDF Public Methods
    Spectrum f(const Vector3f &wo, const Vector3f &wi) const;
    FourierBSDF(const FourierBSDFTable &bsdfTable, TransportMode mode)
        : BxDF(BxDFType(BSDF_REFLECTION | BSDF_TRANSMISSION | BSDF_GLOSSY)),
          bsdfTable(bsdfTable),
          mode(mode) {}
    Spectrum Sample_f(const Vector3f &wo, Vector3f *wi, const Point2f &u,
                      Float *pdf, BxDFType *sampledType) const;
    Float Pdf(const Vector3f &wo, const Vector3f &wi) const;
    std::string ToString() const;

  private:
    // FourierBSDF Private Data
    const FourierBSDFTable &bsdfTable;
    const TransportMode mode;
};



// BSDF Inline Method Definitions
inline int BSDF::NumComponents(BxDFType flags) const {
    int num = 0;
    for (int i = 0; i < nBxDFs; ++i)
        if (bxdfs[i]->MatchesFlags(flags)) ++num;
    return num;
}

}  // namespace pbrt

#endif  // PBRT_CORE_REFLECTION_H
