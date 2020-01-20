
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

#ifndef PBRT_CORE_TEXTURE_H
#define PBRT_CORE_TEXTURE_H

// core/texture.h*
#include "pbrt.h"
#include "spectrum.h"
#include "geometry.h"
#include "transform.h"
#include "memory.h"

namespace pbrt {

// TextureMapping2D/3D 通过 SurfaceInteraction 的参数化坐标 uv([0, 1]^2), 生成纹理坐标 st([0, 1] ^ 2)
// 然后 Texture 根据 st 坐标及其对屏幕空间 xy 坐标的变化率进行采样
// uv -> st -> value

// Texture Declarations
class TextureMapping2D {
  public:
    // TextureMapping2D Interface
    virtual ~TextureMapping2D();
    // It also returns estimates for the change in (s, t) with respect
    // to pixel x and y coordinates in the dstdx and dstdy parameters so that textures that use
    // the mapping can determine the (s, t) sampling rate and filter accordingly.
    // 根据表面交点 si 上的信息生成纹理坐标 st 的同时, 
    // 还需要根据表面参数化 uv 坐标相对屏幕空间 xy 坐标的变化率,
    // 计算纹理空间相对屏幕空间的变化率
    virtual Point2f Map(const SurfaceInteraction &si, Vector2f *dstdx,
                        Vector2f *dstdy) const = 0;
};

class UVMapping2D : public TextureMapping2D {
  public:
    // UVMapping2D Public Methods
    UVMapping2D(Float su = 1, Float sv = 1, Float du = 0, Float dv = 0);
    Point2f Map(const SurfaceInteraction &si, Vector2f *dstdx,
                Vector2f *dstdy) const;

  private:
    const Float su, sv, // scale
                du, dv; // offset
};

class SphericalMapping2D : public TextureMapping2D {
  public:
    // SphericalMapping2D Public Methods
    SphericalMapping2D(const Transform &WorldToTexture)
        : WorldToTexture(WorldToTexture) {}
    Point2f Map(const SurfaceInteraction &si, Vector2f *dstdx,
                Vector2f *dstdy) const;

  private:
    Point2f sphere(const Point3f &P) const;
    const Transform WorldToTexture;
};

class CylindricalMapping2D : public TextureMapping2D {
  public:
    // CylindricalMapping2D Public Methods
    CylindricalMapping2D(const Transform &WorldToTexture)
        : WorldToTexture(WorldToTexture) {}
    Point2f Map(const SurfaceInteraction &si, Vector2f *dstdx,
                Vector2f *dstdy) const;

  private:
    // CylindricalMapping2D Private Methods
    Point2f cylinder(const Point3f &p) const {
        Vector3f vec = Normalize(WorldToTexture(p) - Point3f(0, 0, 0));
        return Point2f((Pi + std::atan2(vec.y, vec.x)) * Inv2Pi, vec.z);
    }
    const Transform WorldToTexture;
};

class PlanarMapping2D : public TextureMapping2D {
  public:
    // PlanarMapping2D Public Methods
    Point2f Map(const SurfaceInteraction &si, Vector2f *dstdx,
                Vector2f *dstdy) const;
    PlanarMapping2D(const Vector3f &vs, const Vector3f &vt, Float ds = 0,
                    Float dt = 0)
        : vs(vs), vt(vt), ds(ds), dt(dt) {}

  private:
    const Vector3f vs, vt;
    const Float ds, dt;
};



class TextureMapping3D {
  public:
    // TextureMapping3D Interface
    virtual ~TextureMapping3D();
    virtual Point3f Map(const SurfaceInteraction &si, Vector3f *dpdx,
                        Vector3f *dpdy) const = 0;
};

class IdentityMapping3D : public TextureMapping3D {
  public:
    // IdentityMapping3D Public Methods
    IdentityMapping3D(const Transform &WorldToTexture)
        : WorldToTexture(WorldToTexture) {}
    Point3f Map(const SurfaceInteraction &si, Vector3f *dpdx,
                Vector3f *dpdy) const;

  private:
    const Transform WorldToTexture; // SurfaceInteraction 是定义在世界空间中的
};



template <typename T>
class Texture {
  public:
    // Texture Interface
    virtual T Evaluate(const SurfaceInteraction &) const = 0;
    virtual ~Texture() {}
};



Float Lanczos(Float, Float tau = 2);
Float Noise(Float x, Float y = .5f, Float z = .5f);
Float Noise(const Point3f &p);
Float FBm(const Point3f &p, const Vector3f &dpdx, const Vector3f &dpdy,
          Float omega, int octaves);
Float Turbulence(const Point3f &p, const Vector3f &dpdx, const Vector3f &dpdy,
                 Float omega, int octaves);

}  // namespace pbrt

#endif  // PBRT_CORE_TEXTURE_H
