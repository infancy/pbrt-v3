
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

#ifndef PBRT_CORE_SHAPE_H
#define PBRT_CORE_SHAPE_H

/*
pbrt hides details about primitives behind a two-level abstraction. 
The Shape class provides access to the raw geometric properties of the primitive, such as its surface area and bounding box, and provides a ray intersection routine. 
The Primitive class encapsulates additional nongeometric information about the primitive, such as its material properties. The rest of the renderer then deals only with the abstract Primitive interface. This chapter will focus on the geometry-only Shape class; the Primitive interface is a key topic of Chapter 4.
各种 Shape 只表达了场景中对象的几何结构, GeometricPrimitive 则包含了 Shape, Material, AreaLight

All shapes are defined in object coordinate space; for example, all spheres are defined in a coordinate system where the center of the sphere is at the origin.
所有 shape 都被定义在对象坐标空间中(因此 sphere 对象无需定义球心的位置, 球心总是在对象坐标中的原点位置)



*/

// core/shape.h*
#include "pbrt.h"
#include "geometry.h"
#include "interaction.h"
#include "memory.h"
#include "transform.h"

namespace pbrt {

// Shape Declarations
class Shape {
  public:
    // Shape Interface
    Shape(const Transform *ObjectToWorld, const Transform *WorldToObject,
          bool reverseOrientation);
    virtual ~Shape();

	// 分别返回对象坐标系和世界坐标系下的包围盒
	// 使用对象的包围盒进行提前相交测试, 可以节省相交计算的开销
    virtual Bounds3f ObjectBound() const = 0;
    virtual Bounds3f WorldBound() const;

	// 判断 ray 是否与 shape 相交, 若相交则继续计算交点上的微分几何属性
	// 细节见 Section 3.1.3 Intersection Tests
  // tHit: 交点 isect.hit 到 ray 的距离, 可用于提前终止相交测试
    virtual bool Intersect(const Ray &ray, Float *tHit,
                           SurfaceInteraction *isect,
                           bool testAlphaTexture = true) const = 0;
	// 只判断 ray 是否与 shape 相交
    virtual bool IntersectP(const Ray &ray,
                            bool testAlphaTexture = true) const {
        return Intersect(ray, nullptr, nullptr, testAlphaTexture);
    }

    virtual Float Area() const = 0;



    // AreaLight 依赖于 Shape, Sample 和 Pdf 主要在采样 AreaLight 时使用
    // P837, Sample a point on the surface of the shape and return the PDF with
    // respect to area on the surface.
    virtual Interaction Sample(const Point2f &u, Float *pdf) const = 0;
    // Shapes almost always sample uniformly by area on their surface
    virtual Float Pdf(const Interaction &) const { return 1 / Area(); }

    // P837, Sample a point on the shape given a reference point ref and
    // return the PDF with respect to solid angle from ref
    virtual Interaction Sample(const Interaction &ref, const Point2f &u,
                               Float *pdf) const;
    virtual Float Pdf(const Interaction &ref, const Vector3f &wi) const;

    // Returns the solid angle subtended by the shape w.r.t. the reference
    // point p, given in world space. Some shapes compute this value in
    // closed-form, while the default implementation uses Monte Carlo
    // integration; the nSamples parameter determines how many samples are
    // used in this case.
    virtual Float SolidAngle(const Point3f &p, int nSamples = 512) const;



    // Shape Public Data
	// 在 ray-intersection test 时需要将 ray 变换到对象坐标系中, 计算...时需要将 shape 变换到世界坐标系中, shape 同时包含了这两个变换
	// because multiple shapes in the scene will frequently have the same transformation applied to them, pbrt keeps a pool of Transforms so that they can be reused and passes pointers to the shared Transforms to the shapes.
	// As such, the Shape destructor does not delete its Transform pointers, leaving the Transform management code to manage that memory instead.
    const Transform *ObjectToWorld, *WorldToObject;

	// Shapes also take a Boolean parameter, reverseOrientation, that indicates whether their surface normal directions should be reversed from the default. 
	// This capability is useful because the orientation of the surface normal is used to determine which side of a shape is “outside.”
	// For example, shapes that emit illumination are emissive only on the side the surface normal lies on. 
	// The value of this parameter is managed via the ReverseOrientation statement in pbrt input files.
    const bool reverseOrientation; // 是否翻转法线朝向, 或者说定义 shape 的外部
    const bool transformSwapsHandedness; // ObjectToWorld 变换是否改变了坐标系的手性(左手 -> 右手)
};

}  // namespace pbrt

#endif  // PBRT_CORE_SHAPE_H
