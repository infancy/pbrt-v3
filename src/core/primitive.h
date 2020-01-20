
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

#ifndef PBRT_CORE_PRIMITIVE_H
#define PBRT_CORE_PRIMITIVE_H

// core/primitive.h*
#include "pbrt.h"
#include "shape.h"
#include "material.h"
#include "medium.h"
#include "transform.h"

namespace pbrt {

// The abstract Primitive base class is the bridge between the geometry processing and shading subsystems of pbrt.

// Primitive Declarations
class Primitive {
  public:
    // Primitive Interface
    virtual ~Primitive();

	// 返回图片在世界空间中的包围盒, 最大的用处是构建空间加速结构
    virtual Bounds3f WorldBound() const = 0;

	// Primitive::Intersect() is responsible for updating Ray::tMax with this value if an intersection is found.
	// 每次求交后都更新 ray.tMax, 避免比该交点更远的求交计算
    virtual bool Intersect(const Ray &r, SurfaceInteraction *) const = 0;
    virtual bool IntersectP(const Ray &r) const = 0;

	// P249
	// GetAreaLight() returns a pointer to the AreaLight that describes the primitive’s emission distribution,
    // if the primitive is itself a light source.If the primitive is not emissive, this method should return nullptr.
	// 如果图元会自发光, 那么在着色计算中也要计入
    virtual const AreaLight *GetAreaLight() const = 0;
	// GetMaterial() returns a pointer to the material instance assigned to the primitive. If
	// nullptris returned, ray intersections with the primitive should be ignored; the primitive
	// only serves to delineate a volume of space for participating media. This method is also
	// used to check if two rays have intersected the same object by comparing their Material pointers.
    virtual const Material *GetMaterial() const = 0;

	// The third material-related method, ComputeScatteringFunctions(), initializes represen-
	// tations of the light-scattering properties of the material at the intersection point on
	// the surface. The BSDF object (introduced in Section 9.1) describes local light-scattering
	// properties at the intersection point.
	// 计算交点处的散射属性(漫反射, 镜面反射...参见 enum BxDFType), 存入 isect->bsdf 中
	// P250: MemoryArena 负责分配内存; TransportMode 记录 ray path 的起点; allowMultipleLobes ...
    virtual void ComputeScatteringFunctions(SurfaceInteraction *isect,
                                            MemoryArena &arena,
                                            TransportMode mode,
                                            bool allowMultipleLobes) const = 0;
};

// The GeometricPrimitive class represents a single shape (e.g., a sphere) in the scene. 
// One GeometricPrimitive is allocated for each shape in the scene description provided by the user. 
// GeometricPrimitive = Shape + Material (+ AreaLight), 分离几何与着色数据
// GeometricPrimitive Declarations
class GeometricPrimitive : public Primitive {
  public:
    // GeometricPrimitive Public Methods
    virtual Bounds3f WorldBound() const;
    virtual bool Intersect(const Ray &r, SurfaceInteraction *isect) const;
    virtual bool IntersectP(const Ray &r) const;
    GeometricPrimitive(const std::shared_ptr<Shape> &shape,
                       const std::shared_ptr<Material> &material,
                       const std::shared_ptr<AreaLight> &areaLight,
                       const MediumInterface &mediumInterface);
    const AreaLight *GetAreaLight() const;
    const Material *GetMaterial() const;
    void ComputeScatteringFunctions(SurfaceInteraction *isect,
                                    MemoryArena &arena, TransportMode mode,
                                    bool allowMultipleLobes) const;

  private:
    // GeometricPrimitive Private Data
    std::shared_ptr<Shape> shape;
    std::shared_ptr<Material> material;
    std::shared_ptr<AreaLight> areaLight;
	// the MediumInterface attribute encodes information about the participating media on the inside and outside of the primitive.
    MediumInterface mediumInterface;
};

// TransformedPrimitive: OBJECT INSTANCING AND ANIMATED PRIMITIVES
// The TransformedPrimitive class handles two more general uses of Shapes in the scene:
// shapes with animated transformation matrices and object instancing, which can greatly
// reduce the memory requirements for scenes that contain many instances of the same
// geometry at different locations (such as the one in Figure 4.1). Implementing each of
// these features essentially requires injecting an additional transformation matrix between
// the Shape’s notion of world space and the actual scene world space. Therefore, both are 
// handled by a single class.

// TransformedPrimitive Declarations
class TransformedPrimitive : public Primitive {
  public:
    // TransformedPrimitive Public Methods
    TransformedPrimitive(std::shared_ptr<Primitive> &primitive,
                         const AnimatedTransform &PrimitiveToWorld);
    bool Intersect(const Ray &r, SurfaceInteraction *in) const;
    bool IntersectP(const Ray &r) const;
    const AreaLight *GetAreaLight() const { return nullptr; }
    const Material *GetMaterial() const { return nullptr; }
    void ComputeScatteringFunctions(SurfaceInteraction *isect,
                                    MemoryArena &arena, TransportMode mode,
                                    bool allowMultipleLobes) const {
		// isect 中记录的是 TransformedPrimitive 持有的 primitive, 所以前者的这三个函数都不会被调用 
        LOG(FATAL) <<
            "TransformedPrimitive::ComputeScatteringFunctions() shouldn't be "
            "called";
    }
    Bounds3f WorldBound() const {
        return PrimitiveToWorld.MotionBounds(primitive->WorldBound());
    }

  private:
    // TransformedPrimitive Private Data
	// The key task of the TransformedPrimitive is to bridge the Primitive interface that it implements 
	// and the Primitive that it holds a pointer to, accounting for the effects of the additional transformation that it holds
	// 可以把这个 primitive 想像成在 TransformedPrimitive World 中, 还需要经过 PrimitiveToWorld 的变换后才能做光线相交等计算
    std::shared_ptr<Primitive> primitive; 
    const AnimatedTransform PrimitiveToWorld;
};

// This chapter also introduces the Aggregate class, which represents a container that can
// hold many Primitives. pbrt uses this class as a base for acceleration structures—data
// structures that help reduce the otherwise O(n) complexity of testing a ray for intersec-
// tion with all n objects in a scene. Most rays will intersect only a few primitives and miss
// the others by a large distance. If an intersection acceleration algorithm can reject whole
// groups of primitives at once, there will be a substantial performance improvement com-
// pared to simply testing each ray against each primitive in turn. One beneﬁt from reusing
// the Primitiveinterface for these acceleration structures is that it makes it easy to support
// hybrid approaches where an accelerator of one type holds accelerators of other types.
// 空间加速结构的抽象接口

// Broadly speaking, there are two main approaches to this problem: spatial subdivision and object subdivision

// Like TransformedPrimitives do, the implementation of the Aggregateintersection methods
// leave the SurfaceInteraction::primitive pointer set to the primitive that the ray actually hit

// Aggregate Declarations
class Aggregate : public Primitive {
  public:
    // Aggregate Public Methods
    const AreaLight *GetAreaLight() const;
    const Material *GetMaterial() const;
    void ComputeScatteringFunctions(SurfaceInteraction *isect,
                                    MemoryArena &arena, TransportMode mode,
                                    bool allowMultipleLobes) const;
};

}  // namespace pbrt

#endif  // PBRT_CORE_PRIMITIVE_H
