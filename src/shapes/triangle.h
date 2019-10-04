
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

#ifndef PBRT_SHAPES_TRIANGLE_H
#define PBRT_SHAPES_TRIANGLE_H

// Triangles have a dual role among the shapes in pbrt: 
// not only are they frequently directly speci?ed in scene description ?les, 
// but other shapes often tessellate themselves into triangle meshes.

// shapes/triangle.h*
#include "shape.h"
#include "stats.h"
#include <map>

namespace pbrt {

// triMeshBytes 记录三角形所有数据所占据的内存大小
STAT_MEMORY_COUNTER("Memory/Triangle meshes", triMeshBytes);

// TriangleMesh 使用顶点属性列表+索引列表的方式, 存储了一个三角形网格的所有数据
// Triangle Declarations
struct TriangleMesh {
    // TriangleMesh Public Methods
    TriangleMesh(const Transform &ObjectToWorld, int nTriangles,
                 const int *vertexIndices, int nVertices, const Point3f *P,
                 const Vector3f *S, const Normal3f *N, const Point2f *uv,
                 const std::shared_ptr<Texture<Float>> &alphaMask,
                 const std::shared_ptr<Texture<Float>> &shadowAlphaMask,
                 const int *faceIndices);

    // TriangleMesh Data
    const int nTriangles, nVertices;
    std::vector<int> vertexIndices; // vertexIndices.size() == 3 * nTriangles

    // 顶点属性(vertex attribute)
    // For the ith triangle, its three vertex positions are 
    // P[vertexIndices[3*i]], P[vertexIndices[3*i+1]], and P[vertexIndices[3*i+2]].
    std::unique_ptr<Point3f[]> p;
    std::unique_ptr<Normal3f[]> n;
    std::unique_ptr<Vector3f[]> s; // 只记录切向量 tangent, 另一个切向量 bitangent 可以由 tangent 和 normal 叉乘得到
    std::unique_ptr<Point2f[]> uv;

    std::shared_ptr<Texture<Float>> alphaMask, shadowAlphaMask;
    std::vector<int> faceIndices; // 面数的索引
};

// The Triangleclass actually implements the Shapeinterface. **It represents a single triangle**
class Triangle : public Shape {
  public:
    // Triangle Public Methods
    Triangle(const Transform *ObjectToWorld, const Transform *WorldToObject,
             bool reverseOrientation, const std::shared_ptr<TriangleMesh> &mesh,
             int triNumber) // triNumber 改为 triIndex 更好理解些
        : Shape(ObjectToWorld, WorldToObject, reverseOrientation), mesh(mesh) {
        v = &mesh->vertexIndices[3 * triNumber]; // 只存储首个顶点索引的地址, 以节省存储空间
        triMeshBytes += sizeof(*this);
        faceIndex = mesh->faceIndices.size() ? mesh->faceIndices[triNumber] : 0;
    }

    Bounds3f ObjectBound() const override;
    Bounds3f WorldBound() const override;

    bool Intersect(const Ray &ray, Float *tHit, SurfaceInteraction *isect,
                   bool testAlphaTexture = true) const;
    bool IntersectP(const Ray &ray, bool testAlphaTexture = true) const;

    Float Area() const;

    using Shape::Sample;  // Bring in the other Sample() overload.
    Interaction Sample(const Point2f &u, Float *pdf) const;

    // Returns the solid angle subtended by the triangle w.r.t. the given
    // reference point p.
    Float SolidAngle(const Point3f &p, int nSamples = 0) const;

  private:
    // Triangle Private Methods
    // 分别获得三个顶点的 UV 坐标
    void GetUVs(Point2f uv[3]) const {
        if (mesh->uv) {
            uv[0] = mesh->uv[v[0]];
            uv[1] = mesh->uv[v[1]];
            uv[2] = mesh->uv[v[2]];
        } else {
            uv[0] = Point2f(0, 0);
            uv[1] = Point2f(1, 0);
            uv[2] = Point2f(1, 1);
        }
    }

    // Triangle Private Data
    std::shared_ptr<TriangleMesh> mesh;
    const int *v;
    int faceIndex;
};

std::vector<std::shared_ptr<Shape>> CreateTriangleMesh(
    const Transform *o2w, const Transform *w2o, bool reverseOrientation,
    int nTriangles, const int *vertexIndices, int nVertices, const Point3f *p,
    const Vector3f *s, const Normal3f *n, const Point2f *uv,
    const std::shared_ptr<Texture<Float>> &alphaTexture,
    const std::shared_ptr<Texture<Float>> &shadowAlphaTexture,
    const int *faceIndices = nullptr);
std::vector<std::shared_ptr<Shape>> CreateTriangleMeshShape(
    const Transform *o2w, const Transform *w2o, bool reverseOrientation,
    const ParamSet &params,
    std::map<std::string, std::shared_ptr<Texture<Float>>> *floatTextures =
        nullptr);

bool WritePlyFile(const std::string &filename, int nTriangles,
                  const int *vertexIndices, int nVertices, const Point3f *P,
                  const Vector3f *S, const Normal3f *N, const Point2f *UV,
                  const int *faceIndices);

}  // namespace pbrt

#endif  // PBRT_SHAPES_TRIANGLE_H
