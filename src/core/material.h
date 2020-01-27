
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

#ifndef PBRT_CORE_MATERIAL_H
#define PBRT_CORE_MATERIAL_H

// core/material.h*
#include "pbrt.h"
#include "memory.h"

namespace pbrt {
	
// indicates whether the ray path that found this intersection point started from the camera or from a light source
// TransportMode Declarations
enum class TransportMode 
{ 
    Radiance,  // 发射量(light -> ... -> camera)
    Importance // 接收量(camera -> ... -> light)
};

// Material Declarations
// P571, by separating these two components and having the Material return a BSDF, pbrt is better able to handle a variety of light transport algorithms.
class Material {
  public:
    // Material Interface
    //  takes a point on a surface and creates a BSDF object (and that represents scattering at the point
    virtual void ComputeScatteringFunctions(
        SurfaceInteraction *si, // contains geometric properties at an intersection point on the surface of a shape.
        MemoryArena &arena,     // allocate memory
        TransportMode mode,     // path starting from the camera or one starting from a light source
        bool allowMultipleLobes // aggregate multiple types of scattering into a single BxDF when such BxDFs are available
    ) const = 0;

    virtual ~Material();
    static void Bump(const std::shared_ptr<Texture<Float>> &d,
                     SurfaceInteraction *si);
};

}  // namespace pbrt

#endif  // PBRT_CORE_MATERIAL_H
