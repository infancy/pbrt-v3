
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

#ifndef PBRT_CAMERAS_PERSPECTIVE_H
#define PBRT_CAMERAS_PERSPECTIVE_H

/*
侧视示意图, lens 一般是圆形的
                                                                 focal
                                                                *|
                                                        *   *    |
                                                 *     *         |
|                                         *       *              |
|            lens         |        *         *                   |
|              |          | *           *                        |
|              |     *    |        *                             |
|             *|          |    *                                 |
|           *  |          |*                                     |
|         *    |      *   |                                      |
|       *      |  *       |                                      |
|     *       *|          |
|   *     *    |          |
| *  *         |          |
|*             |          |
|                         |
|                     view(视平面)
image(像平面)

通过穿过透镜中心的中心光线确定位置，而后生成主光线进行着色
虽然主光线可能未通过视平面的相应像素的位置，但这并没有关系
当薄透镜的半径为 0 时，则透镜相机变成了针孔相机
*/

// cameras/perspective.h*
#include "pbrt.h"
#include "camera.h"
#include "film.h"

namespace pbrt {

// PerspectiveCamera Declarations
class PerspectiveCamera : public ProjectiveCamera {
  public:
    // PerspectiveCamera Public Methods
    PerspectiveCamera(const AnimatedTransform &CameraToWorld,
                      const Bounds2f &screenWindow, Float shutterOpen,
                      Float shutterClose, Float lensRadius, Float focalDistance,
                      Float fov, Film *film, const Medium *medium);
    Float GenerateRay(const CameraSample &sample, Ray *) const;
    Float GenerateRayDifferential(const CameraSample &sample,
                                  RayDifferential *ray) const;
    Spectrum We(const Ray &ray, Point2f *pRaster2 = nullptr) const;
    void Pdf_We(const Ray &ray, Float *pdfPos, Float *pdfDir) const;
    Spectrum Sample_Wi(const Interaction &ref, const Point2f &sample,
                       Vector3f *wi, Float *pdf, Point2f *pRaster,
                       VisibilityTester *vis) const;

  private:
    // PerspectiveCamera Private Data
    Vector3f dxCamera, dyCamera;
    Float A; // Area
};

PerspectiveCamera *CreatePerspectiveCamera(const ParamSet &params,
                                           const AnimatedTransform &cam2world,
                                           Film *film, const Medium *medium);

}  // namespace pbrt

#endif  // PBRT_CAMERAS_PERSPECTIVE_H
