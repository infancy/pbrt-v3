
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

#ifndef PBRT_CORE_CAMERA_H
#define PBRT_CORE_CAMERA_H

// core/camera.h*
#include "pbrt.h"
#include "geometry.h"
#include "transform.h"
#include "film.h"

namespace pbrt {

// Camera Declarations
class Camera {
  public:
    // Camera Interface
    Camera(const AnimatedTransform &CameraToWorld, Float shutterOpen,
           Float shutterClose, Film *film, const Medium *medium);
    virtual ~Camera();

    // 返回世界空间中的 ray, ray 是经过归一化的, 整个系统依赖于这个行为
    // 最后返回的浮点值代表在 film 上的贡献值, 用来模拟真实相机的渐晕效果
    virtual Float GenerateRay(const CameraSample &sample, Ray *ray) const = 0;
    // 在 GenerateRay 生成的主光线的基础上, 往胶片平面的x,y方向上各偏移一像素生成辅助光线
    // 生成的 RayDifferential 结构主要用于纹理反走样
    virtual Float GenerateRayDifferential(const CameraSample &sample,
                                          RayDifferential *rd) const;

    // used in Section16, bidirectional light transport algorithms
    virtual Spectrum We(const Ray &ray, Point2f *pRaster2 = nullptr) const;
    virtual void Pdf_We(const Ray &ray, Float *pdfPos, Float *pdfDir) const;
    virtual Spectrum Sample_Wi(const Interaction &ref, const Point2f &u,
                               Vector3f *wi, Float *pdf, Point2f *pRaster,
                               VisibilityTester *vis) const;

    // Camera Public Data
    // P356
    AnimatedTransform CameraToWorld; // 可以把相机变换到世界空间中
    const Float shutterOpen, shutterClose; // 快门时间, 用于模拟运动模糊
    Film *film;
    const Medium *medium;
};

// Camera 生成光线时用到的所有参数
struct CameraSample {
    Point2f pFilm; // 在 film 上的采样点
    Point2f pLens; // lens 上的随机采样点
    Float time; // 采样这条光线的时间
};

inline std::ostream &operator<<(std::ostream &os, const CameraSample &cs) {
    os << "[ pFilm: " << cs.pFilm << " , pLens: " << cs.pLens <<
        StringPrintf(", time %f ]", cs.time);
    return os;
}

// 投影相机, 按投影方式区分可有正交投影, 透视投影
class ProjectiveCamera : public Camera {
  public:
    // ProjectiveCamera Public Methods
    ProjectiveCamera(const AnimatedTransform &CameraToWorld,
                     const Transform &CameraToScreen, // 派生类根据实现方式, 将坐标从相机空间变换到屏幕空间
                     const Bounds2f &screenWindow, // 并指定屏幕空间中要显示的区域
                     Float shutterOpen, Float shutterClose, 
                     Float lensr, Float focald, 
                     Film *film, const Medium *medium)
        : Camera(CameraToWorld, shutterOpen, shutterClose, film, medium),
          CameraToScreen(CameraToScreen) 
    {
        // Initialize depth of field parameters
        // 用于模拟景深效果
        lensRadius = lensr;
        focalDistance = focald;

        // Compute projective camera transformations

        // Compute projective camera screen transformations

/*
         screen space -> NDC space -> raster space

                                     Y
                                     ^
                                     |
                                     |
                                     |                
         (pMin.x, pMax.y)            |
                        W-------------------------screenWindow.pMax
                        |            |     |      |
                        |            |     |      |
                        |            |     |      |
        ----------------|------------O-----|------|----------------> X
                        |            |     |      |
                        |------------------S      |
                        |            |            |    O(0, 0)
        screenWindow.pMin--------------------------    S(x, y)
                                     |
                                     |
                                     |
                                     |                                     
               
        (pMin.x, -pMax.y)             
                        W--------------------------
                        |            |     |      |
                        |            |     |      |
                        |            |     |      |
                        |------------O-----|------|----------------> X
                        |            |     |      |
                        |------------------S      |
                        |            |            |    O(0,  0)
                        ---------------------------    S(x, -y)
                                     |
                                     |
                                     |
                                     |
                                     V
                                     Y

                   (0, 0)
                        W------------------------------------> X
                        |            |     |      |
                        |            |     |      |
                        |------------O     |      |
                        |                  |      | 
                        |------------------S      |
                        |                         |    O(-pMin.x, pMax.y)        by (0,  0) - (pMin.x, -pMax.y) 
                        |--------------------------
                        |                              S(x - pMin.x, pMax.y - y) by (x, -y) - (pMin.x, -pMax.y) 
                        |
                        |
                        |
                        V
                        Y
                              
                虽然屏幕空间是个三维空间, 但不妨把它当作二维平面
                而 screenWindow 即屏幕空间中的一个窗口, 在这个窗口中的对象会被渲染到胶片上

                1.首先将屏幕空间中的坐标变换到 screenWindow 中
                    1).翻转 y 轴, (x_s, y_s) -> (x_s, -y_s), (pMin.x, pMax.y) -> (pMin.x, -pMax.y)
                    2).将 screenWindow 的左上角作为原点, 计算相对于 screenWindow 的坐标
                2.归一化, 得到 NDC 空间中的坐标 [0, 1]^2
                3.再将 NDC 的坐标按 film 的分辨率进行放缩, 得到光栅化空间中的坐标

                e.g.
                设屏幕空间中 screenWindow 的 (pMin.x, pMax.y) 为 (-1.78, 1), 并将 (0, 0) 转换到光栅化空间中:
                (-1.78, 1) * (1, -1) => (-1.78, -1)
                (0s, 0s)   * (1, -1) => (0s, 0s)

                (0s, 0s) - (-1.78, -1) => (1.78, 1)
                (1.78, 1) / (3.56, 2) => (0.5_nd, 0.5_nd)

                (0.5_nd, 0.5_nd) * (1920, 1080) => (960_r, 540_r)

                P.S. screen -> NDC 的操作, 在 OpenGL/D3D 中被拆分到透视变换和视口变换中去了
        */

        // 这个变换要从下往上理解
        ScreenToRaster =
            Scale(film->fullResolution.x, film->fullResolution.y, 1) *
            Scale(
                1 / (screenWindow.pMax.x - screenWindow.pMin.x), 
                1 / (screenWindow.pMin.y - screenWindow.pMax.y),  1) *
            Translate(Vector3f(-screenWindow.pMin.x, -screenWindow.pMax.y, 0));

        // 从 screen 到 raster, z 坐标不会发生变化, 所以反过来, raster.z = 0, 则 screen.z = 0 => camera.z = near, 也就是光栅化平面的坐标会变换到相机空间中的近平面上
        // 正交相机的近平面在 z = 0 上, 透视相机默认的近平面是 z = 0.01f
        RasterToScreen = Inverse(ScreenToRaster);
        RasterToCamera = Inverse(CameraToScreen) * RasterToScreen;
    }

  protected:
    // ProjectiveCamera Protected Data
    Transform CameraToScreen, RasterToCamera;
    Transform ScreenToRaster, RasterToScreen;
    Float lensRadius, focalDistance;
};

}  // namespace pbrt

#endif  // PBRT_CORE_CAMERA_H
