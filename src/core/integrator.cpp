
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

// core/integrator.cpp*
#include "integrator.h"
#include "scene.h"
#include "interaction.h"
#include "sampling.h"
#include "parallel.h"
#include "film.h"
#include "sampler.h"
#include "integrator.h"
#include "progressreporter.h"
#include "camera.h"
#include "stats.h"

namespace pbrt {

STAT_COUNTER("Integrator/Camera rays traced", nCameraRays);

// Integrator Method Definitions
Integrator::~Integrator() {}

// Integrator Utility Functions
Spectrum UniformSampleAllLights(const Interaction &it, const Scene &scene,
                                MemoryArena &arena, Sampler &sampler,
                                const std::vector<int> &nLightSamples,
                                bool handleMedia) 
{
    ProfilePhase p(Prof::DirectLighting);
    Spectrum L(0.f);
    
	// 计算所有光源经由交点 it 向 it.wo 方向发射的辐射度
    for (size_t j = 0; j < scene.lights.size(); ++j) 
    {
        // Accumulate contribution of _j_th light to _L_
        const std::shared_ptr<Light> &light = scene.lights[j];
        int nSamples = nLightSamples[j];
		
        // 使用 sampler.Get1D/2DArray 前, 要保证已调用 sampler.Request1D/2DArray 生成了采样点
        const Point2f *uLightArray = sampler.Get2DArray(nSamples);
        const Point2f *uScatteringArray = sampler.Get2DArray(nSamples);

		// 如果样本数组被用光时，只进行单次采样
        if (!uLightArray || !uScatteringArray) 
        {
            // Use a single sample for illumination from _light_
            Point2f uLight = sampler.Get2D();
            Point2f uScattering = sampler.Get2D();
            L += EstimateDirect(it, uScattering, *light, uLight, scene, sampler,
                                arena, handleMedia);
        } 
		// 否则根据 nSamples 的大小，在光源上取 n 个采样点，计算平均辐射度
        else 
        {
            // Estimate direct lighting using sample arrays
            Spectrum Ld(0.f);
            for (int k = 0; k < nSamples; ++k)
                Ld += EstimateDirect(it, uScatteringArray[k], *light,
                                     uLightArray[k], scene, sampler, arena,
                                     handleMedia);
            L += Ld / nSamples;
        }
    }
    return L;
}

Spectrum UniformSampleOneLight(const Interaction &it, const Scene &scene,
                               MemoryArena &arena, Sampler &sampler,
                               bool handleMedia, const Distribution1D *lightDistrib) 
{
    ProfilePhase p(Prof::DirectLighting);

    // Randomly choose a single light to sample, _light_
	// 从光源中随机选取一个进行采样
    int nLights = int(scene.lights.size());
    if (nLights == 0) 
        return Spectrum(0.f);

    int lightNum;
    Float lightPdf;
	// 如果光源的功率分布可用，则根据其功率分布（由逆变换算法）选取一个光源，并计算其概率密度
    if (lightDistrib) {
        lightNum = lightDistrib->SampleDiscrete(sampler.Get1D(), &lightPdf);
        if (lightPdf == 0) return Spectrum(0.f);
    } else {
        lightNum = std::min((int)(sampler.Get1D() * nLights), nLights - 1);
        lightPdf = Float(1) / nLights;
    }

    const std::shared_ptr<Light> &light = scene.lights[lightNum];
    Point2f uLight = sampler.Get2D();
    Point2f uScattering = sampler.Get2D();
	// 采样单个光源，并除以其功率的概率密度（power_pdf），得到近似采样所有光源的结果
    return EstimateDirect(it, uScattering, *light, uLight,
                          scene, sampler, arena, handleMedia) / lightPdf;
}

// 计算单个 light 经由 isect 向 wo 方向发射的辐射度
//
// prev   light
// -----  -----
//   ^      ^
//    \    /
//  wo \  / wi
//      \/
//    ------
//    isect
//
//	当 bsdf 呈镜面状态而光源分布较广时从 bsdf 采样较为高效
//	当 bsdf 呈漫反射分布状态而光源较小时从光源采样更高效
//	因而使用多重重要性采样（MIS）分别对 light 和 bsdf 进行采样
//	在光源上采样一点 p 计算 Le，在 bsdf 上采样一方向 wi 计算 Li，最后 Ld = MIS(Le, Li)

Spectrum EstimateDirect(const Interaction &it, const Point2f &uScattering,
                        const Light &light, const Point2f &uLight,
                        const Scene &scene, Sampler &sampler,
                        MemoryArena &arena, bool handleMedia, bool specular) 
{
    BxDFType bsdfFlags =
        specular ? BSDF_ALL : BxDFType(BSDF_ALL & ~BSDF_SPECULAR);
    Spectrum Ld(0.f);

	// Sample light source with multiple importance sampling
	// 从光源部分进行多重重要性采样
    Vector3f wi;
    Float lightPdf = 0, scatteringPdf = 0;
    VisibilityTester visibility;
	// 传入交点，在光源上选取一点，计算选到该点的概率密度，交点到该点的方向 wi、入射辐射度 Li 及可见性
    Spectrum Li = light.Sample_Li(it, uLight, &wi, &lightPdf, &visibility);
    VLOG(2) << "EstimateDirect uLight:" << uLight << " -> Li: " << Li << ", wi: "
            << wi << ", pdf: " << lightPdf;
    if (lightPdf > 0 && !Li.IsBlack())
    {
        // Compute BSDF or phase function's value for light sample
		// 计算 BSDF（如果交点在 surface 中）或相位函数（如果交点在 medium 中）的值
        Spectrum f;
        if (it.IsSurfaceInteraction()) {
            // Evaluate BSDF for light sampling strategy
            const SurfaceInteraction &isect = (const SurfaceInteraction &)it;
			// 计算 f(p, wo, wi) * cos(eta_light_isect) 项
            f = isect.bsdf->f(isect.wo, wi, bsdfFlags) *
                AbsDot(wi, isect.shading.n);
            scatteringPdf = isect.bsdf->Pdf(isect.wo, wi, bsdfFlags);
            VLOG(2) << "  surf f*dot :" << f << ", scatteringPdf: " << scatteringPdf;
        } else {
            // Evaluate phase function for light sampling strategy
            const MediumInteraction &mi = (const MediumInteraction &)it;
            Float p = mi.phase->p(mi.wo, wi);
            f = Spectrum(p);
            scatteringPdf = p;
            VLOG(2) << "  medium p: " << p;
        }
        if (!f.IsBlack()) {
            // Compute effect of visibility for light source sample
			// visibility 负责处理光源上的采样点和交点的可见性问题，当两点不可见时 Li = 0
			// 当需要考虑 medium 时，则调用 Tr(scene, sampler) 进一步计算 Li 的值
            if (handleMedia) {	
                Li *= visibility.Tr(scene, sampler);
                VLOG(2) << "  after Tr, Li: " << Li;
            } else {			
              if (!visibility.Unoccluded(scene)) {
                VLOG(2) << "  shadow ray blocked";
                Li = Spectrum(0.f);
              } else
                VLOG(2) << "  shadow ray unoccluded";
            }

            // Add light's contribution to reflected radiance
            if (!Li.IsBlack()) {
                if (IsDeltaLight(light.flags))	// 如果是 delta 光源，无需使用 MIS
                    Ld += f * Li / lightPdf;	// return f * Li / lightPdf;
                else {
                    Float weight =
                        PowerHeuristic(1, lightPdf, 1, scatteringPdf);	// 使用功率启发式方法计算权重
                    Ld += f * Li * weight / lightPdf;
                }
            }
        }
    }

    // Sample BSDF with multiple importance sampling
	// 从 BSDF 部分进行多重重要性采样
    if (!IsDeltaLight(light.flags)) 
    {		   
        Spectrum f;
        bool sampledSpecular = false;
        if (it.IsSurfaceInteraction()) {
            // Sample scattered direction for surface interactions
            BxDFType sampledType;
            const SurfaceInteraction &isect = (const SurfaceInteraction &)it;
			// 在交点上半球方向随机选取一点，计算相应的 wi、scatteringPdf 和 sampledType
			// 并计算 f(p, wo, wi) 项
            f = isect.bsdf->Sample_f(isect.wo, &wi, uScattering, &scatteringPdf,
                                     bsdfFlags, &sampledType);
			// 计算 cos(eta_light_isect) 项
            f *= AbsDot(wi, isect.shading.n);
            sampledSpecular = (sampledType & BSDF_SPECULAR) != 0;
        } else {
            // Sample scattered direction for medium interactions
            const MediumInteraction &mi = (const MediumInteraction &)it;
            Float p = mi.phase->Sample_p(mi.wo, &wi, uScattering);
            f = Spectrum(p);
            scatteringPdf = p;
        }
        VLOG(2) << "  BSDF / phase sampling f: " << f << ", scatteringPdf: " <<
            scatteringPdf;
        if (!f.IsBlack() && scatteringPdf > 0) {
            // Account for light contributions along sampled direction _wi_
			// 计算沿采样方向 wi 的光照贡献

            Float weight = 1;
            if (!sampledSpecular) {
                lightPdf = light.Pdf_Li(it, wi);	// 计算沿 wi 方向采样到 light 的概率
                if (lightPdf == 0) return Ld;
                weight = PowerHeuristic(1, scatteringPdf, 1, lightPdf);
            }

            // Find intersection and compute transmittance
			// 计算沿 wi 方向是否与（面积）光源相交
            SurfaceInteraction lightIsect;
            Ray ray = it.SpawnRay(wi);
            Spectrum Tr(1.f);
            bool foundSurfaceInteraction =
                handleMedia ? scene.IntersectTr(ray, sampler, &lightIsect, &Tr)
                            : scene.Intersect(ray, &lightIsect);

            // Add light contribution from material sampling
			// 如果 wi 与面积光源相交，则需要进一步处理，否则直接计算 light 向 ray 方向发射的辐射度
            Spectrum Li(0.f);
            if (foundSurfaceInteraction) {
                if (lightIsect.primitive->GetAreaLight() == &light)
                    Li = lightIsect.Le(-wi);
            } else
                Li = light.Le(ray);
            if (!Li.IsBlack()) Ld += f * Li * Tr * weight / scatteringPdf;
        }
    }

    return Ld;
}

std::unique_ptr<Distribution1D> ComputeLightPowerDistribution(
    const Scene &scene) 
{
    if (scene.lights.empty()) 
        return nullptr;

    std::vector<Float> lightPower;
    for (const auto &light : scene.lights)
        lightPower.push_back(light->Power().y());

    return std::unique_ptr<Distribution1D>(
        new Distribution1D(&lightPower[0], lightPower.size()));
}

// SamplerIntegrator Method Definitions
void SamplerIntegrator::Render(const Scene &scene) 
{
    Preprocess(scene, *sampler);

    // Render image tiles in parallel
	// (0,0) _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ 
	//		|		|		|		|		|
	//		| tile	| tile	| tile	| tile	|
	//		|_ _ _ _|_ _ _ _|_ _ _ _|_ _ _ _|		
	//		|		|		|		|		|
	//		| tile	| tile	| tile	| tile	|
	//		|_ _ _ _|_ _ _ _|_ _ _ _|_ _ _ _|
	//		|		|		|		|		|
	//		| tile	| tile	| tile	| tile	|
	//		|_ _ _ _|_ _ _ _|_ _ _ _|_ _ _ _|
	//		|		|		|		|		|
	//		| tile	| tile	| tile	| tile	|
	//		|_ _ _ _|_ _ _ _|_ _ _ _|_ _ _ _|
	//										 (filmResolution.width, filmResolution.height)
	// 将图像平面分成一块块 tile，每个线程每次渲染一块 tile

    // Compute number of tiles, _nTiles_, to use for parallel rendering
	// 计算 tile 的数量（为了简单起见，pbrt 总是使用 16*16 像素大小的 tile）
    const int tileSize = 16;

	Bounds2i sampleBounds = camera->film->GetSampleBounds(); // 考虑默认的三角过滤器时, floatBounds 为 [-2, -2] 到 [1922, 1082]	
    Vector2i sampleExtent = sampleBounds.Diagonal();	// 相当于计算生成图像的分辨率
	// 计算 tile_numbers 时向上取整
    Point2i nTiles((sampleExtent.x + tileSize - 1) / tileSize,
                   (sampleExtent.y + tileSize - 1) / tileSize);

	// 提供一个关于 pbrt 当前进度的直观反馈
	ProgressReporter reporter(nTiles.x * nTiles.y, "Rendering");
    {
		// 每个 tile 由单独的线程执行，每次对 lambda 表达式传入一个该 tile 在 nTiles 中的位置
        ParallelFor2D([&](Point2i tile) {
            // Render section of image corresponding to _tile_

            // Allocate _MemoryArena_ for tile
			// 每个线程使用单独的内存池
            MemoryArena arena;

            // Get sampler instance for tile
			// 每个线程使用单独的采样器
            int seed = tile.y * nTiles.x + tile.x;
			std::unique_ptr<Sampler> tileSampler = sampler->Clone(seed);

            // Compute sample bounds for tile
			// 计算这个 tile 在 image 中覆盖到的像素范围

			// 假设 image 的分辨率为 1920 * 1080，则 nTiles 的大小为 120 * 80（即有 119 * 79 块 tile 等待渲染，不包含右下方的边界）
			// 当某个线程在某次调用 lambda 表达式并接收到一个位置为 (119, 79) 的 tile 时
			// 这个 tile 在 image 中覆盖到的像素范围为（则假设 pMin 为 0）：
			// x0 = 119 * 16 = 1904，y0 = 79 * 16 = 1064，x1 = 1920， y1 = 1080
			// 即这个线程需要渲染从图像坐标 (1904, 1064) 开始，(1920, 1080) 结束的这个矩形区域内的像素（也不包含右边和下边的边界）
            // （如果限制渲染的图像长宽总为 16（或 tileSize）的倍数，可以把这段代码简化不少）
			int x0 = sampleBounds.pMin.x + tile.x * tileSize;
            int x1 = std::min(x0 + tileSize, sampleBounds.pMax.x);	
            int y0 = sampleBounds.pMin.y + tile.y * tileSize;
            int y1 = std::min(y0 + tileSize, sampleBounds.pMax.y);
            Bounds2i tileBounds(Point2i(x0, y0), Point2i(x1, y1));
            LOG(INFO) << "Starting image tile " << tileBounds;

            // Get _FilmTile_ for tile
			// 先将渲染得到的图像存入这个 filmTile 中，待渲染结束后将 filmTile 合并到 film 中
            std::unique_ptr<FilmTile> filmTile =
                camera->film->GetFilmTile(tileBounds);

            // Loop over pixels in tile to render them
            for (Point2i pixel : tileBounds) 
            {	// 遍历该二维包围盒上的每一个像素（Bounds2i 定义了相应的迭代器和 begin()、end() 函数）
                {
                    ProfilePhase pp(Prof::StartPixel);
					// 采样一个新的像素前，先对采样器进行一些设置
                    tileSampler->StartPixel(pixel);		
                }
                
				// 检查该像素是否在 pixelBounds 内（为了照顾过滤器, pixelBounds 可能超出了图像平面）
				// 这可以保持（大多数）采样器中所使用的 RNG 值的一致性，方便重现？？？、调试
                // Do this check after the StartPixel() call; this keeps
                // the usage of RNG values from (most) Samplers that use
                // RNGs consistent, which improves reproducability /
                // debugging.
                if (!InsideExclusive(pixel, pixelBounds))
                    continue;

                do 
                {
                    // Initialize _CameraSample_ for current sample
                    CameraSample cameraSample =
                        tileSampler->GetCameraSample(pixel);

                    // Generate camera ray for current sample
					// 为当前样本生成相机光线
                    RayDifferential ray;
                    Float rayWeight = camera->GenerateRayDifferential(cameraSample, &ray);
                    ray.ScaleDifferentials(1 / std::sqrt((Float)tileSampler->samplesPerPixel));

                    ++nCameraRays;

                    // Evaluate radiance along camera ray
					//计算沿这条光线（-ray.direction)的辐射度
					Spectrum L(0.f);
                    if (rayWeight > 0) L = Li(ray, scene, *tileSampler, arena);

                    // Issue warning if unexpected radiance value returned
					//对得到的辐射度做检查
                    {
					    if (L.HasNaNs()) {
                            LOG(ERROR) << StringPrintf(
                                "Not-a-number radiance value returned "
                                "for pixel (%d, %d), sample %d. Setting to black.",
                                pixel.x, pixel.y,
                                (int)tileSampler->CurrentSampleNumber());
                            L = Spectrum(0.f);
                        } else if (L.y() < -1e-5) {
                            LOG(ERROR) << StringPrintf(
                                "Negative luminance value, %f, returned "
                                "for pixel (%d, %d), sample %d. Setting to black.",
                                L.y(), pixel.x, pixel.y,
                                (int)tileSampler->CurrentSampleNumber());
                            L = Spectrum(0.f);
                        } else if (std::isinf(L.y())) {
                              LOG(ERROR) << StringPrintf(
                                "Infinite luminance value returned "
                                "for pixel (%d, %d), sample %d. Setting to black.",
                                pixel.x, pixel.y,
                                (int)tileSampler->CurrentSampleNumber());
                            L = Spectrum(0.f);
                        }
                        VLOG(1) << "Camera sample: " << cameraSample << " -> ray: " <<
                            ray << " -> L = " << L;
                    }

                    // Add camera ray's contribution to image
                    filmTile->AddSample(cameraSample.pFilm, L, rayWeight);

                    // Free _MemoryArena_ memory from computing image sample
                    // value
                    arena.Reset();
                } 
                while (tileSampler->StartNextSample());
            }
            LOG(INFO) << "Finished image tile " << tileBounds;

            // Merge image tile into _Film_
			// 将 filmTile 合并到 film 中
            camera->film->MergeFilmTile(std::move(filmTile));
            reporter.Update();
        }, nTiles);

        reporter.Done();
    }
    LOG(INFO) << "Rendering finished";

    // Save final image after rendering
	// 渲染结束，保存图片
    camera->film->WriteImage();
}

Spectrum SamplerIntegrator::SpecularReflect(
    const RayDifferential &ray, const SurfaceInteraction &isect,
    const Scene &scene, Sampler &sampler, MemoryArena &arena, int depth) const 
{
    // Compute specular reflection direction _wi_ and BSDF value
    Vector3f wo = isect.wo, wi;
    Float pdf;
    BxDFType type = BxDFType(BSDF_REFLECTION | BSDF_SPECULAR);
	// 在 isect 上半球面选取一方向 wi，并计算相应的概率密度及 f(p, wo, wi)
    Spectrum f = isect.bsdf->Sample_f(wo, &wi, sampler.Get2D(), &pdf, type);

    // Return contribution of specular reflection
	// 计算镜面反射的贡献
    const Normal3f &ns = isect.shading.n;
    if (pdf > 0.f && !f.IsBlack() && AbsDot(wi, ns) != 0.f) 
    {
        // Compute ray differential _rd_ for specular reflection
        RayDifferential rd = isect.SpawnRay(wi);
		// 关于光线微分的代码没细看，暂时不做解释
        if (ray.hasDifferentials) 
        {
            rd.hasDifferentials = true;
            rd.rxOrigin = isect.p + isect.dpdx;
            rd.ryOrigin = isect.p + isect.dpdy;

            // Compute differential reflected directions
            Normal3f dndx = isect.shading.dndu * isect.dudx +
                            isect.shading.dndv * isect.dvdx;
            Normal3f dndy = isect.shading.dndu * isect.dudy +
                            isect.shading.dndv * isect.dvdy;
            Vector3f dwodx = -ray.rxDirection - wo,
                     dwody = -ray.ryDirection - wo;
            Float dDNdx = Dot(dwodx, ns) + Dot(wo, dndx);
            Float dDNdy = Dot(dwody, ns) + Dot(wo, dndy);

            rd.rxDirection =
                wi - dwodx + 2.f * Vector3f(Dot(wo, ns) * dndx + dDNdx * ns);
            rd.ryDirection =
                wi - dwody + 2.f * Vector3f(Dot(wo, ns) * dndy + dDNdy * ns);
        }

        return f * Li(rd, scene, sampler, arena, depth + 1) * AbsDot(wi, ns) /
               pdf;
    } 
    else
        return Spectrum(0.f);
}

Spectrum SamplerIntegrator::SpecularTransmit(
    const RayDifferential &ray, const SurfaceInteraction &isect,
    const Scene &scene, Sampler &sampler, MemoryArena &arena, int depth) const {
    Vector3f wo = isect.wo, wi;
    Float pdf;
    const Point3f &p = isect.p;
    const BSDF &bsdf = *isect.bsdf;
	// 在 isect 下半球面选取一方向 wi，并计算相应的概率密度及 f(p, wo, wi)
    Spectrum f = bsdf.Sample_f(wo, &wi, sampler.Get2D(), &pdf,
                               BxDFType(BSDF_TRANSMISSION | BSDF_SPECULAR));
    Spectrum L = Spectrum(0.f);
    Normal3f ns = isect.shading.n;
    if (pdf > 0.f && !f.IsBlack() && AbsDot(wi, ns) != 0.f) {
        // Compute ray differential _rd_ for specular transmission
        RayDifferential rd = isect.SpawnRay(wi);
        if (ray.hasDifferentials) {
            rd.hasDifferentials = true;
            rd.rxOrigin = p + isect.dpdx;
            rd.ryOrigin = p + isect.dpdy;

            Normal3f dndx = isect.shading.dndu * isect.dudx +
                            isect.shading.dndv * isect.dvdx;
            Normal3f dndy = isect.shading.dndu * isect.dudy +
                            isect.shading.dndv * isect.dvdy;

            // The BSDF stores the IOR of the interior of the object being
            // intersected.  Compute the relative IOR by first out by
            // assuming that the ray is entering the object.
            Float eta = 1 / bsdf.eta;
            if (Dot(wo, ns) < 0) {
                // If the ray isn't entering, then we need to invert the
                // relative IOR and negate the normal and its derivatives.
                eta = 1 / eta;
                ns = -ns;
                dndx = -dndx;
                dndy = -dndy;
            }

            /*
              Notes on the derivation:
              - pbrt computes the refracted ray as: \wi = -\eta \omega_o + [ \eta (\wo \cdot \N) - \cos \theta_t ] \N
                It flips the normal to lie in the same hemisphere as \wo, and then \eta is the relative IOR from
                \wo's medium to \wi's medium.
              - If we denote the term in brackets by \mu, then we have: \wi = -\eta \omega_o + \mu \N
              - Now let's take the partial derivative. (We'll use "d" for \partial in the following for brevity.)
                We get: -\eta d\omega_o / dx + \mu dN/dx + d\mu/dx N.
              - We have the values of all of these except for d\mu/dx (using bits from the derivation of specularly
                reflected ray deifferentials).
              - The first term of d\mu/dx is easy: \eta d(\wo \cdot N)/dx. We already have d(\wo \cdot N)/dx.
              - The second term takes a little more work. We have:
                 \cos \theta_i = \sqrt{1 - \eta^2 (1 - (\wo \cdot N)^2)}.
                 Starting from (\wo \cdot N)^2 and reading outward, we have \cos^2 \theta_o, then \sin^2 \theta_o,
                 then \sin^2 \theta_i (via Snell's law), then \cos^2 \theta_i and then \cos \theta_i.
              - Let's take the partial derivative of the sqrt expression. We get:
                1 / 2 * 1 / \cos \theta_i * d/dx (1 - \eta^2 (1 - (\wo \cdot N)^2)).
              - That partial derivatve is equal to:
                d/dx \eta^2 (\wo \cdot N)^2 = 2 \eta^2 (\wo \cdot N) d/dx (\wo \cdot N).
              - Plugging it in, we have d\mu/dx =
                \eta d(\wo \cdot N)/dx - (\eta^2 (\wo \cdot N) d/dx (\wo \cdot N))/(-\wi \cdot N).
             */
            Vector3f dwodx = -ray.rxDirection - wo,
                     dwody = -ray.ryDirection - wo;
            Float dDNdx = Dot(dwodx, ns) + Dot(wo, dndx);
            Float dDNdy = Dot(dwody, ns) + Dot(wo, dndy);

            Float mu = eta * Dot(wo, ns) - AbsDot(wi, ns);
            Float dmudx =
                (eta - (eta * eta * Dot(wo, ns)) / AbsDot(wi, ns)) * dDNdx;
            Float dmudy =
                (eta - (eta * eta * Dot(wo, ns)) / AbsDot(wi, ns)) * dDNdy;

            rd.rxDirection =
                wi - eta * dwodx + Vector3f(mu * dndx + dmudx * ns);
            rd.ryDirection =
                wi - eta * dwody + Vector3f(mu * dndy + dmudy * ns);
        }
        L = f * Li(rd, scene, sampler, arena, depth + 1) * AbsDot(wi, ns) / pdf;
    }
    return L;
}

}  // namespace pbrt
