
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

#ifndef PBRT_CORE_SAMPLER_H
#define PBRT_CORE_SAMPLER_H

// core/sampler.h*
#include "pbrt.h"
#include "geometry.h"
#include "rng.h"
#include <inttypes.h>

namespace pbrt {

// https://zhuanlan.zhihu.com/p/73943687

// The task of a Sampleris to **generate a sequence** of n-dimensional samples in $[0, 1)^n$ 
// 单个图像采样点需要的维度 n 由具体的光线传输算法来决定， 如 (x, y, t, u, v, [(u0, u1), (u2, u3), (u4, u5), (u6, u7)])

// Sampler Declarations
class Sampler {
  public:
    // Sampler Interface
    virtual ~Sampler();
    // 每个像素点的采样数量
    Sampler(int64_t samplesPerPixel);

    /*
        for (Point2i pixel : tileBounds)
        {
            sampler->StartPixel(pixel);
            // 隐含了 sampler->StartThisSample(), 会对该采样点生成新的采样序列

            do 
            {
                //Point2f p = a(sampler->Get2D()); // get x, y
                //Float t = c(sampler->Get1D()); // get t
                //...
                auto s = sampler->GetCameraSample(); // get x, y, t, u, v
            }
            while (sampler->StartNextSample());
        }
    */

    // 开始渲染这个像素前, 先传入其坐标, 有的采样器根据这个信息可以生成更好的采样点
    virtual void StartPixel(const Point2i &p);
    // 通知采样器对当前像素的下一个采样点进行采样
	virtual bool StartNextSample();

    // 一些积分器(如 SPPM)并不按顺序使用样本, 所以需要提供指定使用某个样本的方法
    // virtual bool SetCurrentPixelSampleIndex(int64_t sampleIndex); 
	virtual bool SetSampleNumber(int64_t sampleNum);
    // int64_t GetCurrentSampleNumberIndex() const { return currentPixelSampleIndex_; }
	int64_t CurrentSampleNumber() const { return currentPixelSampleIndex; }


	// 使用 Get1D()、Get2D()、GetCameraSample() 获得独立的样本 
	virtual Float   Get1D() = 0;
	virtual Point2f Get2D() = 0; //相当于连续调用两次 Get1D, 但有的采样器可以借此生成更好的分布
	CameraSample GetCameraSample(const Point2i &pRaster);


    // 在渲染之前生成所有采样点, 而不是像 Get1D 那样即时生成
	// 这一对函数(Request/Get)和渲染的整个过程是紧密相关的，刚开始接触时可以跳过. 在阅读 DirectLightingIntegrator 的实现时, 可以参考下方的注释来理解
    void Request1DArray(int n);
    void Request2DArray(int n);
    const Float *Get1DArray(int n);
    // 形如 (x, y, t, u, v, [(u0, u1), (u2, u3), (u4, u5), (u6, u7)]) 的采样序列, Get2DArray 返回的是中括号部分的数据
	const Point2f *Get2DArray(int n);

    // 传入预计的每像素采样数量, 返回最佳的采样数量
    // 例如 ZeroTwoSequenceSampler 使用 2^n 可以获得最好的质量
    virtual int RoundCount(int n) const { return n; }

    // 积分器会使用多线程进行计算, 使用 Clone 给每个线程生成一个单独的采样器
    virtual std::unique_ptr<Sampler> Clone(int seed) = 0;

    std::string StateString() const {
      // pixel coord: ({0}, {1}), current pixel sample index: {2}
      return StringPrintf("(%d,%d), sample %" PRId64, currentPixel.x,
                          currentPixel.y, currentPixelSampleIndex);
    }

    // Sampler Public Data
    const int64_t samplesPerPixel;

  protected:
    // Sampler Protected Data
    Point2i currentPixel;
    int64_t currentPixelSampleIndex;

	/*
        sampleArray1D/2D 的生成, 使用和整个渲染流程耦合在一起, 是 Sampler 的接口里比较难理解的地方

        void DirectLightingIntegrator::Preprocess(...)
        {
            //...
            for (int i = 0; i < maxDepth; ++i)
            {
                for (size_t j = 0; j < scene.lights.size(); ++j)
                {
                    sampler.Request2DArray(nLightSamples[j]);
                    sampler.Request2DArray(nLightSamples[j]);
                
                    void Sampler::Request2DArray(int nLightSamples)
                    {
                        //...
                        samples2DArraySizes.push_back(nPixel);
                        sampleArray2D.push_back(std::vector<Point2f>(nLightSamples * samplesPerPixel));
                    }
                }
            }
        }

	    每个子数组的大小为 n * samplesPerPixel，samples1/2DArraySizes 记录的是这个 n 的大小, 其主要作用是在 Get1DArray 时对传入参数进行校验
    */
    std::vector<int> samples1DArraySizes, samples2DArraySizes;	
    std::vector<std::vector<Float>> sampleArray1D;
    std::vector<std::vector<Point2f>> sampleArray2D;

  private:
    /*
        void SamplerIntegrator::Render(const Scene &scene)
        {
            // ...
            for (Point2i pixel : tileBounds)
            {
                // ...
                tileSampler->StartPixel(pixel);	// currentPixelSampleIndex = 0; array1DOffset = array2DOffset = 0;

                do
                {
                    Spectrum L = Li(ray, scene, *tileSampler, arena);

                    Spectrum DirectLightingIntegrator::Li(...)
                    {
                        // ...
                        L += UniformSampleAllLights(isect, scene, arena, sampler, nLightSamples);

                        // 递归追踪镜面项, 会用到 sampleArray2D[nLightSize * depth ... nLightSize * (depth + n)] 中的数据
                    }

                    Spectrum UniformSampleAllLights(...)
                    {
                        //TODO: 麻烦的地方在于每次从 sampleArray2D 取出的是一列数据, 而不是一行
                        for (size_t j = 0; j < scene.lights.size(); ++j)
                        {
                            const Point2f *Sampler::Get2DArray(int n)
                            {
                                 //...
                                 return &sampleArray2D[array2DOffset++][currentPixelSampleIndex * n];
                            }

                            int nLightSamples = nLightSamples[j];
                            const Point2f *uLightArray = sampler.Get2DArray(nLightSamples); // 获得 uLightArray[nLightSamples]
                            const Point2f *uScatteringArray = sampler.Get2DArray(nLightSamples);

                            if (!uLightArray || !uScatteringArray)
                            {
                                // ...
                            }
                            else
                            {
                                // ...

                                for (int k = 0; k < nLightSamples; ++k)
                                {
                                    Ld += EstimateDirect(it, uScatteringArray[k], *light,
                                                         uLightArray[k], scene, sampler, arena, handleMedia);
                                }
                    
                                // ...
                            }
                        }
                    }

                    filmTile->AddSample(cameraSample.pFilm, L, rayWeight);
                }
                while (tileSampler->StartNextSample()); // ++currentPixelSampleIndex; array1DOffset = array2DOffset = 0;
            }
        }

	    array1D/2DOffset 记录当前样本在 sampleArray1/2D 中的索引
    */
    size_t array1DOffset, array2DOffset;
};

// 针对单个像素点生成样本的采样器
class PixelSampler : public Sampler {
  public:
    // PixelSampler Public Methods
    // 每像素的采样数量, 每个采样点需要采样的维度
    PixelSampler(int64_t samplesPerPixel, int nSampledDimensions);

    // 像素采样器的子类会在 StartPixel 中生成 samples1D 和 samples2D 中的采样点(此外还有 sampleArray1D 和 sampleArray2D)

    bool StartNextSample();
    bool SetSampleNumber(int64_t);
    Float Get1D();
    Point2f Get2D();

  protected:
    // PixelSampler Protected Data
    // 根据 nSampledDimensions 提前生成所有供 Get1D/2D 使用的样本
    std::vector<std::vector<Float>> samples1D;
    std::vector<std::vector<Point2f>> samples2D;
    int current1DDimension = 0, current2DDimension = 0;

    // 1.在构造时生成初始的种子
    // 2.当申请超出 nSampledDimensions 范围的样本时, 返回使用 rng 生成的 [0, 1) 内均匀分布的样本
    RNG rng;
};

// 针对整个图像平面生成样本的采样器(如 Halton Sampler)
class GlobalSampler : public Sampler {
  public:
    // GlobalSampler Public Methods
    bool StartNextSample();
    void StartPixel(const Point2i &);
    bool SetSampleNumber(int64_t sampleNum);
    Float Get1D();
    Point2f Get2D();

    GlobalSampler(int64_t samplesPerPixel) : Sampler(samplesPerPixel) {}

    // GlobalSampler 的派生类必须实现这两个接口, GlobalSampler 从这两个接口得到采样数据
    // 参考 P429/Table7.2 理解

    // 获得 currentPixel 的第 sampleNum 个采样的全局索引(Table7.2 的第一列)
    virtual int64_t GetIndexForSample(int64_t sampleNum) const = 0;

    // 输入第一列的全局索引, 返回第三列像素点上的采样值
    // 例如采样坐标 (1.250000, 2.333333) 的采样值是 (0.250000, 0.333333), 根据第二个参数 dimension 来返回 0.250000 还是 0.333333
    virtual Float SampleDimension(int64_t index, int dimension) const = 0;

  private:
    // GlobalSampler Private Data
    int dimension;
    int64_t intervalSampleIndex; // 全局索引

    // PBRT 假设前几个维度生成采样点的质量更好, 所以把他们保留给 (x, y, t, u, v), 接下来的采样点分配给 sampleArray1D/2D, 其余的给 Get1D/2D 使用
    static const int arrayStartDim = 5; // (x, y, t, u, v, [(u0, u1), (u2, u3), (u4, u5), (u6, u7)])
    int arrayEndDim;
};

}  // namespace pbrt

#endif  // PBRT_CORE_SAMPLER_H
