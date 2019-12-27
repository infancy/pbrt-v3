
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

#ifndef PBRT_CORE_LOWDISCREPANCY_H
#define PBRT_CORE_LOWDISCREPANCY_H

// core/lowdiscrepancy.h*
#include "pbrt.h"
#include "rng.h"
#include "sampling.h"
#include "sobolmatrices.h"

namespace pbrt {

// Low Discrepancy Declarations
Float RadicalInverse(int baseIndex, uint64_t a);

static PBRT_CONSTEXPR int PrimeTableSize = 1000; // 重排序表的最大长度
extern const int Primes[PrimeTableSize];
extern const int PrimeSums[PrimeTableSize];
std::vector<uint16_t> ComputeRadicalInversePermutations(RNG &rng);
Float ScrambledRadicalInverse(int baseIndex, uint64_t a, const uint16_t *perm);

inline void Sobol2D(int nSamplesPerPixelSample, int nPixelSamples, Point2f *samples, RNG &rng);

extern uint32_t CMaxMinDist[17][32];

inline uint64_t SobolIntervalToIndex(const uint32_t log2Resolution,
                                     uint64_t sampleNum, const Point2i &p);
inline float SobolSampleFloat(  int64_t index, int dimension, uint32_t scramble = 0);
inline double SobolSampleDouble(int64_t index, int dimension, uint64_t scramble = 0);



#pragma region used by HaltonSampler

// Low Discrepancy Inline Functions
// 反转比特位
inline uint32_t ReverseBits32(uint32_t n) 
{
    n = (n << 16) | (n >> 16);
    n = ((n & 0x00ff00ff) << 8) | ((n & 0xff00ff00) >> 8); // 交换前面16位的两个8位, 交换后面16位的两个8位
    n = ((n & 0x0f0f0f0f) << 4) | ((n & 0xf0f0f0f0) >> 4); // 递归这个过程
    n = ((n & 0x33333333) << 2) | ((n & 0xcccccccc) >> 2);
    n = ((n & 0x55555555) << 1) | ((n & 0xaaaaaaaa) >> 1);
    return n;
}

inline uint64_t ReverseBits64(uint64_t n) {
    uint64_t n0 = ReverseBits32((uint32_t)n);
    uint64_t n1 = ReverseBits32((uint32_t)(n >> 32));
    return (n0 << 32) | n1;
}

// InverseRadicalInverse<2>(94, 7) => 94=0b1011110, 反转比特位, 0b0111101=61
template <int base>
inline uint64_t InverseRadicalInverse(uint64_t inverse, int nDigits) {
    uint64_t index = 0;
    for (int i = 0; i < nDigits; ++i) {
        uint64_t digit = inverse % base;
        inverse /= base;
        index = index * base + digit;
    }
    return index;
}

#pragma endregion used by HaltonSampler



#pragma region used by ZeroTwoSequenceSampler and MaxMinDistSampler

// P459 式 7.10 矩阵-向量乘法的优化实现
// 其中 C[i] 代表矩阵的一列, a 的一个比特位代表了向量的一个元素
// 因为矩阵 C 比较特殊, 每列只有一个比特位为1, 所以这里的乘法被优化成了异或操作
// MultiplyGenerator 只被 SampleGeneratorMatrix 调用
inline uint32_t MultiplyGenerator(const uint32_t *C, uint32_t a) 
{
    uint32_t v = 0;
    for (int i = 0; a != 0; ++i, a >>= 1)
        if (a & 1) v ^= C[i];
    return v;
}

inline Float SampleGeneratorMatrix(const uint32_t *C, uint32_t a,
                                   uint32_t scramble = 0) {
#ifndef PBRT_HAVE_HEX_FP_CONSTANTS
    return std::min((MultiplyGenerator(C, a) ^ scramble) * Float(2.3283064365386963e-10),
                    OneMinusEpsilon);
#else
    // 为了避免 P458 式 7.11 的运算, 这里又有一个优化, 提前翻转了矩阵 C 的每一列(注意式 7.10 的计算)
    // 所以这里只需要移位, 不需要再翻转
    // 同时还用 scramble 对结果进行置乱, 避免不同像素点生成一样的采样值.
    // 参考 P459, 这里的置乱也有优化
    return std::min((MultiplyGenerator(C, a) ^ scramble) * Float(0x1p-32),
                    OneMinusEpsilon);
#endif
}

// https://zhuanlan.zhihu.com/p/20374706

// GrayCodeSample 的参考代码, 没有实际用到
inline uint32_t GrayCode(uint32_t v) { return (v >> 1) ^ v; }

// 相比 MultiplyGenerator 要更高效的生成矩阵的方法, 参考 P460 的 Table7.4 和 P461
inline void GrayCodeSample(const uint32_t *C, uint32_t n, uint32_t scramble,
                           Float *p) 
{
    uint32_t v = scramble;

    for (uint32_t i = 0; i < n; ++i) 
    {
#ifndef PBRT_HAVE_HEX_FP_CONSTANTS
        p[i] = std::min(v * Float(2.3283064365386963e-10) /* 1/2^32 */,
                        OneMinusEpsilon);
#else
        p[i] = std::min(v * Float(0x1p-32) /* 1/2^32 */,
                        OneMinusEpsilon);
#endif
        // C ^ a, 根据格雷码的性质, 每次改变 a 的一个比特位得到新的采样值(根据 i+1 的尾随 0 的数量, 决定变化的比特位)
        // 对比式 7.10, 实际效果是从结果里加上或减去矩阵的一列
        v ^= C[CountTrailingZeros(i + 1)];
    }
}

inline void GrayCodeSample(const uint32_t *C0, const uint32_t *C1, uint32_t n,
                           const Point2i &scramble, Point2f *p) 
{
    uint32_t v[2] = {(uint32_t)scramble.x, (uint32_t)scramble.y};

    for (uint32_t i = 0; i < n; ++i) 
    {
#ifndef PBRT_HAVE_HEX_FP_CONSTANTS
        p[i].x = std::min(v[0] * Float(2.3283064365386963e-10), OneMinusEpsilon);
        p[i].y = std::min(v[1] * Float(2.3283064365386963e-10), OneMinusEpsilon);
#else
        p[i].x = std::min(v[0] * Float(0x1p-32), OneMinusEpsilon);
        p[i].y = std::min(v[1] * Float(0x1p-32), OneMinusEpsilon);
#endif
        v[0] ^= C0[CountTrailingZeros(i + 1)];
        v[1] ^= C1[CountTrailingZeros(i + 1)];
    }
}

// \param nSamplesPerPixelSample 这个采样点上需要的样本数量
// \param nPixelSamples 每个像素的采样点数量
// \param samples 待填充的数组
inline void VanDerCorput(int nSamplesPerPixelSample, int nPixelSamples,
                         Float *samples, RNG &rng) 
{
    // Define _CVanDerCorput_ Generator Matrix
    /*
        这实际上是一个 32x32 的单位矩阵，但因为每个元素都是0或1，可以用一个比特位来表示，
        参考 P458 的式 7.10 的推导结果（省略了中间的合并），这里把**每一列**都压缩到了一个 uint32_t 里 
        因为 P458 页的优化, 这个翻转了矩阵 C 的每一列

        列方向<----------------
                              |
                              |
                              |
                              |
                              |
                              |
                              v
                            行方向
    */
    const uint32_t CVanDerCorput[32] = 
    {
#ifdef PBRT_HAVE_BINARY_CONSTANTS
     // clang-format off
      0b10000000000000000000000000000000,
      0b1000000000000000000000000000000,
      0b100000000000000000000000000000,
      0b10000000000000000000000000000,
      // Remainder of Van Der Corput generator matrix entries
      0b1000000000000000000000000000,
      0b100000000000000000000000000,
      0b10000000000000000000000000,
      0b1000000000000000000000000,
      0b100000000000000000000000,
      0b10000000000000000000000,
      0b1000000000000000000000,
      0b100000000000000000000,
      0b10000000000000000000,
      0b1000000000000000000,
      0b100000000000000000,
      0b10000000000000000,
      0b1000000000000000,
      0b100000000000000,
      0b10000000000000,
      0b1000000000000,
      0b100000000000,
      0b10000000000,
      0b1000000000,
      0b100000000,
      0b10000000,
      0b1000000,
      0b100000,
      0b10000,
      0b1000,
      0b100,
      0b10,
      0b1,
      // clang-format on
#else
        0x80000000, 0x40000000, 0x20000000, 0x10000000, 0x8000000, 0x4000000,
        0x2000000,  0x1000000,  0x800000,   0x400000,   0x200000,  0x100000,
        0x80000,    0x40000,    0x20000,    0x10000,    0x8000,    0x4000,
        0x2000,     0x1000,     0x800,      0x400,      0x200,     0x100,
        0x80,       0x40,       0x20,       0x10,       0x8,       0x4,
        0x2,        0x1
#endif  // PBRT_HAVE_BINARY_CONSTANTS
    };

    int totalSamples = nSamplesPerPixelSample * nPixelSamples;
    uint32_t scramble = rng.UniformUInt32();
    GrayCodeSample(CVanDerCorput, totalSamples, scramble, samples); // 一次生成了所有的采样点

    // Randomly shuffle 1D sample points
    // P464, 打乱每个维度的样本, 去掉不同采样模式间的相关性???
    for (int i = 0; i < nPixelSamples; ++i)
        Shuffle(samples + i * nSamplesPerPixelSample, nSamplesPerPixelSample, 1, rng);

    Shuffle(samples, nPixelSamples, nSamplesPerPixelSample, rng);
}

inline void Sobol2D(int nSamplesPerPixelSample, int nPixelSamples,
                    Point2f *samples, RNG &rng) 
{
    // Define 2D Sobol$'$ generator matrices _CSobol[2]_
    /*static???*/ const uint32_t CSobol[2][32] = 
    {
        {0x80000000, 0x40000000, 0x20000000, 0x10000000, 0x8000000, 0x4000000,
         0x2000000, 0x1000000, 0x800000, 0x400000, 0x200000, 0x100000, 0x80000,
         0x40000, 0x20000, 0x10000, 0x8000, 0x4000, 0x2000, 0x1000, 0x800,
         0x400, 0x200, 0x100, 0x80, 0x40, 0x20, 0x10, 0x8, 0x4, 0x2, 0x1
        },
        {0x80000000, 0xc0000000, 0xa0000000, 0xf0000000, 0x88000000, 0xcc000000,
         0xaa000000, 0xff000000, 0x80800000, 0xc0c00000, 0xa0a00000, 0xf0f00000,
         0x88880000, 0xcccc0000, 0xaaaa0000, 0xffff0000, 0x80008000, 0xc000c000,
         0xa000a000, 0xf000f000, 0x88008800, 0xcc00cc00, 0xaa00aa00, 0xff00ff00,
         0x80808080, 0xc0c0c0c0, 0xa0a0a0a0, 0xf0f0f0f0, 0x88888888, 0xcccccccc,
         0xaaaaaaaa, 0xffffffff
        }
    };

    Point2i scramble;
    scramble[0] = rng.UniformUInt32();
    scramble[1] = rng.UniformUInt32();

    GrayCodeSample(CSobol[0], CSobol[1], nSamplesPerPixelSample * nPixelSamples,
                   scramble, samples);

    for (int i = 0; i < nPixelSamples; ++i)
        Shuffle(samples + i * nSamplesPerPixelSample, nSamplesPerPixelSample, 1, rng);

    Shuffle(samples, nPixelSamples, nSamplesPerPixelSample, rng);
}

#pragma endregion used by ZeroTwoSequenceSampler and MaxMinDistSampler



#pragma region used by SobolSampler

inline uint64_t SobolIntervalToIndex(const uint32_t m, uint64_t frame,
                                     const Point2i &p) 
{
    if (m == 0) return 0;

    const uint32_t m2 = m << 1;
    uint64_t index = uint64_t(frame) << m2;

    uint64_t delta = 0;
    for (int c = 0; frame; frame >>= 1, ++c)
        if (frame & 1)  // Add flipped column m + c + 1.
            delta ^= VdCSobolMatrices[m - 1][c]; // 这个和下面的异或都类似于 MultiplyGenerator

    // flipped b
    uint64_t b = 
        (
            (   (uint64_t)((uint32_t)p.x) << m   ) | ((uint32_t)p.y)
        ) ^ delta;

    for (int c = 0; b; b >>= 1, ++c)
        if (b & 1)  // Add column 2 * m - c.
            index ^= VdCSobolMatricesInv[m - 1][c];

    return index;
}

inline Float SobolSample(int64_t index, int dimension, uint64_t scramble = 0) {
#ifdef PBRT_FLOAT_AS_DOUBLE
    return SobolSampleDouble(index, dimension, scramble);
#else
    return SobolSampleFloat(index, dimension, (uint32_t)scramble);
#endif
}

inline float SobolSampleFloat(int64_t a, int dimension, uint32_t scramble)
{
    CHECK_LT(dimension, NumSobolDimensions) <<
        "Integrator has consumed too many Sobol' dimensions; you "
        "may want to use a Sampler without a dimension limit like "
        "\"02sequence.\"";

    // SobolSample 的实现也类似于 MultiplyGenerator
    uint32_t v = scramble;
    for (int i = dimension * SobolMatrixSize; a != 0; a >>= 1, i++)
        if (a & 1) v ^= SobolMatrices32[i];

#ifndef PBRT_HAVE_HEX_FP_CONSTANTS
    return std::min(v * 2.3283064365386963e-10f /* 1/2^32 */,
                    FloatOneMinusEpsilon);
#else
    return std::min(v * 0x1p-32f /* 1/2^32 */,
                    FloatOneMinusEpsilon);
#endif
}

inline double SobolSampleDouble(int64_t a, int dimension, uint64_t scramble)
{
  CHECK_LT(dimension, NumSobolDimensions) <<
      "Integrator has consumed too many Sobol' dimensions; you "
      "may want to use a Sampler without a dimension limit like "
      "\"02sequence\".";

    uint64_t result = scramble & ~ - (1LL << SobolMatrixSize);
    for (int i = dimension * SobolMatrixSize; a != 0; a >>= 1, i++)
        if (a & 1) result ^= SobolMatrices64[i];

    return std::min(result * (1.0 / (1ULL << SobolMatrixSize)),
                    DoubleOneMinusEpsilon);
}

#pragma endregion used by SobolSampler

}  // namespace pbrt

#endif  // PBRT_CORE_LOWDISCREPANCY_H
