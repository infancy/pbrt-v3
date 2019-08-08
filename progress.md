# 进度

## Code

### src/accelerators

- [ ] bvh （Bounding Volume Hierarchy，层次包围盒）
- [ ] kdtreeaccel (kdtree)

### src/camera

- [ ] environment
- [ ] orthographic	（正视投影相机）
- [ ] perspective （透视投影相机）
- [ ] realistic

### src/core

- [ ] api （pbrt 对外提供的 api）
- [ ] bssrdf
- [ ] camera
- [ ] efloat
- [ ] error
- [ ] fileuti
- [ ] film
- [ ] filter
- [ ] floatfile
- [ ] geometry （包含了 vector3、point3、normal3、bbox、ray 等基础类）
- [ ] imageio
- [x] integrator （积分器的基类，包含了核心的渲染算法）
- [ ] interaction （记录了光线与 surface 或 volume 的交点）
- [ ] interpolation
- [ ] light
- [ ] lightdistrib （光源功率的分布）
- [ ] lowdiscrepancy
- [ ] material
- [ ] medium
- [ ] memory
- [ ] microfacet
- [ ] mipmap.h
- [ ] parallel
- [ ] paramset
- [ ] parser
- [ ] pbrt.h （公共头文件）
- [ ] pbrtlex.cpp
- [ ] pbrtlex.ll
- [ ] pbrtparse
- [ ] pbrtparse.output
- [ ] pbrtparse.y
- [ ] primitive
- [ ] progressreporter （用于报告进度）
- [ ] quaternion
- [ ] reflection （brdf）
- [ ] rng.h （Random Number Generator，随机数生成器）
- [x] sampler 
- [ ] sampling （独立的采样算法）
- [ ] scene
- [ ] shape
- [ ] sobolmatrices
- [ ] spectrum （光谱类，可简单的视为一个 RGBColor）
- [ ] stats
- [ ] stringprint
- [ ] texture
- [ ] transform （矩阵类）

### src/filter

### src/integrators

- [ ] ao      （Ambient Occlusion，环境光遮蔽）
- [ ] bdpt    （Bidirectional Path Tracing，双向路径跟踪）
- [x] dl      （directlighting，直接光照）
- [ ] mlt     （Metropolis Light Transport，梅特波利斯光线传输）
- [x] path    （Path Tracing，路径跟踪）
- [ ] sppm	  （Stochastic Progressive Photon Mapping，随机渐进光子映射）
- [ ] volpath （Volume Path Tracing，体积路径跟踪）
- [ ] whitted （Recursive Ray Tracing，递归光线跟踪）

### src/materials

- [ ] disney
- [ ] fourier
- [ ] glass
- [ ] hair
- [ ] kdsubsurface
- [ ] matte
- [ ] metal
- [ ] mirror
- [ ] mixmat
- [ ] plastic
- [ ] substrate
- [ ] subsurface
- [ ] translucent 
- [ ] uber

### src/samplers

- [ ] halton 
- [ ] maxmin
- [ ] random
- [ ] sobol
- [ ] stratified
- [ ] zerotwosequence

施工中...

## Book

Preface
Further Reading
Preface to the Online Edition

1 Introduction
1.1 Literate Programming
1.2 Photorealistic Rendering and the Ray-Tracing Algorithm
1.3 pbrt: System Overview
1.4 Parallelization of pbrt
1.5 How to Proceed through This Book
1.6 Using and Understanding the Code
1.7 A Brief History of Physically Based Rendering
Further Reading
Exercises

2 Geometry and Transformations
2.1 Coordinate Systems
2.2 Vectors
2.3 Points
2.4 Normals
2.5 Rays
2.6 Bounding Boxes
2.7 Transformations
2.8 Applying Transformations
2.9 Animating Transformations
2.10 Interactions
Further Reading
Exercises

3 Shapes
3.1 Basic Shape Interface
3.2 Spheres
3.3 Cylinders
3.4 Disks
3.5 Other Quadrics
3.6 Triangle Meshes
3.7 Curves
3.8 Subdivision Surfaces
3.9 Managing Rounding Error
Further Reading
Exercises

4 Primitives and Intersection Acceleration
4.1 Primitive Interface and Geometric Primitives
4.2 Aggregates
4.3 Bounding Volume Hierarchies
4.4 Kd-Tree Accelerator
Further Reading
Exercises

5 Color and Radiometry
5.1 Spectral Representation
5.2 The SampledSpectrum Class
5.3 RGBSpectrum Implementation
5.4 Radiometry
5.5 Working with Radiometric Integrals
5.6 Surface Reflection
Further Reading
Exercises

6 Camera Models
6.1 Camera Model
6.2 Projective Camera Models
6.3 Environment Camera
6.4 Realistic Cameras
Further Reading
Exercises

7 Sampling and Reconstruction
7.1 Sampling Theory
7.2 Sampling Interface
7.3 Stratified Sampling
7.4 The Halton Sampler
7.5 (0, 2)-Sequence Sampler
7.6 Maximized Minimal Distance Sampler
7.7 Sobol’ Sampler
7.8 Image Reconstruction
7.9 Film and the Imaging Pipeline
Further Reading
Exercises

8 Reflection Models
8.1 Basic Interface
8.2 Specular Reflection and Transmission
8.3 Lambertian Reflection
8.4 Microfacet Models
8.5 Fresnel Incidence Effects
8.6 Fourier Basis BSDFs
Further Reading
Exercises

9 Materials
9.1 BSDFs
9.2 Material Interface and Implementations
9.3 Bump Mapping
Further Reading
Exercises

10 Texture
10.1 Sampling and Antialiasing
10.2 Texture Coordinate Generation
10.3 Texture Interface and Basic Textures
10.4 Image Texture
10.5 Solid and Procedural Texturing
10.6 Noise
Further Reading
Exercises

11 Volume Scattering
11.1 Volume Scattering Processes
11.2 Phase Functions
11.3 Media
11.4 The BSSRDF
Further Reading
Exercises

12 Light Sources
12.1 Light Emission
12.2 Light Interface
12.3 Point Lights
12.4 Distant Lights
12.5 Area Lights
12.6 Infinite Area Lights
Further Reading
Exercises

13 Monte Carlo Integration
13.1 Background and Probability Review
13.2 The Monte Carlo Estimator
13.3 Sampling Random Variables
13.4 Metropolis Sampling
13.5 Transforming between Distributions
13.6 2D Sampling with Multidimensional Transformations
13.7 Russian Roulette and Splitting
13.8 Careful Sample Placement
13.9 Bias
13.10 Importance Sampling
Further Reading
Exercises

14 Light Transport I: Surface Reflection
14.1 Sampling Reflection Functions
14.2 Sampling Light Sources
14.3 Direct Lighting
14.4 The Light Transport Equation
14.5 Path Tracing
Further Reading
Exercises

15 Light Transport II: Volume Rendering
15.1 The Equation of Transfer
15.2 Sampling Volume Scattering
15.3 Volumetric Light Transport
15.4 Sampling Subsurface Reflection Functions
15.5 Subsurface Scattering Using the Diffusion Equation
Further Reading
Exercises

16 Light Transport III: Bidirectional Methods
16.1 The Path-Space Measurement Equation
16.2 Stochastic Progressive Photon Mapping
16.3 Bidirectional Path Tracing
16.4 Metropolis Light Transport
Further Reading
Exercises

17 Retrospective and The Future
17.1 Design Retrospective
17.2 Alternative Hardware Architectures
17.3 Conclusion
Further Reading

A Utilities
A.1 Main Include File
A.2 Image File Input and Output
A.3 Communicating with the User
A.4 Memory Management
A.5 Mathematical Routines
A.6 Parallelism
A.7 Statistics
Further Reading
Exercises

B Scene Description Interface
B.1 Parameter Sets
B.2 Initialization and Rendering Options
B.3 Scene Definition
B.4 Adding New Object Implementations
Further Reading
Exercises

References
Index of Fragments
Index of Identifiers