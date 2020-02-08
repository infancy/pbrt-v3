<!--
-->

# pbrt-v3 阅读笔记(二): 更多细节

tree_crown.jpg

本篇就真的是笔记了, 记录了前两篇文章未提及的不少阅读细节和笔者理解以作备忘, 并(应该会)持续更新. 文章的各个部分可以说是有相互联系但又各自独立的, 您可以按自己的兴趣分别阅读各个部分.



- Interaction 解耦了 Geometry 和系统的其它部分, 当通过 ray casting 在场景中找到一个交点后, 我们就暂时不需要再关注场景中的 shapes 了.
- BRDF 解耦了 Shading 和 Integrator, 当通过 Interaction 上关联的 Material 计算出对应的 BRDF 后, 同样的, 我们不用再关注材质, 结合 Interaction 上的其它信息, 我们不用再关注 primitives 了.

类似的, sampler, camera, light 等都通过各自的接口隐藏了其细节, 让我们可以把注意力集中在某一点上, 而不必了解整个系统的所有细节
通过接口把细节隐藏起来, 不然没法读下去


# 引论和模块划分

第一章的 1.7 小结和每一章的 Further Reading 都会介绍这一章内容的历史发展, 梳理脉络, 非常有意义

## 资源和引用

[]()

# 第一部分: Geometry

## 几何 & 变换

### Transform

(坐标)变换在实时渲染里应该能更快理解含义

unity.Transform, unreal.Transform

## shape/primitive

### 浮点误差

todo

## 光线加速结构

## 资源和引用

[]()


# 第二部分: spectrum/radiance


## 资源和引用

[]()

# 第三部分: sample/filter/denoise


## 资源和引用

[]()

# 第四部分: Camera/Light

## Camera

## Light

### 区域光源 AreaLight

在 pbrt-v3 里, 如果某个 shape 是 AreaLight 的话, 那它必须一整个都是 AreaLight. 

## 资源和引用

[]()


# 第五部分: Shading

## BRDFs

在 xxx 中我们介绍了 BRDF 的基本概念, 这里将进行拓展

基于经验, 测量, 物理定律

## 材质

## 纹理

### 纹理过滤, 圆锥跟踪, 光线微分


圆锥型传感器 => 光线型传感器

## 资源和引用

[基于物理着色：BRDF](https://zhuanlan.zhihu.com/p/21376124)

[BRDF - Wakapon](http://wiki.nuaj.net/index.php?title=BRDF)

[More DETAILS! PBR的下一个发展在哪里？](https://zhuanlan.zhihu.com/p/32209554)

[]()

[]()

[]()

[]()



# 第六部分: Rendering

## 资源和引用

[]()



# more



# 资源和引用

## pbrt

- code: https://github.com/mmp/pbrt-v3
- book: http://www.pbr-book.org/3ed-2018/contents.html
- errata: https://www.pbrt.org/errata-3ed.html
- courses: https://www.pbrt.org/courses.html

## 中文翻译&笔记

- https://github.com/zq317157782/raiden
- https://book.douban.com/subject/26974215/
- https://zhuanlan.zhihu.com/wayonpbrt
- https://www.zhihu.com/people/lou-jia-jie-95/posts
- https://www.qiujiawei.com/

# detail

- Moving Frostbite to PBR
- https://www.scratchapixel.com/lessons/digital-imaging/colors

## light sources

## light transport

- https://www.scratchapixel.com/lessons/3d-basic-rendering/global-illumination-path-tracing

## monte carlo integration

- https://www.scratchapixel.com/lessons/mathematics-physics-for-computer-graphics/monte-carlo-methods-in-practice/monte-carlo-methods

## ray–object intersections, acceleration structures and more

- 用线性代数知识解决光线和三角形的交点问题: https://www.qiujiawei.com/triangle-intersect/

- 微分几何与渲染(1): https://www.qiujiawei.com/partial-derivatives/

- https://www.qiujiawei.com/bvh-1/
- https://www.qiujiawei.com/bvh-2/

## material and texture

- https://zhuanlan.zhihu.com/p/53086060
- https://belcour.github.io/blog/research/2018/05/05/brdf-realtime-layered.html

## surface scattering

- 基于物理着色: https://zhuanlan.zhihu.com/p/20091064
- 基于物理着色：BRDF https://zhuanlan.zhihu.com/p/21376124

## camera

## sampling and reconstruction

- 低差异序列（一）- 常见序列的定义及性质: https://zhuanlan.zhihu.com/p/20197323
- 低差异序列（二）- 高效实现以及应用: https://zhuanlan.zhihu.com/p/20374706
- http://extremelearning.com.au/unreasonable-effectiveness-of-quasirandom-sequences/
