<!--
-->

# pbrt-v3 阅读笔记(三): 更多细节

tree_crown.jpg

本篇就真的是笔记了, 记录了前两篇文章未提及的不少阅读细节和笔者理解以作备忘, 并(应该会)持续更新. 文章的各个部分可以说是有相互联系但又各自独立的, 您可以按自己的兴趣分别阅读各个部分.



- Interaction 解耦了 Geometry 和系统的其它部分, 当通过 ray casting 在场景中找到一个交点后, 我们就暂时不需要再关注场景中的 shapes 了.
- BRDF 解耦了 Shading 和 Integrator, 当通过 Interaction 上关联的 Material 计算出对应的 BRDF 后, 同样的, 我们不用再关注材质, 结合 Interaction 上的其它信息, 我们不用再关注 primitives 了.

类似的, sampler, camera, light 等都通过各自的接口隐藏了其细节, 让我们可以把注意力集中在某一点上, 而不必了解整个系统的所有细节
通过接口把细节隐藏起来, 不然没法读下去

# Geometry

## 几何 & 变换

### Transform

(坐标)变换在实时渲染里应该能更快理解含义

unity.Transform, unreal.Transform

## shape/primitive

### 浮点误差

todo

## 光线加速结构



# spectrum/radiance



# sample/filter/denoise



# Camera/Light

## Camera

## Light

### 区域光源 AreaLight

在 pbrt-v3 里, 如果某个 shape 是 AreaLight 的话, 那它必须一整个都是 AreaLight. 

# Shading

## BRDFs

## 材质

## 纹理

### 纹理过滤, 圆锥跟踪, 光线微分


圆锥型传感器 => 光线型传感器

# Rendering


