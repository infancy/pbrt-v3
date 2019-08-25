# pbrt-v3 代码&书籍阅读笔记

## 3 shapes

### 3.1 basic shape interface

3.1.5 Sidedness

 The potential for improved performance is reduced when using this technique with ray tracing, however, 
 since it is often necessary to perform the ray–object intersection before determining the surface normal to do the back-facing test. 
 Furthermore, this feature can lead to a physically inconsistent scene description if one-sided objects are not in fact closed. 
 支持背面剔除会导致潜在的性能损失???, 针对非封闭的图元还会带来渲染的不真实性. 
 传统的光栅化渲染流水线计算的往往是局部光照, 当使用全局光照模型时(real time raytracing)也需要开始考虑这个问题吧.