<!--
假设读者已经了解了一点内容, 不然要讲的东西太多了
面向知识本身吧, 不然太难写了
讲到直接光照

前文尽量不要依赖于后文

## 蒙特卡洛积分
-->

# pbrt-v3 阅读笔记(一): 主要知识


## 基础概念

root.jpg

本篇文章说是对 pbr 的阅读笔记, 笔者也想借此用自己的理解来讲讲"基于物理的渲染"是什么一回事, 以及对这本书的阅读建议. 以下内容假设读者有对图形渲染(如 Phong 着色模型, 针孔相机, 点光源), 光线跟踪算法的基本了解, 因此一些基本概念就略过不讲了. 但还未接触过基于物理的渲染(光谱, 辐射度, 亮度, BRDF, 反射方程, 光线传输方程), 

pbrt 里有很多优化, 初学的时候要避开
如果有阅读 <> 过书籍, 课程, 对 pbr 的基本概念有了解, 就更好了

```C++
// core/geometry.h

class Vector2;
class Vector3;

class Point3;
class Point2;

class Normal3;

class Ray;

class Transform;

```

限于篇幅, 这里没有按历史介绍光线投射(ray cast), 递归光线跟踪(recursive ray tracing), 随机光线跟踪(stochastic ray tracing) 等内容, 以上算法可以参考 [Ray Tracing Essentials, Part 1: Basics of Ray Tracing(中文字幕)][RTE1], [全局光照技术进化史1-光线追踪篇][GIBook]

## 1. 经典的 Phong 光照模型(向量, 矩阵, 坐标/空间变换, 相机, 图像)

## 2. 简单的光线跟踪程序(光线跟踪, 光线-对象相交, )

回顾一个简单的光线跟踪实例, 我们会设置一个场景(scene), 场景中包含相机(camera), 若干个光源(lights)和要渲染的模型(models)

图一:

```C++
// 一个典型的例子 raytracing_v1:

class Camera;
class Model;
class Light;

class Scene
{
    Camera camera;
    Model  models[];
    Light  lights[];
};

class Renderer
{
    void render(const Scene& scene);
};
```



## 3. 基于物理的渲染

能量守恒

![](img/wikipedia_raytracing_variant.png)

<font size=2 color="gray">wikipedia, 略有修改</font>

让我们从这个场景开始讲起吧, 可以看到



很多简化, 假设, 在一开始学习时是有帮助的, 就好像我们中学时学的是牛顿三大定律而非相对论, 几何光学而非波动光学一样, 可以

如果把这个程序和现实世界联系在一起, 会发现它有很多不足之处(或者说过于理想): 

- 光源是理想化的一个点, 而现实中最为广泛的是面积(和体积光源); 
- 着色计算中假设模型表面的漫反射是完美漫反射, 光泽反射是完美镜面反射;
- 模型表面的着色计算过于简单, 能表达的材质种类太少; 

在刚开始接触图形渲染时, 这些简化处理是可以帮助隐藏细节, 理解主要原理的. 但是当我们希望渲染出更真实的图像时, 这些简化处理就限制了渲染器的表现能力. 那么很自然的, 如果我们能让渲染器(的各个组件)更真实去模拟现实环境中的效果, 也就是基于实际的物理法则去实现这个渲染器, 也许就能产生更为逼真的效果, 实际上也确实如此. 简单来说, 就是要"基于物理渲染"

那么在实现这样一个 "Physical Based Renderer" 之前, 不妨先来看看现实世界的"实现原理"

总而言之, 现在的这个程序有太多理想的假设和简化, 渲染出的场景还不够真实, 更进一步来说, 是不足以模拟真实世界中人的视觉感知.那不妨先来看看人的视觉是怎么产生的

在进入下一节之前, 让我们对式xx 的描述进行一些修改:


你可以以后再回过头来看这个修改的意义

## 资源和引用

[如何开始学习图形学编程](https://zhuanlan.zhihu.com/p/55518151)([原文](https://erkaman.github.io/posts/beginner_computer_graphics.html))
[系统的学习计算机图形学，有哪些不同阶段的书籍的推荐？](https://www.zhihu.com/question/26720808)
[]()
[]()
[]()


[Ray Tracing Essentials, Part 1: Basics of Ray Tracing(中文字幕)][RTE1]: 通过视频简介光线跟踪算法的历史
[全局光照技术进化史1-光线追踪篇][GIBook]: <全局光照技术>一书中对光线跟踪的历史简介

[Ray Tracing in One Weekend Series](https://raytracing.github.io/): 大概是目前最简洁的光线跟踪教程
[Learn OpenGL](https://learnopengl-cn.github.io/)([原博客](https://learnopengl.com/))


[编译 pbrt-v3]()
[编译 mitsuba]()


[RTE1]: https://developer.nvidia.com/rtx/raytracing-essentials-part1/cn
[GIBook]: https://zhuanlan.zhihu.com/p/24063586

# 场景一: 只包含相机和光源

![](img/wikipedia_raytracing_emitted1.png)

<!--
需要的时候才加上解释
-->
<font size=2 color="gray">图 wikipedia_raytracing_emitted1, 这里请先忽略掉右下的物体</font>

下面我们就先考虑整个场景只包含光源的情况, 这种情况下 Film 上自然也只有光源, 这样就只需要考虑相机和光源的关系了, 即(这里为了统一, 光线是从相机指向光源的, 与现实中相反):

![](img/wikipedia_raytracing_emitted2.png)


虽然图上画了一个比较大的光源, 但为了方便分析, 我们可以把它当作是理想的点光源, 那光源射向相机的, 也就只有一条光线.

## 相机的成像原理


### 1. 针孔相机
在分析人眼之前, 我们不妨先来看看相机的成像原理

当然, 在实际开发时, 我们一般把图像平面 Film 放在 Camera 前面, 这样更方便计算

<details><summary><b>* 拓展 </b></summary><p>

针孔相机确实是可以工作的

</p></details>


### 2. 数码相机

虽然针孔相机也可以工作, 但现实中常用的数码相机, 一般是采用薄透镜模组作为

图像处理器: 光信号转换为电信号

<!--
数码相机的简单理解
薄透镜
DCC, CMOS 传感器的工作原理

上面是一张人眼的示意图, 晶状体, 视网膜...比较, 我们可以直观的看到, 人眼的成像其实类似于数码相机, 视网膜接受到光信号后, 由大脑形成最后的视觉, 也就是各种各样的颜色. 这里有一点非常重要, 那就是颜色是在人的视觉中枢中形成的, 视网膜上接受的只是一段电磁波而已.
-->

### 3. 人的视觉

在简单分析了数码相机的成像原理后, 我们来看看人的视觉是如何工作

人眼的成像原理.

![](img/eye.jpg)

外界的光信号最后会, 最后光线汇聚到眼球后方的视网膜上

![](img/eye_nerve.jpg)

在这之后, 视网膜还会继续通过, 将传送给视觉神经中枢, 在视觉中枢中才会形成最后的视觉. 在这个过程中, 眼球同时承担了收集和采集光信号的作用.

![](img/camera_eye_light.jpg)

<font size=2 color="gray"> 数码相机的传感器应该集成在主板上, 左侧 '3' 代表的导线只是方便示意而已 </font>

我们可以通过上面这张图来比较数码相机和人类的视觉系统, 两者有很大的相似之处, 都通过一个镜头(薄透镜/视网膜)收集光信号, 通过传感器(图像传感器/视网膜)采集信号, 处理器(图像处理器/视觉神经中枢)来处理信号. 当然, 以上都是简略的分析, 还有很多细节未涉及到.

<details><summary><b>* 拓展</b></summary><p>
</p></details>


在这里笔者想强调一个重点, 那就是视觉是在视神经中枢形成的, 我们所看到的五颜六色的世界, 眼球接受到的仅仅只是光信号而已, 换句话说, **颜色是人的主观感觉, 可见光并不是颜色**

<!--
不知道读到这里, 您会不会有一点反应不过来, 
-->
笔者在第一次认识到这一点时, 是有点反应不过来的, 毕竟无论是自己长期的生活经验, 还是中学的光学课程, 都让笔者形成了可见光就是颜色的看法, 忽然要重新认识这两者的关系, 打破以往长期的生活经验, 还是有点不适应的.(多少有点像第一次听到相对论的概念一样)

如果您也有类似困惑的话, 请继续阅读下面的拓展部分

<details><summary><b>* 拓展</b></summary><p>

### 视网膜

这里不妨多举几个例子, 譬如狗只有..., 只能看到黑白的世界; 墨鱼/鸟类的视网膜上有四种感光细胞, 可以感知到更多的光线, 有比人眼更丰富的颜色; 蜜蜂能看到紫外线, 帮助它们寻找花粉...

<!--
**可见光是客观存在的**, 是空间中的一段电磁波. 而大千世界的生物具体感知到了什么, 则取决于它们自己, 具体到人类, 看到了颜色, **颜色是人的主观感受**

不过为了实现基于物理的渲染, 与笔者有着类似观念的读者们, 需要好好重新认识可见光和颜色的关系
-->

</p></details>



## 从光源射向相机的光线

在上一小结, 我们简单介绍了成像原理, 以及光和颜色的关系, 下面

人眼所能看到的光, 是电磁波中 300nm 到 700nm 左右的一段, 也就是可见光是电磁波的一个子集. 那么电磁波具有什么样的属性呢? 频率和波长.

这个时候如果再联系光谱图, 

<!--
首先简单来说, 我们能看见物体, 是因为光线或经由空间中的介质(空气, 水)传播直接进入到了人眼, 或经由物体表面的反射, 折射间接照射到人眼
-->


### 颜色/能量

瞬间的能量分布, 严格一点来说, 就是单位时间的能量(不用一开始详细的讲)

### 辐射度和光度

## 光源

光源暂且没什么好讲的, 就是在不停向外辐射能量

## 回顾场景一

现在再来看场景一, 我们知道了 

<!--
怎么介绍自发光呢, 放场景二讲完区域光源后再说
-->

## small pbrt v1


## 资源和引用

[]()


# 场景二: 加入物体

![](img/scene_direct_lighting.png)

<font size=2 color="gray">图 scene_direct_lighting, 考虑到右上的光源是个点光源, 我们在 Film 上就忽略掉它吧 </font>

我们管这种光线从光源发射后, 经过单个物体反射后才入射相机的情况, 叫"直接光照"(可以理解为光线直接照射到了物体上, 反正就是一种约定的叫法, 对于场景一中光线从光源照射到相机的情况, 有另外的叫法, 下文会将), 而对于光线在物体间经过多次反射后才抵达相机的情况, 则称为间接光照, 间接光照会留待最后才分析.

这种情况下, 我们应该会对光线经过物体表面时的反射比较感兴趣:

![](img/scene_direct_lighting_model.png)

简单联系现实生活, 就会发现不同物体表面的反射有非常大的差别

## 物体表面的材质属性

<details hidden><summary></summary>

#### > BRDF

</details>

### 1. 漫反射

### 2. 镜面反射

### 3. 光泽反射

当光线入射至物体表面时，表面根据其材质特性, 会对光线进行不同种类的反射(我们暂且不考虑折射等现象), 比较理想的两种模型, 是完全漫反射和完全镜面反射, 此外现实中比较常见的则是光泽反射 glossy reflect:

![](img/wikipedia_brdf.png)

这里如果假设反射发生的位置是 p, 入射光线(黑色)的方向是 wi(incoming direction), 那么给定一个反射光线的方向 wo(outgoing direction), 对照图 wikipedia_brdf, 就有以下情况:

- Diffuse: 给定任意出射方向 wo, 出射方向上的辐射度都是同样大的
- Glossy: 辐射度在 wi 的反射向量上最大
- Mirror: 辐射度只分布在 wi 的反射方向上, 其它方向都为 0

我们可以用一个函数来描述这些反射在空间上的分布情况, 即给定位置 p, 入射方向 wi 和出射方向 wo, 这个函数可以告诉我们入射方向的辐射度在出射方向上的分布, 即:

$$
f_{\mathrm{r}}\left(\mathrm{p}, \omega_{\mathrm{o}}, \omega_{\mathrm{i}}\right)=\frac{\mathrm{d} L_{\mathrm{o}}\left(\mathrm{p}, \omega_{\mathrm{o}}\right)}{\mathrm{d} E\left(\mathrm{p}, \omega_{\mathrm{i}}\right)}=\frac{\mathrm{d} L_{\mathrm{o}}\left(\mathrm{p}, \omega_{\mathrm{o}}\right)}{L_{\mathrm{i}}\left(\mathrm{p}, \omega_{\mathrm{i}}\right) \cos \theta_{\mathrm{i}} \mathrm{d} \omega_{\mathrm{i}}} \tag{BRDF}
$$

<!--
这里要讲点微积分吗
找张图来解释这个公式
再放点代码
-->

在图形学中, 这个函数被称为**双向反射分布函数 BRDF(),** 一般记作 f(p, wi, wo), 当给定 BRDF 函数满足能量守恒(即出射能量 <= 入射能量) 和互异性(交换入射和出射方向, 函数返回值不变)时, 我们称这个 BRDF 是基于物理的, 即 Physical based BRDF. 

通过 BRDF, 我们可以描述 lambert 漫反射为 

$$
f_{\mathrm{lambert}}\left(\mathrm{p}, \omega_{\mathrm{o}}, \omega_{\mathrm{i}}\right)= albedo / pi
$$

这里的 albedo, 指的是漫反射表面的反射光与入射光的比率, 也就是表面反射能量的多少, 简称反射率
<!--
关于 albedo, 之后再详细解释吧
理解为, 看作是
-->

而完美镜面反射则要特殊一点, 因为它只在入射光的反射方向上有出射光线, 所以:

$$
f_{\mathrm{specluar}}}\left(\mathrm{p}, \omega_{\mathrm{o}}, \omega_{\mathrm{i}}\right)= 0, if wo != wr, where wr is reflect direction of wi
$$

有了 BRDF 以后, 我们就可以重写式 xx 了

关于材质还有很多内容可以讲述, 但目前介绍的知识足够我们继续讲下去, 有关材质的更多细节就留待之后的文章讲吧

## small pbrt v2

## 资源和引用

[brdf为什么要定义为一个单位是sr-1的量？](https://www.zhihu.com/question/28476602)

[]()

[]()

[]()

[]()



# 场景三 加入更多光源

## 加入更多点光源

```C++
// ...

class Renderer
{
    void render(const Scene& scene);
};
```

~熟悉微积分的同学也许能很快意识到, 我们要计算的直接光照, 其实是一个积分方程. 但考虑到应该有不少人~我们从单个点光源开始拓展, 当在场景中增加光源时, 原来单个光源的计算, 变成了以下的求和:

formula

## 加入面积光源

<details hidden><summary></summary>

#### target: 反射方程

</details>

## 用数值方法求解积分

我们先暂且不去理会上面的反射方程, 先去解决另外一个问题: 如何求解积分?

解析法, 数值法

对于

## 求解反射方程

限于篇幅, 这里暂且先略过蒙特卡洛积分的推导过程, 只从概念上帮助大家来理解
<!--TODO
~~限于篇幅, 这里暂且先略过蒙特卡洛积分的推导过程, 只从概念上帮助大家来理解~~
-->

### 1. 选取光源

<details hidden><summary></summary>

#### > 逆变换算法

</details>



### 2. 在选取概率高的位置放置更多采样点

<details hidden><summary></summary> 

#### > 重要性采样

</details>



### 3. 根据光源和 BRDF 的特征选取采样位置

<details hidden><summary></summary>

#### > 多重重要性采样

</details>



## small pbrt v3

"我无法实现的, 我就无法理解" -- 费曼

这篇文章写了这么长, 都还没提到所谓的"pbrt-v3 阅读笔记", 现在我们可以回到正题, 

来实现一个麻雀版的 pbrt 了

```C++
// 极为简化的 pbrt-v3, 即 small_pbrt_v1:

class Light;
class Camera;
class Model;

class Scene
{
    Camera camera;
    Light  lights[];
    Model  models[];
};

class Renderer
{
    void render(const Scene& scene);
};
```

目前这个程序只有一堆接口, 主要是为了表达概念, ~还无法运行, 在下一篇文章里笔者会尝试去实现一个可以跑起来的程序~没有给出实现, 在下一篇文章中, 我们将结合具体的代码来讲解. 另外对文中的错漏之处, 欢迎各位读者批评指正, 笔者也将持续改进本系列文章.


## 小结



## 资源和引用

[解析解、闭合解、数值解](http://blog.sina.com.cn/s/blog_65838c3c0101e7tg.html)





# * 场景四: 加入更多物体

<font size=4 color="gray">本节介绍了一些稍为复杂的内容, 可以留待以后再看</font>


cornell_box.jpg

这次我们换一个新的场景来分析, 上图是图形学中一个非常经典的场景(wiki), 我们的分析就从这里开始.

<!--
间接光照
-->

## LET

## 路径积分

## 求解路径积分

### 终止路径

<details hidden><summary></summary>

#### > 俄罗斯轮盘赌

</details>

理论上这条路径是可以无限递归下去的, 但在实际实现是我们必然要在一定的步长后终止路径, 但又不能影响到最终结果, 那么要怎么做呢?

按概率终止.


## small pbrt v4

## 资源和引用

[]()


# 总结

限于时间和经验, 很多地方讲的比较简略, 希望之后能慢慢补充

## 资源和引用

An Explanation of the Rendering Equation: 反射方程的视频介绍

<!--
全文引用
-->