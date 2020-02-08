<!--
采样器, BSDF 用现成的就行了吧 

## 全局光照
-->

# pbrt-v3 阅读笔记(二): 全局光照

trunk.jpg

leaves.jpg

在[上一篇文章][overview]中我们介绍了 PBR 中的一些基本概念, 

引用一个 zhihu 上的回答作为本篇的主题, 也是下面要解决的主要问题, 那就是 这本书写的太厚导致读不下去.  


## 求解反射方程

## 全局光照

采样器的接口太复杂
有很多优化, 干扰理解
MIS 理解了就还好

## Windows 平台编译并运行 pbrt-v3

以下使用 cmake 和 visual studio 2017 编译笔者 fork 的 pbrt-v3(编译原作者的 pbrt-v3 步骤也是一样的)

1. 用 cmake 生成 vs2017 的工程文件(*.sln, *.vcproj 等)
2. 打开上一步生成的 pbrt.sln, 编译程序
3. 在设置中添加环境变量, 然后运行程序

笔者的电脑上只安装了 visual studio 2017 和 visual studio 2019, 且无法在 vs2019 上编译 pbrt-v3, 这里演示的是用 vs2017 来 pbrt-v3

[overview]:overview.md