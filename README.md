Path Tracer 4CG
====
特性
----
* 求交与加速结构
  - KD树
  - Fast BVH树（thanks to [Brandon Pelfrey](https://github.com/brandonpelfrey/Fast-BVH)）
* 材质
  - Cook-Torrance
  - FrostedGlass（毛玻璃折射）
  - SmoothGlass（透明玻璃折射）
  - 理想镜面反射
  - 漫反射
* 光源
  - 点光源、平面光源、方向光源
* 体素
  - 三角网格
  - 球体
* OpenMP多线程


编译构建与使用
-----

场景文件
-----
自定义格式处理场景设置，支持球体和三角网格的模型体素
场景文件范例见testscenes/test

渲染范例
----
https://github.com/koscielny/path_tracer_4CG/raw/master/result_image/cornell_box_10.png

https://github.com/koscielny/path_tracer_4CG/raw/master/result_image/cornell_box_12.png

https://github.com/koscielny/path_tracer_4CG/raw/master/result_image/cornell_box_21.png

https://github.com/koscielny/path_tracer_4CG/raw/master/result_image/bunny_29.png


参考文献
-----
Jensen, H. W, et al. "Monte Carlo Ray Tracing Course Notes." of Complex Objects”, Computers and Graphics (2003):215--224.

Kajiya, James T. "The rendering equation." Acm Siggraph Computer Graphics 20.4(1986):143-150.

Pharr, Matt, and G. Humphreys. "Physically Based Rendering, Second Edition: From Theory To Implementation." (2010).

[Importance Sampling of the Phong Reflectance Model - Jason Lawrence from Princeton](http://www.cs.princeton.edu/courses/archive/fall12/cos526/papers/importance.pdf)