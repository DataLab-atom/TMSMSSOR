# 数值计算与计算机应用 2023 : 求解大规模稀疏垂直线性互补问题的二步模基并行算法
这里是关于数值计算与计算机应用 2023 [求解大规模稀疏垂直线性互补问题的二步模基并行算法](https://computmath.cjoe.ac.cn/szjs/CN/10.12288/szjs.s2022-0868#4)的代码

# 实验
## 编译运行-环境搭建
###这里是openmp的版本 还有一份openacc的版本详情请参考`openacc`


Linlux 下使用详情请参见`帮助文档.docx`。

Ubuntu 系统下请将`帮助文档.docx`中第三步的编译指令换成如下指令。
```
gcc Main.c -o main -lm -f openmp
```
Windows 系统下请参考知乎博客：
[Visual Studio 2019 下配置 OpenMP 多线程编程环境](https://zhuanlan.zhihu.com/p/86708660)进行环境搭建，然后可直接使用Visual Studio 2019运行代码。

## 详细参数
在程序成功运行时会输出 
```
out_the_tank //openmp 环境未成功加载
in_the_tank  //openmp 环境成功加载 直接无法运行也是环境未准备就绪的情况
```
注意此项指标是程序是否正常运行的关键如果最终程序以`环境未成功加载`的形式运行，将无法复现论文中的效果。

此外你可以在`Main.c line 1578 - 1587`处更改运行时的超级参数。



# 如果你觉得这对你的工作有帮住，请引用我们的论文
```
@article{郑华:350,
author = {郑华},
author = { 温海斌},
author = { 卢晓平},
title = {求解大规模稀疏垂直线性互补问题的二步模基并行算法},
publisher = {数值计算与计算机应用},
year = {2023},
journal = {数值计算与计算机应用},
volume = {44},
number = {4},
eid = {350},
pages = {350-367},
keywords = {垂直线性互补问题|二步方法|多分裂|并行},
url = {https://computmath.cjoe.ac.cn/szjs/CN/10.12288/szjs.s2022-0868},
doi = {10.12288/szjs.s2022-0868},
}
```
