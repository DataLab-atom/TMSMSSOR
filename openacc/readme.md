# 数值计算与计算机应用 2023 : 求解大规模稀疏垂直线性互补问题的二步模基并行算法
这里是关于数值计算与计算机应用 2023 [求解大规模稀疏垂直线性互补问题的二步模基并行算法](https://computmath.cjoe.ac.cn/szjs/CN/10.12288/szjs.s2022-0868#4)的代码

# 实验
## 编译运行-环境搭建
请在确保你准备好了linux下的openacc运行环境nvidiahpc库以及nvidia运行时库(nvidiatoolkit)。

在确保环境准备就绪得到情况下你可以直接通过```sh run.sh``` 获取分析所需数据

此外我们还准备了一份不是很完备的分析代码详情请参考`plot.ipynb`,以及可以自行通过nvc++编译`main_acc.c`来实现一些个性化的设置


# 如果你觉得这对你的工作有帮助，请引用我们的论文
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
