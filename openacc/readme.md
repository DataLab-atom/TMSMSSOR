# 数值计算与计算机应用 2023 : 求解大规模稀疏垂直线性互补问题的二步模基并行算法
这里是关于数值计算与计算机应用 2023 [求解大规模稀疏垂直线性互补问题的二步模基并行算法](https://computmath.cjoe.ac.cn/szjs/CN/10.12288/szjs.s2022-0868#4)的代码

# 实验
## 编译运行-环境搭建
请在确保你准备好了linux下的openacc运行环境nvidiahpc库以及nvidia运行时库(nvidiatoolkit)。
原编译环境如下
```
nvcc -V
out:
  nvcc: NVIDIA (R) Cuda compiler driver
  Copyright (c) 2005-2022 NVIDIA Corporation
  Built on Wed_Jun__8_16:49:14_PDT_2022
  Cuda compilation tools, release 11.7, V11.7.99
  Build cuda_11.7.r11.7/compiler.31442593_0

nvc++ -V
  out:
  nvc++ 22.7-0 64-bit target on x86-64 Linux -tp haswell
  NVIDIA Compilers and Tools
  Copyright (c) 2022, NVIDIA CORPORATION & AFFILIATES.  All rights reserved.

```

## 更直接的方式
在确保环境准备就绪得到情况下你可以直接通过```sh run.sh``` 获取分析所需数据。
此外我们还准备了一份不是很完备的分析代码详情请参考`plot.ipynb`。

## 更个性化的方式
### 单独运行实例
```
// .MCMPLUS
cd openacc
numworkers = 16
m = 256
runAlgorithmTag = 2
gpu_id = 0
./main_acc $numworkers $m $runAlgorithmTag $gpu_id
//参数列表明细：
numworkers  //使用的线程数
m    //实际运算时方阵A的行数以及列数为 m*m
runAlgorithmTag // 请在 1 和 2 中间选择 分别对应 mcm 以及 mcm plus
gpu_id //当你拥有多块gpu时可以通过此项设置，选择使用哪一块gpu, 如果只有一块gpu 请将此项置为0
```

也可以自行通过nvc++编译`main_acc.c`来实现一些个性化的设置。
### 自行编译指令示例
```
clear & nvc -o main_acc -acc main_acc.c -Minfo
```

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

# 其它
如果还有其它问题可以提交issue，也可以通过我的个人邮箱`wen_hai_bin@163.com`跟我联系。

