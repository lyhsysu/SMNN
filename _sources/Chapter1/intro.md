# 自旋玻璃与空腔方法

- 作者：李宇豪
- 日期：2023年1月10日
- 更新：2023年1月13日
- 许可：<a rel="license" href="http://creativecommons.org/licenses/by-nc-sa/4.0/">CC BY-NC-SA 4.0</a>
  
---
`````{admonition} 摘要
:class: attention
本文首先介绍了使用空腔方法求解配分函数流程：(1) 将哈密顿量在因子图上表示；(2) 分别计算向系统中加入一个功能节点和一个变量节点时系统自由能的增量；(3) 考察功能节点和变量节点在计算自由能时的重叠部分，将各个节点的自由能增量不重复的累加成为系统的总自由能。然后介绍了如何用一套消息传递的方法迭代求解磁化强度，最后介绍了空腔方法在 RNN 中的应用，以及在其他方面的应用的参考文献。
`````
```{admonition} 预备知识：因子图
因子图是一类无向概率图模型, 包括变量节点和因子节点（或称为功能节点）。变量节点和功能节点之间用无向边连接。定义在因子图上的联合概率分布可以表示为各个因子的联乘积。

比如，对于一个服从下式的贝叶斯网络

$$
p(I, D, G, S, L)=P(I) P(D) P(G \mid I, D) P(L \mid G) P(S \mid I)
$$

实际上就是将一个联合概率分布 $p(I, D, G, S, L)$ 表达为多个局部的概率分布 $P(I)$、 $P(D)$、 $P(G \mid I, D)$、 $P(L \mid G)$、 $P(S \mid I)$ 的联乘。将各个概率分布抽象为一个函数，$f_I(I)$、 $f_D(D)$、 $f_G(G, I, D)$、 $f_S(L, G)$、 $f_L(S, I)$，在因子图中可以表示为

<p style="text-align:center"><img src="https://user-images.githubusercontent.com/106574511/211964166-c73b162c-dffa-47f2-9557-b4481bbce559.png" alt="" class="bg-primary" width="400px"></p>

其中，$I$、$D$、$G$ 等圆圈称为“变量节点”，代表函数的变量；$f_I$、$f_D$、$f_G$ 等正方形称为“功能节点”，代表连乘中的一项（一个局部的概率分布或其他函数）。
```

```{admonition} Assignment 1
用空腔方法计算 S-K 模型的自由能和磁化强度，并把结果与数值计算的结果进行比较。

详见 Appendix B
```

```{admonition} Assignment 2
a: 编写一个程序来实现 Soulas 码的编码和解码方案 $\left(K \equiv|\partial a|=3, \mathrm{R} \equiv \frac{N}{M}=0.5\right)$；

b: 然后考虑在特定温度 $\beta_p=\frac{1}{2} \ln \frac{1-p}{p}$ 下解码性能随着翻转概率 $p$。

详见 Appendix C
```

