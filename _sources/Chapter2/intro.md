# 变分平均场理论与信念传播算法

- 作者：李宇豪
- 日期：2023年1月12日
- 更新：2023年1月13日
- 许可：<a rel="license" href="http://creativecommons.org/licenses/by-nc-sa/4.0/">CC BY-NC-SA 4.0</a>

---

`````{admonition} 摘要
:class: attention
摘要
`````

```{admonition} Assignment 1
关于树状图中Bethe近似对应的分布

在本篇Sec.3.1中，我们使用对树状图划分区域的方式得到了Bethe自由能

$$
\begin{aligned}
F_{\text {Bethe }}
= & -\sum_a \sum_{x_a} b_a\left(\boldsymbol{x}_a\right) \ln f_a\left(\boldsymbol{x}_a\right)+\sum_a \sum_{x_a} b_a\left(\boldsymbol{x}_a\right) \ln b_a\left(\boldsymbol{x}_a\right) \\
&-\sum_i\left(d_i-1\right) \sum_{x_i} b_i\left(x_i\right) \ln b_i\left(x_i\right)
\end{aligned}
$$

实际上，根据以上结果，以及变分自由能的定义

$$
F(b)=U(b)-H(b)
$$

其中，$U(b)=\sum_x b(x) E(x) ,\; H(b)=-\sum_x b(x) \ln b(x)$

我们可以直接总结出 $b(x)$ 的形式为

$$
b_{B A}(\boldsymbol{x})=\frac{\prod_a b_a\left(\boldsymbol{x}_a\right)}{\prod_i b_i\left(x_i\right)^{d_i-1}}
$$

文献[2]中提到这是在文献[3]中提出的一个“广为人知”的结论，但是我目前还没有理解文献[3]中的内容，留待以后补充。

关于 $b_{BA}(x)$，有两个小题：(1) 证明 $b_{BA}(x)$ 是归一化的；(2) 证明在无环的树状图中，$b_{BA}(x)$ 是精确的
```

```{admonition} Assignment 2
用Belief Propagation求解SK model中的 $C_{ij}$
```

```{admonition} Assignment 3
用Bethe Approximation求解inverse Ising problem
```

```{admonition} Reference

[1] Yedidia, Jonathan S., William T. Freeman, and Yair Weiss. "[Understanding belief propagation and its generalizations.](https://www.cs.huji.ac.il/course/2005/pmai/tirguls/TR2001-22.pdf)" *Exploring artificial intelligence in the new millennium* 8.236-239 (2003): 0018-9448.

[2] Yedidia, Jonathan S., William T. Freeman, and Yair Weiss. "[Constructing free-energy approximations and generalized belief propagation algorithms.](https://www.cs.princeton.edu/courses/archive/spring06/cos598C/papers/YedidaFreemanWeiss2004.pdf)" *IEEE Transactions on information theory* 51.7 (2005): 2282-2312.

[3] Cowell, Robert. "[Advanced inference in Bayesian networks.](https://link.springer.com/chapter/10.1007/978-94-011-5014-9_2)" *Learning in graphical models*. Springer, Dordrecht, 1998. 27-49.

[4] Heskes, Tom. "[On the uniqueness of loopy belief propagation fixed points.](https://citeseerx.ist.psu.edu/document?repid=rep1&type=pdf&doi=98a0dd903f0eea9fc35bd5f163563252c118ca60)" *Neural Computation* 16.11 (2004): 2379-2413.

[5] Yuille, Alan. "[Belief propagation, mean-field, and bethe approximations.](https://www.cs.jhu.edu/~ayuille/JHUcourses/VisionAsBayesianInference2020/10/BeliefPropagationMFT.pdf)" *Dept. Stat., Univ. California, Los Angeles, Los Angeles, CA, USA, Tech. Rep* 4 (2010).

[6] Yedidia, Jonathan S., William T. Freeman, and Yair Weiss. "[Bethe free energy, Kikuchi approximations, and belief propagation algorithms.](https://merl.com/publications/docs/TR2001-16.pdf)" *Advances in neural information processing systems* 13 (2001): 689.

[7] Eric P. Xing "[Variational Inference: Loopy Belief Propagation](https://www.cs.cmu.edu/~epxing/Class/10708-17/notes-17/10708-scribe-lecture12.pdf)" Probabilistic Graphical Models (10-708, Spring 2017) Lecture notes, School of Computer Science, Carnegie Mellon University

[8] Ricci-Tersenghi, Federico. "[The Bethe approximation for solving the inverse Ising problem: a comparison with other inference methods.](https://arxiv.org/abs/1112.4814)" *Journal of Statistical Mechanics: Theory and Experiment* 2012.08 (2012): P08015.

```