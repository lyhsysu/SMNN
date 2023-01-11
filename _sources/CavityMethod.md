# 空腔方法与消息传递方程

## Pre-因子图

## 一、空腔方法求解配分函数

``````{margin} 
```{admonition} 注 <font color="blue">1</font>
这是因为我们之后将用本章得到的结论对 Sourlas 码进行一个编程模拟。
```
``````

在本节，我们以无限范围模型（Infinite-range model）为例<sup><font color="blue">1</font></sup>介绍如何使用空腔方法求解配分函数。无限范围模型是对 S-K 模型的推广，S-K 模型只考虑二体相互作用，无限范围模型则考虑了多个自旋的相互作用，其哈密顿量为

$$
H(\boldsymbol{\sigma})=-\sum_{a=1}^MJ_a\prod_{i\in\partial a}\sigma_i
$$ (H0)

式 {eq}`H0` 是一个连乘的形式，因此可以用因子图表示为

<p style="text-align:center"><img src="https://user-images.githubusercontent.com/106574511/211736226-918e1c52-2392-4e1a-9ef8-bb247015186d.png" alt="" class="bg-primary" width="400px"></p>

这是一个稀疏连接的网络，可以展开成树状图

<p style="text-align:center"><img src="https://user-images.githubusercontent.com/106574511/211737170-439e7f40-c215-4ef2-a4cf-82a4b1cf9045.png" alt="" class="bg-primary" width="400px"></p>

系统的配分函数为

$$
Z=\mathrm{Tr~} \exp(-\beta H(\boldsymbol{\sigma}))
$$ (Z0)

从物理的角度考虑，配分函数 $Z$ 不是一个广延量，我们更关注的是作为广延量的自由能 $F=-\frac{1}{\beta}\ln Z$。空腔方法的核心是利用因子图将整个网络划分为功能节点和变量节点，分别求出加上一个功能节点和加上一个变量节点（先将这个节点拿掉，再加上，称为“空腔”）时自由能的增量，然后再将各个节点的自由能累加起来，得到总的自由能。接下来三个小节分别是求解功能节点的自由能增量、变量节点的自由能增量、累加成为总的自由能。

### 1.1 功能节点的自由能增量

我们首先考虑将一个功能节点 $a$ 拿走，剩下的网络构成 $\text{old}$ 网络，然后再将 $a$ 加回去，作为 $\text{new}$ 网络。新旧网络的自由能差值即这个功能节点的自由能。新网络的哈密顿量可以写为

$$
H^{\text{new} } = H^{\text{old} } - J_a\prod_{i\in \partial a} \sigma_i
$$ (H1)

因此配分函数可以写成：

$$
Z^{\text{new} } =Z^{\text{old} } \sum_{\{\sigma_i | i \in \partial a \} }\left[  \exp(\beta J_a \prod_{i \in \partial a} \sigma_i ) \sum_{\{\sigma_i | i \notin \partial a \} } \frac{\exp(-\beta H^{\text{old} })}{Z^{\text{old} } } \right]
$$ (Z1)

``````{admonition} 推导细节 a
:class: dropdown
$$
\begin{aligned}
 &Z^{\text{new} } = \sum_{\{\sigma_i\} }\exp \left(-\beta H^{\text{old} }+\beta J_a \prod_{i \in \partial a} \sigma_i \right)\\
 &=Z^{\text{old} } \sum_{\{\sigma_i\} }\left[\frac{\exp(-\beta H^{\text{old} })}{Z^{\text{old} } } \exp(\beta J_a \prod_{i \in \partial a} \sigma_i ) \right] \\
 &=Z^{\text{old} } \sum_{\{\sigma_i | i \in \partial a \} }\sum_{\{\sigma_i | i \notin \partial a \} }\left[\frac{\exp(-\beta H^{\text{old} })}{Z^{\text{old} } } \exp(\beta J_a \prod_{i \in \partial a} \sigma_i ) \right] \\
 &=Z^{\text{old} } \sum_{\{\sigma_i | i \in \partial a \} }\left[  \exp(\beta J_a \prod_{i \in \partial a} \sigma_i ) \sum_{\{\sigma_i | i \notin \partial a \} } \frac{\exp(-\beta H^{\text{old} })}{Z^{\text{old} } } \right]
\end{aligned}
$$

```{attention}
这里将求和号拆开是

$$
 \sum_{\{\sigma_i\} } = \sum_{\{\sigma_i | i \in \partial a \} }\sum_{\{\sigma_i | i \notin \partial a \} }
$$

而非使用加号（+）

$$
 \sum_{\{\sigma_i\} } = \sum_{\{\sigma_i | i \in \partial a \} }+\sum_{\{\sigma_i | i \notin \partial a \} }
$$

这里的“拆开”主要是表示求和顺序的不同，即先对内部求和号求和，再对外部求和号求和，由于二者的指标中没有交叉部分，这实际上是等效于对总体求和的。
```
``````

考虑式 {eq}`Z1` 后面一项 $\sum_{\{\sigma_i | i \notin \partial a \} } \frac{\exp(-\beta H^{\text{old} })}{Z^{\text{old} } } $ 的物理意义。如果求和号中 $\sigma_i$ 的指标 $i$ 是对旧网络中的所有节点，即 $\{\sigma_i|i\notin a\}$，那么这一项的意义是整个旧网络所有节点的联合概率分布。现在求和指标中去掉了边界上的点 $\{\sigma_i|i\in \partial a\}$，其表示的是边界点的边缘概率分布。

**空腔近似（cavity aproximation）**

在这里我们要做一个“空腔假设”，把**空腔的边缘概率分布认为是空腔的概率分布**，即

$$
P_{\text{Cavity} }(\{\sigma_i|i\in\partial a \}) = \sum_{\{\sigma_i | i \notin \partial a \} } \frac{\exp(-\beta H^{\text{old} })}{Z^{\text{old} } } 
$$ (a1)

并且认为**边界上各个节点是相互独立**的。因此近似有

$$
P_{\text{Cavity} }(\{\sigma_i|i\in\partial a \}) \approx \prod_{i\in\partial a} q_{i\to a}(\sigma_i)
$$ (a2)

其中，$q_{i\to a}(\sigma_i)$ 是一个节点的空腔概率分布，我们用 $m_{i \to a}$ 将其参数化为

$$
q_{i\to a}(\sigma) = \frac{1+\sigma_i m_{i\to a} }{2}
$$ (a3)

```{admonition} 推导细节 b
:class: dropdown
对于一个节点来说，其概率分布为 $q_{i\to a}(\sigma)$，因此我们写出（拿掉节点 $a$ 后）节点 $i$ 上自旋变量取值的期望，即物理上的磁化强度 $m$

$$
m_{i\to a}=(+1)q_{i\to a}(\sigma_i = +1) + (-1)q_{i\to a}(\sigma_i = -1)
$$

根据上式，可以反解出

$$
q_{i\to a}(\sigma) = \frac{1+\sigma_i m_{i\to a} }{2}
$$

或者说，只有当 $q_{i\to a}(\sigma)$ 是这种形式时，才满足磁化强度的定义。可以将

$$
q_{i\to a}(\sigma=+1) = \frac{1+m_{i\to a} }{2},\quad q_{i\to a}(\sigma=-1) = \frac{1-m_{i\to a} }{2}
$$

代入磁化强度的定义式验证

$$
m_{i\to a}=(+1)\frac{1+m_{i\to a} }{2} + (-1)\frac{1-m_{i\to a} }{2} = \frac{1+m_{i\to a} }{2}+\frac{m_{i\to a}-1}{2}=m_{i\to a}
$$

这是两条自洽的方程
```

因此，由空腔近似得到

$$
\sum_{\{\sigma_i | i \notin \partial a \} } \frac{\exp(-\beta H^{\text{old} })}{Z^{\text{old} } } \approx \prod_{i\in\partial a}\frac{1+\sigma_i m_{i\to a} }{2}
$$ (a4)

代回式 {eq}`Z1` 得

$$
Z^{\text{new} } =Z^{\text{old} } \sum_{\{\sigma_i | i \in \partial a \} }\left[  \exp(\beta J_a \prod_{i \in \partial a} \sigma_i ) \prod_{i\in\partial a}\frac{1+\sigma_i m_{i\to a} }{2}\right]
$$ (Z2)

``````{margin} 
```{admonition} 注 <font color="blue">2</font>
请注意，从式 {eq}`Z2` 到式 {eq}`Z3` 是精确的，没有做任何近似。
```
``````

通过一个 trick（见 **推导细节 c**），可以将其化简为<sup><font color="blue">2</font></sup>

$$
Z^{\text{new} }=Z^{\text{old} }\cosh(\beta J_a)\left(1+\tanh(\beta J_a)\prod_{i\in\partial a}m_{i\to a} \right)
$$ (Z3)

```{admonition} 推导细节 c
:class: dropdown
考虑将 {eq}`Z2` 中的 $\prod_{i\in\partial a}\frac{1+\sigma_i m_{i\to a} }{2}$ 展开处理：

$$
\begin{aligned}
	\prod_{i\in\partial a}\frac{1+\sigma_i m_{i\to a} }{2} &=\frac{1}{2^N}(1+\sigma_1 m_{1\to a})(1+\sigma_2 m_{2\to a})\cdots (1+\sigma_N m_{N\to a})\\
	&=\frac{1}{2^N}( 1 + \sigma_1 m_{1\to a} + \sigma_2 m_{2\to a} + \cdots \\
	&\;\,\,\qquad\quad + \sigma_1\sigma_2 m_{1\to a}m_{2\to a} + \sigma_1\sigma_3 m_{1\to a}m_{3\to a} + \cdots \\
	&\;\,\,\qquad\quad + \cdots \\
	&\;\,\,\qquad\quad + \prod_{i\in\partial a} \sigma_i \prod_{i\in\partial a} m_{i\to a}) 
\end{aligned}
$$

在大括号外要对 $\{\sigma_i\}$ 做求和，可以确定 $\prod_{i\in\partial a}\sigma_i = \pm 1$。但是对于 $\sigma_1\sigma_2$，$\sigma_2\sigma_4\sigma_5$ 等中间项，可以写作 $\prod_{i\in D}\sigma_i= \prod_{i\in \partial a}\sigma_i / \prod_{i\in \partial a,i \notin D}\sigma_i $，这里 $\prod_{i\in \partial a,i \notin D}\sigma_i$ 的取号是不确定的（有可能是 $1$ 也有可能是 $-1$ ），因此对这个展开式做求和时只有第一项和最后一项会确定下来，其他项都会因为有 $+1$ 和 $-1$ 两种可能值而抵消掉。因此有

$$
\prod_{i\in\partial a}\frac{1+\sigma_i m_{i\to a} }{2} =\frac{1}{2}(1+ \prod_{i\in\partial a} \sigma_i \prod_{i\in\partial a} m_{i\to a})=\frac{1}{2}(1+ \prod_{i\in\partial a} \sigma_i m_{i\to a})
$$

进而得到

$$
\begin{aligned}
	Z^{\text{new} } &= Z^{\text{old} } \sum_{\{ {\prod_{i\in\partial a}\sigma_i} = \pm 1 \} }\left[  \exp(\beta J_a \prod_{i \in \partial a} \sigma_i )\frac{1}{2} (1+ \prod_{i\in\partial a} \sigma_i m_{i\to a})\right]\\
	&=Z^{\text{old} }\left[ \exp(\beta J_a)\frac{1}{2} (1+ \prod_{i\in\partial a}m_{i\to a})+\exp(-\beta J_a)\frac{1}{2} (1- \prod_{i\in\partial a}m_{i\to a})\right]\\
	&=Z^{\text{old} } \left[\frac{1}{2}[\exp(\beta J_a)+\exp(-\beta J_a)]+ \frac{1}{2}[\exp(\beta J_a)-\exp(-\beta J_a)]\prod_{i\in\partial a}m_{i\to a} \right]\\
	&=Z^{\text{old} }\left[ \cosh(\beta J_a)+\sinh(\beta J_a)\prod_{i\in\partial a}m_{i\to a} \right]  \\
	&=Z^{\text{old} }\cosh(\beta J_a)\left(1+\tanh(\beta J_a)\prod_{i\in\partial a}m_{i\to a} \right)
\end{aligned}
$$
```

因此加上一个功能节点后自由能的增量为

$$
\Delta F_a=-\frac{1}{\beta}\ln \frac{Z^{\text{new} } }{Z^{\text{old} } }=\ln \left[\cosh(\beta J_a)\left(1+\tanh(\beta J_a)\prod_{i\in\partial a}m_{i\to a} \right)\right]
$$ (a5)

### 1.2 变量节点的自由能增量

对于变量节点，我们采取与上一小节同样的方法，但是在添加变量节点 $i$ 的同时要将与其邻近的功能节点 $\{b|b\in \partial i\}$ 一起添加进网络中。配分函数的变化写为

$$
Z^{\text{new} } = \sum_{\sigma_i^{\text{old} } }\sum_{\sigma_i}\exp \left(-\beta H^{\text{old} }+\beta \sum_{b\in \partial i} J_b \prod_{j \in \partial b} \sigma_j \right) 
$$ (Z4)

``````{margin} 
```{admonition} 注 <font color="blue">3</font>
为了书写方便，以下写到 $j \in \partial b\backslash i$ 时默认 $b \in \partial i$。
```
``````

这里的空腔就包括了变量节点 $i$ 和与其邻近的功能节点 $\{b|b\in \partial i\}$。将式 {eq}`Z4` 改写为<sup><font color="blue">3</font></sup>

$$
\begin{aligned}
Z^{\text{new} } = &Z^{\text{old} }\sum_{\{\sigma_j | j \in \partial b\backslash i \} } \sum_{\sigma_i}\sum_{\{\sigma_j | j \notin \partial b\backslash i \} }\frac{\exp (-\beta H^{\text{old} })}{Z^{\text{old} } } \exp \left[\beta \sum_{b\in \partial i} J_b \left(\sigma_i \prod_{j \in \partial b\backslash i} \sigma_j\right) \right]
\end{aligned}
$$ (Z5)

```{admonition} 推导细节 d
:class: dropdown
首先将式 {eq}`Z4` 中的空腔分离出来，即，将 $\prod_{j \in \partial b}$ 写成 $\sigma_i \prod_{j \in \partial b\backslash i} \sigma_j$。然后将旧网络分离出来（与 **推导细节a** 相同的操作）

$$
\begin{aligned}
	Z^{\text{new} } &= \sum_{\sigma_i^{\text{old} } }\sum_{\sigma_i}\exp \left(-\beta H^{\text{old} }+\beta \sum_{b\in \partial i} J_b \prod_{j \in \partial b} \sigma_j \right) \\
	& = \sum_{\sigma_i^{\text{old} } }\sum_{\sigma_i}\exp \left[-\beta H^{\text{old} }+\beta \sum_{b\in \partial i} J_b \left(\sigma_i \prod_{j \in \partial b\backslash i} \sigma_j\right) \right]\\
	&=Z^{\text{old} }\sum_{\sigma_i^{\text{old} } }\sum_{\sigma_i}\frac{\exp (-\beta H^{\text{old} })}{Z^{\text{old} } } \exp \left[\beta \sum_{b\in \partial i} J_b \left(\sigma_i \prod_{j \in \partial b\backslash i} \sigma_j\right) \right] 
\end{aligned}
$$

将第一个求和号 $\sum_{\sigma_i^{\text{old} } }$ 中的指标分成属于空腔边缘的节点和其他节点两部分，即

$$
\sum_{\sigma_i^{\text{old} } } = \sum_{\{\sigma_j | j \in \partial b\backslash i\} } \sum_{\{\sigma_j | j \notin \partial b\backslash i \} }
$$

因此

$$
Z^{\text{new} } = Z^{\text{old} }\sum_{\{\sigma_j | j \in \partial b\backslash i \} } \sum_{\sigma_i}\sum_{\{\sigma_j | j \notin \partial b\backslash i \} }\frac{\exp (-\beta H^{\text{old} })}{Z^{\text{old} } } \exp \left[\beta \sum_{b\in \partial i} J_b \left(\sigma_i \prod_{j \in \partial b\backslash i} \sigma_j\right) \right]
$$

考虑到 $\sum_{\{\sigma_j | j \notin \partial b\backslash i\} }$ 这个求和在上式中只对 $\frac{\exp (-\beta H^{\text{old} })}{Z^{\text{old} } }$ 起作用，将其放到 $\sum_{\sigma_i}$ 后面，即可得到式 {eq}`Z5`。
```

``````{margin} 
```{admonition} 注 <font color="blue">4</font>
注意这里多了一个指标 $\prod_{j \in \partial b\backslash i}$，是因为……
```
``````

式 {eq}`Z5` 的 $\sum_{\{\sigma_j | j \notin \partial b\backslash i \} } \frac{\exp(-\beta H^{\text{old} })}{Z^{\text{old} } }$ 部分表示空腔的边缘概率分布，采用与式 {eq}`a1` 相同的空腔近似，即<sup><font color="blue">4</font></sup>

$$
\sum_{\{\sigma_j | j \notin \partial b\backslash i \} } \frac{\exp(-\beta H^{\text{old} })}{Z^{\text{old} } } = P_{\text{cavity} }(\{\sigma_j | j \in \partial b\backslash i \}) \approx \prod_{b \in \partial i}\prod_{j \in \partial b\backslash i} \frac{1+\sigma_jm_{j \to b} }{2}
$$ (a6)

可得

$$
Z^{\text{new} }\approx Z^{\text{old} } \sum_{\{\sigma_j | j \in \partial b\backslash i\} } \sum_{\sigma_i} \prod_{b \in \partial i}\prod_{j \in \partial b\backslash i} q_{j\to b}(\sigma_j) \exp \left[\beta \sum_{b\in \partial i} J_b \left(\sigma_i \prod_{j \in \partial b\backslash i} \sigma_j\right) \right]
$$ (Z6)

再进行化简，得

$$
Z^{\text{new} } \approx Z^{\text{old} } \left( \prod_{b \in \partial i} \Lambda^+_{b\to i} + \prod_{b \in \partial i} \Lambda^-_{b\to i}  \right)
$$ (Z7)

其中，

$$
\Lambda^\pm_{b\to i} \equiv\cosh(\beta J_b)\left( 1\pm \tanh(\beta J_b) \prod_{j \in \partial b\backslash i} m_{j\to b} \right)
$$ (a7)

```{admonition} 推导细节 e
:class: dropdown
待补充
```

自由能的变化量为

$$ 
\Delta F_i = -\frac{1}{\beta}\ln \frac{Z^{\text{new} } }{Z^{\text{old} } }=\ln \left[ \prod_{b \in \partial i} \Lambda^+_{b\to i} + \prod_{b \in \partial i} \Lambda^-_{b\to i} \right]
$$ (a8)

### 1.3 总的自由能

从上面两个小节，我们已经得到了加上一个功能节点时自由能的增量【式 {eq}`a5`】和加上一个变量节点时自由能的增量【式 {eq}`a8`】

求系统的总自由能即相当于求加上各个节点后自由能的增量之和，但是注意在计算变量节点的自由能时把与其临近的功能节点也算入了，因此要减掉重复了 $|\partial a|$ 次的对功能节点的求和，即

$$
F = \sum_{i} \Delta F_i + \sum_{a} \Delta F_a - \sum_{a} |\partial a| \Delta F_a
$$ (a9)

注意到，这里的 $\Lambda^{\pm}_{b\to i}$ 在定义时【式 {eq}`a7`】用到了 $m_{j\to b}$，下一节介绍如何迭代求解。