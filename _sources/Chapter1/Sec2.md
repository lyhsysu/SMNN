# Sec 2. 从空腔方法到消息传递算法

<p style="text-align:center"><img src="https://user-images.githubusercontent.com/106574511/211763400-a7df0856-317b-4693-818a-320fb64c0c0f.png" alt="" class="bg-primary" width="400px"></p>

我们首先定义一个**空腔磁化强度**（cavity magnetization）$m_{i\to a}$，即变量节点 $i$ 和功能节点 $a$ 之间的连接不存在时，节点 $i$ 的磁化强度

$$
m_{i\to a}=\sum_{x_i}x_iP_{i\to a}(x_i)=\sum_{x_i}x_i\frac{\exp\left[-\beta H_{i\to a}(x_i)\right]}{\sum_{x_i}\exp\left[-\beta H_{i\to a}(x_i)\right]}
$$ (a10)

其中，$H_{i \rightarrow a}$ 是指把节点 $a$ 拿掉后节点 $i$ 的哈密顿量，

$$
H_{i \rightarrow a}=H_{\text {cavity } }-\sum_{b \in \partial i \backslash a} J_b \sigma_i \prod_{j \in \partial b \backslash i} \sigma_j
$$ (a11)

$H_{\text {cavity } }$是指空腔内所有节点作为一个系统的哈密顿量，它包含了上图中 $b$、$c$ 等节点与空腔外界的连接，即用红线围起来的边界部分 $\mathcal{B}$。求 $H_{i \rightarrow a}$ 时要从 $H_{\text {cavity } }$ 中减去空腔与 $\mathcal{B}$ 连接的部分。

我们通过空腔近似对式 {eq}`a10` 做一个trick，化简得到

$$
m_{i \rightarrow a}=\frac{\prod_{b \in \partial i \backslash a} \Lambda_{b \rightarrow i}^{+}-\prod_{b \in \partial i \backslash a} \Lambda_{b \rightarrow i}^{-} }{\prod_{b \in \partial i \backslash a} \Lambda_{b \rightarrow i}^{+}+\prod_{b \in \partial i \backslash a} \Lambda_{b \rightarrow i}^{-} }
$$ (a12)

其中，$\Lambda^{\pm}_{b\to i}$ 仍由式 {eq}`a7` 定义。

```{admonition} 推导细节 f
:class: dropdown
首先将式 {eq}`a10` 改写为

$$
m_{i\to a}=\frac{\sum_{\sigma}\sigma_i\exp{-\beta H_{i\to a}(\boldsymbol{\sigma})} }{\sum_{\sigma}\exp{-\beta H_{i\to a}(\boldsymbol{\sigma})} }=\frac{\frac{\sum_\sigma \sigma_i \exp \left(-\beta H_{i \rightarrow a}(\sigma)\right)}{Z_{\text {cavity } } } }{\frac{\sum_\sigma \exp \left(-\beta H_{i \rightarrow a}(\sigma)\right)}{Z_{\text {cavity } } } }
$$ 

考虑分母部分

$$
\begin{aligned}
	&\frac{\sum_\sigma \exp \left(-\beta H_{i \rightarrow a}(\sigma)\right)}{Z_{\text {cavity } } } \\
	&= \frac{\sum_\sigma \exp \left[-\beta \left( H_{\text {cavity } }-\sum_{b \in \partial i \backslash a} J_b \sigma_i \prod_{j \in \partial b \backslash i} \sigma_j \right) \right]}{Z_{\text {cavity } } }\\
	&=\frac{\sum_\sigma \exp(-\beta H_{\text {cavity } })  \exp \left(\beta \sum_{b \in \partial i \backslash a} J_b \sigma_i \prod_{j \in \partial b \backslash i} \sigma_j \right)}{Z_{\text {cavity } } }\\
	&=\sum_\sigma\frac{\exp(-\beta H_{\text {cavity } })}{Z_{\text {cavity } } }\exp \left(\beta \sum_{b \in \partial i \backslash a} J_b \sigma_i \prod_{j \in \partial b \backslash i} \sigma_j \right)
\end{aligned}
$$

注意到 $\frac{\exp(-\beta H_{\text {cavity } })}{Z_{\text {cavity } } }$ 正是空腔的联合概率分布，根据空腔假设，其等于空腔的边缘概率分布

$$
\frac{\exp(-\beta H_{\text {cavity } })}{Z_{\text {cavity } } } = P_{\text {cavity } }(\mathcal{B})
$$

因此分母部分改写为

$$
\sum_{\sigma}P_{\text {cavity } }(\mathcal{B}) \exp \left(\beta \sum_{b \in \partial i \backslash a}  J_b \sigma_i \prod_{j \in \partial b \backslash i} \sigma_j\right)
$$

将对自旋变量的求和 $\sum_{\sigma}$ 分为对节点 $i$ 和对边界 $\mathcal{B}$ 两部分，并将空腔近似

$$
P_{\text {cavity } }(\mathcal{B}) \approx \prod_{b \in \partial i \backslash a} \prod_{j \in \partial b \backslash i} q_{j \rightarrow b}\left(\sigma_j\right)=\prod_{b \in \partial i \backslash a} \prod_{j \in \partial b \backslash i} \frac{1+\sigma_j m_{j\to b} }{2}
$$

代入，得

$$
\begin{aligned}
	&\sum_{\sigma}P_{\text {cavity } }(\mathcal{B}) \exp \left[ \beta\sum_{b \in \partial i \backslash a}  J_b \left( \sigma_i \prod_{j \in \partial b \backslash i}\right) \sigma_j\right]\\
	&=\sum_{\sigma_i} \sum_{\mathcal{B} } P_{\text {cavity } }(\mathcal{B}) \exp \left[ \beta\sum_{b \in \partial i \backslash a}  J_b \left( \sigma_i \prod_{j \in \partial b \backslash i}\right) \sigma_j\right]\\
	&=\quad\sum_{\sigma_i} \sum_{\mathcal{B} } \prod_{b \in \partial i \backslash a} \prod_{j \in \partial b \backslash i} \frac{1+\sigma_j m_{j\to b} }{2} \exp \left[ \beta\sum_{b \in \partial i \backslash a}  J_b \left( \sigma_i \prod_{j \in \partial b \backslash i}\right) \sigma_j\right]\\
	&=\sum_{\sigma_i}\prod_{b \in \partial i \backslash a}\sum_{\mathcal{B} }\prod_{j \in \partial b \backslash i}\frac{1+\sigma_j m_{j\to b} }{2} \exp \left[ \beta\sum_{b \in \partial i \backslash a}  J_b \left( \sigma_i \prod_{j \in \partial b \backslash i}\right) \sigma_j\right]
\end{aligned}
$$

这与式 {eq}`Z6` 完全相同，经过同样的 trick 改写为

$$
\prod_{b \in \partial i \backslash a} \Lambda_{b \rightarrow i}^{+}+\prod_{b \in \partial i \backslash a} \Lambda_{b \rightarrow i}^{-} 
$$

对分子做一个类似（但有一点不同的）trick，可以改写为

$$
\prod_{b \in \partial i \backslash a} \Lambda_{b \rightarrow i}^{+}-\prod_{b \in \partial i \backslash a} \Lambda_{b \rightarrow i}^{-}
$$

<font color="red">具体怎么操作待补充</font>

由此，我们得到了式 {eq}`a12` 

```

``````{margin} 
```{admonition} 注 <font color="blue">6</font>
空腔磁化强度 $m_{i\to a}$ 的含义是将节点 $a$ 移除后，节点（自旋）$i$ 的磁化强度；共轭空腔磁化强度 $\hat{m}_{a\to i}$ 的含义是节点 $i$ 只与节点 $a$ 相互作用时，节点（自旋）$i$ 的磁化强度。
```
关于共轭空腔磁化强度的含义，可以通过以下操作进行理解：

将一个变量节点 $i$ 和与其邻近的一个功能节点 $a$ 一起加入到一个空腔系统中，然后计算在新系统中节点 $i$ 的磁化强度（与 **推导细节 f** 类似的操作），最终可以得到

$$
m_i = \frac{\Lambda^+_{a\to i}-\Lambda^-_{a\to i}}{\Lambda^+_{a\to i}+\Lambda^-_{a\to i}}= \hat{m}_{a\to i}
$$
``````

定义空腔磁化强度 $m_{j\to b}$ 的**共轭空腔磁化强度**<sup><font color="blue">6</font></sup>（conjugate cavity magnetization）为

$$
\hat{m}_{b \rightarrow j} \equiv \tanh \left(\beta J_b\right) \prod_{j \in \partial b \backslash i} m_{j \rightarrow b}
$$ (a13)

则式 {eq}`a7` 可以改写为

$$
\Lambda^\pm_{b\to i} =\cosh(\beta J_b)\left( 1\pm \hat{m}_{b \rightarrow j} \right)
$$ (a14)

因此式 {eq}`a12` 改写为

$$
m_{i \rightarrow a}=\frac{\prod_{b \in \partial i \backslash a}\left(1+\hat{m}_{b \rightarrow i}\right)-\prod_{b \in \partial i \backslash a}\left(1-\hat{m}_{b \rightarrow i}\right)}{\prod_{b \in \partial i \backslash a}\left(1+\hat{m}_{b \rightarrow i}\right)+\prod_{b \in \partial i \backslash a}\left(1-\hat{m}_{b \rightarrow i}\right)}
$$ (a15)

式 {eq}`a13` 与式 {eq}`a15` 的组成一组迭代方程。

与单个节点的空腔概率分布 $q_{i \rightarrow a}$ 对应，我们定义 $p_{a \to i}$ 为节点 $i$ 只与节点 $a$ 相互作用时，节点 $i$ 的概率分布。并且，与式 {eq}`a3` 相同，我们可以用 $\hat{m}_{a \rightarrow i}$ 将 $p_{a \to i}$ 参数化

$$
p_{a \rightarrow i}\left(\sigma_i\right) =\frac{1+\sigma_i \hat{m}_{a \rightarrow i} }{2}
$$ (a16)

接下来我们再用两个局域场来参数化 $q_{i \rightarrow a}$ 和 $p_{a \to i}$。定义“空腔场”$h_{i\to a}$ 为**将节点 $a$ 移除后节点 $i$ 处受到的局域场**，定义其“共轭空腔场”$u_{a\to i}$ 为**节点 $i$ 只与节点 $a$ 相互作用时，节点 $i$ 受到的局域场**。

考虑一个一般的局域场 $h_i$ 对自旋 $\sigma_i$ 作用时，有

$$
H(\sigma_i)=h_i\sigma_i,\quad Z=\sum_{\sigma_i=\pm 1}e^{\beta h_i\sigma_i}=2\cosh(\beta h_i),\quad P(\sigma_i)=\frac{e^{\beta h_i\sigma_i}}{Z}=\frac{e^{\beta h_i\sigma_i}}{2\cosh(\beta h_i)}
$$ (a16_5)

同样，可以写出 $h_{i\to a}$ 和 $u_{a\to i}$ 对自旋 $\sigma_i$（即节点 $i$）作用时的 $q_{i \rightarrow a}$ 和 $p_{a \to i}$：

$$
\begin{aligned}
	q_{i \rightarrow a}\left(\sigma_i\right) & \equiv \frac{\exp \left(\beta h_{i \rightarrow a} \sigma_i\right)}{2 \cosh \beta h_{i \rightarrow a} }\\
	p_{a \rightarrow i}\left(\sigma_i\right) & \equiv \frac{\exp \left(\beta u_{a \rightarrow i} \sigma_i\right)}{2 \cosh \beta u_{a \rightarrow i} }
\end{aligned}
$$ (a17)

结合式 {eq}`a3` {eq}`a16` 可以得到

$$
\begin{aligned}
m_{i \rightarrow a} &=\tanh \beta h_{i \rightarrow a} \\
\hat{m}_{a \rightarrow i} &=\tanh \beta u_{a \rightarrow i}
\end{aligned}
$$ (a18)

再结合式 {eq}`a13` {eq}`a15` 可以得到

$$
\begin{aligned}
&h_{i \rightarrow a}=\frac{1}{\beta}\left(\sum_{b \in \partial i \backslash a} \beta u_{b \rightarrow i}\right) \\
&u_{a \rightarrow i}=\frac{1}{\beta} \tanh ^{-1}\left[\tanh \left(\beta J_a\right) \prod_{j \in \partial a \backslash i} \tanh \left(\beta h_{j \rightarrow a}\right)\right]
\end{aligned}
$$ (a19)

```{admonition} 推导细节 g
:class: dropdown

考虑

$$
q_{i \rightarrow a}\left(\sigma_i\right) \equiv \frac{\exp \left(\beta h_{i \rightarrow a} \sigma_i\right)}{2 \cosh \beta h_{i \rightarrow a} }= \frac{1+\sigma_i m_{i\to a} }{2}
$$

将其变形为

$$
m_{i\to a} = \frac{\exp \left(\beta h_{i \rightarrow a} \sigma_i\right)/\sigma_i}{\cosh \beta h_{i \rightarrow a}}-1 = \frac{\left(2e^{\beta h_{i \rightarrow a}}\sigma_i -e^{\beta h_{i \rightarrow a}} - e^{-\beta h_{i \rightarrow a}}\right) /\sigma }{e^{\beta h_{i \rightarrow a}} + e^{-\beta h_{i \rightarrow a}}}
$$

考虑分子部分 $\frac{2e^{\beta h_{i \rightarrow a}}\sigma_i -e^{\beta h_{i \rightarrow a}} - e^{-\beta h_{i \rightarrow a}}}{\sigma_i}$

当 $\sigma_i=1$ 时：

$$
numerator = \frac{2e^{\beta h_{i \rightarrow a}} -e^{\beta h_{i \rightarrow a}} - e^{-\beta h_{i \rightarrow a}}}{1} = e^{\beta h_{i \rightarrow a}} - e^{-\beta h_{i \rightarrow a}}
$$

当  $\sigma_i=-1$ 时：

$$
numerator = \frac{2e^{-\beta h_{i \rightarrow a}} -e^{\beta h_{i \rightarrow a}} - e^{-\beta h_{i \rightarrow a}}}{-1} = e^{\beta h_{i \rightarrow a}} - e^{-\beta h_{i \rightarrow a}}
$$

即无论 $\sigma_i$ 取何值，恒有

$$
\frac{2e^{\beta h_{i \rightarrow a}}\sigma_i -e^{\beta h_{i \rightarrow a}} - e^{-\beta h_{i \rightarrow a}}}{\sigma_i}\equiv e^{\beta h_{i \rightarrow a}} - e^{-\beta h_{i \rightarrow a}}
$$

因此

$$
\hat{m}_{a \rightarrow i} =\frac{e^{\beta h_{i \rightarrow a}} - e^{-\beta h_{i \rightarrow a}}}{e^{\beta h_{i \rightarrow a}} + e^{-\beta h_{i \rightarrow a}}}=\frac{\sinh(\beta h_{i \rightarrow a})}{\cosh(\beta h_{i \rightarrow a})}=\tanh(\beta h_{i \rightarrow a})
$$

同理，可得

$$
m_{i\to a} =\tanh \beta u_{a \rightarrow i}
$$

这样就得到了式 {eq}`a18`，下面考虑如何得到式 {eq}`a19`。

根据 $\hat{m}_{b \rightarrow j} = \tanh \left(\beta J_b\right) \prod_{j \in \partial b \backslash i} m_{j \rightarrow b}$ 和 $\hat{m}_{b \rightarrow j} =\tanh(\beta h_{j \rightarrow b})$，有

$$
\tanh \beta u_{b \rightarrow j} =\tanh \left(\beta J_b\right) \prod_{j \in \partial b \backslash i} m_{j \rightarrow b}
$$

再将 $m_{j\to b} =\tanh \beta u_{b \rightarrow j}$ 代入，有

$$
\tanh \beta u_{b \rightarrow j} =\tanh \left(\beta J_b\right) \prod_{j \in \partial b \backslash i} \tanh(\beta h_{j \rightarrow b})
$$

因此得到（为了与式 {eq}`a19` 表示一致，已将 $b$ 用 $a$ 替换）

$$
u_{a \rightarrow i}=\frac{1}{\beta} \tanh ^{-1}\left[\tanh \left(\beta J_a\right) \prod_{j \in \partial a \backslash i} \tanh \left(\beta h_{j \rightarrow a}\right)\right]
$$

**那么问题来了，$h_{i\to a}$ 的表达式是怎么推出来的呢**
```

式 {eq}`a19` 是一个迭代方程，通过迭代求解 $u_{b \rightarrow i}$ 后，即可通过下式求解 $m_i$

$$
m_i=\tanh \left(\sum_{b \in \partial i} \beta u_{b \rightarrow i}\right)
$$ (a20)

```{tip}
式 {eq}`a20` 的含义：

$\hat{m}_{a \rightarrow i} =\tanh \beta u_{a \rightarrow i}$ 表示节点 $i$ 只与节点 $a$ 相互作用时节点 $i$ 的磁化强度。那么将 $a$ 对 $i$ 周围所有功能节点 $b\in\partial i$ 求和，即表示节点 $i$ 与周围所有节点都相互作用时的磁化强度。这不再是一个空腔的磁化强度，是没有移除任何节点的完整系统中的磁化强度。
```