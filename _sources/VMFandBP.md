# 变分平均场理论与信念传播算法




## 一、变分方法与变分自由能



## 二、平均场近似



## 三、Bethe 近似
平均场近似要求各个格点之间完全相互独立，这只能在温度比较高的时候实现，对于温度较低的情况，必须要考虑格点之间的关联。本节中的 Bethe 近似就引入了格点间的部分关联。

仍然使用如下的分布

$$
P(\boldsymbol{x})=\frac{e^{-\beta E(x)}}{Z} \qquad E(\boldsymbol{x})=-\sum_a J_a \prod_{i \in \partial a} x_i
$$ (b14)

为了更好地在因子图中表示，我们首先进行一个变量替换，令

$$
f_a(\boldsymbol{x}_a)=e^{J_a \prod_{i \in \partial a} x_i}
$$ (b15)

$f_a(\boldsymbol{x}_a)$ 表示一个相互作用 $a$，其中 $\boldsymbol{x}_a$ 表示参与该相互作用的所有自旋变量。同样令 $\beta=1$，因此式 {eq}`b14` 改写为

$$
\begin{aligned}
& P(\boldsymbol{x})=\frac{1}{Z} \prod_a f_a\left(\boldsymbol{x}_a\right) \\
& E(\boldsymbol{x})=-\sum_a \ln f_a\left(\boldsymbol{x}_a\right) 
\end{aligned}
$$ (b16)

显然，此时分布 $P(x)$ 是因子连乘的形式，可以用因子图表示为

<p style="text-align:center"><img src="https://user-images.githubusercontent.com/106574511/212213084-192b03b7-3597-4594-b90c-b406ea370776.png" alt="" class="bg-primary" width="400px"></p>

对这个因子图进行**区域划分**，每一个区域都要将与其中的功能节点相连的所有变量节点都包含在里面。比如，上图中的两个蓝色线圈出来的是区域，但黄色虚线圈出来的不能算作区域，因为与功能节点A相连的变量节点2没有被包含在其中。

对于一个区域 $R$，我们可以用下式计算区域内部的变分内能和变分熵：

$$
\begin{aligned}
& E_R\left(x_R\right)=-\sum_{a \in R} \ln f_a\left(\boldsymbol{x}_a\right) \\
& U_R\left(b_R\right)=\sum_{\boldsymbol{x}_R} b_R\left(\boldsymbol{x}_R\right) E_R\left(\boldsymbol{x}_R\right) \\
& H_R\left(b_R\right)=-\sum_{\boldsymbol{x}_R} b_R\left(\boldsymbol{x}_R\right) \ln b_R\left(\boldsymbol{x}_R\right)
\end{aligned}
$$ (b17)

然后将各个区域的变分内能和变分熵相加，即

$$
\begin{gathered}
U_{\mathcal{R}}=\sum_{R \in \mathcal{R}} C_R U_R\left(b_R\right)\\
H_{\mathcal{R}}=\sum_{R \in \mathcal{R}} C_R H_R\left(b_R\right)
\end{gathered}
$$ (b18)

其中，$\mathcal{R}$ 是所有区域的集合，$C_R$ 是用于对区域进行计数的变量，满足

$$
\begin{gathered}
\sum_{R \in \mathcal{R}} \mathbb{I}[a \in R] C_R=1, \forall a \\
\sum_{R \in \mathbb{R}} \mathbb{I}[i \in R] C_R=1, \forall i
\end{gathered}
$$ (b19)

$\mathbb{I}[a \in R]$ 表示如果功能节点 $a$ 在区域 $R$ 中，则 $C_R=1$，否则 $C_R=0$；同理，$\mathbb{I}[i \in R]$ 表示如果变量节点 $i$ 在区域 $R$ 中，则 $C_R=1$，否则 $C_R=0$。$\mathbb{I}[a \in R]$ 与 $C_R$ 相乘再求和后等于 $1$ 表示节点 $a$ 只被考虑了一次，因此上式对于 $C_R$ 的约束实际上是要求所有的变量节点和功能节点在求和时只被考虑一次，这是就划分区域的要求。

### 3.1 Bethe 近似与信念传播算法

Bethe 近似划分区域的方式如下图所示

<p style="text-align:center"><img src="https://user-images.githubusercontent.com/106574511/212213770-4f299654-11d5-406d-944d-88eb0e3da531.png" alt="" class="bg-primary" width="400px"></p>

划分的区域有两种，一种称作“大区域”，是包含一个功能节点以及所有与其相连的变量节点的区域；另一种称作“小区域”，只包含一个变量节点。大区域的计数设为 $1$，小区域的计数设为 $1-d_i$，其中 $d_i$ 是变量节点 $i$ 参与的功能节点的数量，即

$$
C_{R_{L}}=1,\quad C_{R_{S}}=1-d_i
$$ (b20)

```{tip}
这种划分方式和计数方式很自然地满足了上面对于 $C_R$ 的两条约束：

对于第一条（针对功能节点的）约束，只有大区域中有功能节点参与了求和，并且对于每个功能节点 $a$，当其包含在某个区域 $R$ 中使得 $\mathbb{I}[a \in R]=1$ 时，都有且仅有这一个区域的计数 $C_{R_L}=1$；

对于第二条（针对变量节点的）约束，每个变量节点 $i$ 在其参与的 $d_i$ 个大区域（功能节点）中使得 $\mathbb{I}[a \in R]=1,\,C_{R_{L}}=1$，并且使得其自身所在的小区域 $\mathbb{I}[a \in R]=1,\,C_{R_{S}}=1-d_i$，所有大区域的求和为 $d_i$，只有一个小区域参与了求和为 $1-d_i$，因此总的求和为 $1$。
```

``````{margin} 
```{admonition} 注 <font color="blue">2</font>
在这里指大区域，小区域内部没有相互作用
```
``````

Bethe 近似划分区域后，**只考虑各个区域<sup><font color="blue">2</font></sup>内部的相互作用，不考虑区域与区域间的相互作用**，实际上是一个只考虑**最近邻相互作用**的模型。

根据以上划分的区域，计算变分自由能

$$
\begin{aligned}
&U_{\text {Bethe }}= -\sum_a \sum_{x_a} b_a\left(\boldsymbol{x}_a\right) \ln f_a\left(\boldsymbol{x}_a\right) \\
&H_{\text {Bethe }}=-\sum_a \sum_{x_a} b_a\left(\boldsymbol{x}_a\right) \ln b_a\left(\boldsymbol{x}_a\right)+\sum_i\left(d_i-1\right) \sum_{x_i} b_i\left(x_i\right) \ln b_i\left(x_i\right) \\
&F_{\text {Bethe }}= -\sum_a \sum_{x_a} b_a\left(\boldsymbol{x}_a\right) \ln f_a\left(\boldsymbol{x}_a\right)+\sum_a \sum_{x_a} b_a\left(\boldsymbol{x}_a\right) \ln b_a\left(\boldsymbol{x}_a\right) \\
&\qquad\qquad\, -\sum_i\left(d_i-1\right) \sum_{x_i} b_i\left(x_i\right) \ln b_i\left(x_i\right)
\end{aligned}
$$ (b21)

``````{admonition} 推导细节 d
:class: dropdown
变分内能
$$
\begin{aligned}
U_{\text {Bethe }}= &\sum_{R \in \mathcal{R}} C_R U_R\left(b_R\right)=\sum_{R \in \mathcal{R}} C_R\sum_{\boldsymbol{x}_R} b_R\left(\boldsymbol{x}_R\right) E_R\left(\boldsymbol{x}_R\right)\\
=& \sum_a\sum_{x_a} b_a\left(\boldsymbol{x}_a\right) \left[-\sum_{a \in R} \ln f_a\left(\boldsymbol{x}_a\right)\right] \\
=&-\sum_a \sum_{x_a} b_a\left(\boldsymbol{x}_a\right) \ln f_a\left(\boldsymbol{x}_a\right) \\
\end{aligned}
$$

```{admonition} 注 
这里由于哈密顿量没有外场，实际上只有大区域贡献了内能，因此第一行中的求和下标 $R\in\mathcal{R}$ 有 $a$ 个，并且每个大区域中只有一个功能节点，第二行中的求和下标 $a\in R$ 实际上只对一个 $a$ 求和，$C_R=C_{R_L}=1$。如果哈密顿量是具有外场的形式，那么求Bethe内能时也需要考虑小区域的作用。
```


变分熵（大区域和小区域分别计算后求和）
$$
\begin{aligned}
H_{\text {Bethe }}= &\sum_{R \in \mathcal{R}} C_R H_R\left(b_R\right)=\sum_{R \in \mathcal{R}} C_R \left[-\sum_{\boldsymbol{x}_R} b_R\left(\boldsymbol{x}_R\right) \ln b_R\left(\boldsymbol{x}_R\right)\right]\\
=&-\sum_aC_{R_L}\sum_{\boldsymbol{x}_a} b_a\left(\boldsymbol{x}_a\right) \ln b_a\left(\boldsymbol{x}_a\right)-\sum_iC_{R_S}\sum_{\boldsymbol{x}_i} b_i\left(\boldsymbol{x}_i\right) \ln b_i\left(\boldsymbol{x}_i\right) \\
=& -\sum_a \sum_{x_a} b_a\left(\boldsymbol{x}_a\right) \ln b_a\left(\boldsymbol{x}_a\right)+\sum_i\left(d_i-1\right) \sum_{x_i} b_i\left(x_i\right) \ln b_i\left(x_i\right) 
\end{aligned}
$$

变分自由能
$$
\begin{aligned}
F_{\text {Bethe }}= & \; U_{\text {Bethe }}-H_{\text {Bethe }}\\
= & -\sum_a \sum_{x_a} b_a\left(\boldsymbol{x}_a\right) \ln f_a\left(\boldsymbol{x}_a\right)+\sum_a \sum_{x_a} b_a\left(\boldsymbol{x}_a\right) \ln b_a\left(\boldsymbol{x}_a\right) \\
&-\sum_i\left(d_i-1\right) \sum_{x_i} b_i\left(x_i\right) \ln b_i\left(x_i\right)
\end{aligned}
$$
``````

使用拉格朗日乘子法对变分自由能取极值时存在如下约束：

$$
\begin{aligned}
& \sum_{x_i} b_i\left(x_i\right)=1, \forall i \\
& \sum_{\boldsymbol{x}_a} b_a\left(\boldsymbol{x}_a\right)=1, \forall a \\
& \sum_{\boldsymbol{x}_a / x_i} b_a\left(\boldsymbol{x}_a\right)=b_i\left(x_i\right), \forall i, a
\end{aligned}
$$ (b22)

前两条是归一化条件，第三条则是所谓的**边际的一致性条件**，即如果一个变量节点 $i$ 参与了功能节点 $a$，那么除了 $i$ 以外的所有参与功能节点 $a$ 的变量节点的分布求和后等于 $i$ 的分布。

构造Lagrangian

$$
\begin{aligned}
\mathcal{L}_{\text {Bethe }}= & -\sum_a \sum_{\boldsymbol{x}_a} b_a\left(\boldsymbol{x}_a\right) \ln f_a\left(\boldsymbol{x}_a\right)+\sum_a \sum_{\boldsymbol{x}_a} b_a\left(\boldsymbol{x}_a\right) \ln b_a\left(\boldsymbol{x}_a\right) \\
&-\sum_i\left(d_i-1\right) \sum_{x_i} b_i\left(x_i\right) \ln b_i\left(x_i\right)\\
&+\sum_i \lambda_i\left(\sum_{x_i} b_i\left(x_i\right)-1\right)+\sum_a \lambda_a\left(\sum_{\boldsymbol{x}_a} b_a\left(\boldsymbol{x}_a\right)-1\right)\\
&+\sum_{i, a} \sum_{x_i} \rho_{i, a}\left(x_i\right)\left(\sum_{\boldsymbol{x}_a / x_i} b_a\left(\boldsymbol{x}_a\right)-b_i\left(x_i\right)\right)
\end{aligned}
$$ (b23)

通过 $\frac{\partial \mathcal{L}_{\text {Bethe }}}{\partial b_a}=0$ 和 $\frac{\partial \mathcal{L}_{\text {Bethe }}}{\partial b_i}=0$ 得

$$
\begin{aligned}
&\hat{b}_i\left(x_i\right) \propto \prod_{a \in \partial i} P_{a \rightarrow i}\left(x_i\right)\\
&\hat{b}_a\left(\boldsymbol{x}_a\right)\propto f_a\left(\boldsymbol{x}_a\right)\prod_{i\in\partial a} P_{i \rightarrow a}\left(x_i\right)
\end{aligned}
$$ (b24)

其中，

$$
\begin{aligned}
&P_{a \rightarrow i}\left(x_i\right)=\sum_{\boldsymbol{x}_j: j \in \partial a \backslash i} f_a\left(\boldsymbol{x}_a\right) \prod_{j \in \partial a \backslash i} P_{j \rightarrow a}\left(x_j\right)\\
&P_{i \rightarrow a}=\frac{1}{Z_{i\to a}}\prod_{b\in i\backslash a}P_{b\to i}
\end{aligned}
$$ (b25)

```{admonition} 推导细节 e
:class: dropdown
再说吧
```

### 3.2 Bethe 近似在树状图中是精确的

如下图所示是一个无环的树状因子图（从任意一个节点出发无法回到自身节点），这种因子图表示的网络只存在最近邻相互作用，因此在这种树状图中 Bethe 近似是精确的。

<p style="text-align:center"><img src="https://user-images.githubusercontent.com/106574511/212237161-ad9281e7-c51a-4056-8950-694274b9e9d1.png" alt="" class="bg-primary" width="400px"></p>

对于该因子图，其联合概率分布可以拆分为

$$
p\left(x_1, x_2, x_3, x_4\right)=\frac{1}{Z} f_A\left(x_1, x_2\right) f_B\left(x_2, x_3, x_4\right) f_C\left(x_4\right)
$$ (c1)

我们考虑第二个变量节点，即下图中红色的节点

<p style="text-align:center"><img src="https://user-images.githubusercontent.com/106574511/212237476-ecca974c-59cb-4865-b670-566491f1da2b.png" alt="" class="bg-primary" width="400px"></p>

根据上面得到的迭代方程

$$
\begin{aligned}
&b_i\left(x_i\right) \propto \prod_{a \in \partial i} P_{a \rightarrow i}\left(x_i\right)\\
&P_{a \rightarrow i}\left(x_i\right)=\sum_{\boldsymbol{x}_j: j \in \partial a \backslash i} f_a\left(\boldsymbol{x}_a\right) \prod_{j \in \partial a \backslash i} P_{j \rightarrow a}\left(x_j\right)
\end{aligned}
$$ (c2)

则

$$
\begin{aligned}
b_i\left(x_2\right) & \propto P_{A \rightarrow 2}\left(x_2\right) P_{B \rightarrow 2}\left(x_2\right) \\
& \propto\left(\sum_{x_1} f_A\left(x_1, x_2\right) P_{1 \rightarrow A}\left(x_1\right)\right)\left(\sum_{x_3, x_4} f_B\left(x_2, x_3, x_4\right) P_{3 \rightarrow B}\left(x_3\right) P_{4 \rightarrow B}\left(x_4\right)\right)
\end{aligned}
$$ (c3)

注意到变量节点 $1$ 和 $3$ 都只参与了一个功能节点，因此 $P_{1 \rightarrow A}$ 和 $P_{3 \rightarrow B}$ 都等于 $1$。变量节点 $4$ 参与了 $B$、$C$ 两个功能节点，根据迭代方程

$$
P_{i \rightarrow a}=\frac{1}{Z_{i\to a}}\prod_{b\in i\backslash a}P_{b\to i}
$$ (c4)

有

$$
P_{4 \rightarrow B}=P_{C\to 4}(x_4)
$$ (c5)

再将 $P_{C\to 4}(x_4)$ 放入迭代方程中，由于功能节点 $C$ 只与 $4$ 这一个变量节点相连，$P_{C\to 4}(x_4)=f_C(x_4)$。因此

$$
\begin{aligned}
b_i\left(x_2\right) & \propto\left(\sum_{x_1} f_A\left(x_1, x_2\right) P_{1 \rightarrow A}\left(x_1\right)\right)\left(\sum_{x_3, x_4} f_B\left(x_2, x_3, x_4\right) P_{3 \rightarrow B}\left(x_3\right) P_{4 \rightarrow B}\left(x_4\right)\right)\\
& \propto\left(\sum_{x_1} f_A\left(x_1, x_2\right)\right)\left(\sum_{x_1, x_4} f_B\left(x_2, x_3, x_4\right) P_{C \rightarrow 4}\left(x_4\right)\right) \\
& \propto\left(\sum_{x_1} f_A\left(x_1, x_2\right)\right)\left(\sum_{x_3, x_4} f_B\left(x_2, x_3, x_4\right) f_C\left(x_4\right)\right) \\
& \propto \sum_{x_1, x_3, x_4} f_A\left(x_1, x_2\right) f_B\left(x_2, x_3, x_4\right) f_C\left(x_4\right)
\end{aligned}
$$ (c6)

注意到 $\sum_{x_1, x_3, x_4} f_A\left(x_1, x_2\right) f_B\left(x_2, x_3, x_4\right) f_C\left(x_4\right)$ 正是因子图中变量节点 $a$ 的边际概率，即我们想要的 $b_i(x_2)$ 的形式。

同样地，对于功能节点 $B$，有

$$
\begin{aligned}
b_a\left(x_B\right) & \propto f_B\left(x_2, x_3, x_4\right) P_{2 \rightarrow B}\left(x_2\right) P_{3 \rightarrow B}\left(x_3\right) P_{4 \rightarrow B}\left(x_4\right) \\
& \propto f_B\left(x_2, x_3, x_4\right)\left(P_{A \rightarrow 2}\left(x_2\right)\right)\left(P_{C \rightarrow 4}\left(x_4\right)\right) \\
& \propto f_B\left(x_2, x_3, x_4\right)\left(\sum_{x_1} f_A\left(x_1, x_2\right)\right)\left(P_{B \rightarrow 2}\left(x_2\right) P_{B \rightarrow 3}\left(x_3\right)\right) \\
& \propto f_B\left(x_2, x_3, x_4\right)\left(\sum_{x_1} f_A\left(x_1, x_2\right)\right)\left(f_C\left(x_4\right)\right) \\
& \propto \sum_{x_1} f_A\left(x_1, x_2\right) f_B\left(x_2, x_3, x_4\right) f_C\left(x_4\right)
\end{aligned}
$$ (c7)

这也是 $b_a\left(x_B\right)$ 的精确的形式。

类似地，我们可以证明，对于这个树状因子图中所有的 $b_i(x_i)$ 和 $b_a(x_a)$ 都是精确的。


### 3.3 Bethe 近似与空腔方法是等价的

### 3.4 Bethe 近似可以退化到平均场近似

## 四、作业

### Assignment 1

### Assignment 2

### Assignment 3