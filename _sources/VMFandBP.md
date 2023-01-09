# 变分平均场理论与信念传播算法

## 一、引言

对于一个真实分布

$$
P(\boldsymbol{x})=\frac{e^{-\beta E(x)}}{Z} \qquad Z=\sum_x e^{-\beta E(x)} \qquad E(\boldsymbol{x})=-\sum_a J_a \prod_{i \in \partial a} x_i
$$

这里采用热力学的写法，将能量写为 $E(x)$。

根据 $Z$ 可以计算Helmholtz自由能 $F_H=-\ln Z$，但是 $Z$ 是 $\mathcal{O}(N^2)$，很难直接计算。本节通过构造一个Gibbs自由能 $F(b)$（其中 $b(x)$ 是真实分布，$F(b)$ 是泛函）并使用变分方法使其逼近真实的Helmholtz自由能 $F_H=-\ln Z$。

**最大熵规则：**  给定分布的约束下，系统的熵越大，信息量越大，越应该是真实的分布。说明如下：

对于一个未知分布 $P_r$，其熵定义为

$$
S=-k\sum_{r}P_r\ln P_r
$$

该分布要满足归一化和能量守恒（系统的总能量为定值）的约束，即

$$
\sum_r P_r=1 \qquad \sum_r E_r P_r=\mu
$$

使用拉格朗日乘子法求解约束极值，构造Lagrangian

$$
L=-k\sum_{r}P_r\ln P_r+\lambda_1\left(\sum_r P_r-1\right)+\lambda_2\left(\sum_r E_r P_r-\mu\right)
$$

取其变分为零

$$
\frac{\delta L}{\delta P_r}=\sum_{r}\left[-k(\ln P_r+1)+\lambda_1+\lambda_2E_r \right]=0
$$

即

$$
P_r=\exp\left( -1+\frac{\lambda_1}{k}+\frac{\lambda_2}{k}E_r\right)
$$

代入 $P_r$ 的归一化条件得

$$
\begin{aligned}
    &\sum_r P_r=\sum_r \exp\left( -1+\frac{\lambda_1}{k}+\frac{\lambda_2}{k}E_r\right)=1\\
    \Longrightarrow &\exp\left( -1+\frac{\lambda_1}{k}\right)\sum_r \exp\left(\frac{\lambda_2}{k}E_r\right)=1\\
    \Longrightarrow &\exp\left( -1+\frac{\lambda_1}{k}\right)=\frac{1}{\sum_r \exp\left(\frac{\lambda_2}{k}E_r\right)}
\end{aligned}
$$

令 $-\beta \equiv \frac{\lambda_2}{k}$，因此 $P_r$ 可写为

$$
P_r=\frac{\exp\left(-\beta E_r\right)}{\sum_r \exp\left(-\beta E_r\right)}
$$

这正是真实分布中 $P(\boldsymbol{x})$ 的表达式。



## 二、变分自由能

对于真实分布 $b(x)$，Gibbs自由能 $F(b)$ 是 $b(x)$ 的泛函，并且设 $\beta=1$。在热力学中，Gibbs自由能写为

$$
F(b)=U(b)-H(b)
$$

其中，$U(b)$ 是内能；$H(b)$ 是熵（而不是哈密顿量），他们分别写为

$$
\begin{aligned}
    &U(b)=\sum_x b(x) E(x) \\
    &H(b)=-\sum_x b(x) \ln b(x)
\end{aligned}
$$

所谓的变分自由能就是指我们基于真实分布 $b(x)$ 构造的Gibbs自由能 $F(b)$，因为是用变分的方法去逼近真实的Helmholtz自由能，所以又称为变分自由能。对Gibbs自由能 $F(b)$ 和Helmholtz自由能 $F_H$ 作差，得

$$
\begin{aligned}
F(b)-F_H &=\sum_{\boldsymbol{x}} b(\boldsymbol{x}) E(\boldsymbol{x})+\sum_{\boldsymbol{x}} b(\boldsymbol{x}) \ln b(\boldsymbol{x})+\ln Z \\
&=\sum_{\boldsymbol{x}} b(\boldsymbol{x})(-\ln Z-\ln P(\boldsymbol{x}))+\sum_{\boldsymbol{x}} b(\boldsymbol{x}) \ln b(\boldsymbol{x})+\ln Z \\
&=\sum_{\boldsymbol{x}} b(\boldsymbol{x}) \ln \frac{b(\boldsymbol{x})}{P(\boldsymbol{x})}=D_{KL}(b \| P)
\end{aligned}
$$

其中，第二步用到了由 $P(\boldsymbol{x})=\frac{e^{- E(\boldsymbol{x})}}{Z}$ 得到的 $E(\boldsymbol{x})=-\ln Z-\ln P(\boldsymbol{x})$。可以看到，二者之差为一个KL散度的形式。

**Kullback-Leibler Divergence：  ** KL散度是用来度量两个概率分布相似度的指标。假设对随机变量 $\xi$ 存在两个概率分布 $P$、$Q$，定义从 $P$ 到 $Q$ 的KL散度为
$$
\mathbb{D}_{\mathrm{KL}}(P \| Q)=\sum_i P(i) \ln \left(\frac{P(i)}{Q(i)}\right) \quad\text{或}\quad \int_{-\infty}^{\infty} p(\mathbf{x}) \ln \left(\frac{p(\mathbf{x})}{q(\mathbf{x})}\right) \mathrm{~d} \mathbf{x}
$$

可以证明，KL散度具有非负性。以离散形式为例：根据Jensen's不等式，

$$
    \sum_{i=1}^n \lambda_i f\left(x_i\right) \leq f\left(\sum_{i=1}^n \lambda_i x_i\right), \quad\sum_{i=1}^n \lambda_i=1
$$

并且有 $\sum_i P(i)=1$、$\sum_i Q(i)=1$

$$
\sum_i P(i) \ln \left(\frac{Q(i)}{P(i)}\right) \leq \ln\left(\sum_i P(i)\cdot \frac{Q(i)}{P(i)} \right)=\ln\left(\sum_i Q(i)\right)=\ln(1)=0
$$

$$
    \mathbb{D}_{\mathrm{KL}}(P \| Q)=-\sum_i P(i) \ln \left(\frac{Q(i)}{P(i)}\right) \geq 0
$$

当且仅当$P=Q$时取等号。

根据KL散度的非负性，可以知道 $F(b)-F_H=D_{KL}(b \| P)\geq 0$，即变分自由能的最小值为Helmholtz自由能，此时分布 $b(x)=P(x)$。

为了求解自由能，我们首先要给定 $b(x)$ 的一个具体形式。接下来两个小节分别给出了 $b(x)$ 的两种近似：Mean-Field Approximation和Bethe Approximation。



## 三、Mean-Field Approximation

Mean-Field Approximation是假设该分布由互相独立的自旋格点 $x_i\in\{+1,-1\}$ 组成，其均值即磁化强度为
$$
m_i=(+1)\times \frac{1+m_{i}(+1)}{2}+(-1)\times \frac{1+m_{i}(-1)}{2}
$$
因此可以用 $m_i$ 将单个格点的分布 $b_i$ 参数化：
$$
b_i(x_i)=\frac{1+m_ix_i}{2}
$$
由于各个格点是相互独立的，系统的总分布就等于各个格点的分布相乘，因此系统的Mean-Field Approximation分布为
$$
b_{MF}=\prod_ib_i(x_i)=\prod_i\frac{1+m_ix_i}{2}
$$
在此近似下，变分内能为
$$
\begin{aligned}
U_{M F} & =\sum_{\boldsymbol{x}} b_{M F}(\boldsymbol{x}) E(\boldsymbol{x}) \\
& =\sum_{\boldsymbol{x}} \prod_i b_i\left(x_i\right)\left(-\sum_a J_a \prod_{i \in \partial a} x_i\right) \\
& = -\sum_a J_a \left[\sum_{\boldsymbol{x}} \prod_i b_i\left(x_i\right)\prod_{i \in \partial a} x_i \right]\\

\end{aligned}
$$
注意到 $x_i$ 是相互独立的，接下来的操作相当于是 $\langle ab \rangle=\langle a \rangle\langle b \rangle$，其中 $a$、$b$ 相互独立
$$
\begin{aligned}
U_{M F} & = -\sum_a J_a \left[\sum_{\boldsymbol{x}} \prod_i b_i\left(x_i\right)\prod_{i \in \partial a} x_i \right]\\
& = -\sum_a J_a \prod_{i \in \partial a}\left[\sum_{\boldsymbol{x}} \prod_i b_i\left(x_i\right) x_i \right] \\
& =-\sum_a J_a \prod_{i \in \partial a} m_i,
\end{aligned}
$$
变分熵为
$$
\begin{aligned}
H_{M F} & =-\sum_x b_{M F}(x) \ln b_{M F}(x) \\
& =-\sum_x \prod_j \frac{1+m_j x_j}{2} \ln \prod_i \frac{1+m_i x_i}{2} \\
& =-\sum_i \sum_x \prod_j \frac{1+m_j x_j}{2} \ln \frac{1+m_i x_i}{2} \\
& =-\sum_i \sum_{x_i} \sum_{x \backslash x_i} \prod_{j(\neq i)} \frac{1+m_j x_j}{2} \frac{1+m_i x_i}{2} \ln \frac{1+m_i x_i}{2} \\
& =-\sum_i \sum_{x_i} \frac{1+m_i x_i}{2} \ln \frac{1+m_i x_i}{2} \\
& =\sum_i S_i,
\end{aligned}
$$
其中，
$$
\sum_{x \backslash x_i} \prod_{j(\neq i)} \frac{1+m_j x_j}{2}=\sum_{x \backslash x_i}b_{MF\backslash x_i}=1
$$
即相当于去掉格点 $i$ 后的系统的总分布，由于各格点相互独立，此时系统的分布仍是归一化的。$S_i$ 是单个自旋变量的熵
$$
S_i=-k_B\sum_ip_i\ln p_i=- \sum_{x_i} \frac{1+m_i x_i}{2} \ln \frac{1+m_i x_i}{2}
$$
式中，已令 $\beta=1$，因此 $k_B=\frac{1}{\beta}=1$，$p_i$ 在此处是Mean-Field Approximation下的 $b_i$。

因此变分自由能为
$$
\begin{aligned}
F_{M F} & =U_{M F}-H_{M F} \\
& =-\sum_a J_a \prod_{i \in \partial a} m_i+\sum_i \sum_{x_i} \frac{1+m_i x_i}{2} \ln \frac{1+m_i x_i}{2}
\end{aligned}
$$
令变分自由能取极值
$$
\begin{aligned}
\frac{\partial F_{M F}}{\partial m_i} & =-\sum_a J_a \prod_{j \in \partial a \backslash i} m_j+\sum_{x_i} \frac{x_i}{2} \ln \frac{1+m_i x_i}{2}+\frac{x_i}{2} \\
& =-\sum_a J_a \prod_{j \in \partial a \backslash i} m_j+\frac{1}{2} \ln \frac{1+m_i}{1-m_i}=0
\end{aligned}
$$
即解
$$
\ln \frac{1+m_i}{1-m_i}=2\sum_a J_a \prod_{j \in \partial a \backslash i} m_j=2A
$$
由
$$
\frac{1+m_i}{1-m_i}=e^{2A}\quad \Longrightarrow\quad m_i=\frac{e^{2A}-1}{e^{2A}+1}=\frac{e^{A}-e^{-A}}{e^{A}+e^{-A}}
$$
即
$$
m_i=\tanh \left(\sum_{a \in \partial i} J_a \prod_{j \in \partial a \backslash i} m_j\right)
$$
这是一个迭代方程，将其迭代到不动点即可解出 $m_i$ 并找到此时对应的真实分布



## 四、Bethe Approximation

### 4.1 Bethe Approximation与信念传播算法

注意到Mean-Field Approximation要求各个格点之间完全相互独立，这只能在温度比较高的时候实现，对于温度较低的情况，必须要考虑格点之间的关联。本节中的Bethe Approximation就引入了格点间的部分关联。

仍然使用如下的分布
$$
P(\boldsymbol{x})=\frac{e^{-\beta E(x)}}{Z} \qquad E(\boldsymbol{x})=-\sum_a J_a \prod_{i \in \partial a} x_i
$$
为了更好地在因子图中表示，我们首先进行一个变量替换，令
$$
f_a(\boldsymbol{x}_a)=e^{J_a \prod_{i \in \partial a} x_i}
$$
$f_a(\boldsymbol{x}_a)$ 表示一个相互作用 $a$，其中 $\boldsymbol{x}_a$ 表示参与该相互作用的所有自旋变量。同样令 $\beta=1$，因此
$$
\begin{aligned}
& P(\boldsymbol{x})=\frac{1}{Z} \prod_a f_a\left(\boldsymbol{x}_a\right) \\
& E(\boldsymbol{x})=-\sum_a \ln f_a\left(\boldsymbol{x}_a\right) 
\end{aligned}
$$
显然，此时分布 $P(x)$ 是因子连乘的形式，用因子图表示为

![](https://images.cnblogs.com/cnblogs_com/blogs/779311/galleries/2257667/o_221226125622_%E8%B4%9D%E8%92%82%E8%BF%91%E4%BC%BC%E7%9A%84%E5%9B%A0%E5%AD%90%E5%9B%BE.png)

这里variable node是表示自旋的变量节点，function node是表示相互作用的功能节点。对这个因子图进行区域划分，每一个区域都要将与其中的功能节点相连的所有变量节点都包含在里面。上图中的黄色虚线圈出来的不能算作一个区域，因为与功能节点A相连的变量节点2没有被包含在其中。

对于一个区域 $R$，我们可以计算其内部的变分内能和变分熵：
$$
\begin{aligned}
& E_R\left(x_R\right)=-\sum_{a \in R} \ln f_a\left(\boldsymbol{x}_a\right) \\
& U_R\left(b_R\right)=\sum_{\boldsymbol{x}_R} b_R\left(\boldsymbol{x}_R\right) E_R\left(\boldsymbol{x}_R\right) \\
& H_R\left(b_R\right)=-\sum_{\boldsymbol{x}_R} b_R\left(\boldsymbol{x}_R\right) \ln b_R\left(\boldsymbol{x}_R\right)
\end{aligned}
$$
 然后将各个区域的变分内能和变分熵相加，即
$$
\begin{gathered}
U_{\mathcal{R}}=\sum_{R \in \mathcal{R}} C_R U_R\left(b_R\right)\\
H_{\mathcal{R}}=\sum_{R \in \mathcal{R}} C_R H_R\left(b_R\right)
\end{gathered}
$$
其中，$\mathcal{R}$ 是所有区域的集合，$C_R$ 是用于对区域进行计数的变量，满足
$$
\begin{gathered}
\sum_{R \in \mathcal{R}} \mathbb{I}[a \in R] C_R=1, \forall a \\
\sum_{R \in \mathbb{R}} \mathbb{I}[i \in R] C_R=1, \forall i
\end{gathered}
$$
$\mathbb{I}[a \in R]$ 表示如果功能节点 $a$ 在区域 $R$ 中，则 $C_R=1$，否则 $C_R=0$；同理，$\mathbb{I}[i \in R]$ 表示如果变量节点 $i$ 在区域 $R$ 中，则 $C_R=1$，否则 $C_R=0$。

$\mathbb{I}[a \in R]$ 与 $C_R$ 相乘再求和后等于 $1$ 表示节点 $a$ 只被考虑了一次，因此上式对于 $C_R$ 的约束实际上是要求所有的变量节点和功能节点在求和时只被考虑一次，这是就划分区域的要求。

Bethe Approximation划分区域的方式如下图所示

![](https://images.cnblogs.com/cnblogs_com/blogs/779311/galleries/2257667/o_221226132342_%E8%B4%9D%E8%92%82%E8%BF%91%E4%BC%BC%E7%9A%84%E5%88%92%E5%88%86%E5%8C%BA%E5%9F%9F%E7%9A%84%E6%96%B9%E5%BC%8F.png)

划分的区域有两种，一种称作“大区域”，是包含一个功能节点以及所有与其相连的变量节点的区域；另一种称作“小区域”，只包含一个变量节点。大区域的计数设为 $1$，小区域的计数设为 $1-d_i$，其中 $d_i$ 是变量节点 $i$ 参与的功能节点的数量，即
$$
C_{R_{L}}=1,\quad C_{R_{S}}=1-d_i
$$
这种划分方式和计数方式很自然地满足了上面对于 $C_R$ 的两条约束：对于第一条（针对功能节点的）约束，只有大区域中有功能节点参与了求和，并且对于每个功能节点 $a$，当其包含在某个区域 $R$ 中使得 $\mathbb{I}[a \in R]=1$ 时，都有且仅有这一个区域的计数 $C_{R_L}=1$；对于第二条（针对变量节点的）约束，每个变量节点 $i$ 在其参与的 $d_i$ 个大区域（功能节点）中使得 $\mathbb{I}[a \in R]=1,\,C_{R_{L}}=1$，并且使得其自身所在的小区域 $\mathbb{I}[a \in R]=1,\,C_{R_{S}}=1-d_i$，所有大区域的求和为 $d_i$，只有一个小区域参与了求和为 $1-d_i$，因此总的求和为 $1$。

Bethe Approximation划分区域后，只考虑各个区域（主要指大区域）内部的相互作用，不考虑区域与区域间的相互作用，实际上是一个只考虑最近邻的相互作用的模型。

根据以上划分的区域，变分内能为
$$
\begin{aligned}
U_{\text {Bethe }}= &\sum_{R \in \mathcal{R}} C_R U_R\left(b_R\right)=\sum_{R \in \mathcal{R}} C_R\sum_{\boldsymbol{x}_R} b_R\left(\boldsymbol{x}_R\right) E_R\left(\boldsymbol{x}_R\right)\\
=& \sum_a\sum_{x_a} b_a\left(\boldsymbol{x}_a\right) \left[-\sum_{a \in R} \ln f_a\left(\boldsymbol{x}_a\right)\right] \\
=&-\sum_a \sum_{x_a} b_a\left(\boldsymbol{x}_a\right) \ln f_a\left(\boldsymbol{x}_a\right) \\

\end{aligned}
$$
这里由于哈密顿量没有外场，实际上只有大区域贡献了内能，因此第一行中的求和下标 $R\in\mathcal{R}$ 有 $a$ 个，并且每个大区域中只有一个功能节点，第二行中的求和下标 $a\in R$ 实际上只对一个 $a$ 求和，$C_R=C_{R_L}=1$。如果哈密顿量是具有外场的形式，那么求Bethe内能时也需要考虑小区域的作用。

变分熵为（大区域和小区域分别计算后求和）
$$
\begin{aligned}
H_{\text {Bethe }}= &\sum_{R \in \mathcal{R}} C_R H_R\left(b_R\right)=\sum_{R \in \mathcal{R}} C_R \left[-\sum_{\boldsymbol{x}_R} b_R\left(\boldsymbol{x}_R\right) \ln b_R\left(\boldsymbol{x}_R\right)\right]\\
=&-\sum_aC_{R_L}\sum_{\boldsymbol{x}_a} b_a\left(\boldsymbol{x}_a\right) \ln b_a\left(\boldsymbol{x}_a\right)-\sum_iC_{R_S}\sum_{\boldsymbol{x}_i} b_i\left(\boldsymbol{x}_i\right) \ln b_i\left(\boldsymbol{x}_i\right) \\
=& -\sum_a \sum_{x_a} b_a\left(\boldsymbol{x}_a\right) \ln b_a\left(\boldsymbol{x}_a\right)+\sum_i\left(d_i-1\right) \sum_{x_i} b_i\left(x_i\right) \ln b_i\left(x_i\right) 
\end{aligned}
$$
因此变分自由能为
$$
\begin{aligned}
F_{\text {Bethe }}= & \; U_{\text {Bethe }}-H_{\text {Bethe }}\\
= & -\sum_a \sum_{x_a} b_a\left(\boldsymbol{x}_a\right) \ln f_a\left(\boldsymbol{x}_a\right)+\sum_a \sum_{x_a} b_a\left(\boldsymbol{x}_a\right) \ln b_a\left(\boldsymbol{x}_a\right) \\
&-\sum_i\left(d_i-1\right) \sum_{x_i} b_i\left(x_i\right) \ln b_i\left(x_i\right)
\end{aligned}
$$
在这里我们对变分自由能取机制时存在一些约束：
$$
\begin{aligned}
& \sum_{x_i} b_i\left(x_i\right)=1, \forall i \\
& \sum_{\boldsymbol{x}_a} b_a\left(\boldsymbol{x}_a\right)=1, \forall a \\
& \sum_{\boldsymbol{x}_a / x_i} b_a\left(\boldsymbol{x}_a\right)=b_i\left(x_i\right), \forall i, a
\end{aligned}
$$
前两条是归一化条件，第三条则是边际的一致性条件，即如果一个变量节点 $i$ 参与了功能节点 $a$，那么参与功能节点 $a$ 的除了 $i$ 以外的所有变量节点的分布求和后等于 $i$ 的分布。

因此用拉格朗日乘子法求条件极值时首先构造Lagrangian
$$
\begin{aligned}
L= & -\sum_a \sum_{\boldsymbol{x}_a} b_a\left(\boldsymbol{x}_a\right) \ln f_a\left(\boldsymbol{x}_a\right)+\sum_a \sum_{\boldsymbol{x}_a} b_a\left(\boldsymbol{x}_a\right) \ln b_a\left(\boldsymbol{x}_a\right) \\
&-\sum_i\left(d_i-1\right) \sum_{x_i} b_i\left(x_i\right) \ln b_i\left(x_i\right)\\
&+\sum_i \lambda_i\left(\sum_{x_i} b_i\left(x_i\right)-1\right)+\sum_a \lambda_a\left(\sum_{\boldsymbol{x}_a} b_a\left(\boldsymbol{x}_a\right)-1\right)\\
&+\sum_{i, a} \sum_{x_i} \rho_{i, a}\left(x_i\right)\left(\sum_{\boldsymbol{x}_a / x_i} b_a\left(\boldsymbol{x}_a\right)-b_i\left(x_i\right)\right)
\end{aligned}
$$
令 $\frac{\partial L}{\partial b_a}=0$，得
$$
\hat{b}_a\left(\boldsymbol{x}_a\right)=f_a\left(\boldsymbol{x}_a\right) \exp \left[\lambda_a-1+\sum_{i} \rho_{i, a}\left(x_i\right)\right]
$$
令 $\frac{\partial L}{\partial b_i}=0$，得
$$
\hat{b}_i\left(x_i\right)=\exp \left[\frac{1}{d_i-1}\left(1-\lambda_i+\sum_{a } \rho_{i, a}\left(x_i\right)\right)\right]
$$
这两个解来自于文献[2]和[5]，我自己还没有求解

我们定义从把 $a$ 节点拿掉后 $i$ 节点传递给 $a$ 的“消息”为（文献[1]、[2]中都提到这个定义，但是这里采用的解释来自文献[4]）
$$
P_{i \rightarrow a}\equiv e^{\rho_{i, a}\left(x_i\right)}=\frac{1}{Z_{i\to a}}\prod_{b\in i\backslash a}P_{b\to i}
$$
因此
$$
\begin{aligned}
\hat{b}_a\left(\boldsymbol{x}_a\right)&=f_a\left(\boldsymbol{x}_a\right) \exp \left[\lambda_a-1+\sum_{i\in\partial a} \rho_{i, a}\left(x_i\right)\right]\\
&=e^{\lambda_a-1}f_a\left(\boldsymbol{x}_a\right)\exp \left[\sum_{i\in\partial a} \ln P_{i \rightarrow a}\left(x_i\right)\right]\\
&=e^{\lambda_a-1}f_a\left(\boldsymbol{x}_a\right)\prod_{i\in\partial a}\exp \left[\ln P_{i \rightarrow a}\left(x_i\right)\right]\\
&=e^{\lambda_a-1}f_a\left(\boldsymbol{x}_a\right)\prod_{i\in\partial a} P_{i \rightarrow a}\left(x_i\right)
\end{aligned}
$$

即
$$
\hat{b}_a\left(\boldsymbol{x}_a\right)\propto f_a\left(\boldsymbol{x}_a\right)\prod_{i\in\partial a} P_{i \rightarrow a}\left(x_i\right)
$$

再定义从节点 $a$ 传出到 $i$ 的“消息”为
$$
P_{a\to i}\equiv\frac{\hat{b}_i\left(x_i\right)}{P_{i\to a}}
$$
据文献[4]所说，对上式进行一些操作后可以将 $\hat{b}_i(x_i)$ 写为
$$
\hat{b}_i\left(x_i\right) \propto \prod_{a \in \partial i} P_{a \rightarrow i}\left(x_i\right)
$$
这一步具体怎么操作还没搞明白。并且，有可能是根据上面 $P_{a\to i }$ 的定义和边际的一致性约束条件可以得到
$$
P_{a \rightarrow i}\left(x_i\right)=\sum_{\boldsymbol{x}_j: j \in \partial a \backslash i} f_a\left(\boldsymbol{x}_a\right) \prod_{j \in \partial a \backslash i} P_{j \rightarrow a}\left(x_j\right)
$$
这一步现在也还不理解。

总之，到这里，我们得到了信念传播的迭代方程。

$\quad $

**本节用到的参考文献：**[1]~[7]

文献[2]中有关于本节的较为详细和深入的内容，为整理本节的笔记提供了很大的帮助 

### 4.2 一个例子：Bethe Approximation在树状因子图中是精确的

对于如下图所示的一个无环的树状因子图（从任意一个节点出发无法回到自身节点），对于这类因子图表示的网络，只存在最近邻相互作用，因此Bethe Approximation是精确的。



![](C:\Users\lyhsy\Pictures\Typora\变分平均场\image-20221227171223167.png)

对于该因子图，其联合概率分布可以拆分为
$$
p\left(x_1, x_2, x_3, x_4\right)=\frac{1}{Z} f_A\left(x_1, x_2\right) f_B\left(x_2, x_3, x_4\right) f_C\left(x_4\right)
$$
我们考虑第二个变量节点，即下图中红色的节点

![](https://images.cnblogs.com/cnblogs_com/blogs/779311/galleries/2257667/o_221227092426_BP%E4%B8%BE%E4%BE%8B%E5%B8%A6%E6%B3%A8%E9%87%8A.png)

根据上面得到的迭代方程
$$
b_i\left(x_i\right) \propto \prod_{a \in \partial i} P_{a \rightarrow i}\left(x_i\right), \quad P_{a \rightarrow i}\left(x_i\right)=\sum_{\boldsymbol{x}_j: j \in \partial a \backslash i} f_a\left(\boldsymbol{x}_a\right) \prod_{j \in \partial a \backslash i} P_{j \rightarrow a}\left(x_j\right)
$$
则
$$
\begin{aligned}
b_i\left(x_2\right) & \propto P_{A \rightarrow 2}\left(x_2\right) P_{B \rightarrow 2}\left(x_2\right) \\
& \propto\left(\sum_{x_1} f_A\left(x_1, x_2\right) P_{1 \rightarrow A}\left(x_1\right)\right)\left(\sum_{x_3, x_4} f_B\left(x_2, x_3, x_4\right) P_{3 \rightarrow B}\left(x_3\right) P_{4 \rightarrow B}\left(x_4\right)\right)
\end{aligned}
$$
注意到变量节点 $1$ 和 $3$ 都只参与了一个功能节点，因此 $P_{1 \rightarrow A}$ 和 $P_{3 \rightarrow B}$ 都等于 $1$。变量节点 $4$ 参与了 $B$、$C$ 两个功能节点，根据迭代方程
$$
P_{i \rightarrow a}=\frac{1}{Z_{i\to a}}\prod_{b\in i\backslash a}P_{b\to i}
$$
有
$$
P_{4 \rightarrow B}=P_{C\to 4}(x_4)
$$
再将 $P_{C\to 4}(x_4)$ 放入迭代方程中，由于功能节点 $C$ 只与 $4$ 这一个变量节点相连，$P_{C\to 4}(x_4)=f_C(x_4)$。因此
$$
\begin{aligned}
b_i\left(x_2\right) & \propto\left(\sum_{x_1} f_A\left(x_1, x_2\right) P_{1 \rightarrow A}\left(x_1\right)\right)\left(\sum_{x_3, x_4} f_B\left(x_2, x_3, x_4\right) P_{3 \rightarrow B}\left(x_3\right) P_{4 \rightarrow B}\left(x_4\right)\right)\\
& \propto\left(\sum_{x_1} f_A\left(x_1, x_2\right)\right)\left(\sum_{x_1, x_4} f_B\left(x_2, x_3, x_4\right) P_{C \rightarrow 4}\left(x_4\right)\right) \\
& \propto\left(\sum_{x_1} f_A\left(x_1, x_2\right)\right)\left(\sum_{x_3, x_4} f_B\left(x_2, x_3, x_4\right) f_C\left(x_4\right)\right) \\
& \propto \sum_{x_1, x_3, x_4} f_A\left(x_1, x_2\right) f_B\left(x_2, x_3, x_4\right) f_C\left(x_4\right)
\end{aligned}
$$
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
$$
这也是 $b_a\left(x_B\right)$ 的精确的形式。

类似地，我们可以证明，对于这个树状因子图中所有的 $b_i(x_i)$ 和 $b_a(x_a)$ 都是精确的。



### 4.3 Bethe Approximation与空腔方法的等价性

接下来我们只考虑空腔概率 $P_{i\to a}$，将
$$
P_{a \rightarrow i}\left(x_i\right)=\sum_{\boldsymbol{x}_j: j \in \partial a \backslash i} f_a\left(\boldsymbol{x}_a\right) \prod_{j \in \partial a \backslash i} P_{j \rightarrow a}\left(x_j\right)
$$
代入
$$
P_{i \rightarrow a}=\frac{1}{Z_{i\to a}}\prod_{b\in i\backslash a}P_{b\to i}
$$
并将 $f_a(x_a)$ 的定义
$$
f_a(\boldsymbol{x}_a)=e^{J_a \prod_{i \in \partial a} x_i}
$$
一并代入，得
$$
P_{i \rightarrow a}\left(x_i\right)=\frac{1}{Z_{i \rightarrow a}} \prod_{b \in \partial i \backslash a} \sum_{x_j: j \in \partial b \backslash i} e^{J_a \prod_{i \in \partial a} x_i} \prod_{j \in \partial b \backslash i} P_{j \rightarrow b}\left(x_j\right)
$$
为了证明Bethe Approximation与空腔方法是等价的，我们先假设 $P_{j\to b}$ 是一个空腔概率，根据[空腔方法](https://www.cnblogs.com/lyhspace/p/17003818.html)第二小节中的内容，若 $x_j\in\{\pm 1\}$，可以将其参数化为
$$
P_{j\to b}(x_j)=\frac{1+m_{j\to b}x_j}{2}
$$
接下来只要能根据该假设得到空腔的磁化强度的迭代方程，即说明该假设是自洽的。

考虑 $x_i=+1$ 时，令
$$
A_b^{+}=\sum_{x_j: j \in \partial b \backslash i} e^{J_b \prod_{j \in \partial b \backslash i} x_j} \prod_{j \in \partial b \backslash i} \frac{1+m_{j \rightarrow b} x_j}{2}
$$
这里与[空腔方法](https://www.cnblogs.com/lyhspace/p/17003818.html)第二小节中对 $\prod_{i \in \partial a} \frac{1+\sigma_i m_{i \rightarrow a}}{2}$ 进行展开处理的操作类似，经过这个trick后可以得到
$$
A_b^{+}=\cosh J_b\left(1+\tanh J_b \prod_{j \in \partial b \backslash i} m_{j \rightarrow b}\right)
$$
类似地，当 $x_i=-1$ 时，可得
$$
A_b^{-}=\cosh J_b\left(1-\tanh J_b \prod_{j \in \partial b \backslash i} m_{j \rightarrow b}\right)
$$
根据空腔磁化强度的定义
$$
m_{i\to a}=\sum_{x_i}x_iP_{i\to a}(x_i)=\sum_{x_i}x_i\frac{\exp\left[-\beta H_{i\to a}(x_i)\right]}{\sum_{x_i}\exp\left[-\beta H_{i\to a}(x_i)\right]}
$$
注意到我们这里定义的 $A_b^\pm$ 与[空腔方法](https://www.cnblogs.com/lyhspace/p/17003818.html)第三小节空腔磁化强度部分的 $\Lambda_{b\to i}^\pm$ 一致，按照相同的步骤可得
$$
\begin{aligned}
m_{i \rightarrow a} 
& =\frac{\prod_{b \in \partial i \backslash a} A_b^{+}-\prod_{b \in \partial i \backslash a} A_b^{-}}{\prod_{b \in \partial i \backslash a} A_b^{+}+\prod_{b \in \partial i \backslash a} A_b^{-}} \\
& =\frac{\prod_{b \in \partial i \backslash a}\left(1+\prod_1 \tanh J_b \prod_{j \in \partial b \backslash i} m_{j \rightarrow b}\right)-\prod_{b \in \partial i \backslash a}\left(1-\tanh J_b \prod_{j \in \partial b \backslash i} m_{j \rightarrow b}\right)}{\prod_{b \in \partial i \backslash a}\left(1+\tanh J_b \prod_{j \in \partial b \backslash i} m_{j \rightarrow b}\right)+\prod_{b \in \partial i \backslash a}\left(1-\tanh J_b \prod_{j \in \partial b \backslash i} m_{j \rightarrow b}\right)}
\end{aligned}
$$
令
$$
\tanh u_{b \rightarrow i}=\tanh J_b \prod_{j \in \partial b \backslash i} m_{j \rightarrow b}
$$
即得[空腔方法](https://www.cnblogs.com/lyhspace/p/17003818.html)第三小节中空腔磁化强度的迭代方程：
$$
\begin{aligned}
&m_{i \rightarrow a}  =\tanh \left(\sum_{b \in \partial i \backslash a} u_{b \rightarrow i}\right) \\
&\tanh u_{b \rightarrow i}  =\tanh J_b \prod_{j \in \partial b \backslash i} m_{j \rightarrow b}
\end{aligned}
$$
从这里可以看到，Bethe Approximation与空腔方法是等价的。



### 4.4 Bethe Approximation与Mean-Field Approximation

Mean-Field Approximation是忽略自旋变量间的所有相互作用，Bethe Approximation考虑了最近邻的相互作用，因此Bethe Approximation是比Mean-Field Approximation更高一级的近似。一般来说，高一级的近似是可以退化到低一级的近似的。本节将展示Bethe Approximation如何退化到Mean-Field Approximation。

上一小节提到，Bethe Approximation与cavity method是等价的。因此我们继续从上一小节空腔的磁化强度出发，将 $\sum_{b \in \partial i \backslash a} u_{b \rightarrow i}$ 改写为
$$
m_{i \rightarrow a}  =\tanh \left(\sum_{b \in \partial i} u_{b \rightarrow i}- u_{a\to i}\right)
$$
根据 $m_i=\tanh\left( \sum_{b \in \partial i} u_{b \rightarrow i}\right)$ （已令 $\beta=1$），可以将 $\sum_{b \in \partial i} u_{b \rightarrow i}$ 反表示为 $\tanh^{-1}(m_i)$。

再利用
$$
\tanh (x+y)=\frac{\tanh x+\tanh y}{1+\tanh x \tanh y}
$$
以及
$$
\tanh u_{a \rightarrow i}  =\tanh J_a \prod_{j \in \partial a \backslash i} m_{j \rightarrow a}
$$
可得
$$
\begin{aligned}
m_{i \rightarrow a} & =\tanh \left(\sum_{b \in \partial i} u_{b \rightarrow i}-u_{a \rightarrow i}\right) \\
& =\tanh \left(\tanh ^{-1}\left(m_i\right)-u_{a \rightarrow i}\right) \\
& =\frac{m_i-\tanh \beta J_a \prod_{j \in \partial a \backslash i} m_{j \rightarrow a}}{1-m_i \tanh \beta J_a \prod_{j \in \partial a \backslash i} m_{j \rightarrow a}},
\end{aligned}
$$
考虑一个 $i$ 和 $j$ 的二体相互作用，上式可以改写为两个迭代式：
$$
\begin{aligned}
& m_{i \rightarrow j}=\frac{m_i-\tanh \beta J_{i j} m_{j \rightarrow i}}{1-m_i \tanh \beta J_{i j} m_{j \rightarrow i}} \\
& m_{j \rightarrow i}=\frac{m_j-\tanh \beta J_{i j} m_{i \rightarrow j}}{1-m_j \tanh \beta J_{i j} m_{i \rightarrow j}}
\end{aligned}
$$
将两个式子互相带入，可分别得到 $m_{i\to j}$ 和 $m_{j\to i}$ 关于 $m_i$、$m_j$ 及 $J_{ij}$ 的表达式，将其记为 $f$，即
$$
\begin{aligned}
m_{i \rightarrow j} & =f\left(m_i, m_j, \tanh \beta J_{i j}\right) \\
m_{j \rightarrow i} & =f\left(m_j, m_i, \tanh \beta J_{i j}\right) 
\end{aligned}
$$
其中，
$$
f(a, b, t) =\frac{1-t^2-\sqrt{\left(1-t^2\right)^2-4 t(a-b t)(b-a t)}}{2 t(b-a t)}
$$
（具体怎么算出来的，等以后再补充）

对于 $m_i=\tanh\left( \sum_{b \in \partial i} u_{b \rightarrow i}\right)$，根据（注意从这里开始不再忽略 $\beta$）
$$
\tanh u_{b \rightarrow i}  =\tanh\beta J_b \prod_{j \in \partial b \backslash i} m_{j \rightarrow b}
$$
在二体相互作用中写为
$$
\tanh u_{j \rightarrow i}  =m_{i \rightarrow j} \tanh \beta J_{ij}
$$
其中，$m_{i \rightarrow j} =f\left(m_i, m_j, \tanh \beta J_{i j}\right)$。因此
$$
u_{j \rightarrow i}  =\tanh^{-1}\Big[ f\left(m_i, m_j, \tanh \beta J_{i j}\right) \tanh \beta J_{ij}\Big]
$$
磁化强度写为
$$
\begin{aligned}
m_i&=\tanh\left( \sum_{j \in \partial i} u_{j \rightarrow i}\right)\\
&=\tanh\left\{\sum_{j \in \partial i}\tanh^{-1}\big[ f\left(m_i, m_j, \tanh \beta J_{i j}\right) \tanh \beta J_{ij}\big]\right\}
\end{aligned}
$$
考虑温度较高时，即 $\beta$ 较小时，根据近似
$$
\begin{aligned}
& \tanh ^{-1} x \approx x, \tanh x \approx x \\
& (1+x)^a \approx 1+a x+\frac{1}{2} a(a-1) x^2
\end{aligned}
$$
有
$$
\begin{aligned}
& \tanh \beta J_{i j}\approx \beta J_{i j} \\
& \tanh^{-1}\big[ f\left(m_i, m_j, \tanh \beta J_{i j}\right) \tanh \beta J_{ij}\big]\approx\beta J_{ij}f\left(m_i, m_j, \beta J_{i j}\right)
\end{aligned}
$$
并且对 $f(a,b,t)$ 做某种复杂的近似操作，可以得到
$$
m_i=\tanh \left[\sum_{j \in \partial i}\left(\beta J_{i j} m_j-\beta^2 J_{i j}^2\left(1-m_j^2\right) m_i\right)\right]
$$
这个方程称为“TAP方程”，会在之后的一章中详细介绍。如果温度足够高，$\beta^2$ 也作为高阶小量忽视掉，则
$$
m_i=\tanh \left(\beta\sum_{j \in \partial i} J_{i j} m_j\right)
$$
这就是我们在第三小节中用Mean-Field Approximation得到的结果。

从上面可以看到，Bethe Approximation在温度较高时可以退化为Mean-Field Approximation，因为高温实际上意味这自旋变量之间的相互作用减弱。

### 4.5 平均场逆伊辛问题







$ \quad$

**4.4、4.5小节参考文献：**[8]



## 作业

### 1. 关于树状图中Bethe Approximation对应的分布的问题

在本篇第四小节，我们使用对树状图划分区域的方式得到了Bethe自由能
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



### 2. 用Belief Propagation求解SK model中的 $C_{ij}$



### 3. 用Bethe Approximation求解inverse Ising problem





## 参考文献

[1] Yedidia, Jonathan S., William T. Freeman, and Yair Weiss. "[Understanding belief propagation and its generalizations.](https://www.cs.huji.ac.il/course/2005/pmai/tirguls/TR2001-22.pdf)" *Exploring artificial intelligence in the new millennium* 8.236-239 (2003): 0018-9448.

[2] Yedidia, Jonathan S., William T. Freeman, and Yair Weiss. "[Constructing free-energy approximations and generalized belief propagation algorithms.](https://www.cs.princeton.edu/courses/archive/spring06/cos598C/papers/YedidaFreemanWeiss2004.pdf)" *IEEE Transactions on information theory* 51.7 (2005): 2282-2312.

[3] Cowell, Robert. "[Advanced inference in Bayesian networks.](https://link.springer.com/chapter/10.1007/978-94-011-5014-9_2)" *Learning in graphical models*. Springer, Dordrecht, 1998. 27-49.

[4] Heskes, Tom. "[On the uniqueness of loopy belief propagation fixed points.](https://citeseerx.ist.psu.edu/document?repid=rep1&type=pdf&doi=98a0dd903f0eea9fc35bd5f163563252c118ca60)" *Neural Computation* 16.11 (2004): 2379-2413.

[5] Yuille, Alan. "[Belief propagation, mean-field, and bethe approximations.](https://www.cs.jhu.edu/~ayuille/JHUcourses/VisionAsBayesianInference2020/10/BeliefPropagationMFT.pdf)" *Dept. Stat., Univ. California, Los Angeles, Los Angeles, CA, USA, Tech. Rep* 4 (2010).

[6] Yedidia, Jonathan S., William T. Freeman, and Yair Weiss. "[Bethe free energy, Kikuchi approximations, and belief propagation algorithms.](https://merl.com/publications/docs/TR2001-16.pdf)" *Advances in neural information processing systems* 13 (2001): 689.

[7] Eric P. Xing "[Variational Inference: Loopy Belief Propagation](https://www.cs.cmu.edu/~epxing/Class/10708-17/notes-17/10708-scribe-lecture12.pdf)" Probabilistic Graphical Models (10-708, Spring 2017) Lecture notes, School of Computer Science, Carnegie Mellon University

[8] Ricci-Tersenghi, Federico. "[The Bethe approximation for solving the inverse Ising problem: a comparison with other inference methods.](https://arxiv.org/abs/1112.4814)" *Journal of Statistical Mechanics: Theory and Experiment* 2012.08 (2012): P08015.



