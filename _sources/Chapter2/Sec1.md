# Sec 1. 变分方法与变分自由能

``````{margin} 
```{admonition} 注 <font color="blue">1</font>
这里采用热力学的写法，将能量写为 $E(x)$。
```
``````

对于一个真实分布<sup><font color="blue">1</font></sup>

$$
P(\boldsymbol{x})=\frac{e^{-\beta E(x)}}{Z} \qquad Z=\sum_x e^{-\beta E(x)} \qquad E(\boldsymbol{x})=-\sum_a J_a \prod_{i \in \partial a} x_i
$$ (b1)

根据 $Z$ 可以计算 Helmholtz 自由能 $F_H=-\ln Z$，但是 $Z$ 是 $\mathcal{O}(N^2)$，很难直接计算。本章通过构造一个Gibbs自由能 $F(b)$（其中 $b(x)$ 是真实分布，$F(b)$ 是泛函）并使用变分方法使其逼近真实的Helmholtz自由能 $F_H=-\ln Z$。

```{prf:theorem} 最大熵原理
:nonumber: 
随机变量的概率分布是很难测定的。给定分布的约束（如期望、方差等可测得值）的情况下，在所有可能的分布中，熵最大的分布代表着信息量最大的分布，是最好的、最接近真实分布的。
```
说明如下：

对于一个未知分布 $P_r$，其熵定义为

$$
S=-k\sum_{r}P_r\ln P_r
$$ (b2)

该分布要满足归一化和能量守恒（系统的总能量为定值）的约束，即

$$
\sum_r P_r=1 \qquad \sum_r E_r P_r=\mu
$$ (b3)

使用拉格朗日乘子法求解条件极值（见 **推导细节 a**），可以得到

$$
P_r=\frac{\exp\left(-\beta E_r\right)}{\sum_r \exp\left(-\beta E_r\right)}
$$ (b4)

这正是真实分布中 $P(\boldsymbol{x})$ 的表达式。

```{admonition} 推导细节 a
:class: dropdown
构造Lagrangian

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
```

对于真实分布 $b(x)$，Gibbs 自由能 $F(b)$ 是 $b(x)$ 的泛函，并且设 $\beta=1$。在热力学中，Gibbs 自由能写为

$$
F(b)=U(b)-H(b)
$$ (b5)

其中，$U(b)$ 是内能；$H(b)$ 是熵（而不是哈密顿量），他们分别写为

$$
\begin{aligned}
    &U(b)=\sum_x b(x) E(x) \\
    &H(b)=-\sum_x b(x) \ln b(x)
\end{aligned}
$$ (b6)

所谓的变分自由能就是指我们基于真实分布 $b(x)$ 构造的这个 Gibbs 自由能 $F(b)$，因为是用变分的方法去逼近真实的 Helmholtz 自由能，所以又称为变分自由能。对 Gibbs 自由能 $F(b)$ 和 Helmholtz 自由能 $F_H$ 作差，得

$$
\begin{aligned}
F(b)-F_H &=\sum_{\boldsymbol{x}} b(\boldsymbol{x}) E(\boldsymbol{x})+\sum_{\boldsymbol{x}} b(\boldsymbol{x}) \ln b(\boldsymbol{x})+\ln Z \\
&=\sum_{\boldsymbol{x}} b(\boldsymbol{x})(-\ln Z-\ln P(\boldsymbol{x}))+\sum_{\boldsymbol{x}} b(\boldsymbol{x}) \ln b(\boldsymbol{x})+\ln Z \\
&=\sum_{\boldsymbol{x}} b(\boldsymbol{x}) \ln \frac{b(\boldsymbol{x})}{P(\boldsymbol{x})}=D_{KL}(b \| P)
\end{aligned}
$$ (b7)

其中，第二步用到了由 $P(\boldsymbol{x})=\frac{e^{- E(\boldsymbol{x})}}{Z}$ 得到的 $E(\boldsymbol{x})=-\ln Z-\ln P(\boldsymbol{x})$。可以看到，二者之差为一个KL散度的形式。

````{prf:definition} Kullback-Leibler 散度
:nonumber: 
KL散度是用来度量两个概率分布相似度的指标。假设对随机变量 $\xi$ 存在两个概率分布 $P$、$Q$，定义从 $P$ 到 $Q$ 的KL散度为

$$
\mathbb{D}_{\mathrm{KL}}(P \| Q)=\sum_i P(i) \ln \left(\frac{P(i)}{Q(i)}\right) \quad\text{或}\quad \int_{-\infty}^{\infty} p(\mathbf{x}) \ln \left(\frac{p(\mathbf{x})}{q(\mathbf{x})}\right) \mathrm{~d} \mathbf{x}
$$

可以证明，KL散度具有非负性。以离散形式为例：根据Jensen's不等式

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
````

根据KL散度的非负性，可以知道 $F(b)-F_H=D_{KL}(b \| P)\geq 0$，即变分自由能的下限为Helmholtz自由能，此时分布 $b(x)=P(x)$。为了求解自由能，我们首先要给定 $b(x)$ 的一个具体形式。接下来两节分别给出了 $b(x)$ 的两种近似：平均场和 Bethe 近似。