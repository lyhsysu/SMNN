# Sec 2. 平均场近似

平均场近似（Mean-Field Approximation）是假设该分布**由互相独立的自旋格点 $x_i\in\{+1,-1\}$ 组成**，其均值（即磁化强度）为

$$
m_i=(+1)\times \frac{1+m_{i}(+1)}{2}+(-1)\times \frac{1+m_{i}(-1)}{2}
$$ (b8)

因此可以用 $m_i$ 将单个格点的分布 $b_i$ 参数化：（这在上一页中提到过）

$$
b_i(x_i)=\frac{1+m_ix_i}{2}
$$ (b9)

由于各个格点是相互独立的，系统的总分布就等于各个格点的分布相乘，因此系统的平均场近似分布为

$$
b_{MF}=\prod_ib_i(x_i)=\prod_i\frac{1+m_ix_i}{2}
$$ (b10)

在此近似下，可以计算变分内能、变分熵和变分自由能（见 **推导细节 b**）：

$$
\begin{aligned}
& U_{M F}= -\sum_a J_a \left[\sum_{\boldsymbol{x}} \prod_i b_i\left(x_i\right)\prod_{i \in \partial a} x_i \right]\\
&H_{M F}=-\sum_i \sum_{x_i} \frac{1+m_i x_i}{2} \ln \frac{1+m_i x_i}{2}\\
&F_{M F} =U_{M F}-H_{M F}=-\sum_a J_a \prod_{i \in \partial a} m_i+\sum_i \sum_{x_i} \frac{1+m_i x_i}{2} \ln \frac{1+m_i x_i}{2}
\end{aligned}
$$ (b11)

```{admonition} 推导细节 b
:class: dropdown
变分内能

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

变分熵

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

其中用到了

$$
\sum_{x \backslash x_i} \prod_{j(\neq i)} \frac{1+m_j x_j}{2}=\sum_{x \backslash x_i}b_{MF\backslash x_i}=1
$$

即去掉格点 $i$ 后的系统的总分布，由于各格点相互独立，此时系统的分布仍是归一化的。$S_i$ 是单个自旋变量的熵

$$
S_i=-k_B\sum_ip_i\ln p_i=- \sum_{x_i} \frac{1+m_i x_i}{2} \ln \frac{1+m_i x_i}{2}
$$

式中，已令 $\beta=1$，因此 $k_B=\frac{1}{\beta}=1$，$p_i$ 在此处是平均场近似下的 $b_i$。

因此变分自由能为

$$
\begin{aligned}
F_{M F} & =U_{M F}-H_{M F} \\
& =-\sum_a J_a \prod_{i \in \partial a} m_i+\sum_i \sum_{x_i} \frac{1+m_i x_i}{2} \ln \frac{1+m_i x_i}{2}
\end{aligned}
$$
```

令变分自由能取极值

$$
\frac{\partial F_{M F}}{\partial m_i} =0
$$ (b12)

解得

$$
m_i=\tanh \left(\sum_{a \in \partial i} J_a \prod_{j \in \partial a \backslash i} m_j\right)
$$ (b13)

```{admonition} 推导细节 c
:class: dropdown

$$
\begin{aligned}
\frac{\partial F_{M F}}{\partial m_i} & =-\sum_a J_a \prod_{j \in \partial a \backslash i} m_j+\sum_{x_i} \frac{x_i}{2} \ln \frac{1+m_i x_i}{2}+\frac{x_i}{2} \\
& =-\sum_a J_a \prod_{j \in \partial a \backslash i} m_j+\frac{1}{2} \ln \frac{1+m_i}{1-m_i}
\end{aligned}
$$

式 {eq}`b13` 即

$$
\sum_a J_a \prod_{j \in \partial a \backslash i} m_j+\frac{1}{2} \ln \frac{1+m_i}{1-m_i} = 0
$$

将其改写为

$$
\ln \frac{1+m_i}{1-m_i}=2\sum_a J_a \prod_{j \in \partial a \backslash i} m_j := 2A
$$

由 $\ln \frac{1+m_i}{1-m_i}=2A$ 得

$$
m_i=\frac{e^{2A}-1}{e^{2A}+1}=\frac{e^{A}-e^{-A}}{e^{A}+e^{-A}}=\tanh A
$$

再代入 $\sum_a J_a \prod_{j \in \partial a \backslash i} m_j = A$ 得

$$
m_i=\tanh \left(\sum_{a \in \partial i} J_a \prod_{j \in \partial a \backslash i} m_j\right)
$$
```

这是一个迭代方程，将其迭代到不动点即可解出 $m_i$ 并找到此时对应的真实分布。