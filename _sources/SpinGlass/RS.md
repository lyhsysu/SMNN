# S-K模型的复本对称解

- 作者：李宇豪
- 日期：2023年1月19日
- 更新：2023年1月19日
- 许可：<a rel="license" href="http://creativecommons。org/licenses/by-nc-sa/4。0/">CC BY-NC-SA 4.0</a>

---

S-K模型的哈密顿量写为

$$
    H=-\sum_{i<j} J_{i j} S_i S_j-h \sum_i S_i
$$ (sg1)

``````{margin} 
```{admonition} 注 <font color="blue">1</font>
注意这里的高斯分布含有归一化因子$1/N$，这是为了保证哈密顿量 $\mathcal{O}(N)$，即是一个广延量。
```
``````

其中，$J_{ij}$ 服从高斯分布 $J_{ij}\sim N\left(\frac{J_0}{N},\frac{J^2}{N}\right)$<sup><font color="blue">1</font></sup>。系统的配分函数为

$$
    Z=\mathrm{Tr~}e^{-\beta H}
$$ (sg2)

其中 $\mathrm{Tr~}:=\sum_{\{S_i=\pm 1\}}$，称为 trace operation。利用配分函数 {eq}`sg2` 可以计算自由能：

$$
    F=-T\log Z=-T\log \mathrm{Tr~}e^{-\beta H}
$$ (sg3)

注意到，$F$ 是 $Z$ 的函数，而 $Z$ 是 $\boldsymbol{J}\equiv\{J_{ij}\}$ 的函数。对于一个确定的系统，$\boldsymbol{J}$ 是从高斯分布中取出来的一个固定的集合（这个过程称为“淬火”，quenched）。为了获得最终的自由能，我们需要将 $F$ 对 $\boldsymbol{J}$ 取平均，也即对 $\log Z$ 取平均（称为构型平均，configurational average）：

$$
    [F]=-T[\log Z]=-T\int \prod_{(ij)}\mathrm{d}J_{ij}P(J_{ij})\log Z
$$ (sg4)

## 1 复本技巧 Replica trick

由于 $\log Z$ 对 $\boldsymbol{J}$ 的依赖关系是十分复杂的，我们无法直接求 $\log Z$。复本技巧（replica trick）是一个可以将对 $\log Z$ 的构型平均转移到对 $Z^n$ 的构型平均的技巧。

```{admonition} replica trick
$$
[\log Z]=\lim_{n\to 0}\frac{[Z^n]-1}{n}    
$$
这个式子的意义是，准备一个系统的 $n$ 份复本，计算其配分函数的乘积的构型平均，然后取极限 $n\to 0$。
```

利用复本技巧，我们可以写出配分函数的复本平均（replica average）：

$$
\left[Z^n\right]=\int \prod_{i<j} d J_{i j} P\left(J_{i j}\right)\mathrm{~~Tr~} \exp \left(\beta \sum_{i<j} J_{i j} \sum_{\alpha=1}^n S_i^\alpha S_j^\alpha+\beta h \sum_{i=1}^N \sum_{\alpha=1}^n S_i^\alpha\right)
$$ (sg5)

其中，$\mathrm{Tr~}=\sum_{\left\{\mathbf{s}^\alpha, \mathbf{s}^\beta, \ldots, \mathbf{s}^n\right\}}$。将 $P\left(J_{i j}\right)$ 的表达式

$$
    P\left(J_{i j}\right)=\frac{1}{J} \sqrt{\frac{N}{2 \pi}} \exp \left\{-\frac{N}{2 J^2}\left(J_{i j}-\frac{J_0}{N}\right)^2\right\}
$$ (sg6)

代入上式，然后进行一个看不懂的积分，得

$$
\left[Z^n\right]=C_1 \mathrm{~~Tr~} \exp \left\{\frac{1}{N} \sum_{i<j}\left(\frac{1}{2} \beta^2 J^2 \sum_{\alpha, \beta} S_i^\alpha S_j^\alpha S_i^\beta S_j^\beta+\beta J_0 \sum_\alpha S_i^\alpha S_j^\alpha\right)+\beta h \sum_{i=1}^N \sum_{\alpha=1}^n S_i^\alpha\right\}
$$ (sg7)

利用如下技巧

$$
    \sum_{\alpha, \beta} S_i^\alpha S_j^\alpha S_i^\beta S_j^\beta=2 \sum_{\alpha<\beta} S_i^\alpha S_j^\alpha S_i^\beta S_j^\beta+\sum_\alpha\left(S_i^\alpha S_j^\alpha\right)^2=2 \sum_{\alpha<\beta} S_i^\alpha S_j^\alpha S_i^\beta S_j^\beta+n
$$ (sg7_5)

式 {eq}`sg7` 的指数部分改写为

$$
\exp \left\{\frac{1}{N} \sum_{i<j}\left(\frac{1}{2} \beta^2 J^2\left(2 \sum_{\alpha<\beta} S_i^\alpha S_j^\alpha S_i^\beta S_j^\beta+n\right)+\beta J_0 \sum_\alpha S_i^\alpha S_j^\alpha\right)+\beta h \sum_{i=1}^N \sum_{\alpha=1}^n S_i^\alpha\right\}
$$ (sg8)

注意到指数上第一项中的求和号中含有与 $i$、$j$ 无关的项 $\frac{1}{2} \beta^2 J^2n$，并且求和 $\sum_{i<j}$ 共有 $N(N-1)/2$ 项，因此将常数项提到前面，改写为

$$
\begin{aligned}
\left[Z^n\right]=&C_1 \exp \left(\frac{(N-1) \beta^2 J^2 n}{4}\right) \\
&\mathrm{Tr~} \exp \left\{\frac{1}{N} \sum_{i<j}\left(\beta^2 J^2 \sum_{\alpha<\beta} S_i^\alpha S_j^\alpha S_i^\beta S_j^\beta+\beta J_0 \sum_\alpha S_i^\alpha S_j^\alpha\right)+\beta h \sum_{i=1}^N \sum_{\alpha=1}^n S_i^\alpha\right\}
\end{aligned}
$$ (sg9)

并且，根据一个技巧 $\left(\sum_iA_i\right)^2=\sum_iA_i^2+\sum_{i\neq j}A_iA_j=\sum_iA_i^2+2\sum_{i< j}A_iA_j$，即

$$
    \sum_{i< j}A_iA_j=\frac{1}{2}\left(\left(\sum_iA_i\right)^2-\sum_iA_i^2\right)
$$ (sg10)

可以分别改写

$$
\begin{aligned}
\frac{1}{N} \sum_{i<j} \beta^2 J^2 \sum_{\alpha<\beta} S_i^\alpha S_j^\alpha S_i^\beta S_j^\beta & =\frac{\beta^2 J^2}{2 N}\left(\sum_{\alpha<\beta}\left(\sum_i S_i^\alpha S_i^\beta\right)^2-\sum_i \sum_{\alpha<\beta}  (S_i^\alpha)^2 (S_i^\beta)^2\right) \\
& =\frac{\beta^2 J^2}{2 N} \sum_{\alpha<\beta}\left(\sum_i S_i^\alpha S_i^\beta\right)^2-\frac{\beta^2 J^2}{2 N} \sum_i \sum_{\alpha<\beta} 1
\end{aligned}
$$ (sg11)

和

$$
\begin{aligned}
\frac{\beta J_0}{N} \sum_{i<j} \sum_\alpha S_i^\alpha S_j^\alpha&=\frac{\beta J_0}{2 N}\left(\sum_\alpha\left(\sum_i S_i^\alpha\right)^2 -\sum_i \sum_{\alpha}  (S_i^\alpha)^2 \right) \\
&=\frac{\beta J_0}{2 N} \sum_\alpha\left(\sum_i S_i^\alpha\right)^2-\frac{\beta J_0}{2 N} \sum_i \sum_\alpha 1
\end{aligned}
$$ (sg12)

忽略掉上两式中的常数项，并且根据大 $N$ 极限下的近似 $\exp \left(\frac{(N-1) \beta^2 J^2 n}{4}\right)\approx \exp \left(\frac{N \beta^2 J^2 n}{4}\right)$，式 {eq}`sg9` 改写为

$$
\begin{aligned}
\left[Z^n\right]&=C_2 \exp \left(\frac{N \beta^2 J^2 n}{4}\right)\\
&=\mathrm{Tr~} \exp \left\{\frac{\beta^2 J^2}{2 N} \sum_{\alpha<\beta}\left(\sum_i S_i^\alpha S_i^\beta\right)^2+\frac{\beta J_0}{2 N} \sum_\alpha\left(\sum_i S_i^\alpha\right)^2+\beta h \sum_{i=1}^N \sum_{\alpha=1}^n S_i^\alpha\right\}
\end{aligned}
$$ (sg13)

再利用如下变换将 $\left(\sum_i S_i^\alpha S_i^\beta\right)^2$ 和 $\left(\sum_i S_i^\alpha\right)^2$ 线性化。

```{admonition} Hubbard-Stratonovich transform
$$
\exp \left(\frac{y^2}{2}\right)=\int_{-\infty}^{\infty} \frac{\mathrm{d} x}{\sqrt{2 \pi}} \exp \left(-\frac{x^2}{2}\right) \exp (x y)
$$
```

令 $x:=\beta J\sqrt{N}q_{\alpha\beta}$，$y:=\beta J\sum_i S_i^\alpha S_i^\beta/\sqrt{N}$，有

$$
\exp \frac{\beta^2 J^2}{2 N}\left(\sum_i S_i^\alpha S_i^\beta\right)^2  =\int_{-\infty}^{\infty} \frac{\mathrm{d} q_{\alpha \beta}}{\sqrt{2 \pi}} \exp \left(-\beta^2 J^2 N \frac{q_{\alpha \beta}^2}{2}+\beta^2 J^2 q_{\alpha \beta} \sum_i S_i^\alpha S_i^\beta\right)
$$ (sg14)

令 $x:=\sqrt{\beta JN}m_{\alpha}$，$y:=\sqrt{\beta J/N}\sum_i S_i^\alpha$，有

$$
\exp \frac{\beta J_0}{2 N}\left(\sum_i S_i^\alpha\right)^2 =\int_{-\infty}^{\infty} \frac{\mathrm{d} m_\alpha}{\sqrt{2 \pi}} \exp \left(-\beta J_0 N m_\alpha^2+\beta J_0 m_\alpha \sum_i S_i^\alpha\right)
$$ (sg15)

其中，忽略了积分号前面的常数（由 $\mathrm{d}x$ 带来的）。式 {eq}`sg13` 改写为

$$
\begin{aligned}
\left[Z^n\right]=&\;C_3 \exp \left(\frac{N \beta^2 J^2 n}{4}\right) \int_{-\infty}^{\infty} \prod_{\alpha<\beta} \mathrm{d}q_{\alpha \beta} \prod_\alpha \mathrm{d}m_\alpha\\
&\mathrm{Tr~} \exp \left(-\frac{\beta^2 J^2 N}{2} \sum_{\alpha<\beta} q_{\alpha \beta}^2-\frac{\beta J_0 N}{2} \sum_\alpha m_\alpha^2\right) \\
&\exp \left(\beta^2 J^2 \sum_{\alpha<\beta} q_{\alpha \beta} \sum_i S_i^\alpha S_i^\beta+\beta J_0 \sum_\alpha m_\alpha \sum_i S_i^\alpha\right) \exp \left(\beta h \sum_{i=1}^N \sum_{\alpha=1}^n S_i^\alpha\right)\\
=&\;C_3 \exp \left(\frac{N \beta^2 J^2 n}{4}\right) \int_{-\infty}^{\infty} \prod_{\alpha<\beta} \mathrm{d}q_{\alpha \beta} \prod_\alpha \mathrm{d}m_\alpha\\
&\;\exp \left(-\frac{\beta^2 J^2 N}{2} \sum_{\alpha<\beta} q_{\alpha \beta}^2-\frac{\beta J_0 N}{2} \sum_\alpha m_\alpha^2\right) \\
&\;\mathrm{Tr~}\exp \left(\beta^2 J^2 \sum_{\alpha<\beta} q_{\alpha \beta} \sum_i S_i^\alpha S_i^\beta+\beta \sum_\alpha\left(J_0 m_\alpha+h\right) \sum_i S_i^\alpha\right)
\end{aligned}
$$ (sg16)

将最后一项指数上的求和号拿下来成为连乘号，即

$$
    \prod_{i=1}^N \mathrm{Tr~}\exp \left(\beta^2 J^2 \sum_{\alpha<\beta} q_{\alpha \beta} S_i^\alpha S_i^\beta+\beta \sum_\alpha\left(J_0 m_\alpha+h\right)  S_i^\alpha\right)
$$ (sg17)

令

$$
    L:=\beta^2 J^2 \sum_{\alpha<\beta} q_{\alpha \beta} S_i^\alpha S_i^\beta+\beta \sum_\alpha\left(J_0 m_\alpha+h\right)  S_i^\alpha
$$ (sg18)

考虑到连乘中的每一项都是对相同的变量 $\left\{\mathbf{s}^\alpha, \mathbf{s}^\beta, \ldots, \mathbf{s}^n\right\}$ 的相同取值 $\{\pm 1\}$ 进行的，与 $i$ 无关，因此

$$
    \prod_{i=1}^N \mathrm{Tr~}\exp \left(L\right)=\Big(\mathrm{Tr~}\exp \left(L\right)\Big)^N=\exp \Big\{ N\log\Big(\mathrm{Tr~}\exp \left(L\right)\Big) \Big\}
$$ (sg19)

式 {eq}`sg16` 改写为

$$
\begin{aligned}
\left[Z^n\right] =&C_4 \exp \left(\frac{N \beta^2 J^2 n}{4}\right) \int_{-\infty}^{\infty} \prod_{\alpha<\beta} d q_{\alpha \beta} \prod_\alpha d m_\alpha \\
& \exp \left\{-\frac{\beta^2 J^2 N}{2} \sum_{\alpha<\beta} q_{\alpha \beta}^2-\frac{\beta J_0 N}{2} \sum_\alpha m_\alpha^2+N\log\Big(\mathrm{Tr~}\exp \left(L\right)\Big)\right\}
\end{aligned}
$$ (sg20)

接下来用鞍点近似法（method of steepest descent）来考虑这个积分

```{admonition} method of steepest descent
考虑如下积分在 $\lambda$ 为正的大实数时的渐进行为

$$
    I(\lambda) =\int_C\mathrm{d}z e^{\lambda g(z)}
$$

可以用下式近似作为积分结果

$$
I(\lambda)=e^{i \theta} \sqrt{\frac{2 \pi}{\lambda\left|g^{\prime \prime}\left(z_0\right)\right|}} e^{\lambda g\left(z_0\right)}[1+\mathcal{O}(1 / \lambda)]
$$

其中，$z_0$ 是函数的鞍点
```

式 {eq}`sg20` 的积分中，最后一项指数上与 $N$ 成正比，在大 $N$极限下，积分结果由指数上的宗量的最大值决定。即

$$
\left[Z^n\right]=C_4 \exp \left\{\frac{N \beta^2 J^2 n}{4}-\frac{\beta^2 J^2 N}{2} \sum_{\alpha<\beta}\left(q_{\alpha \beta}^{\star}\right)^2-\frac{\beta J_0 N}{2} \sum_\alpha\left(m_\alpha^{\star}\right)^2+N \log \mathrm{Tr~} \exp (L)\right\}
$$ (sg21)

其中，$\left\{q_{\alpha \beta}^{\star}, m_\alpha^{\star}\right\}=\arg \max _{\left\{q_{\alpha \beta}, m_\alpha\right\}} E$

$$
\begin{aligned}
E:=-\frac{\beta^2 J^2 N}{2} \sum_{\alpha<\beta} q_{\alpha \beta}^2-\frac{\beta J_0 N}{2} \sum_\alpha m_\alpha^2+N \log \mathrm{Tr~}\exp (L) 
\end{aligned}
$$ (sg22)

令 $\frac{\partial}{\partial q_{\alpha\beta}}E=0$，得

$$
q_{\alpha \beta}^{\star}=\frac{1}{\beta^2 J^2} \frac{\partial}{\partial q_{\alpha \beta}} \log \mathrm{Tr~} \exp (L)=\frac{1}{\beta^2 J^2} \frac{\mathrm{Tr~} \exp (L) \beta^2 J^2}{\mathrm{Tr~}\exp (L)} S^\alpha S^\beta=\frac{\mathrm{Tr~} S^\alpha S^\beta \exp (L)}{\mathrm{Tr~}\exp (L)}
$$ (sg23)

令 $\frac{\partial}{\partial m_{\alpha}}E=0$，得

$$
m_\alpha^{\star}=\frac{1}{\beta J_0} \frac{\partial}{\partial m_\alpha} \log \mathrm{Tr~}  \exp (L)=\frac{1}{\beta J_0}\frac{\mathrm{Tr~}  \exp (L)\beta J_0}{\mathrm{Tr~}  \exp (L)} S^\alpha=\frac{\mathrm{Tr~}  S^\alpha \exp (L)}{\mathrm{Tr~}  \exp (L)}
$$ (sg24)

再将式 {eq}`sg21` 改写为

$$
\left[Z^n\right]=C_4 \exp \left\{N n\left(\frac{\beta^2 J^2}{4}-\frac{\beta^2 J^2}{2 n} \sum_{\alpha<\beta}\left(q_{\alpha \beta}^{\star}\right)^2-\frac{\beta J_0}{2 n} \sum_\alpha\left(m_\alpha^{\star}\right)^2+\frac{1}{n} \log\mathrm{Tr~} \exp (L)\right)\right\}
$$ (sg24_5)

``````{margin} 
```{admonition} 注 <font color="blue">2</font>
这里其实是一个“$0\cdot\infty$”极限未定式，应当讨论无穷的阶数。但是现在先姑且认为$0\cdot\infty=0$
```
``````
考虑 $n\to 0$ 的极限，此时指数上的宗量趋于零<font color="blue">2</font></sup>，利用麦克劳林展开（Maclaurin's series）得

$$
\left[Z^n\right] \approx 1+N n\left\{\frac{\beta^2 J^2}{4}-\frac{\beta^2 J^2}{2 n} \sum_{\alpha<\beta}\left(q_{\alpha \beta}^{\star}\right)^2-\frac{\beta J_0}{2 n} \sum_\alpha\left(m_\alpha^{\star}\right)^2+\frac{1}{n} \log \mathrm{Tr~} \exp (L)\right\}
$$ (sg25)

至此，我们要求的构型平均 $\left[\log Z\right] $可以写为

$$
\begin{aligned}
\left[\log Z\right] & = N\lim _{n \rightarrow 0} \frac{\left[Z^n\right]-1}{N n} \\
& = N\lim _{n \rightarrow 0}\left\{\frac{\beta^2 J^2}{4}-\frac{\beta^2 J^2}{2 n} \sum_{\alpha<\beta}\left(q_{\alpha \beta}^{\star}\right)^2-\frac{\beta J_0}{2 n} \sum_\alpha\left(m_\alpha^{\star}\right)^2+\frac{1}{n} \log \mathrm{Tr~}  \exp (L)\right\}\\
\end{aligned}
$$ (sg26)

## 2 复本对称解

要进一步求解这个问题，就需要考虑参数 $q_{\alpha\beta}$、$m_\alpha$ 对 $\alpha$、$\beta$ 的依赖关系。一个朴素的假设是，各个复本都是等价的，$q_{\alpha\beta}$、$m_\alpha$ 与 $\alpha$、$\beta$ 无关，即 $\forall \alpha, \beta: q_{\alpha \beta}=q, m_\alpha=m$，这就是所谓的复本对称假设。

根据这个假设，式 {eq}`sg26` 可以写为

 $$
\begin{aligned}
[\log Z] & =N\lim _{n \rightarrow 0}\left\{\frac{\beta^2 J^2}{4}-\frac{\beta^2 J^2(n-1)}{4} q^2-\frac{\beta J_0}{2} m^2+\frac{1}{n} \log \mathrm{Tr~} \exp \left(L*\right)\right\} \\
& =N\left\{\frac{\beta^2 J^2}{4}\left(1+q^2\right)-\frac{\beta J_0}{2} m^2+\lim _{n \rightarrow 0} \frac{1}{n} \log \mathrm{Tr~} \exp \left(L*\right)\right\}
\end{aligned}
 $$ (sg27)

其中第一步用到了求和 $\sum_{\alpha<\beta}$ 共有 $n(n-1)/2$ 项，并且

$$
    L*:=L(q_{\alpha \beta}=q, m_\alpha=m)=\beta^2 J^2 q \sum_{\alpha<\beta} S^\alpha S^\beta+\beta\left(J_0 m+h\right) \sum_\alpha S^\alpha
$$ (sg28)

考虑式 {eq}`sg27` 的最后一项

$$
\begin{aligned}
&\frac{1}{n} \log \mathrm{Tr~} \exp \left(L*\right)\\
&=\frac{1}{n} \log \mathrm{Tr~}\exp \left(\beta^2 J^2 q \sum_{\alpha<\beta} S^\alpha S^\beta+\beta\left(J_0 m+h\right) \sum_\alpha S^\alpha\right)\\
&=\frac{1}{n} \log \mathrm{Tr~}\exp \left\{\frac{1}{2}\beta^2 J^2 q \left(\sum_{\alpha} S^\alpha \right)^2-\frac{1}{2}\beta^2 J^2 qn+\beta\left(J_0 m+h\right) \sum_\alpha S^\alpha\right\}\\
&=\frac{1}{n} \log \left\{ \exp\left( -\frac{\beta^2 J^2 qn}{2} \right) \mathrm{Tr~}\exp 
\left[ \frac{1}{2}\beta^2 J^2 q \left(\sum_{\alpha} S^\alpha \right)^2+\beta\left(J_0 m+h\right) \sum_\alpha S^\alpha \right] \right\}
\end{aligned}
$$ (sg29)

对于 $\exp \left[ \frac{1}{2}\beta^2 J^2 q \left(\sum_{\alpha} S^\alpha \right)^2\right]$，再次使用Hubbard-Stratonovich变换对其线性化：令 $x:=\hat{z}\beta J\sqrt{q}$，$y:=\beta J\sqrt{q}\sum_\alpha S^\alpha$，有

 $$
    \exp \left[ \frac{1}{2}\beta^2 J^2 q \left(\sum_{\alpha} S^\alpha \right)^2\right]=\int_C \mathrm{d}\hat{z}\sqrt{\frac{\beta^2 J^2 q}{2\pi}}\exp\left(-\frac{\hat{z}^2}{2}\beta^2 J^2 q\right)\exp\left(\beta^2 J^2 q\hat{z}\sum_{\alpha} S^\alpha\right)
 $$ (sg30)

注意到 

 $$
    \mathrm{d}\hat{z}\sqrt{\frac{\beta^2 J^2 q}{2\pi}}\exp\left(-\frac{\hat{z}^2}{2}\beta^2 J^2 q\right)=\mathrm{d}\hat{z}p(\hat{z}),\quad p(\hat{z}):=\sqrt{\frac{\beta^2 J^2 q}{2\pi}}\exp\left(-\frac{\hat{z}^2}{2}\beta^2 J^2 q\right)
 $$ (sg31)

即 $\hat{z}$ 是一个服从 $p(\hat{z})$ 的高斯变量，将其变换为标准高斯，即令 $z=\sqrt{\beta^2 J^2 q}\hat{z}$，则 $z\sim\mathcal{N}(0,1)$。式 {eq}`sg30` 改写为

 $$
\begin{aligned}
\exp \left[ \frac{1}{2}\beta^2 J^2 q \left(\sum_{\alpha} S^\alpha \right)^2\right]&=\int_C \frac{\mathrm{d}z}{\sqrt{2\pi}} \exp\left(-\frac{z^2}{2}\right)\exp\left(\beta J\sqrt{q}z\sum_{\alpha} S^\alpha\right)\\
&=\int_C Dz\exp\left(\beta J\sqrt{q}z\sum_{\alpha} S^\alpha\right)
\end{aligned}
 $$ (sg32)

其中，$Dz=\mathrm{d}z\exp\left(-z^2/2\right)/\sqrt{2\pi} $ 是标准高斯测度。因此式 {eq}`sg29` 写为

 $$
\begin{aligned}
&\frac{1}{n} \log \mathrm{Tr~} \exp \left(L*\right)\\
&=\frac{1}{n} \log \left\{ \exp\left( -\frac{\beta^2 J^2 qn}{2} \right) \mathrm{Tr~}\exp 
\left[ \frac{1}{2}\beta^2 J^2 q \left(\sum_{\alpha} S^\alpha \right)^2+\beta\left(J_0 m+h\right) \sum_\alpha S^\alpha \right] \right\}\\
&=\frac{1}{n} \log \left\{ \exp\left( -\frac{\beta^2 J^2 qn}{2} \right) \mathrm{Tr~}\int_C Dz\exp\left(\beta J\sqrt{q}z\sum_{\alpha} S^\alpha+\beta\left(J_0 m+h\right) \sum_\alpha S^\alpha\right)\right\}
\end{aligned}
 $$ (sg33)

对于上式最后一项，有

 $$
\begin{aligned}
&\mathrm{Tr~}\int_C Dz\exp\left(\beta J\sqrt{q}z\sum_{\alpha} S^\alpha+\beta\left(J_0 m+h\right) \sum_\alpha S^\alpha\right)\\
&=\int_C Dz\mathrm{~Tr~}\exp\left[\sum_{\alpha} S^\alpha\Big(\beta J\sqrt{q}z+\beta\left(J_0 m+h\right)\Big)\right]\\
&=\int_C Dz \prod_{\alpha=1}^n\mathrm{Tr~}\exp\left[S^\alpha\Big(\beta J\sqrt{q}z+\beta\left(J_0 m+h\right)\Big)\right]\\
&=\int_C Dz \prod_{\alpha=1}^n\mathrm{Tr~}\exp\left[S^\alpha\beta\hat{H}(z)\right]
\end{aligned}
 $$ (sg34)

其中，$\hat{H}(z):=J\sqrt{q}z+\left(J_0 m+h\right)$。注意到上式最后一步连乘中的每一项都是同样对 $\{ S^\alpha =\pm 1 \}$ 求和，因此连乘可以转为指数

 $$
\begin{aligned}
\prod_{\alpha=1}^n\mathrm{Tr~}\exp\left[S^\alpha\beta\hat{H}(z)\right]&=\prod_{\alpha=1}^n \left\{ \exp\Big(-\beta\hat{H}(z)\Big) +\exp\Big(\beta\hat{H}(z)\Big) \right\}\\
&=\left\{ 2\cosh\Big(\beta\hat{H}(z)\Big) \right\}^n=\exp\left\{n\log \left[2\cosh\Big(\beta\hat{H}(z)\Big) \right]\right\}
\end{aligned}
 $$ (sg35)

式 {eq}`sg29` 最终写为
 $$
\begin{aligned}
&\frac{1}{n} \log \mathrm{Tr~} \exp \left(L*\right)\\
&=\frac{1}{n} \log \left\{ \exp\left( -\frac{\beta^2 J^2 qn}{2} \right) \int Dz \exp\left\{n\log \left[2\cosh\Big(\beta\hat{H}(z)\Big) \right] \right\}\right\}\\
&=\frac{1}{n} \log \int Dz~ \exp\left\{n\log \left[2\cosh\Big(\beta\hat{H}(z)\Big)\right]  -\frac{n}{2}\beta^2 J^2 q \right\}
\end{aligned}
 $$ (sg36)

当 $n\to 0$时，利用麦克劳林展开

 $$
\begin{aligned}
\frac{1}{n} \log \mathrm{Tr~} \exp \left(L*\right) &\approx \frac{1}{n} \log \left\{1+n \int D z~ \log\left[ 2\cosh\Big(\beta\hat{H}(z)\Big)\right]-\frac{n}{2} \beta^2 J^2 q\int D z\right\}\\
&\approx \int D z~ \log\left[ 2 \cosh \Big(\beta\hat{H}(z)\Big)\right] -\frac{1}{2} \beta^2 J^2 q
\end{aligned}
 $$ (sg37)

其中，第一步是利用指数函数在 $0$ 附近的展开，第二步是利用对数函数在 $1$ 附近的展开。由此，式 {eq}`sg27` 改写为

 $$
\begin{aligned}
[\log Z] &=N\left\{\frac{\beta^2 J^2}{4}\left(1+q^2\right)-\frac{\beta J_0}{2} m^2+\lim _{n \rightarrow 0} \left\{\int D z~ \log\left[ 2 \cosh \Big(\beta\hat{H}(z)\Big)\right] -\frac{1}{2} \beta^2 J^2 q\right\}\right\}\\
&=N\left\{\frac{\beta^2 J^2}{4}\left(1-q^2\right)-\frac{\beta J_0}{2} m^2+\int D z~ \log\left[ 2 \cosh \Big(\beta\hat{H}(z)\Big)\right]\right\}
\end{aligned}
 $$ (sg38)

根据

 $$
\begin{aligned}
&\frac{\partial}{\partial m}[\log Z]= -\beta J_0m+\int \mathrm{D} z(\tanh \beta \tilde{H}(z)) \cdot \beta J_0=0\\
&\frac{\partial}{\partial q}[\log Z]= \frac{\beta^2 J^2}{2}(q-1)+\int \mathrm{D} z(\tanh \beta \tilde{H}(z)) \cdot \frac{\beta J}{2 \sqrt{q}} z=0
\end{aligned}
 $$ (sg38_5)

得

 $$
\begin{aligned}
&m=\int D z \tanh \beta \tilde{H}(z)\\
&q=1-\int D z \operatorname{sech}^2 \beta \tilde{H}(z)=\int D z \tanh ^2 \beta \tilde{H}(z)
\end{aligned}
 $$ (sg39)