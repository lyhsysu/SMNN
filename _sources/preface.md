<style>
h2 { font: 24px !important; }
h3 { font: 20px !important; }
p { font: 16px !important; }
</style>

# Introduction

<br>

&emsp;&emsp;二级相变又称为连续相变，是一种十分常见的现象，也是近代凝聚态物理和统计物理研究的重要内容。常见的相变有气—液相变、铁磁—顺磁相变、BCS 超导相变等。在这些相变过程中，外界条件（如温度）是连续变化的，但是物质的状态却可以发生突然的改变，导致某些物理量发生不连续的变化或者发散。1933年，朗道（Landau）提出**序参量**（order parameter）的概念，建立了朗道相变理论{cite}`landau1937theory1,landau1937theory2`。根据朗道的理论，二级相变的微观本质是系统的对称性自发破缺，系统的对称性破缺产生了某种”序“，序参量（即描述”序“的参量）就是一个描述系统对称性的变量。序参量为零对应着系统的无序状态，非零对应有序状态，在发生二级相变的临界点，序参量连续地从零变为非零。

&emsp;&emsp;玻璃化是另一种常见的相变，比如高温下熔融的二氧化硅液体经过淬火形成石英玻璃，这是一种宏观上有序、微观上分子排列高度无序的物质结构。对于一般的铁磁体材料，当温度升高时，原子自身的热运动逐渐超过原子之间的相互作用，物质宏观上变为顺磁性；当温度重新降低至一定水平（称为**居里点**）后，原子之间的铁磁相互作用使得所有的磁矩都按同一个方向排列（即**长程有序状态**），物质快速恢复铁磁性。但是实验发现，在 AuFe、CuMn 等一些稀磁合金中，温度下降时，材料的磁化率先缓慢增高，经过一个峰值后再缓慢下降{cite}`shih1931magnetic,kondo1964resistance,cannella1972magnetic`。这是因为这类材料中自旋之间的相互作用是完全随机的，温度下降时，复杂的相互作用产生了阻挫效应（frustration effect），无法形成长程有序状态，各个磁矩被随机地冻结在某个方向，最后呈现无规则的长程无序状态。1970年，P. Anderson 在 B. R. Coles 的建议下提出”自旋玻璃“的概念{cite}`anderson1970localisation`，”玻璃“一词是将这种铁磁材料中自旋排列的无序与化学玻璃中分子位置的无序进行类比。

&emsp;&emsp;1975年，S. Edwards 与 P. Anderson 提出第一个短程相互作用的自旋玻璃模型（E-A 模型），并提出了复本方法（replica method）和复本对称假设（the assumption of replica symmetry）来求解该模型{cite}`edwards1975theory`。这是第一个定性描述自旋玻璃的模型，能够在低温下观察到玻璃相{cite}`nishimori2001statistical`的存在。同年，D. Sherrington 与 S. Kirkpatrick 提出了自旋玻璃的平均场模型（S-K 模型）{cite}`sherrington1975solvable`，将短程相互作用推广到全连接的情形。根据复本对称理论的计算结果，S-K 模型中存在一个临界温度，在临界温度以上序参量为零，临界温度以下序参量大于零，这再一次验证了 Landau 的理论。但是，使用复本对称理论求解出 S-K 模型的零温熵为负值，这是不符合物理直观的，因为 S-K 模型中自旋取向只有两个方向，系统的构型数是可数的，熵不可能为负。1978年，de Almeida 和 Thouless 发现 S-K 模型的复本对称解仅在高温下稳定，在自旋玻璃相中并不稳定，他们找到了复本对称解稳定区间的边界（dAT 线），在 dAT 线之下复本对称性破缺{cite}`de1978stability`。为了尝试打破复本对称性来描述自旋玻璃的低温相，Thouless、Anderson 与 Palmer 提出一种平均场方法（TAP 方法），在自由能中考虑自旋对自身的反弹效应（Onsager 反应项{cite}`onsager1936electric,barker1973monte`），并且仅作无序平均{cite}`thouless1977solution`。但是 TAP 方法也被证明仅在高温下有效{cite}`mezard1987spin`。

&emsp;&emsp;1979年，G. Parisi 提出了复本对称破缺（replica symmetry breaking, RSB）的概念，发展了一套有效的数学方法，并给出了 S-K 模型的一个精确的理论解。随后，M. Mezard、G. Parisi、M.A. Virasoro 等人对 Parisi 解的解释工作揭示了玻璃态低温相的复杂性质，其特征是遍历性破坏（ergodicity breaking）、超度量性（ultrametricity）和非自均性（non-selfaverageness）。

··· ···

1987年，Marc Mézard、Giorgio Parisi 和 Miguel Angel Virasoro 提出了一种新的数学方法用于处理 S-K 模型，被称为空腔方法（cavity method），后来发现，空腔方法对于统计物理中的一些平均场模型，特别是无序系统，具有广泛的适用性，并且相比复本方法更为简单。

··· ···

<br>

**References**

```{bibliography}
:style: unsrt
:filter: docname in docnames
```