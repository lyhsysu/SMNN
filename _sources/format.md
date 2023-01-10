# *技术与排版规范

<br/>

本网站使用 [Jupyter Book](https://jupyterbook.org/en/stable/intro.html) 构建，基于 [GitHub](https://github.com/lyhsysu/lyhsysu) 发布。

本页面第一部分记录了网页更新流程，便于作者进行更新；第二部分是作者基于 Jupyter Book、[MyST Markdown](https://myst-parser.readthedocs.io/en/latest/index.html) 以及 [MathJax](http://docs.mathjax.org/) 的文档制作的一份适用于本网站的文章排版规范，以便全站文章风格统一，同时可供他人参考。

<br/>

## 一、网站更新流程

在图书根目录进行编辑

在图书根目录使用 Anaconda Prompt 命令行

```
jupyter-book build SMNN
```

将 build 完毕后的图书根目录所有文件复制到 GitHub 本地仓库目录

在 GitHub 本地仓库目录使用 Git 命令行

```
ghp-import -n -p -f _build/html
```

网页地址：https://lyhsysu.github.io/SMNN/intro.html

<br/>

## 二、排版规范

### 2.1 通用语法

### 2.2 数学公式

Jupyter Book 使用 [MathJax](http://docs.mathjax.org/) 排版数学公式，在 `_config.yml` 配置文件中，用如下命令将MathJax升级到版本3（默认为版本2）

```
sphinx:
  config:
    mathjax_path: https://cdn.jsdelivr.net/npm/mathjax@3/es5/tex-mml-chtml.js
```

**内联公式**：使用 markdown 的 `$···$` 符号

**公式块**：使用如下代码对公式块进行编号

``````
```{math}
:label: the_label_name
1+1=2
```
``````

使用如下代码对公式编号进行引用

```
式 {eq}`the_label_name` 
```

<div class="admonition note" name="html-admonition" style="background: lightgreen; padding: 10px">
<p class="title">示例</p>
<p>这是一个内联公式 $1+1=2$</p>
下面是一个公式块
```{math}
:label: my_label
1+1=2
```
这是对公式的引用：对于式 {eq}`my_label` 来说
</div>

<br/>

### 2.3 图片

在 `_config.yml` 配置文件中添加 `html_image`

```
parse:
  myst_enable_extensions:
	······
	······
    - html_image
```

在图书根目录新建 `images` 文件夹，用于存放图片。

#### 带题注和编号的图片

使用如下代码插入图片

``````
:::{figure-md} figure_label_name
<img src="figure_address" alt="figure_description" class="bg-primary mb-1" width="200px">

这是图片的题注
:::
``````

使用如下代码引用图片
```
{numref}`figure_label_name`
```
**效果为：**
:::{figure-md} markdown-fig
<img src="https://i.328888.xyz/2023/01/09/05AJE.md.png" alt="斯特恩-盖拉赫实验装置图" class="bg-primary mb-1" width="8cm">

这是图片的题注
:::
这是对图片的引用：对于 {numref}`markdown-fig` 来说

#### 不带题注和编号的图片

直接使用原生 markdown 语法插入图片

```
![图片描述](图片地址)
```

如

![](https://i.328888.xyz/2023/01/09/05AJE.md.png)

<br/>

### 2.4 参考文献

在图书根目录的 `references.bib` 文件中添加参考文献的 BibTex 引用样式

在文中使用 ``{cite}`ref_name` `` 引用参考文献

使用 ``{cite}`ref_name_1,ref_name_2,ref_name_3` ``一次引用多个参考文献 

在文章中使用如下代码添加参考文献目录

``````
```{bibliography}
:style: unsrt
:filter: docname in docnames
```
``````

<div class="admonition note" name="html-admonition" style="background: lightgreen; padding: 10px">
<p class="title">示例</p>
<p>引用一篇文献，比如{cite}`holdgraf_evidence_2014`</p>
<p>引用多篇文献，比如{cite}`holdgraf_rapid_2016,holdgraf_portable_2017,holdgraf_encoding_2017`</p>
<p>生成参考文献目录</p>
```{bibliography}
:style: unsrt
:filter: docname in docnames
```
</div>

<br/>

### 2.5 告示框

在 `_config.yml` 配置文件中添加 `html_admonition`

```
parse:
  myst_enable_extensions:
	······
	······
    - html_admonition
```

<br/>

- 使用如下代码添加一个注释

``````
```{note}
这是一个注释
```
``````

效果为

```{note}
这是一个注释
```

<br/>

- 使用如下代码在大标题下添加摘要

``````
```{admonition} Abstract
这是一个摘要
```
``````

效果为

```{admonition} Abstract
这是一个摘要
```

<br/>

- 使用如下代码添加一个想法/提示

``````
```{admonition} ideal/tips
:class: tip
这是一个想法/提示
```
``````

效果为

```{admonition} ideal/tips
:class: tip
这是一个想法/提示
```