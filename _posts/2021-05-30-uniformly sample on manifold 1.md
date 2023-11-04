---
title: Uniformly Sample on a Manifold - 1.Rationale
layout: post
categories: [Geometry, Probability, Math]
image: /assets/img/SampleOnManifold/cover.jpg
description: "Uniformly Sample on a Manifold - 1.Rationale"
customexcerpt: "A self-explored general method to sample uniformly on a manifold embedded in R^n"
---

# 流形上的均匀采样(一)——理论

* 
{:toc}




在GAMES 101 Homework7 蒙特卡洛路径追踪中，有这么一段代码：

<img src="\assets\img\SampleOnManifold\1.png" alt="1" style="zoom:60%;" />

作用是在三角形光源上均匀采样。我乍一看均匀采样，感觉是一件很平凡(trivial)的事情，不就是给重心坐标的u,v分别用伪随机数生成器生成[0,1]上的均匀点，再带到重心坐标公式里吗？但仔细一看，却发现代码写的不是那么平凡(non-trivial):

代码首先用随机数生成器生成了两个[0,1]上的均匀分布(uniformly distribution)

$$
\xi_1,\xi_2 \sim Uniform[0,1]
$$

接着用如下的方式生成了三角形上的点

$$
\mathbf{p} = (1-\sqrt{\xi_1})\mathbf{v_0}+(\sqrt{\xi_1}(1-\xi_2))\mathbf{v_1}+(\sqrt{\xi_1}\xi_2)\mathbf{v_2}
$$

注意到这里仍然是一个用u,v坐标表示的三角形，$\mathbf{v_0,v_1,v_2}$是三角形的三个顶点。但是u,v的分布并不是均匀分布，而是由均匀分布变换而来

$$
\begin{cases}
u= 1-\sqrt{\xi_1}\\
v= \sqrt{\xi_1}\xi_2
\end{cases}\\
\mathbf{p} = u\mathbf{v_0}+(1-u-v)\mathbf{v_1}+v\mathbf{v_2}
$$

**这是怎么回事？**

本着探究知识的态度，我在网上搜索了一番，很快发现了问题：单纯对参数均匀采样**并不能**给出三角形上的均匀分布（这里特指对于面积的均匀分布，即点落入单位面积中的概率是一样的，符合最简单的蒙特卡洛积分对面积均匀采样的做法）。再另举一个例子，对于用球坐标参数化表示的球面(sphere) $S^2$，

$$
(rsin(\theta)cos(\varphi),rsin(\theta)sin(\varphi),rcos(\theta)), \theta\in[0,\pi/2),\varphi\in[0,2\pi)
$$

如果我们单纯对$\theta,\varphi$在所在区间上均匀采样，得到圆上的点在极点处会分布地比较密集，在赤道附近会分布地比较稀疏，如下图

<img src="\assets\img\SampleOnManifold\2.png" alt="2" style="zoom: 67%;" />

显然不是均匀分布。那么，如果我们按照这样的采样方法计算蒙特卡洛积分，期间除以我们以为的均匀分布的pdf$\frac{1}{4\pi}$，得到的积分的期望应该就是**有偏的(biased)**，理论上是有问题的。

于是我们就得到了一个问题：**该如何在一个曲面上均匀采样呢？进一步的，我们该如何在一些欧氏空间中弯曲的对象，比如说曲面、曲线、特定的空间等，上根据面积/体积/长度均匀采样呢？** 本文将会给出一个较为通用的办法，能够在嵌入$\mathbb{R}^n$的子流形（即前面提到的曲面、曲线、特定的空间等）上均匀采样，并做出一些直观的讲解。由于本人水平有限（就一破工科生，没有什么系统的数学训练），本文不会涉及到太多数学证明，主要是直观的理解；也**不可避免地会出现不严谨的地方**，如有错误，希望各位大佬指正。

那让我们开始，首先是一些预备知识。

## 流形 Manifold

前面我们提到了**流形(manifold)**这个概念，乍一看很高大上，给人一种云里雾里的感觉，实际上很简单。直观上来讲，**流形(manifold)**就是局部具有**欧式空间(Euclidean Space)**性质的空间或者对象。举一个最简单的例子，身处在地球上，我们环顾四周，会觉得地球是平的。要不是近现代科学告诉我们地球是个球，单凭人类自己的感觉，我们大概率会觉得自己踩在一个平面上。所以我们说地球面（假设它是光滑的）局部是个二维的欧氏空间$\mathbb{R}^2$，也就是一个平面——我们称它为**二维流形(2-dimensional manifold)**。再举一个例子，二维平面中有一条光滑曲线，我们把这条曲线放大再放大会发现它局部是一根直线。所以我们说曲线局部是个一维的欧式空间$\mathbb{R}$，它是个**一维流形(1-dimensional manifold)**。从另一个角度来思考，流形也可以看作是一个平直的欧式空间扭曲的结果。就比如球面就是平面“扭一扭”形成的，曲线就是直线“扭一扭”形成的。广义上来讲，一个**d维流形(d-dimensional manifold)**是一个局部具有d维欧氏空间$\mathbb{R}^d$性质的空间，也可以看作是$\mathbb{R}^d$“扭一扭”形成的空间。

<img src="\assets\img\SampleOnManifold\3.png" alt="3" style="zoom:80%;" />

“扭一扭”之后是什么呢？当然是“泡一泡”了！（不舔一添了(笑)）这里我们需要知道的是，流形这个空间，比如说我说的一条曲线，实际上是独立于$\mathbb{R}^n$空间的，它不是在一个什么空间里，它就是它本身，而数学家一般研究的是这个曲线**内蕴的性质(intrinsic property)**。但是它确确实实可以放到$\mathbb{R}^2$里，也可以放到$\mathbb{R}^3$里，也可以放到$\mathbb{R}^n$里，我们把这个放的过程，也就是“泡一泡”的过程，叫做**浸入(immersion)**，而一种特殊的“泡一泡”就是**嵌入(embedding)**，本文要研究的就是嵌入$\mathbb{R}^n$空间的流形，也就是$\mathbb{R}^n$中的曲面/曲线/曲空间...

那么，有了流形这个概念后，我们就可以把$\mathbb{R}^n$中的(比较好的)曲面/曲线/曲空间...统称为$\mathbb{R}^n$中的流形了（其实这才是我的主要目的）。我们的地球面，就是三维欧氏空间$\mathbb{R}^3$中的二维流形；平时做数学题在平面内画的曲线，就是二维欧氏空间$\mathbb{R}^2$中的一维流形，而在三维空间画的曲线，就是$\mathbb{R}^3$中的一维流形；三维空间中的一个球体，是$\mathbb{R}^3$中的一个三维流形(流形也可以不扭曲)；甚至$\mathbb{R}^n$本身也可以当作$\mathbb{R}^n$中的一个n维流形。

至此，我们对流形的介绍就差不多结束了。本文引入流形概念主要是想把n维空间中的曲面/曲线/曲空间...给统一成一个对象，并且用一个方法去研究，要不然每次都说曲面/曲线/曲空间...可能会够呛。当然流形远不止此，抽象元素构成的集合只要满足了特定性质也能叫做流形，比如说一些矩阵的集合，再比如说一些映射的集合...倒不如说数学家研究更多的是这些东西，但是这不重要。本文甚至连流形最基本的**图卡(chart)**、**地图册(atlas)**都没介绍，原因正如之前所说：只是想统一一下本文研究的对象罢了。

## 流形上的均匀分布

既然我们想要研究如何在$\mathbb{R}^n$的子流形上均匀采样，那就得定义什么是均匀，也就要说明什么是流形上的概率分布（如何定义流形上的概率分布才能满足问题的需求）。

让我们首先回忆做这件事的初衷是什么：我们想要通过蒙特卡洛方法对曲面上的光进行积分得到通量。我们知道，在曲面上积分，无非就是对于面积微元$dA$把通过这一面积微元的向量场的通量$F\cdot dA$进行求和。因此，用蒙特卡洛方法进行曲面积分时，自然而然地要对曲面上的面积微元进行采样，最后再除以这样采样的pdf。于是本文要探讨的流形上的概率分布就很自然了：整个$\mathbb{R}^n$的子流形是样本空间，在流形中取每一个区域都有一个相对应的面积(测度)和这块面积相对应的概率(测度)。而均匀分布，就是在流形中取同样面积(测度)的区域，对应的概率都是一样的。采样是取点，于是均匀采样就变成了让采样点落入单位面积中的概率一样。

这里把面积称之为**测度(measure)**是因为面积这个概念是可以推广的。对于曲面，我们说它的表面积是多少；对于曲线，我们说它的长度是多少；对于一个物体，我们说它的体积是多少：这些东西的统称就是**测度(measure)**，准确来说是**勒贝格测度(Lebesgue measure)**。那么在曲面上均匀采样就是让采样点落入单位面积的概率一样；在曲线上均匀采样就是让采样点落入单位长度的概率一样；在空间上均匀采样就是让采样点落入单位体积的概率一样——**在d维流形上均匀采样就是让采样点落入单位(d维)测度的概率相同。**

$\mathbb{R}^n$的子流形$M$上均匀分布的概率密度函数pdf因此很自然的就成为了$f_\mathbf{x}(\mathbf{x}) = \frac{1}{\int_M 1dS}$，其中分母是整个子流形的测度(体积/面积/长度)

## 概率分布的变换 Transformation Between Distributions

于是现在的问题就变成了怎么把流形上的均匀分布和随机数的均匀分布联系起来，这就涉及到概率分布的变换。

---

### 二维概率分布的变换

先考虑一个最基本的问题：定义$\mathbb{R}^2$中的单位圆$$\mathbf{B}_1(0)$$, $\\{(x,y)\|x^2+y^2<1\\}$ 上的均匀分布，pdf 

$$
f_{X,Y}(x,y) = \begin{cases}\frac{1}{\pi}, (x,y)\in \mathbf{B}_1(0)\\0,\text{ otherwise} \end{cases}$$

，如果用极坐标换元

$$
\begin{cases}X=Rcos(\Theta)\\Y=Rsin(\Theta)\end{cases}\\
$$

把这块区域变换到$r,\theta$平面上的$\\{(r,\theta)\|r\in[0,1],\theta\in[0,2\pi]\\}$，概率分布还会是均匀分布吗？我们在概率论中学过，显然不是。同时，我们也学习过这种情况下的，也就是二维随机变量变换后概率分布的变换公式：

---

>考虑一个二维随机变量的变换$$\begin{pmatrix}W\\Z\end{pmatrix} = g(\begin{pmatrix}X\\Y\end{pmatrix})$$，其中$$\begin{pmatrix}X\\Y\end{pmatrix}$$有联合pdf $f_{X,Y}$，$g$是一个从$f_{X,Y}$支集到$\mathbb{R}^2$的一一映射。假设$g$的雅可比矩阵$J_g$几乎处处存在且连续且有非零行列式，则$$\begin{pmatrix}W\\Z\end{pmatrix}$$的联合概率密度函数是
>
>$$
>f_{W,Z}(w,z) = \frac{1}{|det(J_g)|}f_{X,Y}(g^{-1}(\begin{pmatrix}w\\z\end{pmatrix}))
>$$

---

**我们从直观的角度来回顾一下。**

想要求$W,Z$空间中的概率密度函数，就可以考虑求出$W,Z$空间中一点$(w,z)$对应的一个微小面积$dA_{w,z}$和其所对应的概率密度$P_{w,z}$(即$f_{W,Z}(w,z)$，为了表示方便). 由于面积极小，$dA_{w,z}$所对应的概率可以近似表达成$P_{w,z}dA_{w,z}$。 又因为变换前后概率应该是不变的，于是我们就可以回到$X,Y$空间找变换之前的这块面积，再找到这块面积对应的概率就好了。这样可以得到一个等式:

$$
P_{w,z}dA_{w,z} = P_{x,y}dA_{x,y}
$$

所以

$$
P_{w,z} = P_{x,y}\frac{dA_{x,y}}{dA_{w,z}}\\
f_{W,Z}(w,z) = \frac{1}{dA_{w,z}/dA_{x,y}}f_{X,Y}(x,y)
$$

也就是**新概率密度等于旧概率密度除以变换前后的面积之比**。

这个面积之比，就是众所皆知的**雅可比行列式(Jacobian)**。我们知道，一个矩阵的**行列式(determinant)**代表矩阵线性变换前后的面积之比。一个映射在x点的**雅可比矩阵(Jacobi Matrix)**，代表这个映射在x点附近近似的一个局部的线性变换（也就是微分）。因此雅可比矩阵的行列式在x点的值，就代表x附近局部的面积变换比例。

于是之前的公式就很好理解了。$g$是一个从$X,Y$到$W,Z$的映射，它的行列式$det(J_g)$代表点$(x,y)$附近的微小面积$dA_{x,y}$变换到$dA_{w,z}$的面积缩放比例。由于行列式可正可负（正负代表的是面积是否有被“翻转），对面积比例没有影响，我们给它加上绝对值。而变换前的$(w,z)$，就是$$g^{-1}(\begin{pmatrix}w\\z\end{pmatrix})$$

带入上式，就得到了

$$f_{W,Z}(w,z) = \frac{1}{|det(J_g)|}f_{X,Y}(g^{-1}(\begin{pmatrix}w\\z\end{pmatrix}))$$

### 多维概率分布的变换

**这样的做法也可以推广到别的维度。**

就考虑最基本的一维随机变量的变换$Y = g(X)$（要求是一一映射），

1. 首先在$Y$的空间内找一段和$y$对应的微小的长度$dL_y$，和概率密度$P_y$

2. 回到$X$的空间找变换前的$x$和对应的微小长度$dL_x$，和概率密度$P_x$
3. 根据概率相同列出等式$P_ydL_y = P_xdL_x$并变形得到$f_{Y}(y) = \frac{1}{dL_y/dL_x}f_X(x) = \frac{1}{dL_y/dL_x}f_X(g^{-1}(y))$
4. 找到变换前后长度之比

对于一维空间$\mathbb{R}$到一维空间$\mathbb{R}$的映射，变换前后长度之比就是微分的绝对值$\|\frac{dy}{dx\}\| = g'(g^{-1}(y))$，其实就是一维的雅可比行列式。于是就得到了众所皆知的一维随机变量变换的公式

$$
f_{Y}(y) = \frac{1}{g'(g^{-1}(y))}f_X(g^{-1}(y))
$$

那么对于n维随机变量的变换$\vec Y = g(\vec{X})$（要求是一一映射），公式显然就是

$$
f_{\vec{Y}}(\vec{y}) = \frac{1}{|det(J_g)|}f_{\vec{X}}(g^{-1}(\vec{y}))
$$

---

### $\mathbb{R}^n$子流形上概率分布的变换

回到本文想要探讨的问题，对于一个$\mathbb{R}^n$中的子流形，其上有一个均匀的概率分布，我们想让其与$[0,1]$上均匀的概率分布产生联系，就需要做概率分布的变换。但是对于$\mathbb{R}^n$中的流形，说“这块区域的概率”“那块区域的概率”显然是不够明确的，需要有一个明确的方式去描述这些区域，那就是**参数化(parametrization)**。而这个方法就是我们一直在用的，比如说对$\mathbb{R}^3$中二维球面$S^2$（$\mathbb{R}^3$中的二维流形）的描述

$$
\begin{pmatrix}rsin(\theta)cos(\varphi)\\rsin(\theta)sin(\varphi)\\rcos(\theta)\end{pmatrix}, \theta\in(0,\pi),\varphi\in(0,2\pi)
$$

对于$\mathbb{R}^2$中一维球面$S^1$也就是单位圆（$\mathbb{R}^2$中的一维流形）的描述

$$
\begin{pmatrix}rcos(\theta)\\rsin(\theta)\end{pmatrix},\theta\in(0,2\pi)
$$

对空间中三角形的描述（$\mathbb{R}^3$中二维流形），也就是重心坐标的表示方法

$$
\mathbf{p} = u\mathbf{v_0}+(1-u-v)\mathbf{v_1}+v\mathbf{v_2},u>0\wedge v>0\wedge u+v<1
$$

或者是对三维空间中的单位球$\mathbf{B}_1(0)$($\mathbb{R}^3$中的三维流形，这是一个没有扭曲过的流形)的描述

$$
\begin{pmatrix}rsin(\theta)cos(\varphi)\\rsin(\theta)sin(\varphi)\\rcos(\theta)\end{pmatrix}, r\in(0,1),\theta\in(0,\pi),\varphi\in(0,2\pi)
$$

是我们很熟悉的。

而参数化的过程，也可以看成是一个一一映射/变换，把点从参数空间映射到$\mathbb{R}^n$空间的子流形$M$上

**既然有这样的一一映射，自然的我们就可以用上述的方法通过流形上的概率分布得到参数空间的概率分布；之后只需要做参数空间的分布和$[0,1]^n$上的均匀分布之间的转换就行了！**

我们试着把步骤写一下：

对于一个参数化映射$\vec{M} = \gamma(\vec X), x\in \mathbb{R}^d$（假设参数化的是欧氏空间中的d维流形）

1. 首先在参数空间内找一块和$\vec{x}$对应的微小的d维体积$dV_x$，和概率密度$P_x$

2. 到流形$M$上找变换后的$\vec x$和对应的微小的d维体积$dV_m$，和概率密度$P_m$

3. 根据概率相同列出等式
   
   $$
   P_xdV_x = P_mdV_m
   $$
   
   并变形得到
   
   $$
   f_{\vec{X}}(\vec x) = \frac{dV_m}{dV_x}P_m
   $$
   
   由于我们需要在流形上均匀采样，$P_m$就是$\frac{1}{\int_M 1dS}$，其中分母是整个子流形的d维体积(体积/面积/长度)

4. 最后就是找到变换前后d维体积之比

**这里我们会碰到一个问题**

对于$\mathbb{R}^n$中的n维流形，变换前后的d维体积之比是很明显的，就是$\gamma$的雅可比行列式的绝对值$\|det(J_\gamma)\|$。但是对于$\mathbb{R}^n$中的$d$维流形的参数化，$d<n$，我们还能用行列式的绝对值代表变换的比例吗？

**不能！**因为$\gamma$的雅可比矩阵根本就不是一个方阵，而是一个$n\times d$的矩阵（因为$\gamma$是一个$\Omega$到$\mathbb{R}^n$，$\Omega\subseteq \mathbb{R}^d $的映射），**它甚至都没有行列式！**

比如说对三维空间中$S^2$的参数化

$$
\gamma(\theta,\varphi) = \begin{pmatrix}rsin(\theta)cos(\varphi)\\rsin(\theta)sin(\varphi)\\rcos(\theta)\end{pmatrix}, \theta\in(0,\pi),\varphi\in(0,2\pi)
$$

它的雅可比矩阵

$$
J_\gamma = \begin{pmatrix}
\frac{\partial\gamma_x}{\partial\theta} & \frac{\partial\gamma_x}{\partial\varphi}\\
\frac{\partial\gamma_y}{\partial\theta} & \frac{\partial\gamma_y}{\partial\varphi}\\
\frac{\partial\gamma_z}{\partial\theta} & \frac{\partial\gamma_z}{\partial\varphi}\\ \end{pmatrix} = 
\begin{pmatrix}
rcos(\theta)cos(\varphi) & -rsin(\theta)sin(\varphi)\\
rcos(\theta)sin(\varphi) & rsin(\theta)cos(\varphi)\\
-rsin(\theta) & 0\end{pmatrix}
$$

显然不是方阵，没法求雅可比行列式。

难道就没有办法了吗？显然也不是。对于这种情况下的体积比，有一个东西恰好符合，那就是**格拉姆行列式(Gram Determinant)的平方根。**

## 格拉姆行列式 Gram Determinant

### 基本概念

在数学上，我们定义向量$v_1,\cdots,v_n$的**格拉姆矩阵(Gram Matrix)**$G$，它的第$i,j$个元素$G_{i,j} = \left \langle v_i,v_j \right \rangle$，即$v_i,v_j$的内积（点乘）

工程领域中这个矩阵在机器学习中见的比较多，但是事实上它也具有很明显的几何意义。它的行列式，**格拉姆行列式(Gram Matrix)**代表$\mathbb{R}^n$中向量$v_1,\cdots,v_d$张成的**d维平行六面体(d-dimensional parallelepiped)**$P$的d维体积的平方, $n\geq d$。

这里由$v_1,\cdots,v_d$张成的平行六面体$P$定义成

$$
P = \{c_1v_1+\cdots+c_dv_d\;0\leq c_i \leq 1\}
$$

其实是$v_1,\cdots,v_d$张成的那个多边体。叫它六面体是不稳妥的，但是我找不到parallelepiped的别的翻译，只能作罢。

举个例子，三维向量$v_1,v_2$的格拉姆行列式是

$$
\begin{vmatrix}
\left \langle v_1,v_1 \right \rangle & \left \langle v_1,v_2 \right \rangle\\ 
\left \langle v_2,v_1 \right \rangle & \left \langle v_2,v_2 \right \rangle 
\end{vmatrix}
$$

它代表的就是$v_1,v_2$在三维空间中张出的平行四边形的面积（也就是二维体积）的平方（也是$v_1,v_2$叉乘的模）

再比如说，三维向量$v_1,v_2,v_3$的格拉姆行列式是

$$
\begin{vmatrix}
\left\langle v_1,v_1 \right\rangle & \left\langle v_1,v_2 \right\rangle & \left\langle v_1,v_3 \right\rangle\\ 
\left\langle v_2,v_1 \right\rangle & \left\langle v_2,v_2 \right\rangle &\left\langle v_2,v_3 \right\rangle\\
\left\langle v_3,v_1 \right\rangle & \left\langle v_3,v_2 \right\rangle &\left\langle v_3,v_3 \right\rangle
\end{vmatrix}
$$

代表的确实就是$v_1,v_2,v_3$在三维空间中张出的平行六面体的体积的平方

Gram行列式也可以写成这样的形式$det(G) = det(A^TA)$，其中$A = (v_1\|\cdots\|v_d)$，读者可自行验证

于是$v_1,\cdots,v_d$张成的d维parallelepiped的体积就可以表示为$\sqrt{det(A^TA)}$

>#### 验证
>
>不妨用三维向量$v_1,v_2$的格拉姆行列式来验证一下
>
>$$
>\sqrt{det(A^TA)} =
>\sqrt{\begin{vmatrix}
>\left \langle v_1,v_1 \right \rangle & \left \langle v_1,v_2 \right \rangle\\ 
>\left \langle v_2,v_1 \right \rangle & \left \langle v_2,v_2 \right \rangle 
>\end{vmatrix}}
>$$
>
>当$v_1$变成原本的a倍，
>
>$$
>\sqrt{\begin{vmatrix}
>\left \langle av_1,av_1 \right \rangle & \left \langle av_1,v_2 \right \rangle\\ 
>\left \langle v_2,av_1 \right \rangle & \left \langle v_2,v_2 \right \rangle 
>\end{vmatrix}}\\ =
>\sqrt{a^2\begin{vmatrix}
>\left \langle v_1,v_1 \right \rangle & \left \langle v_1,v_2 \right \rangle\\ 
>\left \langle v_2,v_1 \right \rangle & \left \langle v_2,v_2 \right \rangle 
>\end{vmatrix}} = a\sqrt{det(A^TA)}
>$$
>
>平行四边形的面积也变成原本的a倍
>
>当$v_1$向$v_2$方向上偏移一小段$av_2$
>
>$$
>\sqrt{\begin{vmatrix}
>\left \langle v_1+av_2,v_1+av_2 \right \rangle & \left \langle v_1+av_2,v_2 \right \rangle\\ 
>\left \langle v_2,v_1+av_2 \right \rangle & \left \langle v_2,v_2 \right \rangle 
>\end{vmatrix}}\\ =
>\sqrt{\begin{vmatrix}
>\left \langle v_1+av_2,v_1 \right \rangle + a\left \langle v_1+av_2,v_2 \right \rangle
>& \left \langle v_1+av_2,v_2 \right \rangle
>\\ 
>\left \langle v_2,v_1 \right \rangle + a\left \langle v_2,v_2 \right \rangle 
>& \left \langle v_2,v_2 \right \rangle 
>\end{vmatrix}}
>\\
>=\sqrt{\begin{vmatrix}
>\left \langle v_1+av_2,v_1 \right \rangle
>& \left \langle v_1+av_2,v_2 \right \rangle
>\\ 
>\left \langle v_2,v_1 \right \rangle
>& \left \langle v_2,v_2 \right \rangle 
>\end{vmatrix}}
>\\
>=\sqrt{\begin{vmatrix}
>\left \langle v_1,v_1 \right \rangle + a\left \langle v_2,v_1 \right \rangle
>& \left \langle v_1,v_2 \right \rangle + a\left \langle v_2,v_2 \right \rangle
>\\ 
>\left \langle v_2,v_1 \right \rangle
>& \left \langle v_2,v_2 \right \rangle 
>\end{vmatrix}}
>\\
>=\sqrt{\begin{vmatrix}
>\left \langle v_1,v_1 \right \rangle & \left \langle v_1,v_2 \right \rangle\\ 
>\left \langle v_2,v_1 \right \rangle & \left \langle v_2,v_2 \right \rangle 
>\end{vmatrix}}
>\\
>=\sqrt{det(A^TA)}
>$$
>
>平行四边形的面积不变
>
>确确实实符合面积的性质
>
>PS: Gram行列式实际上是$v_1,\cdots,v_n$向量**外积(exterior product)/楔积(wedge product)**的模$\left \|v_1\wedge\cdots\wedge v_n\right \|^2$，涉及到**外代数(exterior algebra）**，是数学家用代数方式描述几何的一种方法。什么行列式，Gram行列式，叉乘等等本质上都是这玩意衍生出来的。解释起来很复杂，这里就不扯开去了。

---

### Gram行列式与体积之比的联系

回想行列式的几何意义，$det(A)$是代表线性变换后面积/体积的缩放比例，也可以代表矩阵$A$列向量张成的面积/体积大小，雅可比矩阵代表局部线性变换，它的行列式$det(J)$代表变换后局部面积的缩放比例；

而Gram行列式开根号$\sqrt{det(A^TA)}$，同样代表$n\times d$矩阵$A = (v_1\|\cdots\|v_d)$列向量张成的面积/体积大小，只不过矩阵$A$可以不是方阵，代表的是一个可能升维度的线性变换(从$\mathbb{R}^d$到$\mathbb{R}^n$)，$\sqrt{det(A^TA)}$代表的就是升维度的线性变换前后d维体积的缩放比例。把一个升维映射$\gamma$的雅可比矩阵$J_\gamma$当成$A$带入式子，我们就得到了一个升维映射变换前后局部的d维体积的缩放比例：

$$
\sqrt{det(J_\gamma^TJ_\gamma)}
$$

特别的，当$d = n$时，该式与雅可比行列式相等。

由此，我们就得到了上文所说的缩放比例，并且在$n = d$时可以替代雅可比行列式。于是，参数空间上的均匀概率分布就有了一个十分通用的表达式：

$$
f_{\vec{X}}(\vec x) = \frac{dV_m}{dV_x}P_m = \frac{\sqrt{det(J_\gamma^T(\vec x)J_\gamma(\vec x))}}{\int_M 1dS},\vec x \in \Omega
$$

$\Omega$是参数所在的空间，即流形$M$的原像。

顺带一提的是，$\mathbb{R}^n$中子流形上的体积积分（即国内高数教材所说的第一类曲线/曲面积分）$\int_M \rho dS$ 可以看成是对流形上一点和该点密度的求和，在参数化前后$dS$的面积比也正是$\sqrt{det(J_\gamma^TJ_\gamma)}$。所以这些曲面/曲线积分事实上也可以概括成这样的通式

$$
\int_M \rho dS = \int_\Omega \rho(\vec{x})\sqrt{det(J_\gamma^TJ_\gamma)}\ dx^n
$$

当$\rho$为1的时候即为流形的体积。读者可以自行验证，该公式适用于所有曲线/曲面/体积的密度积分，第一类曲线/曲面积分只不过是它的特例罢了。

于是之前的概率分布能进一步地写成

$$
f_{\vec{X}}(\vec x) = \frac{\sqrt{det(J_\gamma^TJ_\gamma)}}{\int_\Omega \sqrt{det(J_\gamma^TJ_\gamma)}\ dx^n},\vec x \in \Omega
$$

除此之外为0。即为**最终的通用公式**。

有了这个公式，剩下的事情就很容易了。只需要通过满足$[0,1]^n$上的均匀分布的n个随机变量$\xi_1,\cdots,\xi_n$生成满足$f_{\vec{X}}(\vec x)$联合概率密度函数分布的随机变量$X_1,\cdots,X_n$即可。常用的方法有**逆采样法(Inverse Sampling)**和**拒绝采样法(Reject Sampling)**

## 逆采样 Inverse Sampling 

### 一维情况

逆采样的原理很简单，在概率与统计中也学习过。首先我们考虑随机变量$X$, $F_X$是它的**累积分布函数CDF(Cumulative Distribution Function)**。对$X$施以变换$F_X$, 得到的随机变量$Y = F_X(X)$满足什么样的分布？

由于$F_X$是单调递增的，对于任意一个$y\in [0,1]$，都可以找到$F_X(x_y) = y$，也就是$F_X$有反函数。可以直接证明$Y$的CDF

$$
F_Y(y) = P\{F_X(X)\leq y\} = P\{X\leq F_X^{-1}(y)\} = F_X(F_X^{-1}(y)) = y
$$

因此$f_Y(y) = F_Y'(y) = 1$，$Y$满足$[0,1]$上的均匀分布。

$$
F_X(X) = U \sim Uniform[0,1]
$$

那想要生成满足$X$分布的随机变量，把$F_X^{-1}$施加在$U$上即可

$$
X = F_X^{-1}(U)
$$

正确性也容易验证，$X$的CDF为

$$
F_X(c) = P\{F_X^{-1}(U)\leq c\} = P\{U\leq F_X(c)\} = F_X(x)
$$

### 高维情况

但是上述这只是对于一个随机变量的逆采样。**如果随机变量有多个，而且不一定独立的话，该如何生成满足一定联合概率密度函数的随机变量呢？**PBRT里给出了一种基于逆采样的办法

举一个例子，

如果现在我想要通过$[0,1]$上的均匀分布生成满足$f_{X,Y,Z}(x,y,z)$联合概率密度函数的随机变量$X,Y,Z$

不妨先生成满足$f_X(x)$关系的随机变量$X$

做法就是对联合密度函数积分，得到边缘密度函数$f_X(x) = \int_{-\infty}^{\infty}\int_{-\infty}^{\infty} f_{X,Y,Z}(x,y,z)dzdy$

再求$X$的CDF$F_X$, 最后对均匀分布$U$施加$F_X^{-1}$

到这一步都是常规。接下来，我们开始计算当$X$取某个特定值时，$Y,Z$满足的条件联合概率密度函数$f_{Y,Z\|X}(y,z\|x) = \frac{f_{X,Y,Z}(x,y,z)}{f_X(x)}$

由于当$X$取到一个定值时，我们进入了一个船新的宇宙，一切操作都可以和之前一样而不会出现错误，这是条件概率的性质。于是就可以找到$Y$的条件CDF$F_{Y\|X}(y\|x)$，对于均匀分布$U$施加$F_{Y\|X}^{-1}$就可以得到$Y$

之后就是找到$Z$的条件CDF$F_{Z\|X,Y}(z\|x,y)$再对于均匀分布施加$F_{Z\|X,Y}^{-1}$即可得到$Z$

这个方法从直观上理解也很容易。当$X$取一个特定值时，$Y$满足一个特定的概率分布，于是就能通过这个特定的分布生成$Y$；当$X,Y$取特定值时，$Z$满足一个特定的概率分布，于是也能通过这个特定的分布生成$Z$，以此类推。

而**拒绝采样**的原理也很简单，互联网上有很多资料，这里不加以赘述。

---

**综上所述，在$\mathbb{R}^n$的子流形$M$上根据体积均匀采样有如下几步：**

1. 找到$\mathbb{R}^n$中$M$的一种参数化映射$\gamma$（要求一一映射）

2. 求$M$的长度/面积/体积 ,即${\int_\Omega \sqrt{det(J_\gamma^TJ_\gamma)}\ dx^n}$

3. 通过公式得到参数空间$\Omega$上的概率密度函数

   $f_{\vec{X}}(\vec x) = \frac{\sqrt{det(J_\gamma^TJ_\gamma)}}{\int_\Omega \sqrt{det(J_\gamma^TJ_\gamma)}\ dx^n},\vec x \in \Omega$

4. 通过常用的方法，如逆采样法，拒绝采样法，MCMC等

---

时隔两年的更新：

- 在这里我们只是阐述了一种如何均匀采样的算法。实际上对于任意的分布我们都可以用上述的变换来达成采样的目的。变换系数就是$\sqrt{det(J_\gamma^TJ_\gamma)}$
- MCMC、拒绝采样，或者是resampled impotence sampling (RIS) 等方法都不需要目标分布是归一化 (normalized) 的，所以如果碰到${\int_\Omega \sqrt{det(J_\gamma^TJ_\gamma)}\ dx^n}$不好解析积分时，用上述方法不normalized pdf也是没有关系的。另外一种办法是数值积分，由于${\int_\Omega \sqrt{det(J_\gamma^TJ_\gamma)}\ dx^n}$是预先计算的，所以并不会影响到真正采样时的效率，下一篇文章中也给出了数值积分的例子。

---

下一篇文章会给出一些栗子来具体说明这个方法。