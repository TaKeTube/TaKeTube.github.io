---
title: Uniformly Sample on a Manifold - 2.Examples
layout: post
categories: [Geometry, Probability, Math]
image: /assets/img/SampleOnManifold/cover.jpg
description: "Uniformly Sample on a Manifold - 2.Examples"
customexcerpt: "A self-explored general method to sample uniformly on a manifold embedded in R^n"
---

# 流形上的均匀采样(二)——例子

* 
{:toc}




在上文中，我们给出了一个在嵌入$\mathbb{R}^n$中的子流形上均匀采样的方法。

在本文中，我们会给出一些例子来具体说明这个方法。

## 球面上的均匀采样

$\mathbb{R}^3$中二维球面$S^2$是$\mathbb{R}^3$中的二维流形，也是最常见的例子。它的一种参数化表示如下，即最常见的球坐标表示：

$$
\gamma(\theta,\varphi) = \begin{pmatrix}rsin(\theta)cos(\varphi)\\rsin(\theta)sin(\varphi)\\rcos(\theta)\end{pmatrix}, \theta\in(0,\pi),\varphi\in(0,2\pi)
$$

其中$r$是球的半径。

$S^2$的表面积就是$4\pi r^2$，不需要通过曲面积分计算（读者也可自行通过上述公式计算验证）

接下来不通过直接套用公式的方法，而是用球面的具体例子回顾一遍公式的推导

1. 首先在参数空间$(0,\pi/2)\times(0,2\pi)$内找一块和$(\theta,\varphi)$对应的微小的面积$dS_{\theta,\varphi}$，和概率密度$P_{\theta,\varphi}$

2. 到球面上找变换后的$(\theta,\phi)$和$dS_{\theta,\varphi}$变换后的微小的面积$dS_{S^2}$，和概率密度$P_{S^2}$

3. 因为变换前后概率相同，列出等式
   
   $$
   P_{\theta,\varphi}dS_{\theta,\varphi} = P_{S^2}dS_{S^2}
   $$
   
   并变形得到
   
   $$
   P_{\theta,\varphi} = f_{\Theta,\Phi}(\theta,\varphi) = \frac{dS_{S^2}}{dS_{\theta,\varphi}}P_{S^2} = \frac{1}{4\pi r^2}\frac{dS_{S^2}}{dS_{\theta,\varphi}}
   $$
   
   由于我们需要在球面上均匀采样，$P_{S^2}$就是$\frac{1}{4\pi r^2}$，其中分母就是球的表面积

4. 变换前后面积之比就是$\sqrt{det(J_\gamma^TJ_\gamma)}$，下面开始计算

$\gamma$的雅可比矩阵是

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

所以

$$
J_\gamma^TJ_\gamma = \begin{pmatrix}
r^2 & \\
 & r^2sin^2(\theta)\\
\end{pmatrix}
$$

$$
\sqrt{det(J_\gamma^TJ_\gamma)} = r^2sin(\theta)
$$

（这和我们平时用球坐标推导出的面积微元前的系数也相符合）

因此

$$
f_{\Theta,\Phi}(\theta,\varphi) = \frac{r^2sin(\theta)}{4\pi r^2} = \frac{sin(\theta)}{4\pi}， \theta\in(0,\pi),\varphi\in(0,2\pi)
$$

除此之外等于0。这里由于$\Theta,\Phi$相互独立，可以直接对另一个变量积分得到两个随机变量的边缘密度函数

$$
f_{\Theta}(\theta) = \int_0^{2\pi}\frac{sin(\theta)}{4\pi}d\varphi = \frac{sin(\theta)}{2},\quad \theta\in(0,\pi)\\
f_{\Phi}(\varphi) = \int_0^{\pi}\frac{sin(\theta)}{4\pi}d\theta = \frac{1}{2\pi},\quad \varphi\in(0,2\pi)
$$

两者的CDF​为

$$
F_\Theta(\theta) = \int_0^{\theta} f_{\Theta}(\alpha) d\alpha = \frac{1-cos(\theta)}{2},\quad \theta\in(0,\pi)\\
F_\Phi(\varphi) = \int_0^{\varphi} f_{\Phi}(\alpha) d\alpha = \frac{\varphi}{2\pi},\quad \varphi\in(0,2\pi)
$$

现在用逆采样法采样，直接让均匀分布

$$
\xi_1,\xi_2 \sim Uniform[0,1]
$$

等于CDF即可（相当于求反函数）

$$
\frac{1-cos(\theta)}{2} = \xi_1 \\
\frac{\varphi}{2\pi} = \xi_2
$$

所以

$$
\begin{cases}
\theta= arccos(1-2\xi_1)\\
\varphi= 2\pi\xi_2
\end{cases}\\
$$

带回参数方程得

$$
\begin{cases}
x = 2rcos(2\pi\xi_2)\sqrt{\xi_1(1-\xi_1)}\\
y = 2rsin(2\pi\xi_2)\sqrt{\xi_1(1-\xi_1)}\\
z = r(1-2\xi_1)
\end{cases}\\
$$

实现的时候既可以先求出参数，再带入$x,y,z$；或是直接求出$x,y,z$，取决于计算的快速与否

写一段代码画个图验证一下，左图是直接对参数均匀采样，右图是对面积均匀采样：

![sphere](\assets\img\SampleOnManifold\sphere.png)

效果不错。当然用别的参数化按照步骤也能得到同样的结果，区别只在计算的复杂上。

## 半球面上的均匀采样

$\mathbb{R}^3$中半球面$H^2$是$\mathbb{R}^3$中的二维流形

半球面上的均匀采样是做蒙特卡洛路径积分时最常见的模型，在它上做均匀采样和上述在球上做均匀采样的办法就差了个参数的取值范围。这里不赘述。

$$
\gamma(\theta,\varphi) = \begin{pmatrix}rsin(\theta)cos(\varphi)\\rsin(\theta)sin(\varphi)\\rcos(\theta)\end{pmatrix}, \theta\in(0,\frac{\pi}{2}),\varphi\in(0,2\pi)
$$

半球面的面积是$2\pi r^2$，按照相同步骤可得

$$
f_{\Theta,\Phi}(\theta,\varphi) = \frac{sin(\theta)}{2\pi}， \theta\in(0,\frac{\pi}{2}),\varphi\in(0,2\pi)
$$

$$
f_{\Theta}(\theta) = \int_0^{2\pi}\frac{sin(\theta)}{2\pi}d\varphi = sin(\theta),\quad \theta\in(0,\frac{\pi}{2})\\
f_{\Phi}(\varphi) = \int_0^{\frac{\pi}{2}}\frac{sin(\theta)}{2\pi}d\theta = \frac{1}{2\pi},\quad \varphi\in(0,2\pi)
$$

$$
F_\Theta(\theta) = \int_0^{\theta} f_{\Theta}(\alpha) d\alpha = 1-cos(\theta),\quad \theta\in(0,\frac{\pi}{2})\\
F_\Phi(\varphi) = \int_0^{\varphi} f_{\Phi}(\alpha) d\alpha = \frac{\varphi}{2\pi},\quad \varphi\in(0,2\pi)
$$

令CDF等于均匀分布$\xi_1,\xi_2$

$$
1-cos(\theta) = \xi_1 \\
\frac{\varphi}{2\pi} = \xi_2
$$

可得

$$
\begin{cases}
\theta= arccos(\xi_1)\\
\varphi= 2\pi\xi_2
\end{cases}\\
$$

带回参数方程可得

$$
\begin{cases}
x = rcos(2\pi\xi_2)\sqrt{1-\xi_1^2}\\
y = rsin(2\pi\xi_2)\sqrt{1-\xi_1^2}\\
z = r\xi_1
\end{cases}\\
$$

这基本上就是GAMES101 Homework 7作业框架里的采样方法

<img src="\assets\img\SampleOnManifold\4.png" alt="4" style="zoom:60%;" />

只不过框架里多了一步取绝对值$\|1-2\xi_1\|$，但是这生成的仍然是$[0,1]$上的均匀分布，不知道为什么要多这么一步。

写一段代码画个图验证一下，左图是直接对参数均匀采样，右图是对面积均匀采样：

<img src="\assets\img\SampleOnManifold\hemisphere.png" alt="hemisphere"  />

效果不错。

## 椭圆盘内部的均匀采样

椭圆盘是$\mathbb{R}^2$中的二维流形（是没有扭曲过的空间），在它内部均匀采样就是一个普通的二维随机变量变换的问题，但是仍然可以用上述流程解决（过程实际上是一样的）。一般来说举的例子是圆盘内的随机采样(比如说PBRT中举的)，这里给它拓展成椭圆盘。

椭圆盘的一个参数化$\gamma$是

$$
\gamma(\theta,\varphi) = \begin{pmatrix}kacos(\theta)\\kbsin(\theta)\end{pmatrix}, k\in(0,1),\theta\in(0,2\pi)
$$

它的面积是$\pi ab$，所以PDF在圆盘内部是$\frac{1}{\pi ab}$，外部是0

$\gamma$的雅可比矩阵是

$$
J_\gamma = \begin{pmatrix}
\frac{\partial\gamma_x}{\partial k} & \frac{\partial\gamma_x}{\partial\theta}\\
\frac{\partial\gamma_y}{\partial k} & \frac{\partial\gamma_y}{\partial\theta}\\
\end{pmatrix} = 
\begin{pmatrix}
acos(\theta) & -kasin(\theta)\\
bsin(\theta) & kbcos(\theta)\end{pmatrix}
$$

由于$J_\gamma$是方阵，$\sqrt{det(J_\gamma^TJ_\gamma)}$就是$det(J_\gamma)$（读者可自行验证），所以比例系数是

$$
\begin{vmatrix}
acos(\theta) & -kasin(\theta)\\
bsin(\theta) & kbcos(\theta)
\end{vmatrix}=
kab(cos(\theta)^2+sin(\theta)^2) = kab
$$

因此参数空间的概率密度函数为

$$
f_{K,\Theta}(k,\theta) = \frac{kab}{\pi ab} = \frac{k}{\pi}， k\in (0,1),\theta\in(0,2\pi)
$$

除此之外为0。发现和椭圆长短轴没有关。

这里$K,\Theta$仍然是独立的，可以直接得到两个的PDF

$$
f_{K}(k) = 2k,\quad k\in(0,1)\\
f_\Theta(\theta) = \frac{1}{2\pi},\quad \theta\in(0,2\pi)
$$

CDF为(省略取值范围了)

$$
F_K(k) = k^2\\
F_\Theta(\theta) = \frac{\theta}{2\pi}
$$

令CDF等于均匀分布$\xi_1,\xi_2$

$$
k^2 = \xi_1 \\
\frac{\theta}{2\pi} = \xi_2
$$

可得

$$
\begin{cases}
k = \sqrt{\xi_1}\\
\theta = 2\pi\xi_2
\end{cases}\\
$$

带回参数方程可得

$$
\begin{cases}
x = a\sqrt{\xi_1}cos(2\pi\xi_2)\\
y = b\sqrt{\xi_1}sin(2\pi\xi_2)
\end{cases}\\
$$

当$a = b$时，就是圆盘内部的均匀采样。

写一段代码画个图验证一下，左图是直接对参数均匀采样，右图是对面积均匀采样：

<img src="\assets\img\SampleOnManifold\ellipse.png" alt="ellipse" style="zoom:80%;" />

效果不错。

## 球体内的均匀采样

$\mathbb{R}^3$中的三维球$\mathbf{B}_r(0)$是$\mathbb{R}^3$中的三维流形（是没有扭曲过的空间）。在它内部均匀采样就是一个普通的三维随机变量变换的问题，但是仍然可以用上述流程解决（过程实际上是一样的）。

$\mathbb{R}^3$中以原点为圆心，半径为$a$的三维球$\mathbf{B}_a(0)$的一个参数化$\gamma$是

$$
\gamma(\theta,\varphi) = \begin{pmatrix}rsin(\theta)cos(\varphi)\\rsin(\theta)sin(\varphi)\\rcos(\theta)\end{pmatrix}, r\in(0,a),\theta\in(0,\pi),\varphi\in(0,2\pi)
$$

它的体积是$\frac{4}{3}\pi a^3$，所以PDF在球内部是$\frac{3}{4\pi a^3}$，外部是0

$\gamma$的雅可比矩阵是

$$
J_\gamma = \begin{pmatrix}
\frac{\partial\gamma_x}{\partial r} & \frac{\partial\gamma_x}{\partial\theta} & \frac{\partial\gamma_x}{\partial\varphi}\\
\frac{\partial\gamma_y}{\partial r} &\frac{\partial\gamma_y}{\partial\theta} & \frac{\partial\gamma_y}{\partial\varphi}\\
\frac{\partial\gamma_z}{\partial r} & \frac{\partial\gamma_z}{\partial\theta} & \frac{\partial\gamma_z}{\partial\varphi}\\ \end{pmatrix} = 
\begin{pmatrix}
sin(\theta)cos(\varphi) & rcos(\theta)cos(\varphi) & -rsin(\theta)sin(\varphi)\\
sin(\theta)sin(\varphi) & rcos(\theta)sin(\varphi) & rsin(\theta)cos(\varphi)\\
cos(\theta) & -rsin(\theta) & 0\end{pmatrix}
$$

由于$J_\gamma$是方阵，$\sqrt{det(J_\gamma^TJ_\gamma)}$就是$det(J_\gamma)$，所以比例系数是

$$
\begin{vmatrix}
sin(\theta)cos(\varphi) & rcos(\theta)cos(\varphi) & -rsin(\theta)sin(\varphi)\\
sin(\theta)sin(\varphi) & rcos(\theta)sin(\varphi) & rsin(\theta)cos(\varphi)\\
cos(\theta) & -rsin(\theta) & 0\end{vmatrix} = r^2sin(\theta)
$$

（和我们平时用球坐标推导出的体积微元前的系数也相符合）

因此参数空间的概率密度函数为

$$
f_{R,\Theta, \Phi}(r,\theta, \varphi) = \frac{3r^2sin(\theta)}{4\pi a^3}， r\in (0,a),\theta\in(0,\pi),\varphi\in(0,2\pi)
$$

除此之外为0。这里$R,\Theta,\Phi$还是独立的，可以直接得到三个参数的PDF

$$
f_R(r) = \int_0^{\pi}\int_0^{2\pi} \frac{3r^2sin(\theta)}{4\pi a^3} d\varphi d\theta = \frac{3r^2}{a^3},\quad r\in(0,a)\\
f_{\Theta}(\theta) = \int_0^a\int_0^{2\pi} \frac{3r^2sin(\theta)}{4\pi a^3} d\varphi dr = \frac{sin(\theta)}{2},\quad \theta\in(0,\pi)\\
f_{\Phi}(\varphi) = \int_0^a\int_0^{\pi} \frac{3r^2sin(\theta)}{4\pi a^3} d\theta dr = \frac{1}{2\pi},\quad \varphi\in(0,2\pi)
$$

除此之外都是0。则他们的CDF(省略取值范围)为

$$
F_R(r) = \frac{r^3}{a^3}\\
F_\Theta(\theta) = \frac{1-cos(\theta)}{2}\\
F_\Phi(\varphi) = \frac{\varphi}{2\pi}
$$

令CDF等于$[0,1]$上的均匀分布$\xi_1,\xi_2,\xi_3$

$$
\frac{r^3}{a^3} = \xi_1\\
\frac{1-cos(\theta)}{2} = \xi_2\\
\frac{\varphi}{2\pi} = \xi_3
$$

可得

$$
\begin{cases}
r = a\sqrt[3]{\xi_1}\\
\theta = arccos(1-2\xi_2)\\
\varphi = 2\pi\xi_3
\end{cases}\\
$$

发现$\theta$和$\phi$和球面上采样是一样的，而$r$有了一个立方根的变换。从这里我们也可以推测出如果只是按照参数均匀采样，球心的点会较为密集。

带回参数方程，可得

$$
\begin{cases}
x = 2a\sqrt[3]{\xi_1}\sqrt{\xi_2(1-\xi_2)}cos(2\pi\xi_3)\\
y = 2a\sqrt[3]{\xi_1}\sqrt{\xi_2(1-\xi_2)}sin(2\pi\xi_3)\\
z = a\sqrt[3]{\xi_1}(1-2\xi_2)
\end{cases}\\
$$

写一段代码画个图验证一下，左图是直接对参数均匀采样，右图是对面积均匀采样：

<img src="\assets\img\SampleOnManifold\ball.png" alt="ball" style="zoom:80%;" />

效果很不错，可以很明显地看出按照参数均匀采样带来的不均匀。

## 三角形上的均匀采样

我们来解决本文最开始提出的问题，如何在一个三角形上按照面积均匀采样？

问题中的三角形特指$\mathbb{R}^3$中三角形，是一个$\mathbb{R}^3$中的三维流形。但是事实上，通过接下来的推导，我们可以发现结果适用于任意维度欧氏空间$\mathbb{R}^n$中的三角形，只要$n\geq 2$即可。

首先，找到三角形的参数化$\gamma$，最常见就是重心坐标表示的三角形，即

$$
u\mathbf{v_0}+v\mathbf{v_1}+(1-u-v)\mathbf{v_2},\quad u>0,v>0,u+v<1
$$

其中$\mathbf{v_0},\mathbf{v_1},\mathbf{v_2}$是三角形的三个顶点，他们是$\mathbb{R}^n$中的向量。这可以看作一个映射。

为了方便表示，这里将$\mathbf{v_0},\mathbf{v_1},\mathbf{v_2}$表示为$\vec{A},\vec{B},\vec{C}$

此处，我们把三角形的面积记作$S$，但是不计算它。

之后我们开始计算$\sqrt{det(J_\gamma^TJ_\gamma)}$

这里与之前不同的是，之前的$\gamma$都能明显地把每个坐标值都表示出来，而这里我们甚至不知道向量的维度$n$具体是多少。但是仍然可以通过矩阵分块的办法把雅可比矩阵表示出来。

我们知道雅可比矩阵的第一列就是$\gamma(u,v)$的所有坐标值对$u$求导。$\gamma(u,v)$第$i$列的值是

$$
uA_i+vB_i+(1-u-v)C_i
$$

对$u$求导就是$A_i-C_i$。整体来看，雅可比矩阵的第一列就是

$$
\vec{A}-\vec{C} = \vec{CA}
$$

同理可得雅可比矩阵的第二列就是

$$
\vec{B}-\vec{C} = \vec{CB}
$$

因此，雅可比矩阵就可以这么表示

$$
J_\gamma = \begin{pmatrix} \vec{CA} & \vec{CB}\end{pmatrix}
$$

$J_\gamma^TJ_\gamma$就是

$$
J_\gamma^TJ_\gamma = \begin{pmatrix} \vec{CA}^T \\ \vec{CB}^T\end{pmatrix}\begin{pmatrix} \vec{CA} & \vec{CB}\end{pmatrix}
$$

$J_\gamma^T$中向量横着写，$J_\gamma$中向量竖着写。而将他们乘起来的结果就是向量分别点乘

$$
J_\gamma^TJ_\gamma = 
\begin{pmatrix}
\vec{CA}^T\vec{CA} & \vec{CA}^T\vec{CB}\\
\vec{CB}^T\vec{CA} & \vec{CB}^T\vec{CB}
\end{pmatrix} = 
\begin{pmatrix}
\vec{CA}\cdot\vec{CA} & \vec{CA}\cdot\vec{CB}\\
\vec{CB}\cdot\vec{CA} & \vec{CB}\cdot\vec{CB}
\end{pmatrix}
$$

它的雅可比行列式就是

$$
det(J_\gamma^TJ_\gamma) = (\vec{CA}\cdot\vec{CA})(\vec{CB}\cdot\vec{CB}) - (\vec{CA}\cdot\vec{CB})^2
$$

这里为了方便，将$\vec{CA},\vec{CB},\vec{AB}$对应的三角形的三条边的边长记为$a,b,c$，三条边长对应的角为$A,B,C$，会发现

$$
det(J_\gamma^TJ_\gamma) = ||\vec{CA}||^2||\vec{CB}||^2 - ||\vec{CA}\cdot\vec{CB}||^2\\
=b^2a^2-|a\cdot b\cdot cosC|^2 = a^2b^2sin(C)^2 = 4(\frac{1}{2}absinC)^2 = 4S^2
$$

正是四倍的三角形面积的平方！因此$\sqrt{det(J_\gamma^TJ_\gamma)} = 2S$

带入公式，求得参数空间的概率密度函数为

$$
f_{U,V}(u,v) = \frac{2S}{S} = 2, \quad \quad u>0,v>0,u+v<1
$$

除此之外为0。十分的简单。之后只要求出CDF即可。

但是注意此处$U,V$两个随机变量**并不独立**，我们应该先求得$U$的CDF生成$U$，再通过条件CDF的方法生成$V$。

参数空间（PDF的支集）$\Omega = \\{(u,v)\|u>0,v>0,u+v<1\\}$，是坐标轴和$u+v = 1$包住的三角形区域

<img src="\assets\img\SampleOnManifold\triangleOmega.png" alt="triangleOmega" style="zoom:50%;" />

$U$的边缘概率密度函数为

$$
f_U(u) = \int_0^{1-u}2dv = 2(1-u),\quad u\in(0,1)
$$

CDF为

$$
F_U(u) = 2u-u^2,\quad u\in(0,1)
$$

让$\xi_1$等于CDF，并用求根公式解出$u$

$$
2u-u^2 = \xi_1\\
u^2-2u+\xi_1 = 0\\
u = \frac{2\pm\sqrt{4-4\xi_1}}{2} = 1\pm \sqrt{1-\xi_1}
$$

由于$u$在$(0,1)$上取值，舍去一个根，得

$$
u = 1-\sqrt{1-\xi_1}
$$

又因为$1-\xi_1$同样是一个$[0,1]$上的均匀分布，可以直接将$1-\xi_1$替换成$\xi_1$来节省运算，因此

$$
u = 1-\sqrt{\xi_1}
$$

接着求$v$的条件PDF

$$
f_{V|U}(v|u) = \frac{f_{U,V}(u,v)}{f_{U}(u)} = \frac{2}{2(1-u)} = \frac{1}{1-u}
$$

注意此处$v$的范围是$(0,1-u)$，除此之外都是0或者说不是良定义(well defined)的

因此$v$的条件CDF为，代表当$u$取某个值时$v$的CDF

$$
F_{V|U}(v|u) = \frac{v}{1-u},\quad v\in(0,1-u)
$$

因此在$u$取$1-\sqrt{\xi_1}$时，用还CDF的反函数作用于均匀分布就能得到$u$取$1-\sqrt{\xi_1}$时满足此时概率分布的$V$。另

$$
\frac{v}{1-u} = \xi_2\\
v = \sqrt{\xi_1}\xi_2
$$

综上

$$
\begin{cases}
u= 1-\sqrt{\xi_1}\\
v= \sqrt{\xi_1}\xi_2
\end{cases}\\
$$

带入原参数方程得

$$
\mathbf{p} = (1-\sqrt{\xi_1})\mathbf{v_0}+(\sqrt{\xi_1}\xi_2)\mathbf{v_1}+(\sqrt{\xi_1}(1-\xi_2))\mathbf{v_2}
$$

于是就解决了为什么HW7里是这么写的问题。并且我们从上述推导知道，该公式适用于**任意大于等于二维的欧式空间中的三角形**均匀采样。

写一段代码画个图验证一下，左图是对参数均匀采样，这里由于$u,v$相关，单纯对两者均匀采样会得到一个四边形，因此采用了如下的采样方式

$$\begin{cases}
u= \xi_1\\
v= (1-\xi_1)\xi_2
\end{cases}\\$$

。右图是对面积均匀采样：

<img src="\assets\img\SampleOnManifold\triangle.png" alt="triangle" style="zoom:80%;" />

效果很好。

## 曲线上的均匀采样

一般我们遇见的空间中不自交的光滑曲线基本上是$\mathbb{R}^n$中的一维流形，比如说圆周，比如说贝塞尔曲线。而如果曲线自交的话，就不能称之为流形，因为在自交的那点，无论放的多大，交叉仍然存在，不能和一维欧氏空间建立对应关系（不局部同胚）。但是只要它的性质足够好（比如说交点有限或**可数(countable)**无穷），我们仍然能用之前的方法在它上面均匀采样，因为可数个交点的长度（1维体积/测度）为0，不影响概率。

下面举一个$\mathbb{R}^3$中参数曲线的例子（有一个自交点，不是流形），这个例子的PDF的积分没法正常计算(是椭圆积分)，但是能通过拒绝采样的办法采样。

考虑以下曲线

$$
\gamma(t) = \begin{pmatrix}cos(t)^2\\sin(t)cos(t)\\sin(t)\end{pmatrix}, t\in(0,2\pi)
$$

不难发现，对于1维参数化而言，比例系数

$$
\sqrt{det(J_\gamma^TJ_\gamma)} = \sqrt{\gamma'(t)^T\gamma'(t)} = ||\gamma'(t)||
$$

就是$\gamma$导数的模长，也是曲线长度计算公式里的比例系数。因此

$$
\sqrt{det(J_\gamma^TJ_\gamma)} = \sqrt{(2sin(t)cos(t))^2+(cos(t)^2-sin(t)^2)^2+cos(t)^2}\\
=\sqrt{(cos(t)^2+sin(t)^2)^2+cos(t)^2} = \sqrt{1+cos(t)^2}
$$

$t$的PDF为

$$
f_T(t) = \frac{\sqrt{1+cos(t)^2}}{\int_0^{2\pi}\sqrt{1+cos(t)^2}dt},\quad t\in(0,2\pi)
$$

其中分母的积分，即曲线的长度是第二类椭圆积分，只有数值解，用计算机计算也是方便的，我们记作$L$

$$
L = \sqrt{2}E(t,1/2)\bigg|_0^{2\pi} = 4\sqrt{2}E(\frac{1}{2}) \approx 7.64039557805542403580952
$$

此时即使求的出CDF（写成椭圆积分$E(t,m)$的形式）也很难找出CDF的逆，因此我们采用拒绝采用法对$t$进行采样

$$
f_T(t) = \frac{\sqrt{1+cos(t)^2}}{L},\quad t\in(0,2\pi)
$$

选取了最简单的均匀分布作为建议分布，参数$c = 2\pi max\{f_T(t)\} = 2\pi \sqrt{2}/L$，伪代码如下

```python
while len(X) < sampleNum:
    t = 2*np.pi*np.random.rand()
    if np.sqrt(2)*np.random.rand() < np.sqrt(1+np.cos(t)**2):
        X.append(np.cos(t)**2)
        Y.append(np.cos(t)*np.sin(t))
        Z.append(np.sin(t))
```

得到了如下结果，左图是按参数均匀采样，右图是按长度均匀采样：

<img src="\assets\img\SampleOnManifold\curve.png" alt="curve" style="zoom:80%;" />

效果...不怎么看的出来，甚至有种左边更均匀的错觉。既然画图看不出效果，我们可以直接从做这件事最终的目的入手——**蒙特卡洛积分**。不如来计算两种情况下曲线上的密度积分

$$
\int_M z^2 dL
$$

通过手算和计算机，得到该积分值是

$$
\int_M z^2 dL = \int_\Omega \gamma_z(t)\sqrt{det(J_\gamma^TJ_\gamma)} =  \int_0^{2\pi}sin(t)^2 \sqrt{1+cos(t)^2} dt \\
\approx 3.4960767390561597472864527
$$

而蒙特卡洛积分的公式是

$$
\frac{1}{N}\sum_{i=1}^N\frac{f(\mathbf{v}_i)}{p(\mathbf{v}_i)} = \frac{L}{N}\sum_{i=1}^N \mathbf{v}_{iz}^2
$$

因为我们假设是均匀采样所以$p(\mathbf{v}_i) = \frac{1}{L}$。采样50000个点的结果如下

![curveintegral](\assets\img\SampleOnManifold\curveintegral.png)

很明显地看得出对于参数均匀采样会导致结果错误，而对长度均匀采样蒙特卡洛积分值和真实值相近。

现在来探究一下为什么图像看起来差不多的原因。其实按照上文推导的思路，我们可以直接表示出如果对参数空间均匀采样，流形上的概率分布。因为

$$
f_{\vec{X}}(\vec x) = \sqrt{det(J_\gamma^TJ_\gamma)} P_m
$$

那么流形上的概率分布就是这样

$$
P_m(\gamma(\vec x)) = \frac{f_{\vec{X}}(\vec x)}{\sqrt{det(J_\gamma^TJ_\gamma)}}
$$

因此在本例子中如果$t$是均匀分布$T\sim Uniform[0,2\pi]$，曲线上的概率分布就是

$$
P_m(\gamma(t)) = \frac{1}{2\pi \sqrt{1+cos(t)^2}}
$$

这里暂且仍然用参数$t$表示（所以上式对$t$积分不是1，但是事实上如果考虑曲线积分，结果就是1，符合概率分布的定义）。在参数空间画出图像，如下图，可以发现概率分布接近均匀，所以在视觉上看上去是均匀的。

<img src="\assets\img\SampleOnManifold\curvePDF.png" alt="curvePDF" style="zoom: 50%;" />

把上述值映射到曲线上，并作图，能较为清楚地看到曲线上的概率密度差异

<img src="\assets\img\SampleOnManifold\curve3.png" alt="curve3" style="zoom: 67%;" />

计算一下按照参数采样的错误蒙特卡洛积分的期望

$$
E[\frac{f(\mathbf{v}_i)}{p_{uniform}(\mathbf{v}_i)}] = \int_\gamma \frac{f(\mathbf{v}_i)}{p_{uniform}(\mathbf{v}_i)} p_{actual}(\mathbf{v}_i) dL \\
= \int_0^{2\pi} \frac{sin(t)^2}{1/L}\cdot \frac{1}{2\pi \sqrt{1+cos(t)^2}} \cdot \sqrt{1+cos(t)^2} dt\\
= \frac{L}{2\pi} \int_0^{2\pi} sin(t)^2 dt = \frac{L}{2} \\
\approx 3.82019778903
$$

和按参数采样的数据十分接近。

这样的方法能也够运用到别的例子里。

## 环面上的均匀采样

环面(Torus) $T^n$可以说是拓扑学里比较经典的对象。下面试着对$\mathbb{R}^3$中的二维环面$T^2$（$\mathbb{R}^3$中的二维流形）进行均匀采样。

一个环半径为$R$，管半径为$r$的环面的一个参数化如下

$$
\gamma(t) = \begin{pmatrix}(R+rcos(t))cos(\theta)\\(R+rcos(t))sin(\theta)\\rsin(t)\end{pmatrix},\quad t\in(0,2\pi),\theta\in(0,2\pi)
$$

直接求其雅可比矩阵比较繁琐，我们先求出旋转面的通式，再把运用到环面的特例上。一个由二维参数曲线$$\alpha(t) = \begin{pmatrix}r(t)\\z(t)\end{pmatrix}$$ $t\in I$ 绕z轴旋转得到的旋转面为

$$
\gamma(\theta,t) = \begin{pmatrix}r(t)cos(\theta)\\r(t)sin(\theta)\\z(t)\end{pmatrix},\quad t\in I,\theta\in(0,2\pi)
$$

它的雅可比矩阵为

$$
J_\gamma = \begin{pmatrix}
-r(t)sin(\theta) & r'(t)cos(\theta)\\
r(t)cos(\theta) & r'(t)sin(\theta)\\
0 & z'(t)\end{pmatrix}
$$

$$
det(J_\gamma^TJ_\gamma) = \begin{vmatrix}
r(t)^2 & 0\\
0 & r'(t)^2+z'(t)^2\\
\end{vmatrix} = r(t)^2||\alpha'(t)||^2
$$

因此比例系数$\sqrt{det(J_\gamma^TJ_\gamma)} = \|r(t)\|\|\|\alpha'(t)\|\|$

运用到环面上

$$
\sqrt{det(J_\gamma^TJ_\gamma)} = (R+rcos(t))||(-rsin(t),rcos(t))|| = r(R+rcos(t))
$$

计算出环面的表面积

$$
\int_M 1dS = \int_0^{2\pi}\int_0^{2\pi}r(R+rcos(t))d\theta dt = \\
2\pi r\int_0^{2\pi}R+rcos(t)dt = 4\pi^2rR
$$

求出参数空间的PDF

$$
f_{T,\Theta}(t,\theta) = \frac{r(R+rcos(t))}{4\pi^2rR} = \frac{R+rcos(t)}{4\pi^2R},\quad t\in(0,2\pi),\theta\in(0,2\pi)
$$

由于$\Theta,T$是独立的

先对$T$积分得到$\Theta$的PDF

$$
f_\Theta(\theta) = \int_0^{2\pi} \frac{R+rcos(t)}{4\pi^2R} dt = \frac{1}{2\pi}
$$

求出$\Theta$的CDF

$$
F_\Theta(\theta) = \frac{\theta}{2\pi}
$$

令$\xi_1$等于CDF

$$
\frac{\theta}{2\pi} = \xi_1
$$

可得

$$
\theta = 2\pi \xi_1
$$

之后对$\Theta$求积分得到$T$的边缘PDF

$$
f_T(t) = \frac{R+rcos(t)}{2\pi R},\quad t\in(0,2\pi)
$$

再求出$T$的CDF

$$
F_T(t) = \frac{Rt+rsin(t)}{2\pi R},\quad t\in(0,2\pi)
$$

让其等于$\xi_1$

$$
\frac{Rt+rsin(t)}{2\pi R} = \xi_1
$$

结果发现是超越方程，$t$没有解析解（笑）。当然每次采样的时候用计算机求出方程的数值解不是不行，但是太麻烦且效率低下了。所以这里继续采用一个简单的拒绝采样，建议函数是$(0,2\pi)$上的均匀分布（虽然效率依然底下但是代码好写），参数$c = 2\pi max\{f_T(t)\} = 2\pi\frac{R+r}{2\pi R}$，最终伪代码如下

```python
for i in range(sampleNum):
    theta = 2*np.pi*np.random.rand()
    flag = 0
    while flag != 1:
        t = 2*np.pi*np.random.rand()
        if (R+r)*np.random.rand() < R+r*np.cos(t):
            X.append((R+r*np.cos(t))*np.cos(theta))
            Y.append((R+r*np.cos(t))*np.sin(theta))
            Z.append(r*np.sin(t))
            flag = 1
```

画图结果如下，左图是按参数均匀采样，右图是按长度均匀采样：

<img src="\assets\img\SampleOnManifold\torus.png" alt="torus" style="zoom:60%;" />

视觉上仍然没什么差别。我们继续通过**蒙特卡洛积分**验证。

假设环面的环半径$R = 2$，管半径$r = 1$，考虑以下曲面积分：

$$
\int_M x^2 dS
$$

计算机计算可得

$$
\int_M x^2 dS = \int_\Omega \gamma_x(\theta,t)^2\sqrt{det(J_\gamma^TJ_\gamma)} =\\
\int_0^{2\pi}\int_0^{2\pi}((R+rcos(t))cos(\theta))^2\ r(R+rcos(t))d\theta dt \\ \approx 217.131296824
$$

而蒙特卡洛积分的公式是

$$
\frac{1}{N}\sum_{i=1}^N\frac{f(\mathbf{v}_i)}{p(\mathbf{v}_i)} = \frac{4\pi^2rR}{N}\sum_{i=1}^N \mathbf{v}_{iz} = \frac{8\pi^2}{N}\sum_{i=1}^N \mathbf{v}_{ix}
$$

因为我们假设是均匀采样所以$p(\mathbf{v}_i) = \frac{1}{S} = \frac{1}{4\pi^2rR}$。采样50000个点的结果如下

<img src="\assets\img\SampleOnManifold\torusintegral.png" alt="torusintegral" style="zoom:80%;" />

可以很明显地看出两者的差异。按照面积采样的值和理论值相差不大。

计算一下按照参数均匀采样的错误蒙特卡洛积分的期望：

流形上的概率分布是

$$
P_m(\gamma(\vec x)) = \frac{f_{\vec{X}}(\vec x)}{\sqrt{det(J_\gamma^TJ_\gamma)}} = \frac{1}{4\pi^2r(R+rcos(t))}
$$

因为参数空间是均匀采样的，$f_{\vec{X}}(\vec x) = \frac{1}{4\pi^2}$。接着计算期望

$$
E[\frac{f(\mathbf{v}_i)}{p_{uniform}(\mathbf{v}_i)}] = \int_M \frac{f(\mathbf{v}_i)}{p_{uniform}(\mathbf{v}_i)} p_{actual}(\mathbf{v}_i) dS \\
= \int_0^{2\pi}\int_0^{2\pi} \frac{((R+rcos(t))cos(\theta))^2}{1/4\pi^2rR}\cdot \frac{1}{4\pi^2r(R+rcos(t))} \cdot r(R+rcos(t)) dtd\theta\\
= \int_0^{2\pi}\int_0^{2\pi} Rr((R+rcos(t))cos(\theta))^2 dtd\theta = \frac{L}{2} \\
\approx 177.65287922
$$

和按参数均匀采样的计算结果十分接近。

将概率分布映射到环面上，得到下图

<img src="\assets\img\SampleOnManifold\torus3.png" alt="torus3" style="zoom:67%;" />

可以看到中间采样频率高些，但是整体差异不大，因此采样结果在视觉上比较均匀。

## 函数图像上的均匀采样

对于一个函数 $f: D \rightarrow \mathbb{R}^m, D \subseteq \mathbb{R}^n$, 其**图像(graph)**被定义为如下的点集：

$$
G_f = \{(\textbf{x},f(\textbf{x}));\textbf{x}\in D\} \subseteq \mathbb{R}^n \times \mathbb{R}^m = \mathbb{R}^{n+m}
$$

不难理解，对于一元函数，其图像就是二维平面中的一条曲线；对于二元函数，其图像就是三维空间中的一个曲面。函数图像很自然的就有一个平凡的参数化：$\gamma(\textbf{x}) = (\textbf{x},f(\textbf{x}))$

其雅可比行列式为

$$
J_\gamma = \begin{pmatrix}
\textbf{I}_n\\
\nabla f(\textbf{x})^T\end{pmatrix}
$$

比例系数的平方

$$
det(J_\gamma^TJ_\gamma) = det\big[ 
\begin{pmatrix}\textbf{I}_n & \nabla f(\textbf{x})\end{pmatrix} 
\begin{pmatrix}
\textbf{I}_n\\
\nabla f(\textbf{x})^T\end{pmatrix}
\big]\\
=det(I_n+\nabla f(\textbf{x})\nabla f(\textbf{x})^T)\\
=1+|\nabla f(\textbf{x})|^2
$$

这里的最后一步运用了**matrix determinant lemma**，这里不做赘述。

因此，参数空间的pdf就有这样一个通式

$$
f_{\vec{X}}(\vec x) = \frac{\sqrt{1+|\nabla f(\vec x)|^2}}{\int_\Omega \sqrt{1+|\nabla f(\vec x)|^2}\ dx^n},\vec x \in \Omega
$$

现在找一个栗子验证一下：在抛物面上均匀采样

$$
\{(x,y,x^2+y^2)| (x,y)\in [-L,L]^2\}
$$

这块曲面的面积是

$$
\int_\Omega \sqrt{1+|\nabla f(\vec x)|^2}\ dx^n = \int_{-L}^L\int_{-L}^{L} \sqrt{1+4x^2+4y^2}\ dxdy
$$

参数空间的pdf如下

$$
f_{X,Y}(x,y) = \frac{\sqrt{1+4x^2+4y^2}}{\int_{-L}^L\int_{-L}^{L} \sqrt{1+4x^2+4y^2}\ dxdy}
$$

由于又不是很好求CDF，更别提其逆了，所以用拒绝采用法，建议函数就取$[-L,L]^2$上的均匀分布，参数$c = max\{f_{X,Y}(x,y)\} = \sqrt{1+8L^2}/S$, 其中$S$是曲面面积。在拒绝采用的过程中，$S$可以被约去。最终代码如下:

```python
for i in range(sampleNum):
    x = L * (2 * np.random.rand() - 1)
    y = L * (2 * np.random.rand() - 1)
    if np.sqrt(1+8*L**2)*np.random.rand() <= np.sqrt(1+4*x**2+4*y**2):
        X.append(x)
        Y.append(y)
        Z.append(x**2+y**2)
```

得到这样的结果，左图是按参数均匀采样，右图是按长度均匀采样：

<img src="\assets\img\SampleOnManifold\graph.png" alt="graph" style="zoom:67%;" />

可以看出按参数采样的中心部分要密一些。接下来我们用一个曲面积分验证：

$$
\int_M x^2 dS
$$

计算器计算可得

$$
\int_M x^2 dS = \int_\Omega \gamma_x(\theta,t)^2\sqrt{det(J_\gamma^TJ_\gamma)} =\\
\int_{-L}^{L}\int_{-L}^{L}x^2\sqrt{1+4x^2+4y^2} dx dy \\ \approx  2.85711
$$

蒙特卡洛积分的公式是

$$
\frac{1}{N}\sum_{i=1}^N\frac{f(\mathbf{v}_i)}{p(\mathbf{v}_i)} = \frac{S}{N}\sum_{i=1}^N \mathbf{v}_{ix}
$$

因为我们假设是均匀采样所以$p(\mathbf{v}_i) = \frac{1}{S}$。采样50000个点的结果如下

<img src="\assets\img\SampleOnManifold\graphintegral.png" alt="graphintegral" style="zoom:80%;" />

可以很明显地看出两者的差异。按照面积采样的值和理论值相差不大。

计算一下按照参数均匀采样的错误蒙特卡洛积分的期望：

流形上的概率分布是

$$
P_m(\gamma(\vec x)) = \frac{f_{\vec{X}}(\vec x)}{\sqrt{det(J_\gamma^TJ_\gamma)}} = \frac{1}{4L^2\sqrt{1+4x^2+4y^2}}
$$

因为参数空间是均匀采样的，$f_{\vec{X}}(\vec x) = \frac{1}{4L^2}$。接着计算期望

$$
E[\frac{f(\mathbf{v}_i)}{p_{uniform}(\mathbf{v}_i)}] = \int_M \frac{f(\mathbf{v}_i)}{p_{uniform}(\mathbf{v}_i)} p_{actual}(\mathbf{v}_i) dS \\
= \int_{-L}^{L}\int_{-L}^{L} \frac{x^2}{1/S}\cdot \frac{1}{4L^2\sqrt{1+4x^2+4y^2}} \cdot \sqrt{1+4x^2+4y^2} dxdy\\
= \int_{-L}^{L}\int_{-L}^{L} \frac{Sx^2}{4L^2} dt = \frac{LS}{3} \\
\approx 2.4820855743374532
$$

和按参数均匀采样的计算结果十分接近。

将概率分布映射到曲面上，得到下图

<img src="\assets\img\SampleOnManifold\graph3bar.png" alt="graph3bar" style="zoom: 67%;" />

正如按视觉效果猜想的，中心部分的密度比较高。

## 结

在阅读了部分网络上的文章之后（总的来说相关的资料比较少），我受到了比较多的启发。[鸡哥大佬的文章](https://zhuanlan.zhihu.com/p/49746076)给出了球面上采样的详细推导过程，[PBRT 13.6章](http://www.pbr-book.org/3ed-2018/Monte_Carlo_Integration/2D_Sampling_with_Multidimensional_Transformations.html) 给出了丰富的例子，[Yivan大佬的文章](https://zhuanlan.zhihu.com/p/26052376)提到了流形上概率分布的概念，也是本文想要探讨的概念之一。但是总体来说我看到的资料要么是不够通用，要么没有太详细的解释（我认为这个问题对于大佬来说应该是trivial的）。动用多元微积分和概率论里学到的知识，我算是自己探究出了一个比较通用的在嵌入$\mathbb{R}^n$中的子流形上均匀采样的方法（虽然对于学数学的人来说应该也就适用于toy examples），自行计算验证过以后发现和PBRT里给出的公式吻合。最后我在外网找到了[相关的论文](https://statweb.stanford.edu/~cgates/PERSI/papers/sampling11.pdf)（一般来说这种比较trivial问题肯定很早就有人研究过了，可惜中文互联网上资料还是比较匮乏），大致也验证了我方法的正确性。

题外话就是，这篇文章其实在去年就只剩最后一个栗子没有加上了，可惜我是一个终极拖延症患者，一个不小心半年就过去了。这半年，我写了个类Linux内核，了解了基本的Machine Learning和Deep Learning，学了高阶的信号处理，，，现在放暑假了都，总算是填完了这个坑，可喜可贺，可喜可贺。

