---
layout: post
title: ファインマン物理学 (3)電磁気学 補章 最小作用の原理
date: '2026-11-22 05:15:00 +0900'
categories: feynman
---
<!-- # ファインマン物理学 (3)電磁気学 -->

## 補章 最小作用の原理

### はじめに

光については、Fermatの原理が「光は所要時間の停留する経路を通る」と述べた。力学にも、これと同型の原理がある。**最小作用の原理**である。粒子がある時刻$t_1$に点1を出て、時刻$t_2$に点2へ着く運動を考えると、実際に選ばれる経路は、ある量を停留値(最小とは限らず、変形に対して1次の変化が消える値)にする経路である。その量を**作用**(action)という。

作用は、経路全体に対して定まる積分である。1次元の運動では、位置$x(t)$と速度$\dot{x}(t)$から作られる**ラグランジアン**(Lagrangian)

$$
L=T-U=\frac{1}{2}m\dot{x}^2-U(x)
$$

($T$は運動エネルギー、$U$はポテンシャルエネルギー)を用いて、作用$S$は

$$
S=\int_{t_1}^{t_2}L\bigl(x(t),\dot{x}(t),t\bigr)\,dt
$$

と定義される。問題は、「端点$t_1$、$t_2$における位置を固定したすべての経路のうち、$S$を停留にするものが実際の運動である」という形で述べられる。ここがNewtonの運動方程式との見方の違いで、運動方程式が各瞬間の加速度を局所的に決めるのに対して、作用原理は経路全体を一挙に選ぶ大域的な述べ方になっている。しかも両者は同じ運動を与える(以下で示す)。

※高校時代にBader先生からこの原理を教わったという原典の導入の逸話は、内容の理解には関わらないので省略する。要点は、自然法則が「積分の停留」という形に書けることへの驚きである。

### 試行経路と変分

停留条件を調べる方法は、Fermatの原理で近くの経路と比べたのと同じである。真の経路$x(t)$の近くに、端点を動かさない試行経路

$$
x'(t)=x(t)+\delta x(t),\qquad \delta x(t_1)=\delta x(t_2)=0
$$

を取る。$\delta x(t)$は任意の小さい関数で、これを**変分**という。速度の変分は$\delta\dot{x}=\dfrac{d}{dt}\delta x$である。作用の差を$\delta x$の1次まで展開すると、

$$
\begin{aligned}
S[x']-S[x]
&=\int_{t_1}^{t_2}\left\{L(x+\delta x,\dot{x}+\delta\dot{x},t)-L(x,\dot{x},t)\right\}dt\\
&\approx\int_{t_1}^{t_2}\left\{\frac{\partial L}{\partial x}\delta x+\frac{\partial L}{\partial\dot{x}}\delta\dot{x}\right\}dt
\end{aligned}
$$

となる。近似はTaylor展開の1次の項だけを残したもので、落とした項は$(\delta x)^2$、$(\delta\dot{x})^2$の大きさであり、「停留」とはこの1次の項が消えることだ、というのが定義そのものである。停留条件は

$$
\delta S=\int_{t_1}^{t_2}\left\{\frac{\partial L}{\partial x}\delta x+\frac{\partial L}{\partial\dot{x}}\frac{d}{dt}\delta x\right\}dt=0
$$

である。

### オイラー-ラグランジュ方程式

第2項を部分積分する。$u=\dfrac{\partial L}{\partial\dot{x}}$、$dv=\dfrac{d}{dt}\delta x\,dt$と置くと、

$$
\int_{t_1}^{t_2}\frac{\partial L}{\partial\dot{x}}\frac{d}{dt}\delta x\,dt
=\left[\frac{\partial L}{\partial\dot{x}}\delta x\right]_{t_1}^{t_2}
-\int_{t_1}^{t_2}\frac{d}{dt}\left(\frac{\partial L}{\partial\dot{x}}\right)\delta x\,dt
$$

境界項は、端点で変分を固定した$\delta x(t_1)=\delta x(t_2)=0$によって消える。残りを元の式に代入すると、

$$
\delta S=\int_{t_1}^{t_2}\left\{\frac{\partial L}{\partial x}-\frac{d}{dt}\left(\frac{\partial L}{\partial\dot{x}}\right)\right\}\delta x(t)\,dt=0
$$

$\delta x(t)$は区間の内部では**任意**に選べる関数である。任意の$\delta x$に対して積分が0になるためには、被積分関数の括弧の中身がすべての時刻で0でなければならない(もしある時刻の近くで正なら、その近傍だけ$\delta x$を正に取ることで積分を正にできてしまう)。したがって

$$
\frac{d}{dt}\left(\frac{\partial L}{\partial\dot{x}}\right)-\frac{\partial L}{\partial x}=0
$$

これが**オイラー-ラグランジュ方程式**(Euler–Lagrange equation)である。3次元では$x$、$y$、$z$の各成分について同じ式が成り立つ。

通常の力学で確かめる。$L=\dfrac{1}{2}m\dot{x}^2-U(x)$なら、

$$
\frac{\partial L}{\partial\dot{x}}=m\dot{x},\qquad
\frac{d}{dt}\left(\frac{\partial L}{\partial\dot{x}}\right)=m\ddot{x},\qquad
\frac{\partial L}{\partial x}=-\frac{dU}{dx}
$$

を代入して、

$$
m\ddot{x}=-\frac{dU}{dx}=F
$$

すなわちNewtonの第2法則が戻る。作用原理は新しい運動を予言するのではなく、運動法則の別の(しかし一般化しやすい)定式化である。

### 例：自由粒子と重力

**自由粒子**：$U=0$、$L=\dfrac{1}{2}m\dot{x}^2$なら、オイラー-ラグランジュ方程式は

$$
\frac{d}{dt}(m\dot{x})=0
$$

となり、速度は一定、経路は直線である。作用の言葉では、与えられた時間で2点を結ぶとき、速度のばらつきは運動エネルギーの時間積分を大きくする(2乗の平均は平均の2乗以上)ため、等速直線が停留経路になることに対応する。

**一様な重力**：鉛直方向を$y$として$L=\dfrac{1}{2}m(\dot{x}^2+\dot{y}^2)-mgy$。$x$成分については

$$
\frac{d}{dt}(m\dot{x})=0\quad\Longrightarrow\quad \dot{x}=\text{const.}
$$

$y$成分については$\dfrac{\partial L}{\partial y}=-mg$だから、

$$
\frac{d}{dt}(m\dot{y})+mg=0\quad\Longrightarrow\quad \ddot{y}=-g
$$

2回積分して$y=y_0+v_{y0}t-\dfrac{1}{2}gt^2$、放物線が得られる。原典が強調するのは、経路が放物線になること自体よりも、微分方程式を解かずに「積分の停留」という一段高い視点から運動全体を特徴づけられる点にある。

### 電磁場の中の荷電粒子

作用原理の真価は、力の形が複雑な場合にもラグランジアンさえ書ければ運動方程式が出る点にある。電磁場の中の電荷$q$については、スカラーポテンシャル$\phi$とベクトルポテンシャル$\mathbf{A}$を用いて、

$$
L=\frac{1}{2}m\mathbf{v}^2-q\phi+q\,\mathbf{v}\cdot\mathbf{A}
$$

と取る。成分$x_i$についてオイラー-ラグランジュ方程式を全部書き下す。まず速度で微分すると、

$$
\frac{\partial L}{\partial v_i}=mv_i+qA_i
$$

この組み合わせ$m\mathbf{v}+q\mathbf{A}$が**正準運動量**(canonical momentum)である。時間微分は、$\mathbf{A}$が位置と時刻の両方に依存することに注意して、

$$
\frac{d}{dt}\left(\frac{\partial L}{\partial v_i}\right)
=m\dot{v}_i+q\frac{\partial A_i}{\partial t}+q\sum_j v_j\frac{\partial A_i}{\partial x_j}
$$

位置で微分すると、

$$
\frac{\partial L}{\partial x_i}=-q\frac{\partial\phi}{\partial x_i}+q\sum_j v_j\frac{\partial A_j}{\partial x_i}
$$

オイラー-ラグランジュ方程式$\dfrac{d}{dt}\dfrac{\partial L}{\partial v_i}=\dfrac{\partial L}{\partial x_i}$に代入し、$q\dfrac{\partial A_i}{\partial t}$を左辺に残して整理すると、

$$
m\dot{v}_i=-q\frac{\partial\phi}{\partial x_i}-q\frac{\partial A_i}{\partial t}
+q\sum_j v_j\left(\frac{\partial A_j}{\partial x_i}-\frac{\partial A_i}{\partial x_j}\right)
$$

最後の括弧は磁場で書き直せる。例えば$i=x$の成分では、

$$
\sum_j v_j\left(\frac{\partial A_j}{\partial x}-\frac{\partial A_x}{\partial x_j}\right)
=v_y\left(\frac{\partial A_y}{\partial x}-\frac{\partial A_x}{\partial y}\right)
+v_z\left(\frac{\partial A_z}{\partial x}-\frac{\partial A_x}{\partial z}\right)
=v_y B_z-v_z B_y=(\mathbf{v}\times\mathbf{B})_x
$$

ただし$B_z=\dfrac{\partial A_y}{\partial x}-\dfrac{\partial A_x}{\partial y}$など、$\mathbf{B}=\nabla\times\mathbf{A}$を使った。他の成分も同様である。また$\mathbf{E}=-\nabla\phi-\dfrac{\partial\mathbf{A}}{\partial t}$だから、結局

$$
m\dot{\mathbf{v}}=q\bigl(\mathbf{E}+\mathbf{v}\times\mathbf{B}\bigr)
$$

すなわち**Lorentz力**が得られる。場を$\mathbf{E}$、$\mathbf{B}$ではなくポテンシャルで表すことが、作用原理では本質的である。

### 相対論的な作用と一般化

相対論では、自由粒子のラグランジアンは

$$
L=-mc^2\sqrt{1-\frac{v^2}{c^2}}
$$

である。運動量は速度で微分して、

$$
p=\frac{\partial L}{\partial v}
=-mc^2\cdot\frac{1}{2\sqrt{1-v^2/c^2}}\left(-\frac{2v}{c^2}\right)
=\frac{mv}{\sqrt{1-v^2/c^2}}=\gamma mv
$$

となり、相対論的運動量が出る。非相対論的極限は、$v\ll c$で平方根を展開して確かめられる。

$$
\sqrt{1-\frac{v^2}{c^2}}\approx 1-\frac{v^2}{2c^2}
\quad\Longrightarrow\quad
L\approx -mc^2+\frac{1}{2}mv^2
$$

定数$-mc^2$は運動方程式に影響しないので、通常の$\dfrac{1}{2}mv^2$に戻る。近似の根拠は$v^2/c^2$の2次以降が無視できることである。電磁場がある場合は$-q\phi+q\mathbf{v}\cdot\mathbf{A}$を加えたものが作用になり、上の導出と合わせて相対論的なLorentz力の運動方程式を与える。

一般化はさらに進む。場そのものについても作用が書け、その停留条件としてMaxwell方程式が出る。粒子と場を合わせた作用を停留にすることが、古典電磁気学全体の1つの定式化になっている。

※場の作用の具体形からMaxwell方程式を導く計算は原典では結果の紹介にとどまるので、形だけ補っておく(補足)。場のラグランジアン密度を$\mathcal{L}=\dfrac{\epsilon_0}{2}\bigl(\mathbf{E}^2-c^2\mathbf{B}^2\bigr)-\rho\phi+\mathbf{j}\cdot\mathbf{A}$($\mathbf{E}$、$\mathbf{B}$はポテンシャルで表す)とし、$\phi$と$\mathbf{A}$を独立に変分すると、それぞれGaussの法則とAmpère–Maxwellの法則がオイラー-ラグランジュ方程式として得られる。

量子力学とのつながりも原典が触れる点である。量子力学では、粒子は1つの経路だけでなくすべての経路を「試し」、各経路に作用$S$で決まる位相$e^{iS/\hbar}$を持つ振幅が割り当てられ、その和が到達確率を与える。$S$が$\hbar$に比べて大きい古典的極限では、停留経路の近傍だけ位相が揃って残り、それ以外は打ち消し合う。最小作用の原理は、この意味で量子力学の古典極限として理解できる(Fermatの原理が波動光学の極限だったのと正確に同じ構造である)。
