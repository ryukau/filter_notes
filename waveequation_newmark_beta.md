<style>
body {
  max-width: 704px;
  margin: auto;
  padding: 32px 8px;
}

canvas {
  /* image-rendering: pixelated; */
  display: inline-block;
  border-style: solid;
  border-width: 0.2px;
  border-color: #000000;
}

.numberInputLabel {
  display: inline-block;
  text-align: right;
  width: 160px;
  margin: 0 8px 4px 0;
  padding-right: 8px;
  border-right: solid 3px #606060;
}

.numberInputRange {
  display: inline-block;
  max-width: 220px;
  width: 50%;
  margin: 0 8px 0 0;
}

.numberInputNumber {
  display: inline-block;
  max-width: 100px;
  width: 20%;
  margin: 0;
}

code {
  overflow-x: scroll;
  overflow-y: hidden;
  white-space: pre;
}

.katex {
  font-size: 1.3em !important;
}
</style>

# Newmark-β 法で1次元の波動方程式
## 微分の表記
微分の表記を短く書けるものに変えます。時間の微分に[ニュートンの表記](https://en.wikipedia.org/wiki/Notation_for_differentiation#Newton's_notation)を使います。

$$
\dot{u} = \frac{\partial u}{\partial t}
,\quad
\ddot{u} = \frac{\partial^2 u}{\partial t^2}
,\quad
\ldots
$$

空間方向の微分に[ラグランジュの表記](https://en.wikipedia.org/wiki/Notation_for_differentiation#Lagrange's_notation)を使います。

$$
u' = \frac{\partial u}{\partial x}
,\quad
u'' = \frac{\partial^2 u}{\partial x^2}
,\quad
\ldots
$$

離散化で出てくる差分を $\Delta_t$ や $\Delta_x$ と表記します。以前の記事では $dt$ や $dx$ と書いていましたが、乗算が $a b c dt$ のようになると区切りが分からなくなるので変えました。

$u$ のインデックスを下付き文字で表記します。

$$
\begin{aligned}
u_{n} &= u(x, t + n\,\Delta_t)\\
u_{n+1} &= u(x, t + (n + 1)\,\Delta_t)\\
u_{n, x+1} &= u(x + \Delta_x, t + n\,\Delta_t)\\
u_{n, x-1} &= u(x - \Delta_x, t + (n + 1)\,\Delta_t)\\
\end{aligned}
$$

## 波動方程式
1次元の[波動方程式](https://en.wikipedia.org/wiki/Wave_equation)です。

$$
\ddot{u} = c^2 u''
$$

波動方程式はエネルギーの減衰が考慮されていないので時間が経っても波が減衰しないはずですが [Implicit FDM](https://en.wikipedia.org/wiki/Finite_difference_method#Implicit_method) で離散化すると減衰が起こります。 Schweickart, James, Marschner の ["Animating Elastic Rods with Sound"](https://www.cs.cornell.edu/projects/rodsound/) に習って [Newmark-β 法](https://en.wikipedia.org/wiki/Newmark-beta_method)を使うことでエネルギーの減衰を改善します。

## Newmark-β 法
Newmark-β 法では速度と位置の計算に次の式を使います。

$$
\begin{aligned}
\dot{u}_{n+1} &= \dot{u}_{n}
+ \Delta_t \left( (1 - \gamma) \ddot{u}_{n} + \gamma \ddot{u}_{n+1} \right)\\

u_{n+1} &= u_{n} + \Delta_t \dot{u}_{n} + \Delta_t^2 \left(
  \left( \frac{1}{2} - \beta \right) \ddot{u}_{n} + \beta \ddot{u}_{n+1}
\right)
\end{aligned}
$$

$\beta$ と $\gamma$ が調整できる形になっていますが $\gamma$ については $1/2$ 以外の値にすると精度が下がるそうです。 $\beta$ の値は[Wikipediaの記事](https://en.wikipedia.org/wiki/Newmark-beta_method)で $1/4$ 、 Gavin の ["Numerical Integration in Structural Dynamics"](http://people.duke.edu/~hpgavin/cee541/NumericalIntegration.pdf) で $1/6$ の場合が紹介されていました。この文章では $\gamma$ は $1/2$ で固定して $\beta$ は調整できる形にしています。

## Newmark-β 法による波動方程式の離散化
位置 $u$ 、速度 $\dot{u}$ 、加速度 $\ddot{u}$ からなる連立微分方程式を立てます。位置と速度の式は Newmark-β 法のもので、加速度の式は1次元の波動方程式です。

$$
\begin{aligned}
\ddot{u}_{n+1} &= c^2 u''_{n+1}\\

\dot{u}_{n+1} &= \dot{u}_{n}
+ \frac{\Delta_t}{2} \left( \ddot{u}_{n} + \ddot{u}_{n+1} \right)\\

u_{n+1} &= u_{n} + \Delta_t \dot{u}_{n} + \Delta_t^2 \left(
  \left( \frac{1}{2} - \beta \right) \ddot{u}_{n} + \beta \ddot{u}_{n+1}
\right)
\end{aligned}
$$

$\ddot{u}_{n+1}$ の右辺について $u''_{n+1}$ を有限差分に変形して $u_{n+1}$ の右辺を代入します。

$$
\begin{aligned}
\ddot{u}_{n+1} &= c^2 u''_{n+1}\\
&= \frac{c^2}{\Delta_x^2} (u_{n+1, x-1} - 2 u_{n+1} + u_{n+1, x+1})\\
&= \frac{c^2}{\Delta_x^2} \Biggl(
  (u_{n, x-1} - 2 u_{n} + u_{n, x+1})
  + \Delta_t (
    \dot{u}_{n, x-1} - 2 \dot{u}_{n} + \dot{u}_{n, x+1})\\
  &\qquad\qquad+ \Delta_t^2 \left( \frac{1}{2} - \beta \right) (
    \ddot{u}_{n, x-1} - 2 \ddot{u}_{n} + \ddot{u}_{n, x+1})\\
  &\qquad\qquad+ \Delta_t^2 \beta (
    \ddot{u}_{n+1, x-1} - 2 \ddot{u}_{n+1} + \ddot{u}_{n+1, x+1})
\Biggr)\\
\end{aligned}
$$

$n+1$ の項を左辺に移項します。

$$
\begin{aligned}
& \left(\frac{\Delta_x^2}{c^2} + 2 \Delta_t^2 \beta \right) \ddot{u}_{n+1}
- \Delta_t^2 \beta (
  \ddot{u}_{n+1, x-1} + \ddot{u}_{n+1, x+1})\\
&= (u_{n, x-1} - 2 u_{n} + u_{n, x+1})
+ \Delta_t (
  \dot{u}_{n, x-1} - 2 \dot{u}_{n} + \dot{u}_{n, x+1})\\
&\quad+ \Delta_t^2 \left( \frac{1}{2} - \beta \right) (
  \ddot{u}_{n, x-1} - 2 \ddot{u}_{n} + \ddot{u}_{n, x+1})\\
\end{aligned}
$$

整理します。

$$
\begin{aligned}
C_0 \ddot{u}_{n+1}
+ C_1 (
  \ddot{u}_{n+1, x-1} + \ddot{u}_{n+1, x+1})
&= (u_{n, x-1} - 2 u_{n} + u_{n, x+1})\\
&\quad+ \Delta_t (
  \dot{u}_{n, x-1} - 2 \dot{u}_{n} + \dot{u}_{n, x+1})\\
&\quad+ C_2 (
  \ddot{u}_{n, x-1} - 2 \ddot{u}_{n} + \ddot{u}_{n, x+1})\\
\end{aligned}
$$

境界となる左右の端では式の形が変わります。次の式は左端の式です。 $x+1$ を $x-1$ に変えると右端の式になります。 $b=1$ のとき固定端、 $b=2$ のとき自由端になります。

$$
\begin{aligned}
C_0 \ddot{u}_{n+1} + b\, C_1 \ddot{u}_{n+1, x+1}
&= (b\, u_{n, x+1} - 2 u_{n})\\
&\quad+ \Delta_t (b\, \dot{u}_{n, x+1} - 2 \dot{u}_{n})\\
&\quad+ C_2 (b\, \ddot{u}_{n, x+1} - 2 \ddot{u}_{n})
\end{aligned}
$$

連立方程式を立てます。

$$
\begin{aligned}
\begin{cases}
\mathbf{A} \ddot{\mathbf{u}}_{n+1}
= \mathbf{u}''_n + \Delta_t \dot{\mathbf{u}}''_{n} + C_2 \ddot{\mathbf{u}}''_{n}\\
u_{n+1} = u_{n} + \Delta_t \dot{u}_{n} + C_2 \ddot{u}_{n} - C_1 \ddot{u}_{n+1}\\
\dot{u}_{n+1} = \dot{u}_{n} + C_3 \left( \ddot{u}_{n} + \ddot{u}_{n+1} \right)\\
\end{cases}
\;,&\quad
\mathbf{A} =
\begin{bmatrix}
 C_0 & b_L C_1 & 0 & 0 & \cdots & 0 & 0\\
 C_1 & C_0 & C_1 & 0 & \cdots & 0 & 0\\
 0 & C_1 & C_0 & C_1 & \cdots & 0 & 0\\
& \vdots & & & & \vdots &\\
 0 & 0 & 0 & 0 & \cdots & b_R C_1 & C_0
\end{bmatrix}\\

C_0 = \frac{\Delta_x^2}{c^2} + 2 \Delta_t^2 \beta
,\quad
C_1 = - \Delta_t^2 \beta
,&\quad
C_2 = \frac{\Delta_t^2}{2} + C_1
,\quad
C_3 = \frac{\Delta_t}{2}
\end{aligned}
$$

実装はGitHubの別ページに分けました。リンク先の `Wave1DNewmarkBeta` クラスになります。

[実装を読む (github.com)](https://github.com/ryukau/filter_notes/blob/b7eb2f3c2258cd587bc24f7984b6b8761e3e414f/docs/demo/waveequation_newmark_beta/wave1d.js#L81)

## Newmark-β 法によるばね-ダンパ波動方程式の離散化
波動方程式にばねの項 $k u$ とダンパの項 $a \dot{u}$ を加えることで波の減衰を調整できるようにします。この文章では次の式をばね-ダンパ波動方程式と呼ぶことにします。

$$
\ddot{u} + a \dot{u} + k u = c^2 u''
$$

Newmark-β 法の連立方程式を立てます。

$$
\begin{aligned}
\ddot{u}_{n+1} &= c^2 u''_{n+1} - a \dot{u}_{n+1} + k u_{n+1}\\

\dot{u}_{n+1} &= \dot{u}_{n}
+ \frac{\Delta_t}{2} \left( \ddot{u}_{n} + \ddot{u}_{n+1} \right)\\

u_{n+1} &= u_{n} + \Delta_t \dot{u}_{n} + \Delta_t^2 \left(
  \left( \frac{1}{2} - \beta \right) \ddot{u}_{n} + \beta \ddot{u}_{n+1}
\right)
\end{aligned}
$$

$\ddot{u}_{n+1}$ の右辺について $\dot{u}_{n+1}$ と $u_{n+1}$ を代入します。

$$
\begin{aligned}
\ddot{u}_{n+1}
&= \frac{c^2}{\Delta_x^2} \Biggl(
  (u_{n, x-1} - 2 u_{n} + u_{n, x+1})
  + \Delta_t (
    \dot{u}_{n, x-1} - 2 \dot{u}_{n} + \dot{u}_{n, x+1})\\
  &\qquad\qquad+ \Delta_t^2 \left( \frac{1}{2} - \beta \right) (
    \ddot{u}_{n, x-1} - 2 \ddot{u}_{n} + \ddot{u}_{n, x+1})\\
  &\qquad\qquad+ \Delta_t^2 \beta (
    \ddot{u}_{n+1, x-1} - 2 \ddot{u}_{n+1} + \ddot{u}_{n+1, x+1})
\Biggr)\\
&\quad- a \Biggl(
  \dot{u}_{n} + \frac{\Delta_t}{2} \left( \ddot{u}_{n} + \ddot{u}_{n+1} \right)
\Biggr)\\
&\quad- k \Biggl(
  u_{n} + \Delta_t \dot{u}_{n} + \Delta_t^2
    \left( \frac{1}{2} - \beta \right) \ddot{u}_{n} + \Delta_t^2 \beta \ddot{u}_{n+1}
\Biggr)
\end{aligned}
$$

$n+1$ の項を左辺に移項します。

$$
\begin{aligned}
& (
  1
  + a \frac{\Delta_t}{2}
  + k \Delta_t^2 \beta
  + 2 \frac{c^2}{\Delta_x^2} \Delta_t^2 \beta
) \ddot{u}_{n+1}
- \frac{c^2}{\Delta_x^2} \Delta_t^2 \beta (
    \ddot{u}_{n+1, x-1} + \ddot{u}_{n+1, x+1})\\
&= \frac{c^2}{\Delta_x^2} \Biggl(
  (u_{n, x-1} + u_{n, x+1})
  + \Delta_t (\dot{u}_{n, x-1} + \dot{u}_{n, x+1})
  + \Delta_t^2 \left( \frac{1}{2} - \beta \right) (
    \ddot{u}_{n, x-1} + \ddot{u}_{n, x+1})
\Biggr)\\
&\quad- \Biggl(
  a \frac{\Delta_t}{2}
  + k \Delta_t^2 \left( \frac{1}{2} - \beta \right)
  + 2 \frac{c^2}{\Delta_x^2} \Delta_t^2 \left( \frac{1}{2} - \beta \right)
\Biggr) \ddot{u}_{n}\\
&\quad- \Biggl(
  a
  + k \Delta_t
  + 2 \frac{c^2}{\Delta_x^2} \Delta_t
\Biggr) \dot{u}_{n} \\
&\quad- \Biggl(
  k
  + 2 \frac{c^2}{\Delta_x^2}
\Biggr) u_{n}
\end{aligned}
$$

整理します。

$$
\begin{aligned}
& C_0 \ddot{u}_{n+1}
+ C_1 (
    \ddot{u}_{n+1, x-1} + \ddot{u}_{n+1, x+1})\\
&= C_2 \bigl(
  (u_{n, x-1} + u_{n, x+1})
  + \Delta_t (\dot{u}_{n, x-1} + \dot{u}_{n, x+1})
  + C_3 (\ddot{u}_{n, x-1} + \ddot{u}_{n, x+1})
\bigr)\\
&\quad- C_4 \ddot{u}_{n} - C_5 \dot{u}_{n} - C_6 u_{n}\\
\end{aligned}
$$

境界となる左右の端では式の形が変わります。次の式は左端の式です。 $x+1$ を $x-1$ に変えると右端の式になります。 $b=1$ のとき固定端、 $b=2$ のとき自由端になります。

$$
\begin{aligned}
C_0 \ddot{u}_{n+1} + b\,C_1 \ddot{u}_{n+1, x+1}
&= b\,C_2 \bigl(
  u_{n, x+1}
  + \Delta_t \dot{u}_{n, x+1}
  + C_3 \ddot{u}_{n, x+1}
\bigr)\\
&\quad- C_4 \ddot{u}_{n} - C_5 \dot{u}_{n} - C_6 u_{n}\\
\end{aligned}
$$

連立方程式を立てます。

$$
\begin{aligned}
C_0 &= 1 + a C_7 + C_6 C_8\\
C_1 &= - C_2 C_8\\
C_2 &= c^2 / \Delta_x^2\\
C_3 &= \Delta_t^2 \left( 1 / 2 - \beta \right)\\
C_4 &= C_3 C_6 + a C_7\\
C_5 &= a + \Delta_t C_6\\
C_6 &= k + 2 C_2\\
C_7 &= \Delta_t / 2\\
C_8 &= \Delta_t^2 \beta\\
\end{aligned}
\qquad
\begin{aligned}
&\quad \begin{cases}
\mathbf{A} \ddot{\mathbf{u}}_{n+1} = \mathbf{b}\\
\dot{u}_{n+1} = \dot{u}_{n} + C_7 \left( \ddot{u}_{n} + \ddot{u}_{n+1} \right)\\
u_{n+1} = u_{n} + \Delta_t \dot{u}_{n} + C_3 \ddot{u}_{n} + C_8 \ddot{u}_{n+1}\\
\end{cases}\\
\\
\mathbf{A} &=
\begin{bmatrix}
 C_0 & b_L C_1 & 0 & 0 & \cdots & 0 & 0\\
 C_1 & C_0 & C_1 & 0 & \cdots & 0 & 0\\
 0 & C_1 & C_0 & C_1 & \cdots & 0 & 0\\
& \vdots & & & & \vdots &\\
 0 & 0 & 0 & 0 & \cdots & b_R C_1 & C_0
\end{bmatrix}\\
\\
\mathbf{b} &= C_2 \bigl(
  (\mathbf{u}_{n, x-1} + \mathbf{u}_{n, x+1})
  + \Delta_t (\dot{\mathbf{u}}_{n, x-1} + \dot{\mathbf{u}}_{n, x+1})\\
  &\qquad+ C_3 (\ddot{\mathbf{u}}_{n, x-1} + \ddot{\mathbf{u}}_{n, x+1})
\bigr)
- C_4 \ddot{\mathbf{u}}_{n} - C_5 \dot{\mathbf{u}}_{n} - C_6 \mathbf{u}_{n}\\
\end{aligned}
$$

実装はGitHubの別ページに分けました。リンク先の `DampedWave1DNewmarkBeta` クラスになります。

[実装を読む (github.com)](https://github.com/ryukau/filter_notes/blob/b7eb2f3c2258cd587bc24f7984b6b8761e3e414f/docs/demo/waveequation_newmark_beta/wave1d.js#L207)

## デモ
キャンバスをクリックすると波が起きます。

$\Delta_t = 1/60$ で固定しています。

ソルバに[ガウス-ザイデル法](https://en.wikipedia.org/wiki/Gauss%E2%80%93Seidel_method)を使っています。 $c \Delta_t / \Delta_x$ が大きいとソルバの収束が遅くなって発散します。ガウス-ザイデル法の反復を増やすことで発散を抑えています。

$\beta < 1/4$ のとき、他のパラメータの値によっては発散します。

<script src="demo/waveequation_newmark_beta/vec2.js"></script>
<script src="demo/waveequation_newmark_beta/canvas.js"></script>
<script src="demo/waveequation_newmark_beta/wave1d.js"></script>

## 参考サイト
- [Energy drift - Wikipedia](https://en.wikipedia.org/wiki/Energy_drift)
- [pde - Energy Conservation - Computational Science Stack Exchange](https://scicomp.stackexchange.com/questions/7202/energy-conservation)
