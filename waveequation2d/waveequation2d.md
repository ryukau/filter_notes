# 2次元のばね-ダンパ波動方程式
2次元のばね-ダンパ[波動方程式](https://en.wikipedia.org/wiki/Wave_equation)を実装します。

## 表記について
時間のインデックスの表記を $n$ から $t$ に変えます。 $n$ は次元の数に使うことにします。

## ラプラシアンの離散化
$n$ 次元の[ラプラシアン](https://en.wikipedia.org/wiki/Laplace_operator)を[有限差分](https://en.wikipedia.org/wiki/Finite_difference)の形で離散化します。

$$
\begin{aligned}
\nabla^{2} u_{t} &= \frac{1}{\Delta_x^{2}} \left(
  \mathtt{neighbor}(u_{t}) - 2n\,u_{t}
\right)\\
\mathtt{neighbor}(u_{t}) &= \sum_{i\,\in \mathtt{axis}(n)} (u_{t,i-1} + u_{t,i+1})
,\qquad
\mathtt{axis}(n) = \{x, y, z, \dots \}
\end{aligned}
$$

$\mathtt{neighbor}(u_{t})$ は[近傍](https://en.wikipedia.org/wiki/Von_Neumann_neighborhood)の総和です。 $\mathtt{axis}(n)$ は $n$ 次元の各軸の集合を表しているつもりですが、フォーマルな書き方がよく分からないのでごまかして書いています。

## 離散化
後で次元を変えられる形で離散化します。 $n$ 次元のばね-ダンパ波動方程式です。

$$
\ddot{u} + a \dot{u} + k u = c^{2} \nabla^{2} u
$$


[Newmark-β法](https://en.wikipedia.org/wiki/Newmark-beta_method)の連立方程式を立てます。

$$
\begin{aligned}
\ddot{u}_{t+1} &= c^{2} \nabla^{2} u_{t+1} - a \dot{u}_{t+1} + k u_{t+1}\\

\dot{u}_{t+1} &= \dot{u}_{t}
+ \frac{\Delta_t}{2} \left( \ddot{u}_{t} + \ddot{u}_{t+1} \right)\\

u_{t+1} &= u_{t} + \Delta_t \dot{u}_{t} + \Delta_t^{2} \left(
  \left( \frac{1}{2} - \beta \right) \ddot{u}_{t} + \beta \ddot{u}_{t+1}
\right)
\end{aligned}
$$

$\ddot{u}_{t+1}$ の右辺について $\dot{u}_{t+1}$ と $u_{t+1}$ を代入します。

$$
\begin{aligned}
\ddot{u}_{t+1}
&= c^{2} \nabla^{2} \left(
  u_{t} + \Delta_t \dot{u}_{t} + \Delta_t^{2} \left(
    \left( \frac{1}{2} - \beta \right) \ddot{u}_{t} + \beta \ddot{u}_{t+1}
  \right)
\right) - a \dot{u}_{t+1} + k u_{t+1}\\

&= \frac{c^{2}}{\Delta_x^{2}} \Biggl(
  ( \mathtt{neighbor}(u_{t}) - 2n\,u_{t})
  + \Delta_t (
    \mathtt{neighbor}(\dot{u}_{t}) - 2n\,\dot{u}_{t})\\
  &\qquad\qquad+ \Delta_t^{2} \left( \frac{1}{2} - \beta \right) (
    \mathtt{neighbor}(\ddot{u}_{t}) - 2n\,\ddot{u}_{t})\\
  &\qquad\qquad+ \Delta_t^{2} \beta (
    \mathtt{neighbor}(\ddot{u}_{t+1}) - 2n\,\ddot{u}_{t+1})
\Biggr)\\
&\quad- a \Biggl(
  \dot{u}_{t} + \frac{\Delta_t}{2} \left( \ddot{u}_{t} + \ddot{u}_{t+1} \right)
\Biggr)\\
&\quad- k \Biggl(
  u_{t} + \Delta_t \dot{u}_{t} + \Delta_t^{2}
    \left( \frac{1}{2} - \beta \right) \ddot{u}_{t} + \Delta_t^{2} \beta \ddot{u}_{t+1}
\Biggr)
\end{aligned}
$$

$t+1$ の項を左辺に移項します。

$$
\begin{aligned}
& (
  1
  + a \frac{\Delta_t}{2}
  + k \Delta_t^{2} \beta
  + 2n \frac{c^{2}}{\Delta_x^{2}} \Delta_t^{2} \beta
) \ddot{u}_{t+1}
- \frac{c^{2}}{\Delta_x^{2}} \Delta_t^{2} \beta \,\mathtt{neighbor}(\ddot{u}_{t+1})\\
&= \frac{c^{2}}{\Delta_x^{2}} \Biggl(
  \mathtt{neighbor}(u_{t})
  + \Delta_t \,\mathtt{neighbor}(\dot{u}_{t})
  + \Delta_t^{2} \left( \frac{1}{2} - \beta \right) \mathtt{neighbor}(\ddot{u}_{t})
\Biggr)\\
&\quad- \Biggl(
  a \frac{\Delta_t}{2}
  + k \Delta_t^{2} \left( \frac{1}{2} - \beta \right)
  + 2n \frac{c^{2}}{\Delta_x^{2}} \Delta_t^{2} \left( \frac{1}{2} - \beta \right)
\Biggr) \ddot{u}_{t}\\
&\quad- \Biggl(
  a
  + k \Delta_t
  + 2n \frac{c^{2}}{\Delta_x^{2}} \Delta_t
\Biggr) \dot{u}_{t} \\
&\quad- \Biggl(
  k
  + 2n \frac{c^{2}}{\Delta_x^{2}}
\Biggr) u_{t}
\end{aligned}
$$

整理します。

$$
\begin{aligned}
& C_0 \ddot{u}_{t+1}
+ C_1 \,\mathtt{neighbor}(\ddot{u}_{t+1})\\
&= C_2 \bigl(
  \mathtt{neighbor}(u_{t})
  + \Delta_t \,\mathtt{neighbor}(\dot{u}_{t})
  + C_3 \,\mathtt{neighbor}(\ddot{u}_{t})
\bigr)\\
&\quad- C_4 \ddot{u}_{t} - C_5 \dot{u}_{t} - C_6 u_{t}\\
\end{aligned}
$$

連立方程式を立てます。[1次元のとき](../waveequation_newmark_beta/waveequation_newmark_beta.html)と比べると、定数については $C_6$ だけが変わっています。

$$
\begin{aligned}
C_0 &= 1 + a C_7 + C_6 C_8\\
C_1 &= - C_2 C_8\\
C_2 &= c^{2} / \Delta_x^{2}\\
C_3 &= \Delta_t^{2} \left( 1 / 2 - \beta \right)\\
C_4 &= C_3 C_6 + a C_7\\
C_5 &= a + \Delta_t C_6\\
C_6 &= k + 2n\,C_2\\
C_7 &= \Delta_t / 2\\
C_8 &= \Delta_t^{2} \beta\\
\end{aligned}
\quad
\begin{aligned}
&\quad \begin{cases}
\mathbf{A} \ddot{\mathbf{u}}_{t+1} = \mathbf{b}\\
\dot{u}_{t+1} = \dot{u}_{t} + C_7 \left( \ddot{u}_{t} + \ddot{u}_{t+1} \right)\\
u_{t+1} = u_{t} + \Delta_t \dot{u}_{t} + C_3 \ddot{u}_{t} + C_8 \ddot{u}_{t+1}\\
\end{cases}\\
\\
\mathbf{b} &= C_2 \bigl(
  \mathtt{neighbor}(\mathbf{u}_{t})
  + \Delta_t \mathtt{neighbor}(\dot{\mathbf{u}}_{t})
  + C_3 \mathtt{neighbor}(\ddot{\mathbf{u}}_{t})
\bigr)\\
&\qquad- C_4 \ddot{\mathbf{u}}_{t} - C_5 \dot{\mathbf{u}}_{t} - C_6 \mathbf{u}_{t}\\
\end{aligned}
$$

### $\mathbf{A}$ と $\mathbf{u}_{t}$ の組み立て
$\mathbf{A}$ と $\mathbf{u}_{t}$ の組み立て方について2次元の例を作って見ていきます。 $4 \times 4$ の $u$ を考えます。

<figure>
<img src="img/2d_lattice.svg" alt="Image of gibbs phenomenon." style="width: 200px;"/>
</figure>

そのまま2次元配列として実装できますが、ソルバで解くときだけ1次元配列にする必要があります。次の図では横書きの文章を書くように、左から右に進み、右端に辿りついたら次の行の左端からまた進むというようにインデックスを振っています。

<figure>
<img src="img/2d_lattice_mod_index.svg" alt="Image of gibbs phenomenon." style="width: 360px;"/>
</figure>

```javascript
for (var x = 0; x < xLength; ++x) {
  for (var y = 0; y < yLength; ++y) {
    u_array[x + y * xLength] = u[x][y]
  }
}
```

これで $\mathbf{u}_{t}$ が組み立てられました。

$\mathbf{A}$ を組み立てます。整理した式の左辺にある関数 $\mathtt{neighbor}$ を展開して2次元の形にします。このときインデックスが配列の端から外に出ないよう例外処理が必要です。

$$
\begin{aligned}
C_0 \ddot{u}_{t+1} + C_1 \,\mathtt{neighbor}(\ddot{u}_{t+1})
&=
C_0 \ddot{u}_{t+1}
+ C_1 \sum_{i\,\in \mathtt{axis}(n)} (
  \ddot{u}_{t,i-1} + \ddot{u}_{t,i+1})\\
&=
C_0 \ddot{u}_{t+1}
+ C_1 (\ddot{u}_{t,x-1} + \ddot{u}_{t,x+1})
+ C_1 (\ddot{u}_{t,y-1} + \ddot{u}_{t,y+1})
\\
&=
C_0 \ddot{u}_{t+1}
+ C_1 \mathtt{sum\_x}(\ddot{u}_{t+1})
+ C_1 \mathtt{sum\_y}(\ddot{u}_{t+1})\\
\end{aligned}
$$

$$
\begin{aligned}
\mathtt{sum\_x}(u_{t})
&=
\begin{cases}
b u_{t,x - 1}                    &,\,x=0\\
b u_{t,x + 1}                    &,\,x=n_x - 1\\
u_{t,x-1} + u_{t,x+1} &,\,\text{otherwise}
\end{cases}\\

\mathtt{sum\_y}(u_{t})
&=
\begin{cases}
b u_{t,y-1}                      &,\,y=0\\
b u_{t,y+1}                      &,\,y=n_y - 1\\
u_{t,y-1} + u_{t,y+1} &,\,\text{otherwise}
\end{cases}\\
\end{aligned}
$$

$n_x$ と $n_y$ はそれぞれ $x$ 軸と $y$ 軸のインデックスの数です。

境界条件は $b$ の値で変更できます。 $b=1$ のとき固定端、 $b=2$ のとき自由端になります。

$4 \times 4$ のときの $\mathbf{A}$ を組み立てます。

$$
\mathbf{A}
=
\begin{bmatrix}
C_0 & b C_1 & 0 & 0 & b C_1 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0\\
C_1 & C_0 & C_1 & 0 & 0 & b C_1 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0\\
0 & C_1 & C_0 & C_1 & 0 & 0 & b C_1 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0\\
0 & 0 & b C_1 & C_0 & 0 & 0 & 0 & b C_1 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0\\

C_1 & 0 & 0 & 0 & C_0 & b C_1 & 0 & 0 & C_1 & 0 & 0 & 0 & 0 & 0 & 0 & 0\\
0 & C_1 & 0 & 0 & C_1 & C_0 & C_1 & 0 & 0 & C_1 & 0 & 0 & 0 & 0 & 0 & 0\\
0 & 0 & C_1 & 0 & 0 & C_1 & C_0 & C_1 & 0 & 0 & C_1 & 0 & 0 & 0 & 0 & 0\\
0 & 0 & 0 & C_1 & 0 & 0 & b C_1 & C_0 & 0 & 0 & 0 & C_1 & 0 & 0 & 0 & 0\\

0 & 0 & 0 & 0 & C_1 & 0 & 0 & 0 & C_0 & b C_1 & 0 & 0 & C_1 & 0 & 0 & 0\\
0 & 0 & 0 & 0 & 0 & C_1 & 0 & 0 & C_1 & C_0 & C_1 & 0 & 0 & C_1 & 0 & 0\\
0 & 0 & 0 & 0 & 0 & 0 & C_1 & 0 & 0 & C_1 & C_0 & C_1 & 0 & 0 & C_1 & 0\\
0 & 0 & 0 & 0 & 0 & 0 & 0 & C_1 & 0 & 0 & b C_1 & C_0 & 0 & 0 & 0 & C_1\\

0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & b C_1 & 0 & 0 & 0 & C_0 & b C_1 & 0 & 0\\
0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & b C_1 & 0 & 0 & C_1 & C_0 & C_1 & 0\\
0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & b C_1 & 0 & 0 & C_1 & C_0 & C_1\\
0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & b C_1 & 0 & 0 & b C_1 & C_0\\
\end{bmatrix}
$$

何となく規則性が見えます。

$$
\mathbf{B} =
\begin{bmatrix}
  C_0 & b C_1 & 0 & 0\\C_1 & C_0 & C_1 & 0\\0 & C_1 & C_0 & C_1\\0 & 0 & b C_1 & C_0\\
\end{bmatrix}
,\quad
\mathbf{D} =
\begin{bmatrix}
  C_1 & 0 & 0 & 0\\0 & C_1 & 0 & 0\\0 & 0 & C_1 & 0\\0 & 0 & 0 & C_1\\
\end{bmatrix}
,\quad
\mathbf{0} =
\begin{bmatrix}
  0 & 0 & 0 & 0\\0 & 0 & 0 & 0\\0 & 0 & 0 & 0\\0 & 0 & 0 & 0\\
\end{bmatrix}
$$

$$
\mathbf{A}
=
\begin{bmatrix}
\mathbf{B} & b\mathbf{D} & \mathbf{0} & \mathbf{0}\\
\mathbf{D} & \mathbf{B} & \mathbf{D} & \mathbf{0}\\
\mathbf{0} & \mathbf{D} & \mathbf{B} & \mathbf{D}\\
\mathbf{0} & \mathbf{0} & b\mathbf{D} & \mathbf{B}\\
\end{bmatrix}
$$

格子の大きさが $n_x \times n_y$ のとき、 $\mathbf{B},\,\mathbf{D},\,\mathbf{0}$ は $n_x \times n_x$ の正方行列、 $\mathbf{A}$ は $n_y \times n_y$ の正方行列になります。

## デモと実装
ソルバに[ガウス-ザイデル法](https://en.wikipedia.org/wiki/Gauss%E2%80%93Seidel_method)を使っています。計算が重たいので反復回数の上限を8にしています。

キャンバスをクリックすると波が起きます。$u_{t}$ の値が大きいと明るく、小さいと暗く表示されます。赤は $u_{t} > 1.0$ 、緑は $u_{t} < -1.0$ を表しています。

<script src="demo/vec2.js"></script>
<script src="demo/canvas.js"></script>
<script src="demo/wave2d.js"></script>

実装はややこしくなっています。JavaScriptでは参照が直接書けないので `array = [{value: value}, ...]` というように `Object` に値を格納することで配列の操作を楽にしています。

[実装を読む (github.io)](https://github.com/ryukau/filter_notes/blob/master/waveequation2d/demo/wave2d.js)

## その他
### ガウス-ザイデル法による歪みとインデックスの振り方
ガウス-ザイデル法の反復回数が少ないときに波が歪みます。この歪みは1次元配列への変換時に蛇行するようにインデックスを振ることで改善できますが、端のノードについての例外処理が増えます。

<figure>
<img src="img/2d_lattice_snake_index.svg" alt="Image of gibbs phenomenon." style="width: 360px;"/>
</figure>

```javascript
for (var x = 0; x < xLength; ++x) {
  for (var y = 0; y < yLength; ++y) {
    var index = y % 2 == 0
      ? x + y * xLength
      : (y + 1) * xLength - x - 1
    u_array[index] = u[x][y]
  }
}
```

### ソルバの発散を防ぐ
2回の反復を1セットにして、1回目はインデックスを昇順、2回目は降順にたどって計算することで、反復回数が少ないときでもソルバによる発散を防げます。

### $u_{t}$ が $5 \times 3$ のときの $\mathbf{A}$
$\mathbf{u}_{t}$ にインデックスを振ります。

<figure>
<img src="img/2d_lattice_mod_index_5_3.svg" alt="Image of gibbs phenomenon." style="width: 432px;"/>
</figure>

$\mathbf{A}$ を組み立てます。

$$
\mathbf{A}
=
\begin{bmatrix}
C_0 & b C_1 & 0 & 0 & 0 & b C_1 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0\\
C_1 & C_0 & C_1 & 0 & 0 & 0 & b C_1 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0\\
0 & C_1 & C_0 & C_1 & 0 & 0 & 0 & b C_1 & 0 & 0 & 0 & 0 & 0 & 0 & 0\\
0 & 0 & C_1 & C_0 & C_1 & 0 & 0 & 0 & b C_1 & 0 & 0 & 0 & 0 & 0 & 0\\
0 & 0 & 0 & b C_1 & C_0 & 0 & 0 & 0 & 0 & b C_1 & 0 & 0 & 0 & 0 & 0\\

C_1 & 0 & 0 & 0 & 0 & C_0 & b C_1 & 0 & 0 & 0 & C_1 & 0 & 0 & 0 & 0\\
0 & C_1 & 0 & 0 & 0 & C_1 & C_0 & C_1 & 0 & 0 & 0 & C_1 & 0 & 0 & 0\\
0 & 0 & C_1 & 0 & 0 & 0 & C_1 & C_0 & C_1 & 0 & 0 & 0 & C_1 & 0 & 0\\
0 & 0 & 0 & C_1 & 0 & 0 & 0 & C_1 & C_0 & C_1 & 0 & 0 & 0 & C_1 & 0\\
0 & 0 & 0 & 0 & C_1 & 0 & 0 & 0 & b C_1 & C_0 & 0 & 0 & 0 & 0 & C_1\\

0 & 0 & 0 & 0 & 0 & b C_1 & 0 & 0 & 0 & 0 & C_0 & b C_1 & 0 & 0 & 0\\
0 & 0 & 0 & 0 & 0 & 0 & b C_1 & 0 & 0 & 0 & C_1 & C_0 & C_1 & 0 & 0\\
0 & 0 & 0 & 0 & 0 & 0 & 0 & b C_1 & 0 & 0 & 0 & C_1 & C_0 & C_1 & 0\\
0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & b C_1 & 0 & 0 & 0 & C_1 & C_0 & C_1\\
0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & b C_1 & 0 & 0 & 0 & b C_1 & C_0\\
\end{bmatrix}
$$

まとめます。

$$
\mathbf{B} =
\begin{bmatrix}
  C_0 & b C_1 & 0 & 0 & 0\\
  C_1 & C_0 & C_1 & 0 & 0\\
  0 & C_1 & C_0 & C_1 & 0\\
  0 & 0 & C_1 & C_0 & C_1\\
  0 & 0 & 0 & b C_1 & C_0\\
\end{bmatrix}
,\quad
\mathbf{D} =
\begin{bmatrix}
  C_1 & 0 & 0 & 0 & 0\\
  0 & C_1 & 0 & 0 & 0\\
  0 & 0 & C_1 & 0 & 0\\
  0 & 0 & 0 & C_1 & 0\\
  0 & 0 & 0 & 0 & C_1\\
\end{bmatrix}
,\quad
\mathbf{0} =
\begin{bmatrix}
  0 & 0 & 0 & 0 & 0\\
  0 & 0 & 0 & 0 & 0\\
  0 & 0 & 0 & 0 & 0\\
  0 & 0 & 0 & 0 & 0\\
  0 & 0 & 0 & 0 & 0\\
\end{bmatrix}
$$

$$
\mathbf{A}
=
\begin{bmatrix}
\mathbf{B} & b\mathbf{D} & \mathbf{0}\\
\mathbf{D} & \mathbf{B} & \mathbf{D}\\
\mathbf{0} & b\mathbf{D} & \mathbf{B}\\
\end{bmatrix}
$$

### $n$ 次元の $\mathbf{A}$
横書き順インデックスを使うときの $n$ 次元の $\mathbf{A}$ を組み立てます。 $4 \times 4 \times 4 \times \dots$ の場合を考えてみます。


$$
\begin{aligned}
  \begin{aligned}
    \mathbf{B} &=
    \begin{bmatrix}
      C_0 & b C_1 & 0 & 0\\
      C_1 & C_0 & C_1 & 0\\
      0 & C_1 & C_0 & C_1\\
      0 & 0 & b C_1 & C_0\\
    \end{bmatrix}
    \\
    \mathbf{B}^{2} &=
    \begin{bmatrix}
      \mathbf{B} & b \mathbf{D} & \mathbf{0} & \mathbf{0}\\
      \mathbf{D} & \mathbf{B} & \mathbf{D} & \mathbf{0}\\
      \mathbf{0} & \mathbf{D} & \mathbf{B} & \mathbf{D}\\
      \mathbf{0} & \mathbf{0} & b \mathbf{D} & \mathbf{B}\\
    \end{bmatrix}
    \\
    \mathbf{B}^{3} &=
    \begin{bmatrix}
      \mathbf{B}^{2} & b \mathbf{D}^{2} & \mathbf{0}^{2} & \mathbf{0}^{2}\\
      \mathbf{D}^{2} & \mathbf{B}^{2} & \mathbf{D}^{2} & \mathbf{0}^{2}\\
      \mathbf{0}^{2} & \mathbf{D}^{2} & \mathbf{B}^{2} & \mathbf{D}^{2}\\
      \mathbf{0}^{2} & \mathbf{0}^{2} & b \mathbf{D}^{2} & \mathbf{B}^{2}\\
    \end{bmatrix}
    \\
    &\hspace{6.5em}\vdots
  \end{aligned}
  \quad
  \begin{aligned}
    \mathbf{D} &=
    \begin{bmatrix}
      C_1 & 0 & 0 & 0\\
      0 & C_1 & 0 & 0\\
      0 & 0 & C_1 & 0\\
      0 & 0 & 0 & C_1\\
    \end{bmatrix}
    \\
    \mathbf{D}^{2} &=
    \begin{bmatrix}
      \mathbf{D} & \mathbf{0} & \mathbf{0} & \mathbf{0}\\
      \mathbf{0} & \mathbf{D} & \mathbf{0} & \mathbf{0}\\
      \mathbf{0} & \mathbf{0} & \mathbf{D} & \mathbf{0}\\
      \mathbf{0} & \mathbf{0} & \mathbf{0} & \mathbf{D}\\
    \end{bmatrix}
    \\
    \mathbf{D}^{3} &=
    \begin{bmatrix}
      \mathbf{D}^{2} & \mathbf{0}^{2} & \mathbf{0}^{2} & \mathbf{0}^{2}\\
      \mathbf{0}^{2} & \mathbf{D}^{2} & \mathbf{0}^{2} & \mathbf{0}^{2}\\
      \mathbf{0}^{2} & \mathbf{0}^{2} & \mathbf{D}^{2} & \mathbf{0}^{2}\\
      \mathbf{0}^{2} & \mathbf{0}^{2} & \mathbf{0}^{2} & \mathbf{D}^{2}\\
    \end{bmatrix}
    \\
    &\hspace{6em}\vdots
  \end{aligned}
  \quad
  \begin{aligned}
    \mathbf{0} &=
    \begin{bmatrix}
      0 & 0 & 0 & 0\\
      0 & 0 & 0 & 0\\
      0 & 0 & 0 & 0\\
      0 & 0 & 0 & 0\\
    \end{bmatrix}
    \\
    \mathbf{0}^{2} &=
    \begin{bmatrix}
      \mathbf{0} & \mathbf{0} & \mathbf{0} & \mathbf{0}\\
      \mathbf{0} & \mathbf{0} & \mathbf{0} & \mathbf{0}\\
      \mathbf{0} & \mathbf{0} & \mathbf{0} & \mathbf{0}\\
      \mathbf{0} & \mathbf{0} & \mathbf{0} & \mathbf{0}\\
    \end{bmatrix}
    \\
    \mathbf{0}^{3} &=
    \begin{bmatrix}
      \mathbf{0}^{2} & \mathbf{0}^{2} & \mathbf{0}^{2} & \mathbf{0}^{2}\\
      \mathbf{0}^{2} & \mathbf{0}^{2} & \mathbf{0}^{2} & \mathbf{0}^{2}\\
      \mathbf{0}^{2} & \mathbf{0}^{2} & \mathbf{0}^{2} & \mathbf{0}^{2}\\
      \mathbf{0}^{2} & \mathbf{0}^{2} & \mathbf{0}^{2} & \mathbf{0}^{2}\\
    \end{bmatrix}
    \\
    &\hspace{5.5em}\vdots
  \end{aligned}
\end{aligned}
$$

$$
\mathbf{A}
=
\begin{bmatrix}
\mathbf{B}^{n} & b\mathbf{D}^{n} & \mathbf{0}^{n} & \mathbf{0}^{n}\\
\mathbf{D}^{n} & \mathbf{B}^{n} & \mathbf{D}^{n} & \mathbf{0}^{n}\\
\mathbf{0}^{n} & \mathbf{D}^{n} & \mathbf{B}^{n} & \mathbf{D}^{n}\\
\mathbf{0}^{n} & \mathbf{0}^{n} & b\mathbf{D}^{n} & \mathbf{B}^{n}\\
\end{bmatrix}
$$
