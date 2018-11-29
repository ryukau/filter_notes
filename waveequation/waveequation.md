# 波のシミュレーション
[波動方程式](https://en.wikipedia.org/wiki/Wave_equation)をコンピュータで計算できる形に変形します。

## 波動方程式
波は波動方程式で表すことができます。

$$
\frac{\partial^2 u(\vec{x},\,t)}{\partial t^2} = c^2 \nabla^2 u(\vec{x},\,t)
$$

$u(\vec{x}, t)$ は、[場 ](https://en.wikipedia.org/wiki/Scalar_field)$u$ の地点 $\vec{x}$ 、時点 $t$ における値です。

$c$ は波の伝わる速度です。

## 1次元の波
1次元の波について考えます。

波の次元が変わると[ラプラシアン](https://en.wikipedia.org/wiki/Laplace_operator) $\nabla^2$ の展開形が変わります。 $n$ 次元のラプラシアンは次のように書けます。

$$
\nabla^2u(\vec{x},\,t) = \sum_{i = 0}^{n - 1} \frac{\partial^2 u(\vec{x},\,t)}{\partial x_i^2}
$$

$n=1$ としたときの波動方程式です。

$$
\frac{\partial^2 u(x,\,t)}{\partial t^2} = c^2 \frac{\partial^2 u(x,\,t)}{\partial x^2}
$$

## 有限差分法
[有限差分法](https://en.wikipedia.org/wiki/Finite_difference_method)は微分方程式を計算できる形にする方法の一つです。

[SymPy](http://www.sympy.org/en/index.html)の [`as_finite_difference`](http://docs.sympy.org/latest/tutorial/calculus.html?highlight=as_finite_difference#finite-differences) を使って1次元の波動方程式を[有限差分](https://en.wikipedia.org/wiki/Finite_difference)の形に変えます。

```python
from sympy import *

c, t, u, x, dx, dt = symbols('c t u x dx dt')

u_xt = u(x, t)
u_dx2 = u_xt.diff(x, x).as_finite_difference(dx)
u_dt2 = u_xt.diff(t, t).as_finite_difference(dt)

eq = Eq(u_dt2, c**2 * u_dx2)
ans = solve(eq, u(x, t + dt))

latex(expand(ans[0]))
```

出力された式です。

$$
\begin{aligned}
u(x,t + dt) = &- \frac{2 c^{2}}{dx^{2}} dt^{2} u{\left (x,t \right )} + \frac{c^{2} dt^{2}}{dx^{2}} u{\left (- dx + x,t \right )} + \frac{c^{2} dt^{2}}{dx^{2}} u{\left (dx + x,t \right )}\\
&+ 2 u{\left (x,t \right )} - u{\left (x,- dt + t \right )}
\end{aligned}
$$

式を整理します。

$$
\begin{aligned}
u(x,t + dt) &= \alpha \left(u{\left (- dx + x,t \right )} + u{\left (dx + x,t \right )}\right)
+ \beta u{\left (x,t \right )} - u{\left (x,- dt + t \right )}
\\
\alpha &= c^{2} \frac{dt^{2}}{dx^{2}}\\
\beta &= 2 \left(1 - \alpha \right)
\end{aligned}
$$

計算できる形になりました。実装します。

```javascript
// u(x, t) -> wave[t][x]
var wave = []
for (var i = 0; i < 3; ++i) {
  wave.push(new Array(512).fill(0))
}

var c = 4       // 波の速度 : m/s
var dt = 1 / 60 // 1ステップの時間 : 秒
var dx = 0.1    // wave の要素間の距離 : メートル

var alpha = (c * dt / dx)**2
var beta = 2 * (1 - alpha)

function step() {
  wave.unshift(wave.pop())

  // 固定端。wave[t][0] と wave [t][last] は 0 で固定。
  var last = wave[0].length - 1
  for (var x = 1; x < last; ++x) {
    wave[0][x] = alpha * (wave[1][x + 1] + wave[1][x - 1])
      + beta * wave[1][x]
      - wave[2][x]
  }
}
```

ここでは $u$ を配列 `wave` として実装しています。配列 `wave` のインデックスは $dt$ や $dx$ の係数として考えることができます。例えば $u(x + 3dx, t - 2dt)$ なら `wave[t - 2][x + 3]` と書きかえられます。

上のコードのままでは波に力を加えていないので動きません。適当に動かしてみてください。

```javascript
function move() {
  // 波を適当に動かす。
  wave[0][Math.floor(wave[0].length / 2)] = 0.01 * Math.sin(Date.now() * 1e-3)
}

function draw() {
  // canvas などに描画。ここでは省略。
}

function animate() {
  draw()
  step()
  move()
  requestAnimationFrame(animate)
}
```

## 境界条件
固定端、自由端、輪について見ていきます。

固定端のときは端を定数で固定します。

```javascript
function step() {
  wave.unshift(wave.pop())

  // wave[0][0] と wave[0][last] に定数を入れておく。
  var last = wave[0].length - 1
  for (var x = 1; x < last; ++x) {
    wave[0][x] = alpha * (wave[1][x + 1] + wave[1][x - 1])
      + beta * wave[1][x]
      - wave[2][x]
  }
}
```

自由端のときは、端の値に、端の一つ内側の値を使います。

```javascript
function step() {
  wave.unshift(wave.pop())

  var last = wave[0].length - 1

  // 何度も繰り返すので for の中に if を入れることを避ける。
  wave[0][0] = alpha * (wave[1][1] + wave[1][1])
    + beta * wave[1][0]
    - wave[2][0]

  wave[0][last] = alpha * (wave[1][last - 1] + wave[1][last - 1])
    + beta * wave[1][last]
    - wave[2][last]

  for (var x = 1; x < last; ++x) {
    wave[0][x] = alpha * (wave[1][x + 1] + wave[1][x - 1])
      + beta * wave[1][x]
      - wave[2][x]
  }
}
```

輪のときは端と端をつなぎます。

```javascript
function step() {
  wave.unshift(wave.pop())

  var last = wave[0].length - 1

  wave[0][0] = alpha * (wave[1][1] + wave[1][last])
    + beta * wave[1][0]
    - wave[2][0]

  wave[0][last] = alpha * (wave[1][0] + wave[1][last - 1])
    + beta * wave[1][last]
    - wave[2][last]

  for (var x = 1; x < last; ++x) {
    wave[0][x] = alpha * (wave[1][x + 1] + wave[1][x - 1])
      + beta * wave[1][x]
      - wave[2][x]
  }
}
```

## デモ
キャンバスをクリックすると左端の固定端を動かして波を起こします。

1フレーム毎にシミュレーションを4ステップ進めています。

<script src="demo/vec2.js"></script>
<script src="demo/canvas.js"></script>
<script src="demo/wave1d.js"></script>

## 問題
### 波が減衰しない
時間が経っても波が減衰しません。最初に出てきた波動方程式にはエネルギーの減衰についての項が無いためです。

波動方程式にばねとダンパの項を加えた式が使えます。

$$
\frac{\partial^2 u}{\partial t^2} + a \frac{\partial u}{\partial t} + ku = c^2 \nabla^2 u
$$

 [A causal and fractional all-frequency wave equation for lossy media](http://heim.ifi.uio.no/~sverre/papers/2011_HolmNasholm-fractZener-JournAcoustSocAm.pdf) の式(10)でエネルギーの減衰について考慮された波動方程式が紹介されています。

$$
\nabla^2 u
- \frac{1}{c_0^2} \frac{\partial^2 u}{\partial t^2}
+ \tau_{\sigma}^{\alpha} \frac{\partial^{\alpha}}{\partial t^{\alpha}} \nabla^2 u
- \frac{\tau_{\varepsilon}^{\beta}}{c_0^2} \frac{\partial^{\beta + 2} u}{\partial t^{\beta + 2}}
= 0
$$

厳密なシミュレーションでなければ、各ステップで適当な減衰係数をかけることで雰囲気は出ます。

```javascript
var attenuation = 0.996 // 減衰係数。 0 < attenuation < 1

function step() {
  wave.unshift(wave.pop())

  // 固定端。
  var last = wave[0].length - 1
  for (var x = 1; x < last; ++x) {
    wave[0][x] = attenuation * (
      alpha * (wave[1][x + 1] + wave[1][x - 1])
      + beta * wave[1][x]
      - wave[2][x]
    )
  }
}
```

### ギブス現象
今回のデモでは無視できないほどの[ギブス現象](https://en.wikipedia.org/wiki/Gibbs_phenomenon)が起こっています。ギブス現象の低減については [On the Gibbs phenomenon and its resolution](http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.386.7325&rep=rep1&type=pdf) でいくつか手法が紹介されています。

### 発散
$c$ 、 $dt$ 、 $dx$ の値によってはシミュレーションが発散します。今回のシミュレーションは explicit な形になっているので implicit な形にすることが考えられます。

## 参考サイト
- [Numerical partial differential equations - Wikipedia](https://en.wikipedia.org/wiki/Numerical_partial_differential_equations)
- [Finite difference coefficient - Wikipedia](https://en.wikipedia.org/wiki/Finite_difference_coefficient)
- [Acoustic attenuation - Wikipedia](https://en.wikipedia.org/wiki/Acoustic_attenuation)
- [Physically Based Modeling](http://www.cs.cmu.edu/~baraff/sigcourse/index.html)
