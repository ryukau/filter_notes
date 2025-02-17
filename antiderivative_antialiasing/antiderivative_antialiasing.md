# 逆微分による歪みのアンチエイリアシング
逆微分による歪みのアンチエイリアシング (ADAA: antiderivative antialiasing) を実装します。

ADAA はいくつかのバリエーションがあります。ここでは以下の 3 つを扱います。

- 連続領域で FIR フィルタを畳み込む手法。 (Parker らによるオリジナル)
- 有限差分に基づく手法。 (Bilbao ら)
- 連続領域での IIR フィルタを近似する手法。 (La Pastina ら)

この文章では主に Bilbao らによる有限差分に基づく手法を扱います。 Parker らによる手法については、連続領域で畳み込む FIR フィルタを三角窓からコサイン窓に変えるとどうなるか調べました。 La Pastina らによる手法は著者らによる実装例の紹介のみです。

基礎となるアイデアについては Parker らの論文がわかりやすいです。[参考文献](#参考文献)の節にリンクを張っています。

以下は「有限差分に基づく手法のレシピ」の内容の実装例へのリンクです。

- [Python 3 による実装 (github.com)](https://github.com/ryukau/filter_notes/blob/master/antiderivative_antialiasing/adaa.py)
- [C++ による実装 (github.com)](https://github.com/ryukau/filter_notes/blob/master/antiderivative_antialiasing/cpp/adaa.cpp)
- [JavaScript による実装 (github.com)](https://github.com/ryukau/UhhyouWebSynthesizers/blob/main/common/dsp/saturator.js)

## 有限差分に基づく手法のレシピ
Bilbao, Esqueda, Parker, Välimäki による "[Antiderivative antialiasing for memoryless nonlinearities](https://drive.google.com/file/d/1SaqbMpxitC8QECkF3OfzHu7cDCnmpzY7/view)" に基づく歪みのアンチエイリアシングを実装します。ここでの歪みはハードクリッピングや $\tanh$ のように入出力を単純にマッピングするタイプの歪みです。

0 次の ADAA 、つまりアンチエイリアシングが無い歪みを以下のように表すことにします。

$$
y^{(0)}_n = {F}^{(0)} ({x}_{n}).
$$

$x$ は入力、 $y$ は出力、 $F$ は歪みの関数です。下付き文字はサンプル数で表された相対的な時間で、 $n$ は現在時刻、 $n-1$ は 1 サンプル前の値、という意味です。上付き文字の $(a)$ は $a$ 次の ADAA であることを表しています。以降では式の見通しをよくするために、 ${F}^{(a)} ({x}_{n}) = {F}^{(a)}_{n}$ と表記します。

$F$ の中身については後述します。

### 1 次の ADAA
1 次の ADAA は以下の式で計算できます。

$$
y^{(1)}_n = \frac{F^{(1)}_{n} - {F}^{(1)}_{n - 1} }{{x}_{n} - {x}_{n - 1}}.
$$

${F}^{(1)}$ は ${F}^{(0)}$ を逆微分した関数です。

上の式は ${x}_{n} - {x}_{n - 1}$ が 0 に近いと、浮動小数点数による計算結果が inf や NaN となる問題があります。この問題のことを ADAA 関連の論文では ill-condition と呼んでいます。 1 次の ADAA では ill-condition となるときは以下の近似式を代わりに使います。

$$
y^{(1)}_n \approx {F}^{(0)} \left( \frac{x_n + x_{n-1}}{2} \right).
$$

以下は Python 3 による 1 次の ADAA の実装例です。

```python
import numpy as np

def adaa1(x, f0, f1):
    """
    `x` : 入力の 1 次元配列。
    `f0`: 歪みの関数。式中の F_0 。
    `f1`: `f0` を逆微分した関数。式中の F_1 。
    """
    tolerance = 1 / 2**24
    y = np.zeros_like(x)
    x1 = 0
    s1 = 0
    for n in range(len(x)):
        x0 = x[n]
        s0 = f1(x0)
        if (x1 == 0 and s1 == 0) or np.abs(x0 - x1) < tolerance:
            y[n] = f0((x0 + x1) / 2)
        else:
            y[n] = (s0 - s1) / (x0 - x1)
        s1 = s0
        x1 = x0
    return y
```

`np.abs(x0 - x1) < tolerance` が ill-condition を避けるための分岐条件です。また、 $F$ の中身によってはポップノイズが出たため `x1 == 0 and s1 == 0` という条件も追加しています。

計算量には振れ幅があります。条件が良ければ 1 サンプルあたり `f1` を 1 回計算するだけで済みます。 Ill-condition に入るときは 1 サンプル当たりで `f1` を 1 回、 `f0` を 1 回、計算する必要があります。 $F$ の計算が重たいときはノイズのサンプル & ホールド (S&H) などをテスト信号として使い、リアルタイムの締め切りに十分間に合うかを確認したほうがよさそうです。

### 2 次の ADAA
2 次の ADAA は以下の式で計算できます。

$$
y^{(2)}_n =
\frac{2}{{x}_{n} - {x}_{n - 2}} \left(
    \frac{{F}^{(2)}_{n} - {F}^{(2)}_{n - 1} }{{x}_{n} - {x}_{n - 1}}
  - \frac{{F}^{(2)}_{n - 1} - {F}^{(2)}_{n - 2}}{{x}_{n - 1} - {x}_{n - 2}}
\right)
$$

${F}^{(2)}$ は ${F}^{(1)}$ を逆微分した関数です。

Ill-condition のときは以下の近似式を使います。

$$
\begin{aligned}
y^{(2)}_n &\approx \frac{2}{\Delta_{n}} \left( {F}^{(1)}(\bar{x}_{n}) + \frac{{F}^{(2)}_{n-1} - {F}^{(2)}(\bar{x}_{n})}{\Delta_{n}} \right), \\
\bar{x}_{n} &= \frac{x_{n} + x_{n-2}}{2}, \quad \Delta_{n} = \bar{x}_{n} - x_{n-1}.
\end{aligned}
$$

この ill-condition の式は $|\Delta_{n}| < \epsilon$ のときに、さらに ill-condition となります。このときはさらに分岐を設けて以下の近似式を使います。

$$
y^{(2)}_{n} \approx {F}^{(0)} \left( \frac{\bar{x}_{n} + x_{n-1}}{2} \right).
$$

以下は Python 3 による 2 次の ADAA の実装例です。

```python
def adaa2(x, f0, f1, f2):
    """
    `x` : 入力の 1 次元配列。
    `f0`: 歪みの関数。式中の F_0 。
    `f1`: `f0` を逆微分した関数。式中の F_1 。
    `f2`: `f1` を逆微分した関数。式中の F_2 。
    """
    tolerance = 1 / 2**24
    y = np.zeros_like(x)
    x1 = 0
    x2 = 0
    s1 = 0
    for n in range(len(x)):
        x0 = x[n]

        f2_x1 = f2(x1)
        s0 = (
            f1((x0 + x1) / 2)
            if np.abs(x0 - x1) < tolerance
            else (f2(x0) - f2_x1) / (x0 - x1)
        )

        if x1 == 0 and x2 == 0:
            y[n] = f0((x0 + 2 * x1 + x2) / 4)
        elif np.abs(x0 - x2) < tolerance:
            x_bar = (x0 + x2) / 2
            delta = x_bar - x1
            if np.abs(delta) < tolerance:
                y[n] = f0((x_bar + x1) / 2)
            else:
                y[n] = (2 / delta) * (f1(x_bar) + (f2_x1 - f2(x_bar)) / delta) # 🤔
        else:
            y[n] = 2 * (s0 - s1) / (x0 - x2)
        s1 = s0
        x2 = x1
        x1 = x0
    return y
```

$F$ の中身によってはポップノイズが出たので、 `x1 == 0 and x2 == 0` の分岐を加えています。この分岐での計算は `y[n] = f0(x0 / 4)` と簡略化できます。

計算量の振れ幅はややこしい組み合わせがあります。上の実装では、最悪の場合 🤔 に分岐して `f2` が 3 回、 `f1` が 1 回、計算されます。

### 歪みの関数 $F$ について
歪みの関数 $F$ は状態を持たない非線形性 (memoryless nonlinearities) に限られます。また、逆微分可能でなければ ADAA を適用できません。大まかにはプログラミング言語の数学ライブラリに入っているような関数であれば、ほぼ使えます。

区間関数 (piecewise function) も使えますが、逆微分した関数が区間のつなぎ目でつながるように積分定数を決めてやる必要があります。区間関数は積分定数によって式が煩雑になりがちです。

以降では $J^n$ で $n$ 回の逆微分を表します。また、以下の数学特殊関数が現れます。

- $\operatorname{Li}_2$: dilogarithm
- $\operatorname{Li}_3$: trilogarithm
- $\Gamma(a, x)$: upper incomplete gamma function
- $\operatorname{Si}$: sine integral
- $\operatorname{Ci}$: cosine integral
- ${}_3 F_3$: hypergeometric function

### hardclip
分岐の境界となる $x = \pm 1$ のときに、分岐の両側の式が同じ値になる必要があります。そのため積分定数は $J^1$ で $1/2$ 、 $J^2$ で $1/6$ となります。

$$
\begin{aligned}
J^0 f_{\mathrm{hardclip}}(x) &= \begin{cases}
  x & |x| < 1 \\
  \mathrm{sgn}(x) & \text{otherwise}
\end{cases} \\
J^1 f_{\mathrm{hardclip}}(x) &= \begin{cases}
  x^2/2 & |x| < 1 \\
  \mathrm{sgn}(x) \cdot x - 1/2 & \text{otherwise}
\end{cases} \\
J^2 f_{\mathrm{hardclip}}(x) &= \begin{cases}
  x^3/6 & |x| < 1 \\
  \mathrm{sgn}(x) \cdot (x^2/2 + 1/6) - x/2 & \text{otherwise}
\end{cases}\\
\end{aligned}
$$

<details>
<summary>詳細</summary>

SymPy で解く。

```python
import sympy

x = sympy.Symbol("x", real=True)

J0_case0 = x
J0_case1 = sympy.sign(x)

J1_case0 = sympy.integrate(J0_case0, x)
J1_case1 = sympy.integrate(J0_case1, x)
J1_case1 += J1_case0.subs(x, 1) - J1_case1.subs(x, 1)  # 積分定数

J2_case0 = sympy.integrate(J1_case0, x)
J2_case1 = sympy.integrate(J1_case1, x)
J2_case1 += J2_case0.subs(x, 1) - J2_case1.subs(x, 1)  # 積分定数

for expr in [J0_case0, J0_case1, J1_case0, J1_case1, J2_case0, J2_case1]:
    print(f"---\n")
    sympy.pprint(expr)
    print()
```

実装例。 `a if c else b` で分岐を 1 行に押し込めるが、 C 言語などの `c ? a : b` へと移植するときに間違えやすいので避けている。

```python
import numpy as np

def hardclipJ0(x):
    return np.clip(x, -1, 1)

def hardclipJ1(x):
    absed = np.abs(x)
    if absed < 1:
        return x * x / 2
    return absed - 1 / 2

def hardclipJ2(x):
    if np.abs(x) < 1:
        return x * x * x / 6
    return (x * x / 2 + 1 / 6) * np.sign(x) - (x / 2)
```

</details>

### halfrect
hardclip と同様ですが、積分定数はすべて 0 です。

$$
\begin{aligned}
J^0 f_{\mathrm{halfrect}}(x) &= \begin{cases}
  0 & x \leq 0 \\
  x & \text{otherwise}
\end{cases} \\
J^1 f_{\mathrm{halfrect}}(x) &= \begin{cases}
  0 & x \leq 0 \\
  x^2/2 & \text{otherwise}
\end{cases} \\
J^2 f_{\mathrm{halfrect}}(x) &= \begin{cases}
  0 & x \leq 0 \\
  x^3/6 & \text{otherwise}
\end{cases}\\
\end{aligned}
$$

### power
出力が実数かつ、負の値が出てほしいので、実装では $J^0 f = \mathrm{sgn}(x) \cdot |x|^β$ としています。また $\beta < 0$ のときの動作は検証していません。

$$
\begin{aligned}
J^0 f_{\mathrm{power}}(x) &= x^β, \quad x \geq 0\\
J^1 f_{\mathrm{power}}(x) &= \begin{cases}
  \dfrac{x^{β + 1}}{β + 1} & \beta \neq -1\\
  \log(x) & \text{otherwise}
\end{cases} \\
J^2 f_{\mathrm{power}}(x) &= \begin{cases}
  \dfrac{x^{β + 2}}{β^2 + 3 β + 2} & \beta \neq -1\\
  x (\log(x) - 1) & \text{otherwise}
\end{cases}\\
\end{aligned}
$$

<details>
<summary>詳細</summary>

Maxima で解きます。

```maxima
/* Maxima */
J0: x^|β|;
J1: ratsimp(integrate(J0, x));
J2: ratsimp(integrate(J1, x));
```

以下は出力です。条件をいくつか聞かれたので答えています。

```
(J0)	x^β
"Is "β" equal to "-1"?"no;
(J1)	x^(β+1)/(β+1)
"Is "β+1" equal to "-1"?"no;
(J2)	x^(β+2)/(β^2+3*β+2)
```

</details>

### softclip2
詳細については以下のリンク先を参照してください。

- [limiter - ソフトクリップ](https://ryukau.github.io/filter_notes/limiter/limiter.html#%E3%82%BD%E3%83%95%E3%83%88%E3%82%AF%E3%83%AA%E3%83%83%E3%83%97)

積分定数によって煩雑な形になっています。

$$
\begin{aligned}
J^0 f_{\mathrm{softclip2}}(x) &= \begin{cases}
  x & \text{if}\ |x| < a_1
    && \text{(linear region)}\\
  h + \mathrm{sgn}(x) \dfrac{0.25 (a_2 - |x|)^2}{a_1 - h}  & \text{if}\ a_1 \leq|x| < a_2
    && \text{(2nd order region)}\\
  h & \text{if}\ a_2 \leq |x|
    && \text{(clipping region)}\\
\end{cases}
\\
a_1 &= rh, \qquad a_2 = 2h - a_1.
\\\\
J^1 f_{\mathrm{softclip2}}(x) &= \begin{cases}
\displaystyle
  \frac{x^2}{2}
  & \text{if}\ |x| < a_1
\\[1em]
\displaystyle
  a_1 \left( \frac{a_1}{2} - h \right)
  + h x
  + \frac{(a_1 - a_2)^3 - (x - a_2)^3}{12 (h - a_1)}
  & \text{if}\ a_1 \leq|x| < a_2
\\[1em]
\displaystyle
  a_1 \left(\frac{a_1}{2} - h\right)
  + h x
  + \frac{\left(a_1 - a_2\right)^3}{12 (h - a_1)}
  & \text{if}\ a_2 \leq |x|
  \end{cases}
\\
J^2 f_{\mathrm{softclip2}}(x) &= \begin{cases}
\displaystyle
  \frac{x^3}{6}
  & \text{if}\ |x| < a_1
\\[1em]
\displaystyle
  \frac{a_1^2 (3 x - 2 a_1)}{6}
  + \frac{h (x - a_1)^2}{2}
  - \frac{
    (x - a_1)^2 \left(
      - 4 a_2 (a_1 - a_2 + x)
      + 2 (a_1 - a_2)^2
      + (a_1 + x)^2
    \right)
  }{48 (h - a_1)}
  & \text{if}\ a_1 \leq|x| < a_2
\\[1em]
\displaystyle
  \frac{a_1^2 (3 x - 2 a_1)}{6}
  + \frac{h (x - a_1)^2}{2}
  + \frac{(a_1 - a_2)^3 \left(- 3 a_1 - a_2 + 4 x\right)}{48 (h - a_1)}
& \text{if}\ a_2 \leq |x|
\end{cases}
\\
\end{aligned}
$$

<details>
<summary>詳細</summary>

SymPy で解きます。

```python
import sympy

x, h, a_1, a_2 = sympy.symbols("x, h, a_1, a_2", real=True)

# J:       |x| < a1,
# K: a1 <= |x| < a2,
# L: a2 <= |x|.
J0 = x
K0 = h + (a_2 - x) ** 2 / (4 * (a_1 - h))
L0 = h

J1 = sympy.integrate(J0, x)
K1 = sympy.integrate(K0, x)
L1 = sympy.integrate(L0, x)

K1 += J1.subs(x, a_1) - K1.subs(x, a_1)
L1 += K1.subs(x, a_2) - L1.subs(x, a_2)

J2 = sympy.integrate(J1, x)
K2 = sympy.integrate(K1, x)
L2 = sympy.integrate(L1, x)

K2 += J2.subs(x, a_1) - K2.subs(x, a_1)
K2 += J2.subs(x, a_2) - K2.subs(x, a_2)

for cases in [[J0, K0, L0], [J1, K1, L1], [J2, K2, L2]]:
    print("---")
    for expr in cases:
        expr = sympy.simplify(expr)
        sympy.pprint(expr)
        print()
```

</details>

### softclipN
以下は softclipN の仕様です。

<figure>
<img src="img/clipping.svg" alt="Input-output plot of arbitrary order polynomial softclip." style="padding-bottom: 12px;"/>
</figure>

式の簡略化のために $\beta \geq 0$ とします。負の値を使うときは、 $J^1$ のときに $\beta = -1$ 、 $J^2$ のときに $\beta=-2$ で分岐が現れますが、ここでは $\beta > 0$ として無視しています。因数分解 (factorization) は SymPy の出力をもとに手で行ったので最適ではありません。

$$
\begin{aligned}
J^0 f_{\mathrm{softclipN}}(x) &= \begin{cases}
x & x \leq r_c\\
C + A (x_c - x)^\beta & r_c < x < x_c \\
C + A (x_c - x_s)^\beta + S (x - x_s) & x_c \leq x \\
\end{cases}
\\
x_c &= r_c + \beta (C - r_c),\quad
A = \frac{r_c - C}{(x_c - r_c)^\beta},\quad
x_s = x_c - \left( -\frac{S}{A \beta} \right)^{1/(\beta - 1)}.
\\\\
J^1 f_{\mathrm{softclipN}}(x) &= \begin{cases}
\dfrac{x^2}{2}
& x \leq r_c \\[1em]
\displaystyle
\frac{A \left(\left(x_{c} - r_{c}\right)^{1 + β} - \left(x_{c} - x\right)^{1 + β}\right)}{1 + β}
+ C \left(x - r_{c}\right)
+ \frac{r_{c}^{2}}{2}
& r_c < x < x_c \\[1em]
\displaystyle
\frac{A \left(x_{c} - r_{c}\right)^{1 + β}}{1 + β}
+ C \left(x_{c} - r_{c}\right)
+ \frac{S \left(x^{2} - x_{c}^{2}\right)}{2}
+ \frac{r_{c}^{2}}{2}
+ \left(x - x_{c}\right) \left(A \left(x_{c} - x_{s}\right)^{β} + C - S x_{s}\right)
& x_c \leq x \\[1em]
\end{cases}
\\
J^2 f_{\mathrm{softclipN}}(x) &= \begin{cases}
\dfrac{x^3}{6}
& x \leq r_c \\[1em]
\displaystyle
A \left(
  \frac{\left(x_{c} - x\right)^{2 + β} - \left(x_{c} - r_{c}\right)^{2 + β}}{\left(1 + β\right) \left(2 + β\right)} + \frac{\left(x_{c} - r_{c}\right)^{1 + β} \left(x - r_{c}\right)}{1 + β}
\right)
+ \frac{C \left(x - r_{c}\right)^{2}}{2}
+ r_{c}^{2} \left(\frac{x}{2} - \frac{r_{c}}{3}\right)
& r_c < x < x_c \\[1em]
\displaystyle
\frac{A \left(x_{c} - r_{c}\right)^{2 + β} \left(1 - \frac{1}{2 + β}\right)}{1 + β}
+ \frac{C \left(x_{c} - r_{c}\right)^{2}}{2}
+ \frac{S \left(x^{3} - x_{c}^{3}\right)}{6}
+ r_{c}^{2} \left(\frac{x_{c}}{2} - \frac{r_{c}}{3}\right)
+ \left(x - x_{c}\right) \left(
  A \left(\frac{\left(x_{c} - r_{c}\right)^{1 + β}}{1 + β} - x_{c} \left(x_{c} - x_{s}\right)^{β}\right)
  - C r_{c}
  + S x_{c} \left(x_{s} - \frac{x_{c}}{2}\right)
  + \frac{r_{c}^{2}}{2}
  + \frac{x + x_{c}}{2} \left(A \left(x_{c} - x_{s}\right)^{β} + C - S x_{s}\right)
\right)
& x_c \leq x \\[1em]
\end{cases}
\end{aligned}
$$

<details>
<summary>詳細</summary>

SymPy で解きます。

```python
import sympy

x, A, C, r_c, x_c, x_s, S = sympy.symbols("x, A, C,r_c, x_c, x_s, S", real=True)

# `positive=True` で β = -1, β = -2 のときの分岐を除去している。
β = sympy.Symbol("β", real=True, positive=True)

# J:        |x| <= r_c,
# K: r_c <  |x| <  x_c,
# L: x_c <= |x|.
J0 = x
K0 = C + A * (x_c - x) ** β
L0 = C + A * (x_c - x_s) ** β + S * (x - x_s)

J1 = sympy.integrate(J0, x)
K1 = sympy.integrate(K0, x)
L1 = sympy.integrate(L0, x)

K1 += J1.subs(x, r_c) - K1.subs(x, r_c)
L1 += K1.subs(x, x_c) - L1.subs(x, x_c)

J2 = sympy.integrate(J1, x)
K2 = sympy.integrate(K1, x)
L2 = sympy.integrate(L1, x)

K2 += J2.subs(x, r_c) - K2.subs(x, r_c)
L2 += K2.subs(x, x_c) - L2.subs(x, x_c)

for cases in [[J0, K0, L0], [J1, K1, L1], [J2, K2, L2]]:
    print("---")
    for expr in cases:
        expr = sympy.simplify(expr)
        sympy.pprint(expr)
        print()
```

以下は C++ による $J^0 f_{\mathrm{softclipN}}$ のリファレンス実装です。

```c++
float processJ0(float x0)
{
  // 数式との対応: C = clipY, R = ratio, A = scale, β = order, S: slope.
  float absed = std::fabs(x0);

  float rc = C * R;
  if (absed <= rc) return x0;

  float xc = rc + beta * (C - rc);
  float A = (rc - C) / std::pow(xc - rc, beta);
  float xs = xc - std::pow(-S / (A * beta), float(1) / (beta - float(1)));
  return absed < xs
    ? std::copysign(C + A * std::pow(xc - absed, beta), x0)
    : std::copysign(S * (absed - xs) + C + A * std::pow(xc - xs, beta), x0);
}
```

</details>

### tanh

$$
\begin{aligned}
J^0 f_{\mathrm{tanh}}(x) &= \tanh{\left(x \right)}\\
J^1 f_{\mathrm{tanh}}(x) &= \log{\left( \cosh{(x)}\right) }\\
J^2 f_{\mathrm{tanh}}(x) &= x \log \left( \frac{\cosh(x)}{e^{2 x} + 1} \right) + \frac{x^2 - \mathrm{Li}_2(-e^{2 x})}{2}\\
\\
\text{where} \quad \mathrm{Li}_2(z) &= \int_1^z \frac{\log(t)}{1 - t} dt.
\end{aligned}
$$

$\mathrm{Li}_2$ は [dilogarithm](https://mathworld.wolfram.com/Dilogarithm.html) と呼ばれる関数です。

ライブラリによっては $\mathrm{Li}_2$ の代わりに [Spence's function](https://mathworld.wolfram.com/SpencesFunction.html) のみが提供されていることがあります。 Spence's function を $\operatorname{Sp}$ とすると、 $\mathrm{Sp}(x) = -\mathrm{Li}_2 (-x)$ となります。例えば Maxima で [`-li[2](-z)`](https://flex.phys.tohoku.ac.jp/texi/maxima/maxima_14.html#IDX463) となるとき、 SciPy では [`spence(z)`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.spence.html) となります。 Spence's function は [cephes](https://netlib.org/cephes/) で実装されているので、 $\mathrm{Li}_2$ よりもよく移植されているのかと思います。

<details>
<summary>詳細</summary>

Maxima で解きます。

```maxima
/* Maxima */
J0: tanh(x);
J1: ratsimp(integrate(J0, x));
J2: ratsimp(integrate(J1, x));
```

以下は出力です。

```
(J0)	tanh(x)
(J1)	log(cosh(x))
(J2)	(2*x*log(cosh(x))-2*x*log(%e^(2*x)+1)-li[2](-%e^(2*x))+x^2)/2
```

</details>

### atan

$$
\begin{aligned}
J^0 f_{\mathrm{atan}}(x) &= \arctan(x)\\
J^1 f_{\mathrm{atan}}(x) &= x \arctan(x) - \frac{\log(x^{2} + 1)}{2}\\
J^2 f_{\mathrm{atan}}(x) &= \frac{x - x \log(x^{2}+ 1) - (1 - x^{2})\arctan(x)}{2}\\
\end{aligned}
$$

<details>
<summary>詳細</summary>

Maxima で解きます。

```maxima
/* Maxima */
J0: atan(x);
J1: ratsimp(integrate(J0, x));
J2: ratsimp(integrate(J1, x));
```

以下は出力です。

```
(J0)	atan(x)
(J1)	-((log(x^2+1)-2*x*atan(x))/2)
(J2)	-((x*log(x^2+1)+(1-x^2)*atan(x)-x)/2)
```

</details>

### algebraic
$$
\begin{aligned}
J^0 f_{\mathrm{algebraic}}(x) &= \frac{x}{z + 1}\\
J^1 f_{\mathrm{algebraic}}(x) &= z - w\\
J^2 f_{\mathrm{algebraic}}(x) &= \mathrm{sgn}(x) \cdot \left(
  z \left(\frac{z}{2} - w + 1 \right) - w
\right)\\
\\
\text{where} \quad z &= |x|, \quad w = \log(z + 1).
\end{aligned}
$$

$w$ の計算には `log1p` が使えます。

<details>
<summary>詳細</summary>

SymPy で解きます。上の式は複素数を無視しているので SymPy の解とは一致していません。

```python
import sympy
x = sympy.Symbol("x", real=True)
J0 = x / (1 + sympy.Abs(x))
J1 = sympy.integrate(J0, x)
J2 = sympy.integrate(J1, x)
```

</details>

### softplus

$$
\begin{aligned}
J^0 f_{\mathrm{softplus}}(x) &= \log(e^x + 1) \\
J^1 f_{\mathrm{softplus}}(x) &= - \mathrm{Li}_2 (-e^x) \\
J^2 f_{\mathrm{softplus}}(x) &= - \mathrm{Li}_3 (-e^x) \\
\end{aligned}
$$

<details>
<summary>詳細</summary>

Wolfram Alpha で解きます。

```mathematica
J0 = log[exp[x] + 1];
J1 = antiderivative[log[exp[x] + 1], x];
J2 = antiderivative[-polylog[2, -exp[x]], x];
```

$\mathrm{Li}_3$ はあまり使われない関数なのか SciPy では実装されていません。 C++ か Fortran であれば [Expander/polylogarithm](https://github.com/Expander/polylogarithm) が使えます。 Python 3 では [`mpmath.polylog`](https://www.mpmath.org/doc/current/functions/zeta.html#polylog) が使えますが、任意精度なので計算に時間がかかります。

```python
import mpmath
import numpy as np
import scipy.special as special

def softplusJ0(x): return np.log(np.exp(x) + 1)
def softplusJ1(x): return -special.spence(np.exp(x))
def softplusJ2(x): return -float(mpmath.polylog(3, -np.exp(x))) # 任意精度なので遅い。
```

</details>

### swish
$$
\begin{aligned}
J^0 f_{\mathrm{swish}}(x) &= \frac{x}{e^{-x β} + 1}\\
J^1 f_{\mathrm{swish}}(x) &= \frac{x \beta \log{\left( e^{x \beta} + 1 \right) } + \mathrm{Li}_2 \left( -e^{x \beta} \right) }{\beta^2}\\
J^2 f_{\mathrm{swish}}(x) &= \frac{2 \mathrm{Li}_3 \left( -e^{x \beta}\right) -x \beta  \mathrm{Li}_2 \left( -e^{x \beta} \right) }{\beta^{3}}\\
\end{aligned}
$$

Softplus と同様に $\mathrm{Li}$ が現れます。

<details>
<summary>詳細</summary>

Maxima で解きます。

```maxima
/* Maxima */
J0: x / (1 + exp(-β * x));
J1: ratsimp(integrate(J0, x));
J2: ratsimp(integrate(J1, x));
```

以下は出力です。

```
(J0)	x/(%e^(-(x*β))+1)
(J1)	(x*β*log(%e^(x*β)+1)+li[2](-%e^(x*β)))/β^2
(J2)	(2*li[3](-%e^(x*β))-x*β*li[2](-%e^(x*β)))/β^3
```

以下は検証のために Wolfram Alpha で解いたときのクエリです。一行ずつコピーして貼り付けました。 Maxima と同様の解が得られます。

```mathematica
J0 = x / (1 + exp(-β * x));
J1 = antiderivative[x / (1 + exp(-β * x)), x];
J2 = antiderivative[(x*β*log(exp(x*β)+1)+polylog[2, -exp(x*β)])/β^2, x];
```

</details>

### exppoly

$$
\begin{aligned}
J^0 f_{\mathrm{exppoly}}(x) &= x^{β} e^{- z}\\
J^1 f_{\mathrm{exppoly}}(x) &= -\Gamma(1+\beta, z) + C_0\\
J^2 f_{\mathrm{exppoly}}(x) &= \Gamma(2+\beta, z) - z \Gamma(1+\beta, z) + C_0 z + C_1\\
\\
\text{where} \quad
z &= |x|, \quad
C_0 = \Gamma(1 + \beta, 0), \quad
C_1 = \Gamma(2 + \beta, 0).
\end{aligned}
$$

$\Gamma(s, x)$ は [upper incomplete gamma function](https://en.wikipedia.org/wiki/Incomplete_gamma_function) です。 $C_0,\,C_1$ は積分定数です。

<details>
<summary>詳細</summary>

SymPy で解きます。

```python
import sympy
x = sympy.Symbol("x", real=True)
β = sympy.Symbol("β", real=True)
J0 = x**β * sympy.exp(-sympy.Abs(x))
J1 = sympy.integrate(J0, x)
J2 = sympy.integrate(J1, x)
```

以下は出力です。

```
J0: integrate[x^β * exp(-x), x]
J1: integrate[-Gamma[1 + β, x], x]
J2: -(x Gamma[1 + β, x]) + Gamma[2 + β, x]
```

</details>

### sinalgexp
$$
\begin{aligned}
J^0 f_{\mathrm{sinalgexp}}(x) &= \sin(\pi (1 - e^{-x}))\\
J^1 f_{\mathrm{sinalgexp}}(x) &= -\operatorname{Si}(\pi e^{-x}) + C_0\\
J^2 f_{\mathrm{sinalgexp}}(x) &= \frac{1}{2} π e^{-x} ({}_3 F_3(1, 1, 1;2, 2, 2;-i e^{-x} π) + {}_3 F_3(1, 1, 1;2, 2, 2;i e^{-x} π)) + C
\end{aligned}
$$

${}_p F_q$ は [generalized hypergeometric function](https://reference.wolfram.com/language/ref/HypergeometricPFQ.html) です。 ${}_p F_q$ は Boost に実装があります。

- [Hypergeometric pFq - 1.87.0](https://www.boost.org/doc/libs/1_87_0/libs/math/doc/html/math_toolkit/hypergeometric/hypergeometric_pfq.html)

<details>
<summary>詳細</summary>

Wolfram Alpha で解きます。

```mathematica
J0 = sin((e^(x) - 1) / e^(x));
J1 = integrate[sin((e^(x) - 1) / e^(x)), x];
J2 = integrate[-SinIntegral[Pi/E^x], x];
```

</details>

### sinatanexp
$$
\begin{aligned}
J^0 f_{\mathrm{sinatanexp}}(x) &= \sin(2 \arctan(e^x - 1))\\
J^1 f_{\mathrm{sinatanexp}}(x) &= \frac{\log((e^x - 1)^2 + 1)}{2} + \arctan(e^x - 1) - x\\
\end{aligned}
$$

$J^2$ は複素解のみしか得られなかったです。テキスト形式の解を以下の詳細欄に掲載しています。

<details>
<summary>詳細</summary>

Maxima で解きます。

```maxima
J0: sin(2 * atan(%e^x - 1));
J1: integrate(J0, x);
J2: integrate(J1, x); /* 解けない */
```

`J2` が解けなかったので Wolfam Alpha で解きます。

```mathematica
J0 = sin(2 * atan(e^x - 1));
J1 = Integrate[sin(2 * atan(e^x - 1)), x];
J2 = Integrate[(-2 x - (1 - I) ArcTan[1 - E^x] + (1 + I) ArcTan[E^x/(2 - E^x)] + Log[2 - 2 E^x + E^(2 x)])/2, x];
```

Wolfam Alpha からの出力です。 `I` は複素単位です。 `J1` については Maxima の解のほうが簡潔な形をしています。

```
J1 = (
  -2 * x
  - (1 - I) * ArcTan(1 - E^x)
  + (1 + I) * ArcTan(E^x/(2 - E^x))
  + Log((E^x - 1)^2 + 1)
)/2;

J2 = (
  -x^2
  + Log(E^x) * Log((E^x - 1)^2 + 1)
  - (1 - I) * x * ArcTan(1 - E^x)
  + (1 + I) * x * ArcTan(E^x/(2 - E^x))
  - (Log(E^x) - I) * x * Log(1 - (1/2 - I/2) * E^x)
  - (Log(E^x) + I) * x * Log(1 - (1/2 + I/2) * E^x)
  - (1 - I) * PolyLog(2, (1/2 - I/2) * E^x)
  - (1 + I) * PolyLog(2, (1/2 + I/2) * E^x)
)/2
```

</details>

### cosdecay

$$
\begin{aligned}
J^0 f_{\mathrm{cosdecay}}(x) &= \frac{1 - \cos(x)}{x}\\
J^1 f_{\mathrm{cosdecay}}(x) &= \log(x) - \operatorname{Ci}(x)\\
J^2 f_{\mathrm{cosdecay}}(x) &= \sin(x) + x (\log(x) - \operatorname{Ci}(x) - 1)\\
\end{aligned}
$$

<details>
<summary>詳細</summary>

SymPy で解きます。

```python
import sympy
x = sympy.Symbol("x", real=True)
J0 = (1 - sympy.cos(x)) / x
J1 = sympy.integrate(J0, x)
J2 = sympy.integrate(J1, x)
```

</details>

### log1p

$$
\begin{aligned}
J^0 f_{\mathrm{log1p}}(x) &= \mathrm{sgn}(z) \log (z + 1) \\
J^1 f_{\mathrm{log1p}}(x) &= (z + 1) \log(z + 1) - z \\
J^2 f_{\mathrm{log1p}}(x) &= \frac{\mathrm{sgn}(z)}{4} \left( 2 (z + 1)^2 \log(z + 1) - 3 z^2 - 2 z \right) \\\\
\text{where} \quad z &= |x|.
\end{aligned}
$$

<details>
<summary>詳細</summary>

Maxima で解きます。

```maxima
J0: log(1+x);
J1: ratsimp(integrate(J0, x));
J2: ratsimp(integrate(J1, x));
```

出力です。

```
(J0)	log(x+1)
(J1)	(x+1)*log(x+1)-x-1
(J2)	((2*x^2+4*x+2)*log(x+1)-3*x^2-6*x)/4
```

`J1` の積分定数を 1 とすると `(x+1)*log(x+1)-x` となります。この式を逆微分すると上で掲載している式が出ます。積分定数を変えたのは、 Maxima の解の形だと `J1` が 0 の周りで不連続になっていたからです。

</details>

### sinrunge
[Runge 関数](https://en.wikipedia.org/wiki/Runge%27s_phenomenon)に sin を乗算した式です。

$$
J^0 f_{\mathrm{sinrunge}}(x) = \frac{\sin(x)}{25 x^2 + 1}.
$$

逆微分については複素解しか得られなかったので、以下の詳細欄にテキスト形式でのみ掲載しています。

<details>
<summary>詳細</summary>

Wolfram Alpha を使います。 Maxima 5.47.0 と SymPy 1.13.3 では解けなかったです。

```mathematica
J1 = integrate[sin(x)/(25*x^2 + 1), x];
J2 = integrate[(CosIntegral[I/5 - x] Sinh[1/5] + CosIntegral[I/5 + x] Sinh[1/5] + I Cosh[1/5] (SinIntegral[I/5 - x] + SinIntegral[I/5 + x]))/10, x];
```

複素解です。実部を取り出せば ADAA に利用できます。

```
J1 = (
  +     sinh(1/5) * (Ci(I/5 - x) + Ci(I/5 + x))
  + I * cosh(1/5) * (Si(I/5 - x) + Si(I/5 + x))
) / 10;

J2 = (
  + (-1 + exp(2/5)) * ((-I + 5   * x) * Ci(I/5 - x) + (+I + 5   * x) * Ci(I/5 + x))
  + (+1 + exp(2/5)) * ((+1 + 5*I * x) * Si(I/5 - x) + (-1 + 5*I * x) * Si(I/5 + x))
)/(100 * exp(1/5));
```

</details>

### 解けなかった逆微分
以下は exppoly と似たような形の関数です。 Upper incomplete gamma function の移植を避けたかったので調べました。

$$
\begin{aligned}
& \int \sin(\pi \tanh(e^x - 1)) dx\\
& \int \sin(\pi \operatorname{erf}(x)) dx\\
& \int \sin \left(\pi \left( \frac{2}{(e^x + 1)} - 1 \right) \right) dx\\
& \int \arctan(\sin(\pi e^{-x})) dx\\
& \int \frac{\sin(\pi e^{-x})}{\sin(\pi e^{-x}) + 1} dx\\
\end{aligned}
$$

### 3 次以上の ADAA
3 次以上の ADAA についても調べたのですが、 ill-condition の扱いが手間だったので途中で諦めました。任意の次数の ADAA の式だけであれば以下にリンクしたコードで計算できます。

- [filter_notes/antiderivative_antialiasing/solve_bilbao.py at master · ryukau/filter_notes · GitHub](https://github.com/ryukau/filter_notes/blob/master/antiderivative_antialiasing/solve_bilbao.py)

以下の詳細の欄の内容は未検証です。

<details>
<summary>詳細</summary>

Bilbao らの手法は Parker らの手法と 1 次のときは同じ。離散系の計算式を有限差分の式から求めている点が Parker らの手法とは異なる。

$$
\begin{aligned}
y^{(1)}_n &= \frac{{F}^{(1)}_{n} - {F}^{(1)}_{n - 1} }{{x}_{n} - {x}_{n - 1}}
\\
y^{(2)}_n &=
\frac{2}{{x}_{n} - {x}_{n - 2}} \left(
    \frac{{F}^{(2)}_{n} - {F}^{(2)}_{n - 1} }{{x}_{n} - {x}_{n - 1}}
  - \frac{{F}^{(2)}_{n - 1} - {F}^{(2)}_{n - 2}}{{x}_{n - 1} - {x}_{n - 2}}
\right)
\\
y^{(3)}_n &=
\frac{2}{{x}_{n - 1} - {x}_{n - 2}} \left(
  \frac{1}{{x}_{n} - {x}_{n - 2}} \left(
    \frac{{F}^{(3)}_{n} - {F}^{(3)}_{n - 1} }{{x}_{n} - {x}_{n - 1}}
  - \frac{{F}^{(3)}_{n - 1} - {F}^{(3)}_{n - 2}}{{x}_{n - 1} - {x}_{n - 2}}
  \right)
- \frac{1}{{x}_{n - 1} - {x}_{n - 3}} \left(
    \frac{{F}_{n - 1} - {F}_{n - 2}}{{x}_{n - 1} - {x}_{n - 2}}
  - \frac{{F}_{n - 2} - {F}_{n - 3}}{{x}_{n - 2} - {x}_{n - 3}}
  \right)
\right)
\\
y^{(4)}_n &=
\frac{4}{{x}_{n - 1} - {x}_{n - 3}} \left(
  \frac{1}{{x}_{n - 1} - {x}_{n - 2}} \left(
    \frac{1}{{x}_{n} - {x}_{n - 2}} \left(
      \frac{{F}^{(4)}_{n} - {F}^{(4)}_{n - 1} }{{x}_{n} - {x}_{n - 1}}
    - \frac{{F}^{(4)}_{n - 1} - {F}^{(4)}_{n - 2}}{{x}_{n - 1} - {x}_{n - 2}}
    \right)
  - \frac{1}{{x}_{n - 1} - {x}_{n - 3}} \left(
      \frac{{F}^{(4)}_{n - 1} - {F}^{(4)}_{n - 2}}{{x}_{n - 1} - {x}_{n - 2}}
    - \frac{{F}^{(4)}_{n - 2} - {F}^{(4)}_{n - 3}}{{x}_{n - 2} - {x}_{n - 3}}
    \right)
  \right)
- \frac{1}{{x}_{n - 2} - {x}_{n - 3}} \left(
    \frac{1}{{x}_{n - 1} - {x}_{n - 3}} \left(
      \frac{{F}^{(4)}_{n - 1} - {F}^{(4)}_{n - 2}}{{x}_{n - 1} - {x}_{n - 2}}
    - \frac{{F}^{(4)}_{n - 2} - {F}^{(4)}_{n - 3}}{{x}_{n - 2} - {x}_{n - 3}}
    \right)
  - \frac{1}{{x}_{n - 2} - {x}_{n - 4}} \left(
      \frac{{F}^{(4)}_{n - 2} - {F}^{(4)}_{n - 3}}{{x}_{n - 2} - {x}_{n - 3}}
    - \frac{{F}^{(4)}_{n - 3} - {F}^{(4)}_{n - 4}}{{x}_{n - 3} - {x}_{n - 4}}
    \right)
  \right)
\right)
\end{aligned}
$$

記号を定義して整理。

$$
q^{(i)}_a = \frac{{F}^{(i)}_{n - a} - {F}^{(i)}_{n - a - 1}}{{x}_{n - a} - {x}_{n - a - 1}}
,\quad
r_{a,b} = {x}_{n - a} - {x}_{n - b}.
$$

$$
\begin{aligned}
y^{(1)}_n &= q^{(1)}_0
\\
y^{(2)}_n &= \frac{2}{r_{0,2}} \left( q^{(2)}_0 - q^{(2)}_1 \right)
\\
y^{(3)}_n &=
\frac{2}{r_{1,2}} \left(
  \frac{1}{r_{0,2}} \left( q^{(3)}_0 - q^{(3)}_1 \right)
- \frac{1}{r_{1,3}} \left( q^{(3)}_1 - q^{(3)}_2 \right)
\right)
\\
y^{(4)}_n &=
\frac{4}{r_{1,3}} \left(
  \frac{1}{r_{1,2}} \left(
    \frac{1}{r_{0,2}} \left( q^{(4)}_0 - q^{(4)}_1 \right)
  - \frac{1}{r_{1,3}} \left( q^{(4)}_1 - q^{(4)}_2 \right)
  \right)
- \frac{1}{r_{2,3}} \left(
    \frac{1}{r_{1,3}} \left( q^{(4)}_1 - q^{(4)}_2 \right)
  - \frac{1}{r_{2,4}} \left( q^{(4)}_2 - q^{(4)}_3 \right)
  \right)
\right)
\end{aligned}
$$

0 除算のフォールバックを無視するなら、以下のように効率よく計算できる。 1 サンプルの計算では、次数によらず `s1` についてのみ $F$ の計算を行えばいい。

```
s1 = (F4(x[n]) - F4(x[n - 1])) / r(0, 1);
s2 = (s1 - y1) / r(0, 2);
s3 = (s2 - y2) / r(1, 2);
s4 = (s3 - y3) / r(1, 3);
output = 4 * s4;
y1 = s1
y2 = s2
y3 = s3
```

`s` の分母の `r` は以下の規則で増える。

- `s` のインデックスが奇数から偶数になるとき、後の引数に +1 。
- `s` のインデックスが偶数から奇数になるとき、前の引数に +1 。

つまり `[r(0, 1), r(0, 2), r(1, 2), r(1, 3), r(2, 3), r(2, 4), r(3, 4), r(3, 5), ...]` 。

0 除算を考慮すると効率のいい計算方法は使えない。 `F` の計算回数は次数に応じて 2 のべき乗のオーダーで増えていくので、実用上は 2 次までと考えていい。

フォールバックの計算を手抜きするなら以下の式が使える。[パスカルの三角形](https://en.wikipedia.org/wiki/Pascal's_triangle)に基づく。

$$
\tilde{y}^{(d)}_n = F^{(0)} \left(
  \frac{1}{2^d}
  \sum_{i=1}^d \left[ (-1)^{i+1} \binom{k}{i} x_{n-i} \right]
\right).
$$

</details>

## 連続領域で FIR フィルタを畳み込む手法の改変
[Parker らの手法](https://dafx16.vutbr.cz/dafxpapers/20-DAFx-16_paper_41-PN.pdf)について、連続領域で畳み込む FIR フィルタを三角窓からコサイン窓に変えるとどうなるか試しました。以下はコサイン窓の定義です。

$$
h_{\mathrm{cos}}(t) = \begin{cases}
\dfrac{1 - \cos(\pi t)}{2}, & 0 \leq t < 2 \\
0, & \text{otherwise}
\end{cases}
$$

Parker らの式 12 と同様に展開します。

$$
\begin{aligned}
\tilde{y}(n)
&=
\int_{-\infty}^{\infty} h_{\mathrm{cos}}(u) y(n - u) du
\\&=
\int_0^2 \frac{1 - \cos(\pi u)}{2} y(n - u) du
\\&=
  \int_0^1 \frac{1 - \cos(\pi \tau)}{2} f(x_n + \tau(x_{n-1} - x_n)) d \tau
+ \int_0^1 \frac{1 - \cos(\pi (\tau + 1))}{2} f(x_{n-1} + \tau(x_{n-2} - x_{n-1})) d \tau
\\&=
  \int_0^1 \frac{1 - \cos(\pi \tau)}{2} f(x_n + \tau(x_{n-1} - x_n)) d \tau
+ \int_0^1 \frac{1 + \cos(\pi \tau)}{2} f(x_{n-1} + \tau(x_{n-2} - x_{n-1})) d \tau
\end{aligned}
$$

ここで integration by substitution が使えないので、 Parker らの式 13 は使えません。 $f$ を直接代入して解きます。

### Ill-condition
Ill-condition に対応するため Parker らの式 39 と式 41 を使います。 $M_1$ はフォールバックの際に $\tau$ に代入する値です。

$$
\begin{aligned}
\int_0^1 y(\tau)\bar{h}(\tau) d\tau &= (x_n + M_1 (x_{n-1} - x_n)) \int_0^1 \bar{h}(\tau) d\tau \\
M_1 &= \frac{\int_0^1 \tau \bar{h}(\tau) d\tau}{\int_0^1 \bar{h}(\tau) d\tau}.
\end{aligned}
$$

$\bar{h}$ に $h_{\mathrm{cos}}$ を代入して解きます。 $f$ は半波整流などの非線形関数です。畳み込みの前半と後半で 2 つの式が現れます。

$$
\begin{aligned}
\int_0^1 h_{\mathrm{cos}}(\tau) d\tau
&= \int_0^1 \frac{1 \mp \cos(\pi \tau)}{2} d\tau
= \frac{1}{2},
\\
M_1
&=\frac{\int_0^1 \tau h_{\mathrm{cos}}(\tau) d\tau}{\int_0^1 h_{\mathrm{cos}}(\tau) d\tau}
= \frac{
  \displaystyle \int_0^1 \tau \frac{1 \mp \cos(\pi \tau)}{2} d\tau
}{
  1/2
}
= \frac{1}{2} \pm \frac{2}{\pi^{2}},
\\
\int_0^1 y(\tau)\bar{h}(\tau) d\tau
&= \begin{cases}
  \displaystyle
  \frac{1}{2} f \left( x_n + \left( \frac{1}{2} + \frac{2}{\pi^{2}} \right) (x_{n-1} - x_n) \right)
  \\
  \displaystyle
  \frac{1}{2} f \left( x_{n-1} + \left( \frac{1}{2} - \frac{2}{\pi^{2}} \right) (x_{n-2} - x_{n-1}) \right)
\end{cases}
\end{aligned}
$$

```python
import sympy
tau = sympy.symbols("τ", real=True)

def M1(h):
    numer = sympy.integrate(tau * h, (tau, 0, 1))
    denom = sympy.integrate(h, (tau, 0, 1))
    print(sympy.latex(numer / denom))
    print(sympy.latex(denom))

M1((1 - sympy.cos(sympy.pi * tau)) / 2)
M1((1 + sympy.cos(sympy.pi * tau)) / 2)
```

### tanh
以下の式を解く必要があるのですが、 Maxima, SymPy, Wolfram Alpha では解けなかったです。

$$
\int_0^1 \frac{1 - \cos(\pi \tau)}{2} \tanh(x_n + \tau(x_{n-1} - x_n)) d \tau
$$

```maxima
expr: (1 - cos(%pi * tau)) / 2 * tanh(x0 + tau * (x1 - x0));
integrate(expr, tau, 0, 1);
```

```python
import sympy
tau, x0, x1 = sympy.symbols("τ, x_0 x_1")
expr = (1 - sympy.cos(sympy.pi * tau)) / 2 * sympy.tanh(x0 + tau * (x1 - x0))
result = sympy.integrate(expr, (tau, 0, 1))
```

### 半波整流
畳み込みの前半だけを扱います。後半は $x_{n+1} \to x_{n+2}$ 、 $x_n \to x_{n+1}$ と置き換わり、 $h_{\mathrm{cos}}$ の展開が $\dfrac{1 - \cos(\pi \tau)}{2}$ から $\dfrac{1 + \cos(\pi \tau)}{2}$ へと変わります。

以下は半波整流の定義です。

$$
f_{\mathrm{HalfRect}}(x) = \begin{cases}
  0, & x < 0 \\
  x, & 0 \leq x
\end{cases}
$$

入力信号 $x$ の値に応じて 4 つの分岐があります。上から順にケース 0 から 3 と番号を振っておきます。

$$
\begin{aligned}
\int_0^1 h_{\mathrm{cos}}(\tau) f_{\mathrm{HalfRect}}(x_n + \tau(x_{n-1} - x_n)) d \tau
&= \begin{cases}
0
, & x_n < 0,    \ x_{n-1} < 0  && \text{(Case 0)}  \\
\displaystyle \int_p^1 h_{\mathrm{cos}}(\tau) (x_n + \tau(x_{n-1} - x_n)) d \tau
, & x_n < 0,    \ x_{n-1} \geq 0 && \text{(Case 1)} \\
\displaystyle \int_0^p h_{\mathrm{cos}}(\tau) (x_n + \tau(x_{n-1} - x_n)) d \tau
, & x_n \geq 0, \ x_{n-1} < 0 && \text{(Case 2)} \\
\displaystyle \int_0^1 h_{\mathrm{cos}}(\tau) (x_n + \tau(x_{n-1} - x_n)) d \tau
, & x_n \geq 0, \ x_{n-1} \geq 0 && \text{(Case 3)} \\
\end{cases}
\\
\text{where}
\quad p &= \frac{-x_n}{x_{n-1} - x_n}.
\end{aligned}
$$

以下は解を計算する SymPy のコードです。手探りで行ったため Term 1 と Term 2 のコードが大幅に重複していますが、整理しても再利用する用途がないのでそのまま掲載しています。

```python
def solveCosine_HalfRect():
    tau, x0, x1, x2 = sympy.symbols("τ, x_0, x_1, x_2", real=True)

    print("--- Term 1")
    expr = (1 - sympy.cos(sympy.pi * tau)) / 2 * (x0 + tau * (x1 - x0))
    p = -x0 / (x1 - x0)
    cases = [
        0,
        sympy.integrate(expr, [tau, p, 1]),
        sympy.integrate(expr, [tau, 0, p]),
        sympy.integrate(expr, [tau, 0, 1]),
    ]
    for index, expr in enumerate(cases):
        expr = sympy.simplify(expr)
        print(f"-- Case {index}")
        print(expr, end="\n\n")

    print("--- Term 2")
    expr = (1 + sympy.cos(sympy.pi * tau)) / 2 * (x1 + tau * (x2 - x1))
    p = -x1 / (x2 - x1)
    cases = [
        0,
        sympy.integrate(expr, [tau, p, 1]),
        sympy.integrate(expr, [tau, 0, p]),
        sympy.integrate(expr, [tau, 0, 1]),
    ]
    for index, expr in enumerate(cases):
        expr = sympy.simplify(expr)
        print(f"-- Case {index}")
        print(expr, end="\n\n")
```

以下は ADAA の実装へのリンクです。動作確認のために SymPy からの出力をそのままコピーしただけなので計算効率が悪いです。手での式の整理はコストパフォーマンスが悪いと判断して中断しました。

- https://github.com/ryukau/filter_notes/blob/e4608cfa031213dcf794204dcc5ada389946bbc3/antiderivative_antialiasing/adaa_parkermod.py#L106-L192

### ハードクリップ
以下はハードクリップの定義です。$|x| \leq 1$ は、 $-1 \leq x \leq 1$ と等価です。以降の表記を簡略化するために絶対値による表記を使っています。

$$
f_{\mathrm{HardClip}}(x) = \begin{cases}
  -1, &  x < -1 \\
   x, & |x| \leq 1 \\
   1, &  1 < x \\
\end{cases}
$$

9 つの分岐があります。上から順にケース 0 から 8 まで番号を振っておきます。 $(-1),\,(+1),\,(L)$ の表記は後で確認しやすいようにしているだけで、関数の引数ではありません。

$$
\begin{aligned}
\int_0^1 h_{\mathrm{cos}}(\tau) f_{\mathrm{HardClip}}(x_n + \tau(x_{n-1} - x_n)) d \tau
&= \begin{cases}
-1/2
, & x_n < -1,     \ x_{n-1} < -1     && \text{(Case 0)} \\
\displaystyle \int_0^a h_{\mathrm{cos}}(\tau) (-1) d \tau
+ \displaystyle \int_a^1 h_{\mathrm{cos}}(\tau) (L) d \tau
, & x_n < -1,     \ |x_{n-1}| \leq 1 && \text{(Case 1)} \\
\displaystyle \int_0^a h_{\mathrm{cos}}(\tau) (-1) d \tau
+ \displaystyle \int_a^b h_{\mathrm{cos}}(\tau) (L) d \tau
+ \displaystyle \int_b^1 h_{\mathrm{cos}}(\tau) (+1) d \tau
, & x_n < -1,     \ x_{n-1} > 1      && \text{(Case 2)} \\
\displaystyle \int_0^a h_{\mathrm{cos}}(\tau) (L) d \tau
+ \displaystyle \int_a^1 h_{\mathrm{cos}}(\tau) (-1) d \tau
, & |x_n| \leq 1, \ x_{n-1} < -1     && \text{(Case 3)} \\
\displaystyle \int_0^1 h_{\mathrm{cos}}(\tau) (L) d \tau
, & |x_n| \leq 1, \ |x_{n-1}| \leq 1 && \text{(Case 4)} \\
\displaystyle \int_0^b h_{\mathrm{cos}}(\tau) (L) d \tau
+ \displaystyle \int_b^1 h_{\mathrm{cos}}(\tau) (+1) d \tau
, & |x_n| \leq 1, \ x_{n-1} > 1      && \text{(Case 5)} \\
\displaystyle \int_0^b h_{\mathrm{cos}}(\tau) (+1) d \tau
+ \displaystyle \int_b^a h_{\mathrm{cos}}(\tau) (L) d \tau
+ \displaystyle \int_a^1 h_{\mathrm{cos}}(\tau) (-1) d \tau
, & x_n > 1,      \ x_{n-1} < -1     && \text{(Case 6)} \\
\displaystyle \int_0^b h_{\mathrm{cos}}(\tau) (+1) d \tau
+ \displaystyle \int_b^1 h_{\mathrm{cos}}(\tau) (L) d \tau
, & x_n > 1,      \ |x_{n-1}| \leq 1 && \text{(Case 7)} \\
1/2
, & x_n > 1,      \ x_{n-1} > 1      && \text{(Case 8)} \\
\end{cases}
\\
\text{where} \quad
L &= x_n + \tau(x_{n-1} - x_n), \quad
a = \frac{-x_n - 1}{x_{n-1} - x_n}, \quad
b = \frac{-x_n + 1}{x_{n-1} - x_n}.
\end{aligned}
$$

以下は解を計算する SymPy のコードです。手探りで行ったため Term 1 と Term 2 のコードが大幅に重複していますが、整理しても再利用する用途がないのでそのまま掲載しています。

```python
import sympy

def solveCosine_Hardclip():
    tau, x0, x1, x2 = sympy.symbols("τ, x_0, x_1, x_2", real=True)

    print("--- Term 1")
    h_cos = (1 - sympy.cos(sympy.pi * tau)) / 2
    exprN = -h_cos
    exprP = h_cos
    exprL = h_cos * (x0 + tau * (x1 - x0))
    a = (-x0 - 1) / (x1 - x0)
    b = (-x0 + 1) / (x1 - x0)
    J = sympy.integrate
    cases = [
        -sympy.Rational(1, 2),
        J(exprN, [tau, 0, a]) + J(exprL, [tau, a, 1]),
        J(exprN, [tau, 0, a]) + J(exprL, [tau, a, b]) + J(exprP, [tau, b, 1]),
        J(exprL, [tau, 0, a]) + J(exprN, [tau, a, 1]),
        J(exprL, [tau, 0, 1]),
        J(exprL, [tau, 0, b]) + J(exprP, [tau, b, 1]),
        J(exprP, [tau, 0, b]) + J(exprL, [tau, b, a]) + J(exprN, [tau, a, 1]),
        J(exprP, [tau, 0, b]) + J(exprL, [tau, b, 1]),
        sympy.Rational(1, 2),
    ]
    for index, expr in enumerate(cases):
        print(f"-- Case {index}")
        expr = sympy.expand(expr)
        expr = sympy.cancel(expr)
        expr = sympy.trigsimp(expr)
        print(expr, end="\n\n")

    print("--- Term 2")
    h_cos = (1 + sympy.cos(sympy.pi * tau)) / 2
    exprN = -h_cos
    exprP = h_cos
    exprL = h_cos * (x1 + tau * (x2 - x1))
    a = (-x1 - 1) / (x2 - x1)
    b = (-x1 + 1) / (x2 - x1)
    J = sympy.integrate
    cases = [
        -sympy.Rational(1, 2),
        J(exprN, [tau, 0, a]) + J(exprL, [tau, a, 1]),
        J(exprN, [tau, 0, a]) + J(exprL, [tau, a, b]) + J(exprP, [tau, b, 1]),
        J(exprL, [tau, 0, a]) + J(exprN, [tau, a, 1]),
        J(exprL, [tau, 0, 1]),
        J(exprL, [tau, 0, b]) + J(exprP, [tau, b, 1]),
        J(exprP, [tau, 0, b]) + J(exprL, [tau, b, a]) + J(exprN, [tau, a, 1]),
        J(exprP, [tau, 0, b]) + J(exprL, [tau, b, 1]),
        sympy.Rational(1, 2),
    ]
    for index, expr in enumerate(cases):
        print(f"-- Case {index}")
        expr = sympy.expand(expr)
        expr = sympy.cancel(expr)
        expr = sympy.trigsimp(expr)
        print(expr, end="\n\n")
```

以下は ADAA の実装へのリンクです。動作確認のために SymPy からの出力をそのままコピーしただけなので計算効率が悪いです。手での式の整理はコストパフォーマンスが悪いと判断して中断しました。

- https://github.com/ryukau/filter_notes/blob/e4608cfa031213dcf794204dcc5ada389946bbc3/antiderivative_antialiasing/adaa_parkermod.py#L195-L403

## 連続領域での IIR フィルタを近似する手法
La Pastina らによってハードクリップのときの実装例が Matlab のコードとして提供されています。以下のリンクから入手できます。

- [Companion page for “Arbitrary-Order IIR Antiderivative Antialiasing”](http://www.dangelo.audio/dafx2021-aaiir.html)

以下は Python 3 へ移植したコードへのリンクです。他の手法の論文通りの実装と混ざっています。 La Pastina らに関連する実装は関数名に `aaiir` を含んでいます。

- [filter_notes/antiderivative_antialiasing/adaa_othermethods.py at master · ryukau/filter_notes · GitHub](https://github.com/ryukau/filter_notes/blob/master/antiderivative_antialiasing/adaa_othermethods.py)

## その他
### C++ で使える数学特殊関数のライブラリ
#### Boost
定番ですが、大きいライブラリなので気軽には使えない印象があります。

- [Chapter 8. Special Functions - 1.87.0](https://www.boost.org/doc/libs/1_87_0/libs/math/doc/html/special.html)

#### Cephes
Cephes は C 言語で書かれた数学関数のライブラリです。実装が簡潔なので比較的手軽に移植できます。ガンマ関数の実装 `gamma` と `lgam` にグローバル変数 `sgngam` が使われているので注意してください。

- [Cephes - netlib.org](https://netlib.org/cephes/)
- [Cephes - www.moshier.net](http://www.moshier.net/#Cephes)

Cephes の Windows 上での ビルドはやや独特で、 cl.exe の `/c` オプションを指定してデバッグビルドしないと動きません。 `/c` オプションを指定しない、あるいはリリースビルドを行うと、 cl.exe が提供する数学ライブラリとリンクされるので `sin` や `exp` といった関数名が競合してしまいます。

CMake でビルドするときは以下のリポジトリの `cephes` ディレクトリが参考になります。

- [GitHub - google-deepmind/torch-cephes: Cephes Mathematical Functions library wrapped for Torch](https://github.com/google-deepmind/torch-cephes)

#### Polylogarithm
Polylogarithm ($\operatorname{Li}_n$) については以下の C++ のライブラリが使えます。 [Clausen function](https://en.wikipedia.org/wiki/Clausen_function) ($\operatorname{Cl}_n$) と Glaisher–Clausen function ($\operatorname{Sl}_n$) の実装もあります。

- [GitHub - Expander/polylogarithm: Implementation of polylogarithms in C/C++/Fortran](https://github.com/Expander/polylogarithm)

### JavaScript での数学特殊関数
ブラウザ上の JavaScript で数学特殊関数を使いたいときは [emscripten](https://emscripten.org/) を使うと楽です。 `emcc` のオプション `-sUSE_BOOST_HEADERS=1` を指定すれば Boost と簡単にリンクできます。

以下は C++ で書かれたバインディングのコード例です。

```c++
// test.cpp
#include <boost/math/special_functions/gamma.hpp>
#include <boost/version.hpp>
#include <emscripten/bind.h>
#include <string>

using namespace emscripten;

// メンテナンス時にバージョンが参照できると便利。
std::string boost_version() { return BOOST_LIB_VERSION; }

EMSCRIPTEN_BINDINGS(BoostMath)
{
  function("boost_version", &boost_version);

  // ガンマ関数のバインディング例。 JavaScript にあわせて double を指定。
  function("gamma_p", &boost::math::gamma_p<double, double>);
  function("gamma_q", &boost::math::gamma_q<double, double>);
  function("tgamma_lower", &boost::math::tgamma_lower<double, double>);
  function("tgamma", &boost::math::tgamma<double, double>);
}
```

上のコードを `text.cpp` として保存したあとに、以下のコマンドで JavaScript ファイルへとコンパイルします。

```ps1
# PowerShell
emcc                    `
  -lembind              `
  -O3                   `
  -sWASM=0              `
  -sUSE_BOOST_HEADERS=1 `
  -sMODULARIZE=1        `
  -sEXPORT_ES6=1        `
  -sPOLYFILL=1          `
  -o test.js            `
  ../test.cpp
```

`emcc` の使い方については以下のページを参照してください。

- [Emscripten Compiler Settings — Emscripten 4.0.0-git (dev) documentation](https://emscripten.org/docs/tools_reference/settings_reference.html#use-boost-headers)
- [Emscripten Compiler Frontend (emcc) — Emscripten 4.0.0-git (dev) documentation](https://emscripten.org/docs/tools_reference/emcc.html)

上のコマンドを実行すると `test.js` が作成されます。以下はバインディングの使用例です。

```html
<!doctype html>
<html>
<script type="module">
  import BoostMath from "./build/somemath.js";

  async function main() {
    const bmath = await BoostMath();

    console.log(`boost_version: ${bmath.boost_version()}`);
    console.log(`gamma_p result: ${bmath.gamma_p(2, 0.5)}`);
    console.log(`gamma_q result: ${bmath.gamma_q(2, 0.5)}`);
    console.log(`tgamma_lower result: ${bmath.tgamma_lower(2, 0.5)}`);
    console.log(`tgamma result: ${bmath.tgamma(2, 0.5)}`);
  }
  main();
</script>

</html>
```

あとは上の HTML を `test.html` などと保存して、 `test.js` と同じディレクトリに配置したあと、 `python3 -m http.server` などを使えばテストできます。 `test.html` をローカルで動作させる方法は以下のリンク先を参考にしてください。

- [How do you set up a local testing server? - Learn web development | MDN](https://developer.mozilla.org/en-US/docs/Learn_web_development/Howto/Tools_and_setup/set_up_a_local_testing_server)

## 参考文献
- Parker, J. D., Zavalishin, V., & Le Bivic, E. (2016, September). [Reducing the aliasing of nonlinear waveshaping using continuous-time convolution](https://dafx16.vutbr.cz/dafxpapers/20-DAFx-16_paper_41-PN.pdf). In Proc. Int. Conf. Digital Audio Effects (DAFx-16), Brno, Czech Republic (pp. 137-144).
- Bilbao, S., Esqueda, F., Parker, J. D., & Välimäki, V. (2017). [Antiderivative antialiasing for memoryless nonlinearities](https://drive.google.com/file/d/1SaqbMpxitC8QECkF3OfzHu7cDCnmpzY7/view). IEEE Signal Processing Letters, 24(7), 1049-1053.
- La Pastina, P. P., D'Angelo, S., & Gabrielli, L. (2021, September). [Arbitrary-order IIR antiderivative antialiasing](https://www.researchgate.net/profile/Stefano-Dangelo/publication/354574545_Arbitrary-Order_IIR_Antiderivative_Antialiasing/links/61408f9c578238365b0981c7/Arbitrary-Order-IIR-Antiderivative-Antialiasing.pdf). In 2021 24th International Conference on Digital Audio Effects (DAFx) (pp. 9-16). IEEE.
- [GitHub - jatinchowdhury18/ADAA: Experiments with Antiderivative Antialiasing](https://github.com/jatinchowdhury18/ADAA?tab=readme-ov-file)
