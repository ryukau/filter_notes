# ステレオパンニング
[Loudness Concepts & Panning Laws](http://www.cs.cmu.edu/~music/icm-online/readings/panlaws/) で紹介されているモノラル → ステレオのパンニング法 (panning law) を基にして、ステレオ → ステレオのパンニングの計算式をいくつか作ります。

## モノラル → ステレオのパンニング
以下はモノラル → ステレオパンニングのブロック線図です。

<figure>
<img src="img/mono_panning_block_diagram.svg" alt="Block diagram of mono to stereo panning." style="padding-bottom: 12px;"/>
</figure>

以下の記号を使います。

- $p$: ユーザが指定したパンの値。 $[0, 1]$ の範囲。 0 で左いっぱい、 0.5 でセンター、 1 で右いっぱい。
- $L$: 左チャンネルのゲイン。
- $R$: 右チャンネルのゲイン。

入力信号を $x$ 、出力信号を $y_L, y_R$ とすると、モノラル → ステレオのパンニングは以下の式で計算できます。

$$
\begin{aligned}
y_L &= L \cdot x\\
y_R &= R \cdot x\\
\end{aligned}
$$

$L^2 + R^2$ で音のパワーを計算できます。パワーは、スピーカから人間の耳に至るまでの反響や位相のずれなどを加味した、大まかな音量を表しています。

### 線形パンニング
線形パンニングは $p$ の値をそのまま使って直線的にゲインを変えるパンニングです。利点はパンを振った後でも左右のチャンネルを足し合わせると元のモノラル信号に戻ることです。欠点は中央にパンを振ると左右にパンを振ったときよりも音が小さく聞こえることです。

線形パンニングの計算式です。

$$
\begin{aligned}
L &= 1 - p\\
R &= p
\end{aligned}
$$

コードにします。 `pan == 0.0` で左いっぱい、 `pan == 0.5` でセンター、 `pan == 1.0` で右いっぱいです。戻り値は `(左チャンネルのゲイン, 右チャンネルのゲイン)` となっています。

```python
def panMonoLinear(pan):
    return (1 - pan, pan)
```

以下の図は線形パンニングのゲイン特性のプロットです。左から順に各チャンネルのゲイン、左右を足し合わせたときのゲイン、パワーを表しています。パンが 0.5 のときパワーが下がっているので再生環境によっては音が小さくなることが確認できます。

<figure>
<img src="img/panMonoLinear.svg" alt="Plot of linear panning gain curve." style="padding-bottom: 12px;"/>
</figure>

### 等パワーパンニング
等パワーパンニングは音の大きさがパンによらず均一に聞こえるようにしたパンニングです。欠点は左右のチャンネルを足し合わせたときに最大で 3 dB 振幅が大きくなることです。

線形パンニングの計算式です。

$$
\begin{aligned}
\theta &= \frac{\pi}{2} p\\
L &= \cos(\theta)\\
R &= \sin(\theta)
\end{aligned}
$$

$L^2 + R^2 = \cos^2(\theta) + \sin^2(\theta)$ となるので三角関数の [Pythagorean identity](https://en.wikipedia.org/wiki/List_of_trigonometric_identities#Pythagorean_identities) よりパワーは常に $1$ です。

コードにします。

```python
def panMonoConstantPower(pan):
    theta = (0.5 * np.pi) * pan
    return (np.cos(theta), np.sin(theta))
```

以下の図は等パワーパンニングのゲイン特性のプロットです。パンが 0.5 のとき $L + R$ が上がっていることが確認できます。

<figure>
<img src="img/panMonoConstantPower.svg" alt="Plot of constant power panning gain curve." style="padding-bottom: 12px;"/>
</figure>

### -4.5 dB パンニング
-4.5 dB パンニングは線形パンニングと等パワーパンニングの中間的なパンニングです。左右のスピーカの間隔が狭いテレビなどでは等パワーパンニングよりも自然に聞こえることがあるそうです。

-4.5 dB パンニングの計算式です。

$$
\begin{aligned}
\theta &= \frac{\pi}{2} p\\
L &= \sqrt{(1 - p)\cos(\theta)}\\
R &= \sqrt{p \sin(\theta)}
\end{aligned}
$$

コードにします。

```python
def panMonoIntermediatePower(pan):
    theta = (0.5 * np.pi) * pan
    return (
        np.sqrt((1 - pan) * np.cos(theta)),
        np.sqrt(pan * np.sin(theta)),
    )
```

以下の図は -4.5 dB パンニングのゲイン特性のプロットです。 $L + R$ は等パワーパンニングの半分、パワーは線形パンニングの半分の値になっています。

<figure>
<img src="img/panMonoIntermediatePower.svg" alt="Plot of -4.5 dB panning gain curve." style="padding-bottom: 12px;"/>
</figure>

## ステレオ → ステレオのパンニング
以下はステレオ → ステレオパンニングのブロック線図です。

<figure>
<img src="img/stereo_panning_block_diagram.svg" alt="Block diagram of stereo to stereo panning." style="padding-bottom: 12px;"/>
</figure>

ステレオ → ステレオパンニングでは、モノラル → ステレオパンニングのゲイン $L, R$ に加えて、入力と出力のチャンネル数を掛け合わせた 4 つのゲインを計算します。それぞれのゲインに記号を割り当てます。

- $G_{LL}$: 左 → 左のゲイン
- $G_{RL}$: 右 → 左のゲイン
- $G_{LR}$: 左 → 右のゲイン
- $G_{RR}$: 右 → 右のゲイン

入力信号を $x$ 、出力信号を $y$ とします。チャンネルは下付き文字で $x_L$ や $y_R$ のように表すことにします。このとき出力信号は以下の式で計算できます。中点 $\cdot$ は乗算を表しています。

$$
\begin{aligned}
y_L &= L \cdot (G_{LL} \cdot x_L + G_{RL} \cdot x_R)\\
y_R &= R \cdot (G_{LR} \cdot x_L + G_{RR} \cdot x_R)\\
\end{aligned}
$$

左右で同じ信号が入力されたときにモノラル → ステレオパンニングと同じ振る舞いをさせたいので、 $L, R$ をかける前のゲインの和が 1 になるようにします。

$$
\begin{aligned}
G_{LL} + G_{RL} &= 1\\
G_{LR} + G_{RR} &= 1\\
\end{aligned}
$$

すると左右で同じ信号 $x$ が入力されたとき、出力信号の式を以下のように変形できます。

$$
\begin{aligned}
y_L = L \cdot x\\
y_R = R \cdot x\\
\end{aligned}
$$

モノラル → ステレオパンニングの計算式と同じになっています。ゲイン $L, R$ はモノラル → ステレオパンニングの式がそのまま使えます。

### 直線的なフェード
ステレオ → ステレオパンニングで導入された 4 つのゲイン $G_{LL}, G_{LR}, G_{RR}, G_{RL}$ をどう決めればいいのか見ていきます。

出力の計算式を再掲します。

$$
\begin{aligned}
y_L &= L \cdot (G_{LL} \cdot x_L + G_{RL} \cdot x_R)
,\qquad G_{LL} + G_{RL} = 1\\
y_R &= R \cdot (G_{LR} \cdot x_L + G_{RR} \cdot x_R)
,\qquad G_{LR} + G_{RR} = 1\\
\end{aligned}
$$

和を 1 以下にする条件から $G_{RL}$ と $G_{LR}$ は以下のように計算できます。

$$
\begin{aligned}
G_{RL} &= 1 - G_{LL}\\
G_{LR} &= 1 - G_{RR}\\
\end{aligned}
$$

つまり $G_{LL}$ と $G_{RR}$ のゲインだけを決めればいいわけです。

また、パンの位置 $p$ によって $G_{LL}$ に以下の条件をつけます。

- パンが中央のとき ($p = 0.5$) 、入力をバイパスしたいので $G_{LL}$ は 1 。
- パンが右に振れているとき ($p > 0.5$) 、左出力 $y_L$ に右入力 $x_R$ を混ぜたくないので、 $G_{LL}$ は 1 。
- パンが左いっぱいに振れたとき ($p = 0$) 、左出力 $y_L$ で左右の入力を均等に混ぜたいので、 $G_{LL}$ は 0.5 。

同様に $G_{RR}$ に以下の条件をつけます。

- パンが中央のとき ($p = 0.5$) 、入力をバイパスしたいので $G_{RR}$ は 1 。
- パンが左に振れているとき ($p < 0.5$) 、右出力 $y_R$ に左入力 $x_L$ を混ぜたくないので、 $G_{RR}$ は 1 。
- パンが右いっぱいに振れたとき ($p = 1$) 、右出力 $y_R$ で左右の入力を均等に混ぜたいので、 $G_{RR}$ は 0.5 。

上の条件を直線でつなぐと直線的なフェードになります。以下は直線的なフェードの図です。

<figure>
<img src="img/stereo_panning_gain.svg" alt="Example plot of stereo to stereo panning gain." style="padding-bottom: 12px;"/>
</figure>

式にします。

$$
\begin{aligned}
G_{LL} &= \begin{cases}
  0.5 + p & \text{if} \enspace p < 0.5\\
  1           & \text{if} \enspace p \geq 0.5\\
\end{cases}\\
G_{RR} &= \begin{cases}
  1          & \text{if} \enspace p \leq 0.5\\
  1.5 - p & \text{if} \enspace p > 0.5\\
\end{cases}\\
G_{RL} &= 1 - G_{LL}\\
G_{LR} &= 1 - G_{RR}\\
\end{aligned}
$$

$G_{LL}$ と $G_{RR}$ は 0.5 を境に線対称になっているので $G_{RR}(p) = G_{LL}(1 - p)$ と書くこともできます。以降では同様の対称性があるものとして $G_{LL}$ だけを求めます。

コードにします。

```python
def getStereoGain(gainFunc, param):
    LL = gainFunc(pan, param)
    RR = gainFunc(1 - pan, param)
    RL = 1 - LL
    LR = 1 - RR
    return (LL, RL, LR, RR)

def panStereoLinear(pan, param=None):
    def gain(pan, unused):
        return np.where(pan >= 0.5, 1, 0.5 + pan)

    return getStereoGain(gain, param)
```

`np.where(a, b, c)`  は C++ の `a ? b : c` に訳せます。 `getStereoGain` は後で繰り返し使うので分けています。

以下は `panStereoLinear` の使用例です。

```python
# source は [左チャンネル, 右チャンネル] の 2 次元配列。

LL, RL, LR, RR = panStereoLinear(pan)
gainL, gainR = panMonoConstantPower(pan)

sigL = gainL * (LL * source[0] + RL * source[1])
sigR = gainR * (LR * source[0] + RR * source[1])
```

以下の図は直線的なフェード曲線のプロットです。

<figure>
<img src="img/panStereoLinear.svg" alt="Plot of linear curve for stereo to stereo panning." style="padding-bottom: 12px;"/>
</figure>

### 部分的に 2 次曲線を使うフェード
以下の図のような $G_{LL}$ の式を作ります。黒い部分は直線、オレンジの部分は 2 次曲線です。

<figure>
<img src="img/stereo_panning_partial_2nd_order.svg" alt="Plot of G_LL curve with partial 2nd order region." style="padding-bottom: 12px;"/>
</figure>

上の図の $G_{LL}$ を $p$ について微分すると以下のような形になります。

<figure>
<img src="img/stereo_panning_partial_2nd_order_diff1.svg" alt="Plot of differentiation of G_LL respect to p." style="padding-bottom: 12px;"/>
</figure>

$\dfrac{dG_{LL}}{dp}$ の図の斜線の領域の面積は 0.5 です。これは横軸 $p$ の 0 から 0.5 の区間での $G_{LL}$ の増分です。台形の面積の式から以下の等式が立ちます。

$$
\frac{(a + 0.5)}{2} h = 0.5
$$

$h$ について解きます。

$$
h = \frac{1}{a + 0.5}
$$

これで $\dfrac{dG_{LL}}{dp}$ が計算できます。

$$
\frac{dG_{LL}}{dp} = \begin{cases}
  h                             & \text{if} \enspace p \leq a\\
  \dfrac{h\,(p - 0.5)}{a - 0.5} & \text{if} \enspace a < p < 0.5\\
  0                             & \text{if} \enspace 0.5 \leq p\\
\end{cases}
$$

$p$ について積分すると $dG_{LL}$ の計算式が得られます。

$$
\begin{aligned}
G_{LL}(p) &= \begin{cases}
  hp + 0.5
    & \text{if} \enspace p \leq a
    & \text{(linear region)}\\
  \dfrac{h(p - 0.5)^2}{2a - 1} + 1
    & \text{if} \enspace a < p < 0.5
    & \text{(2nd order region)}\\
  1
    & \text{if} \enspace 0.5 \leq p
    & \text{(constant region)}\\
\end{cases}\\
h &= \frac{1}{a + 0.5}
\end{aligned}
$$

2 次領域の式は以下の Maxima のコードの出力を整形したものです。 `-C + 1` が積分定数です。

```maxima
G: integrate(h*(p-0.5)/(a-0.5), p);
C: subst(0.5, p, G);
factorout(G - C + 1, p);
```

コードにします。 `param` は値が `[0.0, 1.0]` で正規化されたパラメータで、デフォルト値は適当に決めています。

```python
def panStereoPartial2ndOrder(pan, param=0.8):
    """`param` is in [0.0, 1.0]."""
    def gain(p, a):
        h = 1 / (a + 0.5)
        return np.where(
            p <= a,
            h * p + 0.5,
            np.where(
                p >= 0.5,
                1,
                h * (p - 0.5)**2 / (2 * a - 1) + 1,
            ),
        )

    return getStereoGain(gain, 0.5 * param)
```

以下の図は部分的に 2 次曲線を使ったフェード曲線のプロットです。 $a = 0.4$ です。

<figure>
<img src="img/panStereoPartial2ndOrder.svg" alt="Plot of partial 2nd order curve for stereo to stereo panning." style="padding-bottom: 12px;"/>
</figure>

### 部分的に $\sin$ を使うフェード
部分的に 2 次曲線を使うフェードの 2 次曲線を $\sin$ に置き換えたフェードを作ります。

以下のような $\dfrac{dG_{LL}}{dp}$ の曲線を使います。オレンジの線の部分には $\cos$ を使います。

<figure>
<img src="img/stereo_panning_partial_sin_diff1.svg" alt="Plot of differentiation of G_LL respect to p." style="padding-bottom: 12px;"/>
</figure>

$\dfrac{dG_{LL}}{dp}$ の式を立てます。

$$
\begin{aligned}
\frac{dG_{LL}}{dp} &= \begin{cases}
  h
    & \text{if} \enspace p \leq a\\
  \displaystyle \frac{h}{2} \left( \cos \left( \frac{p - a}{0.5 - a} \pi \right) + 1 \right)
    & \text{if} \enspace a < p < 0.5\\
  0
    & \text{if} \enspace 0.5 \leq p\\
\end{cases}\\
h &= \frac{2} {2a + 1}
\end{aligned}
$$

$h$ は以下の Maxima のコードで求めました。

```maxima
S: integrate(h/2 * (cos((p - a) / (1/2 - a) * %pi) + 1), p, a, 1/2);
A: ratsimp(a * h + S);
solve(1/2 = A, h);
```

$\dfrac{dG_{LL}}{dp}$ の式を積分します。

$$
\begin{aligned}
G_{LL}(p) &= \begin{cases}
  hp + 0.5
    & \text{if} \enspace p \leq a
    & \text{(linear region)}\\
  \displaystyle
  0.5 h \left(
    \frac{0.5 - a}{\pi} \sin \left( \frac{p - a}{0.5 - a} \pi \right)
    + p
    - 0.5
  \right) + 1
    & \text{if} \enspace a < p < 0.5
    & \text{(sin region)}\\
  1
    & \text{if} \enspace 0.5 \leq p
    & \text{(constant region)}\\
\end{cases}\\
h &= \frac{2} {2a + 1}
\end{aligned}
$$

$\sin$ 領域の式は以下の Maxima のコードの `expr` を整形したものです。

```maxima
S: integrate(h/2 * (cos((p - a) / (1/2 - a) * %pi) + 1), p);
C: subst(1/2, p, S);
expr: S - C + 1;
```

コードにします。

```python
def panStereoPartialSin(pan, param=0.8):
    """`param` is in [0.0, 1.0]."""
    def gain(p, a):
        h = 2 / (2 * a + 1)
        return np.where(
            p <= a,
            h * p + 0.5,
            np.where(
                p >= 0.5,
                1,
                0.5 * h * ((0.5 - a) / np.pi * np.sin(
                    (p - a) / (0.5 - a) * np.pi) + p - 0.5) + 1,
            ),
        )

    return getStereoGain(gain, 0.5 * param)
```

以下の図は部分的に $\sin$ を使ったフェード曲線のプロットです。 $a = 0.4$ です。

<figure>
<img src="img/panStereoPartialSin.svg" alt="Plot of partial sin curve for stereo to stereo panning." style="padding-bottom: 12px;"/>
</figure>

### 円形フェード
円の式です。

$$
x^2 + y^2 = r^2
$$

- [Circle - Wikipedia](https://en.wikipedia.org/wiki/Circle)

$r$ を 1 として $y$ について解きます。値が正になる解だけを使います。

$$
y = \sqrt{1 - x^2}\\
$$

$x$ は範囲 $(0, 1]$ のパラメータにして、フェードに使う弧の長さを変えられるようにします。 $x = 0$ のときは $G_{LL}$ の式で 0 除算が起こるので範囲に含めません。

<figure>
<img src="img/half_circle.svg" alt="Image of the interval of arc which corresponds to parameter x." style="padding-bottom: 12px;"/>
</figure>

$G_{LL}$ の計算式です。

$$
G_{LL}(p) = \begin{cases}
  0.5 + 0.5 \dfrac{\sqrt{1 - (x - 2xp)^2} - y}{1 - y}
    & \text{if} \enspace p \leq 0.5\\
  1
    & \text{if} \enspace 0.5 \leq p\\
\end{cases}\\
$$

コードにします。

```python
def panStereoCircle(pan, param=1):
    """`param` is in (0.0, 1.0]."""
    def gain(p, x):
        y = np.sqrt(1 - x * x)
        return np.where(
            p <= 0.5,
            0.5 + 0.5 * np.sqrt(1 - (x - 2 * x * p)**2) - y / (1 - y),
            1,
        )

    return getStereoGain(gain, param)
```

以下の図は円の式を使ったフェード曲線のプロットです。 $x = 1$ です。

<figure>
<img src="img/panStereoCircle.svg" alt="Plot of half circle curve for stereo to stereo panning." style="padding-bottom: 12px;"/>
</figure>

### $n$ 次曲線フェード
$G_{LL}$ の式を立てます。 $n$ はフェード曲線を変えるパラメータで、曲線の次数です。

$$
G_{LL}(p) = \begin{cases}
  1 - 0.5(1 - 2p)^n
    & \text{if} \enspace p \leq 0.5\\
  1
    & \text{if} \enspace 0.5 \leq p\\
\end{cases}\\
$$

コードにします。

```python
def panStereoPoly(pan, param=0.25):
    """`param` is in [0.0, 1.0]."""
    def gain(p, n):
        return np.where(
            p <= 0.5,
            1 - 0.5 * np.power(1 - 2 * p, n),
            1,
        )

    order = (1 + 4 * param)
    return getStereoGain(gain, order)
```

以下の図は $n$ 次曲線を使ったフェード曲線のプロットです。 $n = 2$ です。

<figure>
<img src="img/panStereoPoly.svg" alt="Plot of polynomial curve for stereo to stereo panning." style="padding-bottom: 12px;"/>
</figure>

### $\sin$ を使うフェード
$G_{LL}$ の式を立てます。

$$
G_{LL}(p) = \begin{cases}
  0.5 + 0.5\sin(\pi p)
    & \text{if} \enspace p \leq 0.5\\
  1
    & \text{if} \enspace 0.5 \leq p\\
\end{cases}\\
$$

コードにします。

```python
def panStereoSin(pan, param=None):
    def gain(p, unused):
        return np.where(
            p <= 0.5,
            0.5 + 0.5 * np.sin(np.pi * p),
            1,
        )

    return getStereoGain(gain, param)
```

以下の図は $\sin$ を使ったフェード曲線のプロットです。

<figure>
<img src="img/panStereoSin.svg" alt="Plot of sin curve for stereo to stereo panning." style="padding-bottom: 12px;"/>
</figure>

### S 字のフェード
$G_{LL}$ の式を立てます。 $n$ はフェード曲線を変えるパラメータで、曲線の次数です。

$$
G_{LL}(p) = \begin{cases}
  0.75 - 0.25\cos^n(2 \pi p)
    & \text{if} \enspace p \leq 0.5\\
  1
    & \text{if} \enspace 0.5 \leq p\\
\end{cases}\\
$$

コードにします。 $n$ の範囲は適当に $[1, 4]$ としています。

```python
def panStereoSCurve(pan, param=0):
    def gain(p, n):
        return np.where(
            p <= 0.5,
            0.75 - 0.25 * np.cos(2 * np.pi * p)**n,
            1,
        )

    order = 1 + 3 * param
    return getStereoGain(gain, order)
```

以下の図は S 字のフェード曲線のプロットです。 $n = 1$ です。

<figure>
<img src="img/panStereoSCurve.svg" alt="Plot of S-shaped curve for stereo to stereo panning." style="padding-bottom: 12px;"/>
</figure>

### 漏れのあるフェード (Softplus)
$G_{LL}$ の式を立てます。 $k$ はフェード曲線を変えるパラメータです。範囲は $[1, 100]$ くらいで指数スケールを使うと良さそうです。

$$
\begin{aligned}
G_{LL}(p) &= 0.5 + 0.5 \, \frac{v_{\max} - \mathrm{softplus}(k\,(0.5-p))}{v_{\max} - v_{\min}}\\
\mathrm{softplus}(x) &= \ln (1 + \exp(x))\\
v_{\min} &= \mathrm{softplus}(-0.5 k)\\
v_{\max} &= \mathrm{softplus}(0.5 k)\\
\end{aligned}
$$

コードにします。

```python
def panStereoSoftplus(pan, param=0.5):
    """`param` is in [0.0, 1.0]."""
    def softplus(x):
        return np.log(1 + np.exp(x))

    def gain(p, k):
        v_min = softplus(-0.5 * k)
        v_max = softplus(0.5 * k)
        return 0.5 + 0.5 * (v_max - softplus(k * (0.5 - p))) / (v_max - v_min)

    saturation = np.power(10, 2 * param)
    return getStereoGain(gain, saturation)
```

以下の図は softplus による漏れのあるフェード曲線のプロットです。 $k = 10$ です。

<figure>
<img src="img/panStereoSoftplus.svg" alt="Plot of leaky curve using softplus function for stereo to stereo panning." style="padding-bottom: 12px;"/>
</figure>

- [Rectifier (neural networks) - Wikipedia](https://en.wikipedia.org/wiki/Rectifier_(neural_networks)#Softplus)

### 波打つフェード (Sinc)
$G_{LL}$ の式を立てます。

$$
G_{LL}(p) = \begin{cases}
  0.5 + 0.5 \, \mathrm{sinc}(k (1 - 2p))
    & \text{if} \enspace p \leq 0.5\\
  1
    & \text{if} \enspace 0.5 \leq p\\
\end{cases}\\
$$

コードにします。

```python
def panStereoSinc(pan, param=0.5):
    """`param` is in [0.0, 1.0]."""
    def gain(p, k):
        return np.where(
            p <= 0.5,
            0.5 + 0.5 * np.sinc(k * (1 - 2 * p)),
            1,
        )

    zeroCross = 1 + np.floor(16 * param)
    return getStereoGain(gain, zeroCross)
```

以下の図は波打つフェード曲線のプロットです。 $k = 9$ です。

<figure>
<img src="img/panStereoSinc.svg" alt="Plot of bouncing curve using sinc function for stereo to stereo panning." style="padding-bottom: 12px;"/>
</figure>

- [Sinc function - Wikipedia](https://en.wikipedia.org/wiki/Sinc_function)

## 音のサンプル
データ量を減らすため圧縮に [opus](https://opus-codec.org/) を使っています。

### モノラル → ステレオパンニング
<figure>
  <figcaption>線形パンニング</figcaption>
  <audio controls>
    <source src="snd/opus/panMonoLinear.opus" type="audio/ogg; codecs=opus">
  </audio>
</figure>

<figure>
  <figcaption>等パワーパンニング</figcaption>
  <audio controls>
    <source src="snd/opus/panMonoConstantPower.opus" type="audio/ogg; codecs=opus">
  </audio>
</figure>

<figure>
  <figcaption>-4.5 dB パンニング</figcaption>
  <audio controls>
    <source src="snd/opus/panMonoIntermediatePower.opus" type="audio/ogg; codecs=opus">
  </audio>
</figure>

### ステレオ → ステレオパンニング
#### 線形パンニング
<figure>
  <figcaption>直線的なフェード</figcaption>
  <audio controls>
    <source src="snd/opus/panStereoLinear_panMonoLinear.opus" type="audio/ogg; codecs=opus">
  </audio>
</figure>

<figure>
  <figcaption>部分的に 2 次曲線を使うフェード, a = 0.4</figcaption>
  <audio controls>
    <source src="snd/opus/panStereoPartial2ndOrder_panMonoLinear.opus" type="audio/ogg; codecs=opus">
  </audio>
</figure>

<figure>
  <figcaption>部分的に sin を使うフェード, a = 0.4</figcaption>
  <audio controls>
    <source src="snd/opus/panStereoPartialSin_panMonoLinear.opus" type="audio/ogg; codecs=opus">
  </audio>
</figure>

<figure>
  <figcaption>円形フェード, x = 1</figcaption>
  <audio controls>
    <source src="snd/opus/panStereoCircle_panMonoLinear.opus" type="audio/ogg; codecs=opus">
  </audio>
</figure>

<figure>
  <figcaption>n 次曲線フェード, n = 2</figcaption>
  <audio controls>
    <source src="snd/opus/panStereoPoly_panMonoLinear.opus" type="audio/ogg; codecs=opus">
  </audio>
</figure>

<figure>
  <figcaption>sin を使うフェード</figcaption>
  <audio controls>
    <source src="snd/opus/panStereoSin_panMonoLinear.opus" type="audio/ogg; codecs=opus">
  </audio>
</figure>

<figure>
  <figcaption>S 字のフェード, n = 1</figcaption>
  <audio controls>
    <source src="snd/opus/panStereoSCurve_panMonoLinear.opus" type="audio/ogg; codecs=opus">
  </audio>
</figure>

<figure>
  <figcaption>漏れのあるフェード (Softplus), k = 10</figcaption>
  <audio controls>
    <source src="snd/opus/panStereoSoftplus_panMonoLinear.opus" type="audio/ogg; codecs=opus">
  </audio>
</figure>

<figure>
  <figcaption>波打つフェード (Sinc), k = 9</figcaption>
  <audio controls>
    <source src="snd/opus/panStereoSinc_panMonoLinear.opus" type="audio/ogg; codecs=opus">
  </audio>
</figure>

#### 等パワーパンニング
<figure>
  <figcaption>直線的なフェード</figcaption>
  <audio controls>
    <source src="snd/opus/panStereoLinear_panMonoConstantPower.opus" type="audio/ogg; codecs=opus">
  </audio>
</figure>

<figure>
  <figcaption>部分的に 2 次曲線を使うフェード, a = 0.4</figcaption>
  <audio controls>
    <source src="snd/opus/panStereoPartial2ndOrder_panMonoConstantPower.opus" type="audio/ogg; codecs=opus">
  </audio>
</figure>

<figure>
  <figcaption>部分的に sin を使うフェード, a = 0.4</figcaption>
  <audio controls>
    <source src="snd/opus/panStereoPartialSin_panMonoConstantPower.opus" type="audio/ogg; codecs=opus">
  </audio>
</figure>

<figure>
  <figcaption>円形フェード, x = 1</figcaption>
  <audio controls>
    <source src="snd/opus/panStereoCircle_panMonoConstantPower.opus" type="audio/ogg; codecs=opus">
  </audio>
</figure>

<figure>
  <figcaption>n 次曲線フェード, n = 2</figcaption>
  <audio controls>
    <source src="snd/opus/panStereoPoly_panMonoConstantPower.opus" type="audio/ogg; codecs=opus">
  </audio>
</figure>

<figure>
  <figcaption>sin を使うフェード</figcaption>
  <audio controls>
    <source src="snd/opus/panStereoSin_panMonoConstantPower.opus" type="audio/ogg; codecs=opus">
  </audio>
</figure>

<figure>
  <figcaption>S 字のフェード, n = 1</figcaption>
  <audio controls>
    <source src="snd/opus/panStereoSCurve_panMonoConstantPower.opus" type="audio/ogg; codecs=opus">
  </audio>
</figure>

<figure>
  <figcaption>漏れのあるフェード (Softplus), k = 10</figcaption>
  <audio controls>
    <source src="snd/opus/panStereoSoftplus_panMonoConstantPower.opus" type="audio/ogg; codecs=opus">
  </audio>
</figure>

<figure>
  <figcaption>波打つフェード (Sinc), k = 9</figcaption>
  <audio controls>
    <source src="snd/opus/panStereoSinc_panMonoConstantPower.opus" type="audio/ogg; codecs=opus">
  </audio>
</figure>

#### -4.5 dB パンニング
<figure>
  <figcaption>直線的なフェード</figcaption>
  <audio controls>
    <source src="snd/opus/panStereoLinear_panMonoIntermediatePower.opus" type="audio/ogg; codecs=opus">
  </audio>
</figure>

<figure>
  <figcaption>部分的に 2 次曲線を使うフェード, a = 0.4</figcaption>
  <audio controls>
    <source src="snd/opus/panStereoPartial2ndOrder_panMonoIntermediatePower.opus" type="audio/ogg; codecs=opus">
  </audio>
</figure>

<figure>
  <figcaption>部分的に sin を使うフェード, a = 0.4</figcaption>
  <audio controls>
    <source src="snd/opus/panStereoPartialSin_panMonoIntermediatePower.opus" type="audio/ogg; codecs=opus">
  </audio>
</figure>

<figure>
  <figcaption>円形フェード, x = 1</figcaption>
  <audio controls>
    <source src="snd/opus/panStereoCircle_panMonoIntermediatePower.opus" type="audio/ogg; codecs=opus">
  </audio>
</figure>

<figure>
  <figcaption>n 次曲線フェード, n = 2</figcaption>
  <audio controls>
    <source src="snd/opus/panStereoPoly_panMonoIntermediatePower.opus" type="audio/ogg; codecs=opus">
  </audio>
</figure>

<figure>
  <figcaption>sin を使うフェード</figcaption>
  <audio controls>
    <source src="snd/opus/panStereoSin_panMonoIntermediatePower.opus" type="audio/ogg; codecs=opus">
  </audio>
</figure>

<figure>
  <figcaption>S 字のフェード, n = 1</figcaption>
  <audio controls>
    <source src="snd/opus/panStereoSCurve_panMonoIntermediatePower.opus" type="audio/ogg; codecs=opus">
  </audio>
</figure>

<figure>
  <figcaption>漏れのあるフェード (Softplus), k = 10</figcaption>
  <audio controls>
    <source src="snd/opus/panStereoSoftplus_panMonoIntermediatePower.opus" type="audio/ogg; codecs=opus">
  </audio>
</figure>

<figure>
  <figcaption>波打つフェード (Sinc), k = 9</figcaption>
  <audio controls>
    <source src="snd/opus/panStereoSinc_panMonoIntermediatePower.opus" type="audio/ogg; codecs=opus">
  </audio>
</figure>

## 参考文献
- [Loudness Concepts & Panning Laws](http://www.cs.cmu.edu/~music/icm-online/readings/panlaws/)

## 変更点
- 2021/01/13
  - 記号の間違いを修正。
  - 文章の整理。
  - Softplus と sinc 関数のリンクを追加
- 2021/01/18
  - 文章の整理。
