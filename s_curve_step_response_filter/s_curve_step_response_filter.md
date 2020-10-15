# ステップ応答が S 字を描くフィルタ
リミッタのエンベロープのスムーシングに使うためにステップ応答が S 字を描くフィルタについて調べました。

次の 3 つのフィルタを比較します。

- 移動平均フィルタの重ね掛け
- Bessel フィルタ
- Thiran ローパスフィルタ

## 移動平均フィルタの重ね掛け
移動平均フィルタは [Musicdsp.org の Lookahead Limiter](https://www.musicdsp.org/en/latest/Effects/274-lookahead-limiter.html) のページで紹介されている手法を使えば効率よく計算できます。この手法では重ね掛けを 1 回行ってフィルタ係数を矩形窓から三角窓の形に変えています。

### 移動平均フィルタの重ね掛け
直列につないだ FIR フィルタは、フィルタ係数を畳み込むことで 1 つのフィルタにまとめられます。矩形窓を畳み込んで、重ね掛けしたときのフィルタ係数がどうなるのか見てみます。

```python
# Python3
import matplotlib.pyplot as plt
import numpy as np
import scipy.signal as signal

def average(fir):
    denom = np.sum(fir)
    if denom == 0:
        return fir
    return fir / denom

def movingAverage(nTaps, integrate=0):
    length = nTaps // 2**integrate
    taps1 = average(np.ones(length))
    for _ in range(integrate):
        taps2 = average(np.ones(length + 1)) # フィルタのタップ数をそろえるために +1 。
        taps1 = average(signal.convolve(taps1, taps2))
        length *= 2
    return taps1

def plotFilter(nPlot=4, nTaps = 512):
    cmap = plt.get_cmap("viridis")
    for order in range(nPlot):
        taps = movingAverage(nTaps, order)
        plt.plot(taps, lw=1, color=cmap(order / nPlot), label=f"{order}")
    plt.grid()
    plt.legend()
    plt.show()

plotFilter()
```

出力です。 1 回畳み込んだときに三角窓になっています。 2 回目以降の畳み込みでは[区間ごとに変わる多項式 (piecewise polynomial) になる](https://www.cs.cmu.edu/afs/cs/academic/class/15462-f11/www/homework/hw2_p2_sol.pdf)ようです。

![Convoluted moveing average filter kernel shape.](img/MovingAverageKernel.png)

### ステップ応答
移動平均フィルタを畳み込んだフィルタのステップ応答を確認します。

```python
# Python3
def plotStep(nPlot=4, nTaps = 512):
    sig = np.hstack((np.zeros(nTaps), np.ones(2 * nTaps)))
    cmap = plt.get_cmap("viridis")
    for order in range(nPlot):
        taps = movingAverage(nTaps, order)
        out = signal.lfilter(taps, 1, sig)
        plt.plot(out, lw=1, color=cmap(order / nPlot), label=f"{order}")
    plt.axvline(2 * nTaps, lw=1, color="black", alpha=0.5, ls="--")
    plt.plot(sig, lw=1, color="black", label="Input")
    plt.grid()
    plt.legend()
    plt.show()

plotStep()
```

出力です。 1 回以上畳み込んだフィルタはステップ応答が S 字になっています。

![Moving average and its convolution filters step responses.](img/MovingAverageStepResponse.png)

### 実装
Musicdsp.org の[Lookahead Limiter](https://www.musicdsp.org/en/latest/Effects/274-lookahead-limiter.html) の記事で紹介されていた手法を使って実装します。

整数サンプルの遅延を加えるディレイを用意します。

```c++
// C++
#include <algorithm>
#include <vector>

template<typename Sample> class Delay {
public:
  std::vector<Sample> buffer;
  int32_t wptr = 0;
  int32_t rptr = 0;

  Delay(Sample size = 65536) : buffer(size) {}

  void resize(uint32_t size) {
    buffer.resize(size);
    wptr = 0;
    rptr = 0;
  }

  void reset() { std::fill(buf.begin(), buf.end(), 0); }

  void setFrames(uint32_t delayFrames)
  {
    if (delayFrames >= buf.size()) delayFrames = bufEnd;
    rptr = wptr - int32_t(delayFrames);
    if (rptr < 0) rptr += int32_t(buf.size());
  }

  Sample process(Sample input)
  {
    wptr = (wptr + 1) % int32_t(buf.size());
    buf[wptr] = input;

    rptr = (rptr + 1) % int32_t(buf.size());
    return buf[rptr];
  }
};
```

状態変数 `sum` を用意して、各サンプルでディレイへの入力を加算、ディレイからの出力を減算します。 `sum` をフィルタのタップ数で割った値がフィルタの出力です。

```c++
// C++
template<typename Sample> class MovingAverage {
public:
  Sample sum = 0;
  Delay<Sample> delay;

  void reset() {
    sum = 0;
    delay.reset();
  }

  void setSize(uint32_t size) {
    delay.setFrames(size);
  }

  Sample process(Sample input) {
    sum += input - delay.process(input);
    return sum / delay.buffer.size();
  }
};
```

フィルタのタップ数を $n$ とします。畳み込みを直接計算すると 1 サンプルあたりの計算量は $O(n)$ です。 ディレイを使う方法ならタップ数によらず定数回の計算で畳み込みができるので 1 サンプルあたりの計算量が $O(1)$ に減ります。

#### 補足
Musicdsp.org のページでは、以下のコードのように入力で `1.0 - input` 、出力で `1.0 - output` を計算することで浮動小数点数による計算誤差を減らす手順がありますが、意味があるのかどうかわからないので実装していません。はっきりとは書いていませんが `input` の範囲は `[0.0, 1.0]` のようです。

```c++
// C++
template<typename Sample> class MovingAverage {
public:
  // ...

  Sample process(Sample input) {
    input = Sample(1) - input;
    sum += input - delay.process(input);
    return Sample(1) - sum / delay.buffer.size();
  }
};
```

## Bessel フィルタ
Bessel フィルタは連続系で群遅延が maximally flat になることが特長の IIR フィルタです。離散化にバイリニア変換などを使うと周波数が歪むので特長が損なわれます。

### ステップ応答
SciPy の `signal.bessel()` を試します。 SciPy 1.5.2 では `signal.bessel()` はバイリニア変換を使っています。

```python
# Python3
import matplotlib.pyplot as plt
import numpy as np
import scipy.signal as signal

order = 4
delay = 512

sos = signal.bessel(order, 2 / (np.pi * delay), norm="delay", output="sos")

sig = np.hstack((np.zeros(delay), np.ones(delay * 2)))
out = signal.sosfilt(sos, sig)

plt.axvline(2 * delay, lw=1, color="black", alpha=0.5, ls="--")
plt.plot(sig, lw=1, color="black", label="Input")
plt.plot(out, lw=1, color="red", label="Bessel")
plt.grid()
plt.legend()
plt.show()
```

出力です。ステップ応答は S 字になりますがオーバーシュートが出ます。

![Bessel filter step response.](img/BesselStepResponse.png)

### 2 次セクションに分割したときのフィルタ係数
カットオフ周波数を自由に変更するために、以下の Python3 のコードと、 Neil Robertson さんによる [Design IIR Filters Using Cascaded Biquads](https://www.dsprelated.com/showarticle/1137.php) で紹介されていた手法を組み合わせて C++ に移植します。

```python
# Python3
import numpy as np
import scipy.signal as signal

def bessel(order, delay):
    z, p, k = signal.besselap(order, norm="delay")
    z, p, k = signal.lp2lp_zpk(z, p, k, delay)
    z, p, k, _ = signal.cont2discrete([z, p, k], 1, method="bilinear")
    return signal.zpk2sos(z, p, k)
```

コードの `z, p, k` は `zero, pole, gain` の略です。

#### アナログプロトタイプ
Bessel フィルタは極だけでゼロは一つもありません。このようなフィルタは全極フィルタ (all-pole filter) と呼ばれています。極は伝達関数の分母の多項式の根を求めることで計算できます。

今回は root-finding アルゴリズムを評価、実装する手間を省きたかったので `signal.besselap()` から得られた極をそのまま定数として使うことにしました。

```python
# Python3
import numpy as np
import scipy.signal as signal
z, p, k = signal.besselap(4, norm="delay") # 4 次。
print(p.tolist())
```

整形した出力です。

```
[
  (-2.1037893971796273+2.6574180418567526j), # p1
  (-2.8962106028203722+0.8672341289345038j), # p2
  (-2.8962106028203722-0.8672341289345038j), # p2 の複素共役
  (-2.1037893971796273-2.6574180418567526j)  # p1 の複素共役
]
```

#### カットオフ周波数の設定
`signal.lp2lp_zpk()` はローパスのアナログプロトタイプのカットオフ周波数を設定する関数です。以下はエラー処理などを簡略化した SciPy 1.5.2 での実装の写しです。

```python
# Python3
def lp2lp_zpk(z, p, k, wo=1.0):
    degree = len(p) - len(z)
    z_lp = wo * z
    p_lp = wo * p
    k_lp = k * wo**degree
    return z_lp, p_lp, k_lp
```

Bessel フィルタはゼロが無く、ゲインは別の方法で計算するので `z_lp` と `k_lp` は C++ のコードでは計算していません。

#### バイリニア変換
アナログプロトタイプからの離散化にはバイリニア変換を使います。次のバイリニア変換の式の $s$ に連続系の極を代入することで、離散フィルタの極に変換できます。

$$
z = \frac{2 f_s + s}{2 f_s - s}
$$

`signal.besselap()` から得られた極は $f_s = 1$ とすると正しく変換できました。

$N$ 次の全極フィルタの伝達関数は次の形になります。 $\hat{p}$ はアナログフィルタの極です。

$$
H(z) = \frac{1}{\displaystyle \prod_{i=1}^{N} (s - \hat{p}_i)}
$$

$s$ について解いたバイリニア変換の式です。

$$
s = 2 f_s \frac{1 - z^{-1}}{1 + z^{-1}}
$$

$s$ について解いたバイリニア変換の式を全極フィルタの伝達関数に代入してまとめると次の形になります。 $p$ は離散フィルタの極、 $K$ はゲインです。

$$
H(z) = K \frac{(z + 1)^N}{\displaystyle \prod_{i=1}^{N} (z - p_i)}
$$

#### 2 次セクションの構築
次数が偶数の Bessel フィルタの極は複素共役の組 (conjugated pairs) になります。例えば 4 次なら $(p_1, p_2, p_2^*, p_1^*)$ という 4 つの極が出てきます。 $^*$ は複素共役の演算子です。 2 次セクションを作るときは複素共役の組となる極を使うことで計算を簡略化できます。

$$
\begin{aligned}
H(z)
  &= K \prod_{i=1}^{N / 2} \frac{(z + 1)^2}{(z - p_i)(z - p_i^*)} \\
  &= K \prod_{i=1}^{N / 2} \frac{1 + 2 z^{-1} + z^{-2}}{
    1 - 2\,\mathrm{Re}(p_i) z^{-1} + |p_i|^2 z^{-2}
  }
\end{aligned}
$$

一般的な 2次 IIR フィルタの伝達関数は次のように書けます。

$$
H(z) = \frac{b_0 + b_1 z^{-1} + b_2 z^{-2}}{a_0 + a_1 z^{-1} + a_2 z^{-2}}
$$

$z = e^{j \omega}$ と代入して $\omega = 0$ のときに $H(z) = 1$ となるようなゲイン $K$ は次の式で計算できます。 $z = e^{j \omega}$ と代入するのは振幅応答を求めるときの計算手法です。

$$
\begin{aligned}
H(z) = 1 &= K \frac{b_0 + b_1 + b_2}{a_0 + a_1 + a_2}\\
K &= \frac{a_0 + a_1 + a_2}{b_0 + b_1 + b_2}
\end{aligned}
$$

Bessel フィルタでは $a_0 = 1,\ b_0 + b_1 + b_2 = 4$ で固定なので $K = \dfrac{1 + a_1 + a_2}{4}$ です。

#### C++ での実装
C++17 のフラグを立ててコンパイルしています。

```c++
// C++
#include <array>
#include <complex>

template<typename Sample> struct Bessel4 {
  constexpr static uint8_t halfDegree = 4 / 2;

  // 複素共役の極は無くても計算できる。
  constexpr static std::array<std::complex<Sample>, halfDegree> analogPole{{
    {Sample(-2.1037893971796273), Sample(+2.6574180418567526)},
    {Sample(-2.8962106028203722), Sample(+0.8672341289345038)},
  }};

  std::array<Sample, 2> x0{};
  std::array<Sample, 2> x1{};
  std::array<Sample, 2> x2{};
  std::array<Sample, 2> y0{};
  std::array<Sample, 2> y1{};
  std::array<Sample, 2> y2{};
  std::array<std::array<Sample, 2>, halfDegree> co; // フィルタ係数 a1 と a2 。
  Sample gain = 1;

  Bessel4()
  {
    for (auto &coef : co) coef.fill(0);
  }

  void reset(Sample value = 0)
  {
    x0.fill(value);
    x1.fill(value);
    x2.fill(value);
    y0.fill(value);
    y1.fill(value);
    y2.fill(value);
  }

  void setDelay(Sample delayInSamples)
  {
    auto wo = Sample(0.5) / delayInSamples; // 遅延サンプル数を周波数に変換。

    gain = Sample(1);
    for (uint8_t i = 0; i < co.size(); ++i) {
      std::complex<Sample> pole = wo * analogPole[i]; // カットオフ周波数の適用。
      pole = (Sample(2) + pole) / (Sample(2) - pole); // バイリニア変換。
      co[i][0] = Sample(-2) * pole.real();
      co[i][1] = std::norm(pole);
      gain *= (Sample(1) + co[i][0] + co[i][1]) / Sample(4);
    }
  }

  Sample process(Sample input)
  {
    x0[0] = input;
    // x0[1] = y0[0];
    for (uint8_t i = 1; i < halfDegree; ++i) x0[i] = y0[i - 1];

    for (uint8_t i = 0; i < halfDegree; ++i) {
      y0[0]
        = Sample(1) * x0[i]
        + Sample(2) * x1[i]
        + Sample(1) * x2[i]
        - co[i][0] * y1[i]
        - co[i][1] * y2[i];

      x2[i] = x1[i];
      x1[i] = x0[i];
      y2[i] = y1[i];
      y1[i] = y0[i];
    }

    return gain * y0[1];
  }
};
```

## Thiran ローパスフィルタ
Thiran ローパスフィルタは群遅延が maximally flat になるように設計されたフィルタです。 Maximally flat はフィルタの専門用語で、通過域の特性ができるだけ平坦になるという意味です。

Thiran ローパスフィルタの伝達関数です。

$$
H(z) = \frac{2n!}{n!} \,
\frac{1}{\displaystyle \prod_{i = n + 1}^{2n} (2 \tau + i)} \,
\frac{1}{\displaystyle
  \sum_{k=0}^n \left(
    (-1)^k \binom{n}{k}
    \prod_{i=0}^n \dfrac{2 \tau + i}{2 \tau + k + i}
  \right)
  z^{-k}
}
$$

- $n$: フィルタの次数。
- $\tau$: 所望の遅延 (desired delay) 。

伝達関数をそのまま実装すると発散することがあります。分母の因数分解は無理そうなので、 2 次セクションに分割するなら多項式の根を探す (root-finding) アルゴリズムを使って数値的に計算するしかなさそうです。ここでは SciPy の `tf2sos` を使っています。 C++ では [Boost の root-finding](https://www.boost.org/doc/libs/1_74_0/libs/math/doc/html/root_finding.html) が使えそうです。

- TODO: 係数の計算とテスト


### ステップ応答
ステップ応答を確認します。

```python
# Python3
import numpy as np
import scipy.signal as signal
import scipy.special as special
import matplotlib.pyplot as plt

def thiranLowpass(n, τ):
    """
    n: Order of filter.
    τ: Desired delay in samples.
    """
    a = np.empty(n + 1)
    for k in range(n + 1):
        i = np.arange(n + 1)
        a[k] = (-1)**k * special.binom(n, k) * np.prod((2 * τ + i) / (2 * τ + k + i))
    gain_numer = special.factorial(2 * n, exact=True) / special.factorial(n, exact=True)
    gain_denom = np.prod(2 * τ + np.arange(n + 1, 2 * n + 1))
    b = np.zeros_like(a)
    b[0] = gain_numer / gain_denom
    return (b, a)

order = 4
delay = 512

b, a = thiranLowpass(order, delay / 2)
sos = signal.tf2sos(b, a)

sig = np.hstack((np.zeros(delay), np.ones(delay * 2)))
out = signal.sosfilt(sos, sig)

plt.axvline(2 * delay, lw=1, color="black", alpha=0.5, ls="--")
plt.plot(sig, lw=1, color="black", label="Input")
plt.plot(out, lw=1, color="red", label="Thiran")
plt.grid()
plt.legend()
plt.show()
```

出力です。オーバーシュートが出ています。バイリニア変換した Bessel フィルタと同じに見えます。

![Thiran lowpass filter step response.](img/ThiranStepResponse.png)

ゲインの値を定義通りに計算すると `delay` の値によってはステップ応答が 1 から少しずれた値に収束することがあります。 1 に収束させるときは `signal.freqz()` の値からゲインを計算できます。

```python
# Python3
def thiranLowpassFixedGain(n, τ):
    a = np.empty(n + 1)
    for k in range(n + 1):
        i = np.arange(n + 1)
        a[k] = (-1)**k * special.binom(n, k) * np.prod((2 * τ + i) / (2 * τ + k + i))
    freq, resp = signal.freqz(1, a, worN=1)
    b = np.zeros_like(a)
    b[0] = 1 / resp.real
    return (b, a)
```

## 比較
ここまでに紹介した 3 つのフィルタのステップ応答を比較します。

- 三角窓の FIR フィルタ (移動平均を 1 回畳み込み)
- バイリニア変換した Bessel フィルタ
- Thiran ローパスフィルタ

Bessel フィルタと Thiran ローパスフィルタのステップ応答はほとんど同じで重なっています。

![Comparison of step responase of triangular, Bessel and Thiran filter.](img/StepResponseComparison.png)

離散系で S 字のステップ応答が欲しいときは、オーバーシュートがなく、計算も速くて簡単な移動平均の重ね掛けを使えばよさそうです。

バイリニア変換した Bessel フィルタは、オーバーシュートがあり、移動平均の重ね掛けと比べると計算が重いですが、メモリの使用量が少ないという利点があります。移動平均の重ね掛けは遷移時間と同じサンプル数のメモリが必要です。

Thiran ローパスフィルタはバイリニア変換した Bessel フィルタとほとんど同じですが、 2 次セクションに分割できないので使いどころはなさそうです。

## 参考文献
- [Lookahead Limiter — Musicdsp.org documentation](https://www.musicdsp.org/en/latest/Effects/274-lookahead-limiter.html)
- [hw2_p2_sol.pdf - 矩形窓と三角窓の畳み込み (www.cs.cmu.edu)](https://www.cs.cmu.edu/afs/cs/academic/class/15462-f11/www/homework/hw2_p2_sol.pdf)
- [Bessel filter - Wikipedia](https://en.wikipedia.org/wiki/Bessel_filter)
- [Design IIR Butterworth Filters Using 12 Lines of Code - Neil Robertson](https://www.dsprelated.com/showarticle/1119.php)
- [Design IIR Filters Using Cascaded Biquads - Neil Robertson](https://www.dsprelated.com/showarticle/1137.php)
- [Thiran, J. P. (1971). Recursive digital filters with maximally flat group delay.](https://pdfs.semanticscholar.org/cd84/b18a454fd04a7bfcd1995ddb8c4c6927cc87.pdf)
