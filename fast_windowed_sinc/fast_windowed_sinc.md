# Windowed Sinc フィルタの高速な計算
再帰的に sin を計算する方法を使った windowed sinc フィルタの高速な計算方法について紹介します。

ここではカットオフ周波数と群遅延が 1 サンプルごとに変更される応用を想定しています。以下は主な応用です。

- クロスオーバー周波数が可変のバンド分割
- ディレイ時間変更時のアンチエイリアシング
- 楽器のサンプラーの補間

ここで紹介する高速な計算方法を使うと計算精度は下がります。フィルタ係数を固定できるときには使わないでください。

検証に使ったすべてのコードは以下のリンクから閲覧できます。

- [記事の内容の検証に使ったコード (github.com)](https://github.com/ryukau/filter_notes/tree/master/fast_windowed_sinc)

## sinc 関数
通常、 sinc 関数は以下のような形で定義されています。

$$
\mathrm{sinc}(x) = \frac{\sin(\pi x)}{\pi x}.
$$

ここではローパスのカットオフ周波数 $f_c$ を考慮した以下の式を扱います。

$$
\mathrm{sinc}(x, f_c) = \frac{\sin(2 \pi f_c x)}{\pi x}.
$$

$f_c$ の単位は rad / 2π で、値の範囲が [0.0, 0.5] となるように正規化されてます。

## 正確な実装
以下は矩形窓を用いた windowed sinc フィルタを設計する Python 3 のコードです。

```python
import numpy as np

def modifiedSinc(x, cutoff):
    if x == 0:
        return 2 * cutoff
    return np.sin(np.pi * 2 * cutoff * x) / (np.pi * x)

def lowpassFir(length, cutoff, fractionSample):
    mid = fractionSample - (length // 2 + length % 2)
    fir = np.zeros(length)
    for i in range(length):
        x = i + mid
        fir[i] = modifiedSinc(x, cutoff)
    return fir
```

以下は `lowpassFir` のパラメータの意味です。

- `length` : FIR フィルタのタップ数。
- `cutoff` : 正規化されたカットオフ周波数。単位は rad / 2π 。範囲は [0, 0.5] 。
- `fractionSample` : サンプル数であらわされた小数点以下の群遅延の量。範囲は [0, 1] 。

`fractionSample` は向きがあります。ここでは配列 `fir` のインデックスをさかのぼる方向に固定しています。

## 高速な実装
以下は biquad オシレータと呼ばれる再帰的に sin を計算する方法を使った高速な実装です。見やすさのために Python で書いていますが、特に CPython では `for` が遅いので高速化の恩恵は薄いです。 C++ による実装は「[ディレイのアンチエイリアシング](#ディレイのアンチエイリアシング)」を参照してください。

```python
def lowpassFirBiquadPede(length: int, cutoff: float, fractionSample: float):
    mid = fractionSample - length // 2 - length % 2

    omega = 2 * np.pi * cutoff
    phi = mid * omega
    k = 2 * np.cos(omega)
    u1 = np.sin(phi - omega)
    u2 = np.sin(phi - 2 * omega)

    fir = np.zeros(length)
    for i in range(length):
        u0 = k * u1 - u2
        u2 = u1
        u1 = u0

        x = i + mid
        if abs(x) < 0.1:
            q = np.pi * cutoff * x
            q *= q
            fir[i] = (2 / 3) * cutoff * (15 - 7 * q) / (5 + q)
        else:
            fir[i] = u0 / (np.pi * x)
    return fir
```

再帰的な sin の計算に biquad オシレータを使っています。

`fractionSample` が 0 または 1 に近いときに `x == 0` の周りで値がおかしくなるので、 `if abs(x) < 0.1` の分岐を設けて [Padé approximant](https://en.wikipedia.org/wiki/Pad%C3%A9_approximant) での計算に切り替えています。分岐での近似式は以下の [Maxima](https://maxima.sourceforge.io/) のコードから出力されたものです。

```maxima
sinc: sin(2 * %pi * f_c * x) / (%pi * x);
pade(taylor(sinc, x, 0, 4), 2, 2);
```

以下は整理した出力です。

$$
\mathrm{sinc} (x, f_c) \approx \frac{2}{3} f_c \frac{15 - 7 (\pi f_c x)^2}{5 + (\pi f_c x)^2},
\enspace \text{only when} \ x\ \text{is close to 0}.
$$

Padé approximant への分岐 `if abs(x) < 0.1` のしきい値 `0.1` は適当に決めた値です。適切な分岐先は `cutoff` に応じて変わります。 `cutoff` が低くなるときはしきい値を上げて、 Padé approximant によって計算する範囲を広くしたほうが誤差が減る傾向があります。逆に `cutoff` が正規化されたナイキスト周波数 `0.5` に近いときは、しきい値を下げて `0.001` あたりに設定すると誤差が減る傾向があります。 `cutoff` に応じた適切なしきい値を調べるコードを以下のリンク先の `findBranchingThreshold` に掲載しています。

- [二分探索でしきい値を調べる C++ のコード (github.com)](https://github.com/ryukau/filter_notes/blob/a08bb919d2c2c07e7e21561fd7e3ee12d4c055b2/fast_windowed_sinc/fast_sinc_error/error.cpp#L435-L523)

### 誤差と実装のバリエーション
以下は誤差が変わる定性的な要素の一覧です。

- コンパイラの最適化オプション
- 再帰的な sin の計算方法
- 初期位相の計算方法

`/fp:fast` や `-ffast-math` のような浮動小数点数の計算順序を入れ替えるコンパイラの最適化オプションを使うと、以下で検討するような計算方法の変更では誤差の低減が見込みにくくなります。したがって、以下、この節では `/fp:fast` や `-ffast-math` を使わないことを前提としています。

再帰的な sin の計算方法はいくつか種類があり、別記事の「 [sin, cos を反復的に計算するアルゴリズムのレシピ](../recursive_sine/recursive_sine.html)」にレシピを載せています。簡単に調べたところでは coupled form と呼ばれる形を使うと biquad よりは誤差が減ります。

再帰的な計算を行うので、初期位相の細かいずれによって誤差が変わります。例えば、上で紹介した biquad オシレータを使った高速な実装は `u1` と `u2` の計算において以下の形を使うことで、再帰的な sin の計算精度を少し良くできます。

```python
# 変形した k, u1, u2 の初期値設定。
cω = np.cos(omega)
sω = np.sin(omega)
cφ = np.cos(phi)
sφ = np.sin(phi)
k = 2 * cω
u1 = sφ * cω - cφ * sω    # ~= sin(phi - omega), 加法定理。
c2ω = 2 * cω * cω - 1     # ~= cos(2 * omega), 倍角公式。
s2ω = 2 * sω * cω         # ~= sin(2 * omega), 倍角公式。
u2 = sφ * c2ω - cφ * s2ω  # ~= sin(phi - 2 * omega), 加法定理。
```

## Cosine-sum 窓
再帰的な sin の計算を行うオシレータを 1 つ増やすことで、 [cosine-sum 窓](https://en.wikipedia.org/wiki/Window_function#Cosine-sum_windows)と呼ばれる種類の窓関数を高速にかけることができます。ここでは cosine-sum 窓の 1 つである Blackman-Harris 窓について具体的な実装を紹介します。

以下は Blackman-Harris 窓の計算式です。

$$
\begin{aligned}
w[n]&=a_0 - a_1 \cos \left ( \frac{2 \pi n}{N} \right)+ a_2 \cos \left ( \frac{4 \pi n}{N} \right)- a_3 \cos \left ( \frac{6 \pi n}{N} \right),\\
a_0&=0.35875,\quad a_1=0.48829,\quad a_2=0.14128,\quad a_3=0.01168.
\end{aligned}
$$

- $n$: インデックス。
- $N$: 窓長。

この式は $\cos$ の位相が 2 倍、 3 倍となっているので、以下のように [Chebyshev 多項式](https://en.wikipedia.org/wiki/Chebyshev_polynomials) $T_k$ を使って書くことができます。 $T_k$ の中身は [Wikipedia の記事](https://en.wikipedia.org/wiki/Chebyshev_polynomials#First_kind)に掲載されています。

$$
\begin{aligned}
w[n] &= a_0 - a_1 T_1(c) + a_2 T_2(c) - a_3 T_3(c)\\
&\begin{aligned}
T_1(c) &= c\\
T_2(c) &= 2 c^2 - 1\\
T_3(c) &= 4 c^3 - 3c\\
\end{aligned}\qquad
c = \cos \left ( \frac{2 \pi n}{N} \right).
\end{aligned}
$$

三角関数の計算を $c$ の 1 回のみに減らすことができました。 Chebyshev 多項式を Maxima で展開します。

```maxima
T2: 2*c^2 - 1;
T3: 4*c^3 - 3*c;
w: 0.35875 - 0.48829 * c + 0.14128 * T2 - 0.01168 * T3;
ratsimp(w);
```

以下は整形した出力です。これで高速に計算できる形になりました。

$$
\begin{aligned}
w[n] &= 0.21747 - 0.45325 c + 0.28256 c^2 - 0.04672 c^3,\\
c    &= \cos \left ( \frac{2 \pi n}{N} \right).
\end{aligned}
$$

以下は Python 3 による Blackman-Harris 窓のみの実装です。ここでは windowed sinc フィルタの設計に使うので、左右対称な窓関数を計算しています。

```python
def blackmanHarris(length: int):
    isEven = 1 - length % 2

    ω = 2 * np.pi / float(length - isEven)
    φ = np.pi / 2
    k = 2 * np.cos(ω)
    u1 = np.sin(φ - ω)
    u2 = np.sin(φ - 2 * ω)

    window = np.zeros(length)
    for i in range(length):
        u0 = k * u1 - u2
        u2 = u1
        u1 = u0
        window[i] = 0.21747 + u0 * (-0.45325 + u0 * (0.28256 + u0 * -0.04672))
    return window
```

以下は windowed sinc フィルタの設計と同時に Blackman-Harris 窓をかける実装です。上の `blackmanHarris` と、下の `lowpassBlackmanHarrisBiquad` では、窓関数の計算に使うオシレータの `omega, u1, u2` の初期化が異なっています。 `blackmanHarris` は定義通りに計算しているので窓の両端が 0 になりますが、これだと以下で紹介するディレイのアンチエイリアシングへの応用で `length <= 2` のときに音が止まってしまうので、下では変えています。

```python
def lowpassBlackmanHarrisBiquad(length: int, cutoff: float, fractionSample: float):
    isEven = 1 - length % 2

    mid = fractionSample - (length // 2 + length % 2)

    o1_omega = 2 * np.pi * cutoff
    o1_phi = mid * o1_omega
    o1_k = 2 * np.cos(o1_omega)
    o1_u1 = np.sin(o1_phi - o1_omega)
    o1_u2 = np.sin(o1_phi - 2 * o1_omega)

    o2_omega = 2 * np.pi / float(length + isEven)
    o2_phi = np.pi / 2
    o2_k = 2 * np.cos(o2_omega)
    o2_u1 = np.sin(o2_phi - (1 - isEven) * o2_omega)
    o2_u2 = np.sin(o2_phi - (2 - isEven) * o2_omega)

    fir = np.zeros(length)
    for i in range(length):
        o1_u0 = o1_k * o1_u1 - o1_u2
        o1_u2 = o1_u1
        o1_u1 = o1_u0

        o2_u0 = o2_k * o2_u1 - o2_u2
        o2_u2 = o2_u1
        o2_u1 = o2_u0

        x = i + mid
        if abs(x) < 0.1:
            q = np.pi * cutoff * x
            q *= q
            fir[i] = (2 / 3) * cutoff * (15 - 7 * q) / (5 + q)
        else:
            fir[i] = o1_u0 / (np.pi * x)

        window = 0.21747 + o2_u0 * (-0.45325 + o2_u0 * (0.28256 + o2_u0 * -0.04672))
        fir[i] *= window
    return fir
```

### 他の cosine-sum 窓
他の cosine-sum 窓について、 Chebyshev 多項式を展開して [Horner's method](https://en.wikipedia.org/wiki/Horner%27s_method) の形にした式を掲載します。

SymPy で式変形を行います。

```python
import sympy
from sympy.polys.orthopolys import chebyshevt_poly
from sympy.polys.polyfuncs import horner

cosineSumWindowCoefficients = {
    "blackman": [7938 / 18608, 9240 / 18608, 1430 / 18608],
    "nuttall": [0.355768, 0.487396, 0.144232, 0.012604],
    "blackmannuttall": [0.3635819, 0.4891775, 0.1365995, 0.0106411],
    "blackmanharris": [0.35875, 0.48829, 0.14128, 0.01168],
    "flattop": [0.21557895, 0.41663158, 0.277263158, 0.083578947, 0.006947368],
}

u0 = sympy.symbols("u0")
for key, a_ in cosineSumWindowCoefficients.items():
    expr = 0
    for i in range(len(a_)):
        expr += (-1) ** i * a_[i] * chebyshevt_poly(i, x=u0)
    print(f"{key} = {horner(expr)}")
```

以下は出力です。検証を楽にするために整理せずに掲載しています。 `u0` は上に掲載している `blackmanharrisBiquad` で使われている変数と同じで $\cos(2\pi n /N)$ が入ります。 $n$ はインデックス、 $N$ はサンプル数で表された窓長です。

```python
blackman = u0*(0.153697334479794*u0 - 0.496560619088564) + 0.349742046431642
nuttall = u0*(u0*(0.288464 - 0.050416*u0) - 0.449584) + 0.211536
blackmanharris = u0*(u0*(0.28256 - 0.04672*u0) - 0.45325) + 0.21747
blackmannuttall = u0*(u0*(0.273199 - 0.0425644*u0) - 0.4572542) + 0.2269824
flattop = u0*(u0*(u0*(0.055578944*u0 - 0.334315788) + 0.498947372) - 0.165894739) - 0.05473684
```

## ディレイのアンチエイリアシング
高速な windowed sinc フィルタの計算を応用して、ディレイの時間変更時に生じるエイリアシングを低減します。

### 仕組み
大まかなアイデアとしては 1 サンプルごとに変換レートが変わるダウンサンプラーを実装します。

ここでは音の信号を扱うので 0 Hz からナイキスト周波数まで周波数成分が詰め込まれていることにします。

ディレイ時間の変化は 2 つに場合分けできます。

1. ディレイ時間の伸長。 1 倍速よりも遅い読み取り。低いほうへのピッチシフト。
2. ディレイ時間の短縮。 1 倍速よりも速い読み取り。高いほうへのピッチシフト。

1 のときは通常の分数ディレイフィルタ (fractional delay fitler) による補間で十分です。分数ディレイフィルタとは、カットオフがナイキスト周波数に固定されており、群遅延の大きさを与えると設計できるフィルタのことです。

2 のときはピッチシフトによってナイキスト周波数を超える周波数成分が現れます。これを抑えるために、カットオフを変更できる分数ディレイフィルタが必要です。また、ディレイ時間が変わるということは読み取り時刻の逆走もあり得ます。

以上より、補間に使うフィルタには以下が求められます。

- 時間の逆走への対応
- 非整数サンプルの補間
- カットオフ周波数が可変

素朴に畳み込む windowed sinc フィルタはこれらの要件に対応できます。素朴に畳み込むということは、読み取り時刻の周りのサンプルだけを使って計算できるので逆走に対応できます。また windowed sinc フィルタは補間位置 (群遅延) とカットオフ周波数を容易に変更できます。

### カットオフ周波数の設定
ローパスフィルタのカットオフ周波数は、ディレイ時間の変化によるピッチシフトの量から設定できます。以降ではピッチシフトの量のことを単にピッチと呼びます。

ピッチは以下の式で計算できます。

$$
p = d_1 - d_0 + 1.
$$

- $p$ : ピッチ。比率で表された信号の読み取り速度。
- $d_0$ : 現在のディレイ時間。単位はサンプル数。
- $d_1$ : 1 サンプル前のディレイ時間。単位はサンプル数。

ここでのディレイ時間は、ある時点からの相対的な時間です。 $d_0$ と $d_1$ はディレイ時間の起点が 1 サンプルずれているので注意してください。この定義であればピッチ $p$ の値をそのまま使って $p$ 倍速という表現ができます。例えば:

- $d_1$ が 100 、 $d_0$ が 99 なら、ピッチは 2 倍速。
- $d_1$ が 100 、 $d_0$ が 100 なら、ピッチは 1 倍速。
- $d_1$ が 100 、 $d_0$ が 101 なら、ピッチは 0 倍速。
- $d_1$ が 100 、 $d_0$ が 102 なら、ピッチは -1 倍速。

以下は記号の定義を示した図です。

<figure>
<img src="./img/delaytime_def.svg" alt="A plot to show the definition of delay times d_0 and d_1." style="padding-bottom: 12px;"/>
</figure>

アンチエイリアシングの実装ではピッチの絶対値をカットオフ周波数として使います。

$$
f_c = \begin{cases}
f_s / 2      & \text{if}\ |p| \leq 1 \\
f_s 2^{-|p|} & \text{otherwise.} \\
\end{cases}
$$

- $f_c$: カットオフ周波数。
- $f_s$: サンプリング周波数。

フィルタ設計時のサンプリング周波数を 1 にしてしまえば $f_s$ の乗算を省略できます。

### 短いディレイ時間への対応
Windowed sinc フィルタの長さを固定すると、安定してノイズの少ない出力が得られますが、フィルタの長さの半分よりも短いディレイ時間が指定できなくなります。

短いディレイ時間を設定したいときは、与えられたディレイ時間に応じてフィルタの長さを短くする仕組みが必要です。このとき、 1 サンプル単位でフィルタの長さを変えることも可能ですが、試した範囲では偶数長と奇数長のフィルタを混ぜるとノイズが増えました。よって、以下で紹介している実装ではフィルタの長さを偶数長に固定しています。

### フィードバックコムフィルタ向きの窓関数
フィードバックコムフィルタへの応用では振幅特性が 0 dB を超えていると、フィードバックゲインによっては発振します。高速に計算できる窓関数の向き、不向きを以下に示します。

- 利用可: 三角, Blackman-Harris, Blackman-Nuttall, Nuttall, Flat top.
- 不向き: 矩形, Blackman, Hamming, Hann, Lanczos, Welch.

三角窓は振幅特性がほぼ確実に 0 dB を下回るので安定を優先するときは選択肢になります。他の利用可の窓は振幅特性の通過域がほぼフラットになり、音への応用であればどれも似たような特性です。

### 実装
以下は高速な windowed sinc フィルタの計算を応用した、アンチエイリアシングされたディレイの実装です。 Blackman-Harris 窓を使っています。 Windowed sinc フィルタの長さは偶数に固定しています。 CPU 負荷を減らしたいときは `maxTap` を適当な偶数に減らしてください。

```c++
template<typename Sample, int maxTap = 256> class DelayBlackmanHarrisBiquadSine {
private:
  static_assert(maxTap > 0 && maxTap % 2 == 0);

  static constexpr int minTimeSample = maxTap / 2 - 1;

  Sample maxTime = 0;
  Sample prevTime = 0;
  int wptr = 0;
  std::vector<Sample> buf{maxTap, Sample(0)};

public:
  void setup(size_t maxTimeSample)
  {
    maxTime = Sample(maxTimeSample);
    buf.resize(std::max(size_t(maxTap), maxTimeSample + maxTap / 2 + 1));
  }

  void reset()
  {
    prevTime = 0;
    wptr = 0;
    std::fill(buf.begin(), buf.end(), Sample(0));
  }

  Sample process(Sample input, Sample timeInSample)
  {
    const int size = int(buf.size());

    // Write to buffer.
    if (++wptr >= size) wptr = 0;
    buf[wptr] = input;

    // Shorten FIR filter to some even number length depending on `timeInSample`.
    const int localTap = std::clamp(2 * int(timeInSample), int(2), maxTap);
    const int halfTap = localTap / 2;
    const Sample clamped = std::clamp(timeInSample, Sample(halfTap - 1), maxTime);

    // Set cutoff frequency.
    const Sample timeDiff = std::abs(prevTime - clamped + Sample(1));
    prevTime = clamped;
    Sample cutoff = timeDiff <= Sample(1) ? Sample(0.5) : std::exp2(-timeDiff);

    // Early exit for bypass (0 sample delay) case.
    if (timeInSample <= 0) return input * Sample(2) * cutoff;

    // Setup oscillator 1 for windowed sinc lowpass.
    const int timeInt = int(clamped);
    Sample fraction = clamped - Sample(timeInt);
    const Sample mid = fraction - halfTap;

    constexpr Sample pi = std::numbers::pi_v<Sample>;
    const Sample o1_omega = Sample(2) * pi * cutoff;
    const Sample o1_phi = mid * o1_omega;
    const Sample o1_k = Sample(2) * std::cos(o1_omega);
    Sample o1_u1 = std::sin(o1_phi - o1_omega);
    Sample o1_u2 = std::sin(o1_phi - Sample(2) * o1_omega);

    // Setup oscillator 2 for cosine-sum window function.
    const Sample o2_omega = Sample(2) * pi / Sample(localTap + 1);
    const Sample o2_phi = pi / Sample(2);
    const Sample o2_k = Sample(2) * std::cos(o2_omega);
    Sample o2_u1 = Sample(1); // == sin(o2_phi).
    Sample o2_u2 = std::sin(o2_phi - o2_omega);

    // The rest is convolution.
    int rptr = wptr - timeInt - halfTap;
    if (rptr < 0) rptr += size;

    Sample sum = 0;
    for (int i = 0; i < localTap; ++i) {
      const Sample o1_u0 = o1_k * o1_u1 - o1_u2;
      o1_u2 = o1_u1;
      o1_u1 = o1_u0;

      const Sample o2_u0 = o2_k * o2_u1 - o2_u2;
      o2_u2 = o2_u1;
      o2_u1 = o2_u0;

      const Sample window = Sample(0.21747)
        + o2_u0
          * (Sample(-0.45325) + o2_u0 * (Sample(0.28256) + o2_u0 * Sample(-0.04672)));

      const Sample x = Sample(i) + mid;
      Sample sinc;
      if (std::abs(x) <= Sample(0.1)) {
        // Pade approximant of modified sinc when x is near 0.
        Sample q = pi * cutoff * x;
        q *= q;
        sinc = Sample(2) / Sample(3) * cutoff * (Sample(15) - Sample(7) * q)
          / (Sample(5) + q);
      } else {
        sinc = o1_u0 / pi / x;
      }

      sum += window * sinc * buf[rptr];
      if (++rptr >= size) rptr = 0;
    }
    return sum;
  }
};
```

以下は簡単な使用例です。

```c++
auto processDelay(const std::vector<double> &in, const std::vector<double> &timeInSample) {
  DelayBlackmanHarrisBiquadSine<double> delay;
  delay.setup(48000);

  std::vector<double> out{in.size(), double(0)};
  for (size_t i = 0; i < in.size(); ++i) out[i] = delay.process(in[i], timeInSample[i]);
  return out;
}
```

以下は他の補間方法とエイリアシングノイズを比較したスペクトログラムのプロットです。小さいときは Ctrl + マウスホイールなどで拡大してください。

<figure>
<img src="./img/aa_comparison_on_delaytime_mod.svg" alt="Comparison of aliasing noise produced on delay time change. 6 different interpolation types are tested, that are no interpolation (integer), linear interpolation, order 3 Lagrange interpolation, rectangular window sinc interpolation, triangle window sinc interpolation, and Blackman-Harris window sinc interpolation." style="padding-bottom: 12px;"/>
</figure>

入力信号は 4000 Hz ののこぎり波で、ディレイ時間は sin で変調しています。大まかには縦方向に沿って明るい部分が多いほど、エイリアシングノイズが多いと言えます。

左上が補間なし、中央上が線形補間、右上が 3 次ラグランジュ補間で、この 3 つはエイリアシングが出ていることがはっきりとわかります。下の 3 つは高速な windowed sinc フィルタの計算を応用した実装の結果で、左下は矩形窓、中央下は三角窓、右下は Blackman-Harris 窓 (上記の C++ 実装) です。矩形窓は横線がうねっている部分でのエイリアシングがぼんやりと見えますが、他の 2 つの窓とは異なり、縦の明るい線が見られません。縦の明るい線はフィルタの長さが切り替わるときに生じるポップノイズです。したがって、常にディレイ時間が変わり続けるときは矩形窓にも利点があると考えられますが、前述のようにフィードバックコムフィルタへの応用では通過域のリップルによるフィードバックゲインの制限がかかります。

以下のコードを実行した結果からプロットしました。

- [アンチエイリアシングを施したディレイの C++ 実装 (github.com)](https://github.com/ryukau/filter_notes/blob/master/fast_windowed_sinc/delay/delay.cpp)

### フィルタの長さが変わるときのノイズの低減
以下のようにオシレータの初期設定を変更するとフィルタの長さが変わるときのノイズを減らせます。

```c++
// Setup oscillator 2 for cosine-sum window function.
const Sample o2_omega = Sample(2) * pi / Sample(maxTap + 1);
const Sample o2_phi = pi / Sample(2) + o2_omega * Sample(maxTap / 2 - halfTap);
const Sample o2_k = Sample(2) * std::cos(o2_omega);
Sample o2_u1 = std::sin(o2_phi);
Sample o2_u2 = std::sin(o2_phi - o2_omega);
```

この変更を行うとフィルタの長さが変わるときに窓関数の中央付近だけを切り出すようになります。以下は窓関数の切り詰めを示した図です。

<figure>
<img src="./img/truncation_of_window.svg" alt="Plot of a truncation of Blackman-Harris window for smoother change of FIR length." style="padding-bottom: 12px;"/>
</figure>

以下はフィルタの長さが変わるときのノイズの低減の有無によるエイリアシングの違いを示したスペクトログラムです。上が低減なし、下が低減ありです。縦線が減っているのでポップノイズは減っていますが、エイリアシングは増えています。

<figure>
<img src="./img/blackmanharris_full_vs_smooth.svg" alt="2 spectrograms to show the effect of the full window vs truncated window when reducing the FIR length." style="padding-bottom: 12px;"/>
</figure>

## 参考サイト
- [Window function - Wikipedia](https://en.wikipedia.org/wiki/Window_function)
- [fft - When to use symmetric vs asymmetric (periodic) window functions? - Signal Processing Stack Exchange](https://dsp.stackexchange.com/questions/95448/when-to-use-symmetric-vs-asymmetric-periodic-window-functions)

## 変更点
- 2025/11/22
  - 「フィルタの長さが変わるときのノイズの低減」の節を追加。
