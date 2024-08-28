# 適応ノッチフィルタの実装
Ishibashi らによる "[DSP Implementation of Adaptive Notch Filerts With Overflow Avoidance in Fixed-Point Arithmetic](http://www.apsipa.org/proceedings/2018/pdfs/0001355.pdf)" に基づいた適応ノッチフィルタを実装します。

## CPZ-ANF
CPZ-ANF (adaptive notch filter with constrained poles and zeros) は Ishibashi らの論文で紹介されていた適応ノッチフィルタの名前です。

以下は CPZ の伝達関数です。

$$
\begin{align}
H_N(z) &= \frac{1 + a[n] z^{-1} + z^{-2}}{1 + \rho a[n] z^{-1} + \rho^2 z^{-2}}, \\
a[n] &= -2 \cos(\omega_0[n]).
\end{align}
$$

以下はフィルタ係数 $a$ を更新する式です。

$$
\begin{align}
a[n + 1] &= a[n] - 2 \mu y[n] s[n],\\
s[n] &= z^{-1} - z^{-1} \rho H_N(z). \label{cpz-anf_s}
\end{align}
$$

記号の一覧です。

- $n$ : サンプル数で表された時間。
- $x$ : 入力信号。
- $y$ : 出力信号。
- $\omega_0$ : ノッチの角周波数。
- $\rho$ : ノッチの幅。範囲 $[0, 1)$ 。 1 に近づくほどノッチが狭くなる。
- $\mu$ : 更新の速さ。 $1/f_s$ あたりの正の実数。あまり大きくすると発散する。

Python 3 で実装します。 [Direct-form II](https://ccrma.stanford.edu/~jos/fp/Direct_Form_II.html) での実装が指定されています。

```python
def adaptiveNotchCpz2(x, rho=0.99, mu=1, initialGuess=0.5):
    out = np.zeros_like(x)
    a = -2 * np.cos(2 * np.pi * initialGuess)
    v1 = 0
    v2 = 0
    for i in range(len(x)):
        a1 = rho * a
        a2 = rho * rho

        x0 = x[i]
        v0 = x0 - a1 * v1 - a2 * v2
        y0 = v0 + a * v1 + v2
        out[i] = y0  # ゲインは 0 dB を超えることがあるので注意。

        s0 = (1 - rho) * v0 - rho * (1 - rho) * v2
        a = np.clip(a - 2 * mu * y0 * s0, -2, 2)

        v2 = v1
        v1 = v0
    return out
```

### ゲインの正規化
ノッチフィルタ $H_N$ は $a \geq 0$ のときに周波数 0 のゲインが最大、 $a < 0$ のときにナイキスト周波数のゲインが最大となります。したがって $a \geq 0$ のときに $|H_N(e^{j0})|$ 、$a < 0$ のときに $|H_N(e^{j\pi})|$ を計算すればゲインの最大値が得られます。

$$
\begin{aligned}
\text{Gain at frequency 0 is:} &&
H_N(e^{j0}) &= H_N(1) = \frac{1 + a[n] + 1}{1 + \rho a[n] + \rho^2},
\\
\text{Gain at Nyquist frequency is:} &&
H_N(e^{j\pi}) &= H_N(-1) = \frac{1 - a[n] + 1}{1 - \rho a[n] + \rho^2},
\\
\end{aligned}
$$

ゲインの最大値を 0 dB に正規化するときは、上のコードの `out[i] = y0` の部分を以下に変更します。

```python
if a >= 0:
    out[i] = y0 * (1 + a1 + a2) / (2 + a)
else:
    out[i] = y0 * (1 - a1 + a2) / (2 - a)
```

### 分析
以下は `y0` と `s0` の位相特性です。

<figure>
<img src="img/cpz-anf_phase_diff.svg" alt="Plot of phase difference of output y and gradient s in CPZ-ANF." style="padding-bottom: 12px;"/>
</figure>

$-\pi/2$ と $-3\pi/2$ で横に黒い横線を超えると、入力がサイン波のときに `y0 * s0` の時間平均をとったときの直流成分の符号が変わります。つまり、カットオフ周波数の上下で `y0 * s0` の符号が変わるので、 `a` がターゲット周波数に向かって適応します。

より詳しい分析を [AM 適応ノッチフィルタの分析の節](#分析-1)に掲載しています。 $a$ の計算が異なりますが、位相差についての挙動は同じです。

## AM 適応ノッチフィルタ
以下のように `s0` を `x0` に置き変えても、なぜか適応ノッチフィルタとして動きます。

```python
def adaptiveNotchAM(x, rho=0.99, mu=1, initialGuess=0.5):
    out = np.zeros_like(x)
    a = -2 * np.cos(2 * np.pi * initialGuess)
    x1 = 0
    x2 = 0
    y1 = 0
    y2 = 0
    for i in range(len(x)):
        x0 = x[i]
        y0 = x0 + a * x1 + x2 - rho * a * y1 - rho * rho * y2
        out[i] = y0

        a = np.clip(a - 2 * mu * y0 * x0, -2, 2)

        x2 = x1
        x1 = x0
        y2 = y1
        y1 = y0
    return out
```

上の実装は [Direct-form I](https://ccrma.stanford.edu/~jos/fp/Direct_Form_I.html) の形を使っています。 `a` の計算で `y0 * x0` という振幅変調 (AM) の形が出てきているので AM 適応ノッチフィルタと呼んでいます。

### 分析
サイン波が AM 適応ノッチフィルタに入力されたときの挙動について分析します。

式をたてます。 $x$ と $y$ は $\sin$ でもいいのですが、符号をそろえて楽をするために $\cos$ を使っています。

$$
\begin{aligned}
x(t) &= \cos(\omega t). && \text{Input signal.} \\
y(t) &= A \cos(\omega t + \phi). && \text{Output signal.} \\
s(t)
&= x(t) y(t) = \frac{A}{2} \Big(
    \cos(2 \omega t + \phi)
  + \cos(-\phi)
\Big). && \text{Velocity of filter coefficient}\ a. \\
a(t) &= -2 \cos(\omega_a(t)) = a(t - \delta_t) - 2 \mu s(t),\\
\omega_a &\in [0, \pi).
\end{aligned}
$$

$\phi$ は周波数 $\omega$ でのノッチフィルタの出力の位相差 $\angle H_N(e^{j\omega})$ 、 $A$ は周波数 $\omega$ でのノッチフィルタのゲインです。

$$
\begin{aligned}
\phi &= \angle H_N(e^{j\omega}) \\
A &= |H_N(e^{j\omega})| \\
\end{aligned}
$$

以下はノッチフィルタの位相特性です。

<figure>
<img src="img/notch_phase_response.svg" alt="Plot of phase response of notch filter." style="padding-bottom: 12px;"/>
</figure>

ここで以下の分析ができます。 $\langle s \rangle$ は $s$ の[時間方向の算術平均](https://math.stackexchange.com/questions/535700/bar-mean-vs-bracket-mean)です。

- $a$ はカットオフ周波数が低いと -2 、高いと +2 へと近づく。
- 適応が完了したとき、ノッチによってサイン波が消えるので `y0 = 0` となり、 `y0 * x1` も 0 。
- 適応が完了していないとき：
  - ターゲットよりもカットオフが低いと $y$ の位相が進む。
  - ターゲットよりもカットオフが高いと $y$ の位相が遅れる。
  - $s$ は式に含まれる $- \cos(-\phi)$ より、
    - 位相のずれが 0 なら常に正の値になる。 $s \geq 0$ 。
    - 位相のずれが $\pi$ なら常に負の値になる。 $s \leq 0$ 。
    - 位相のずれが $0 > \phi > -\pi/2$ あるいは $-3\pi/2 > \phi > -2\pi$ の範囲なら、正の値の直流が乗る。 $\langle s \rangle > 0$ 。
    - 位相のずれが $-\pi/2 > \phi > -3\pi/2$ の範囲なら、負の値の直流が乗る。 $\langle s \rangle < 0$ 。

この分析を踏まえてノッチフィルタの位相特性を見直すと、カットオフ周波数の周りでは符号が反転がうまく働いて適応に成功しそうです。しかし、ターゲット周波数が高くなると位相特性が $-3\pi/2$ を下回ってしまうので適応に失敗しそうです。また $\rho$ の値があまりにも 1 に近いとカットオフ周波数の上側で符号が反転する範囲があまりにも狭くなるので、これも失敗しそうです。

ここで肝となるのが、 $s$ が完全に正あるいは負の値となることは稀ということです。つまり $a$ が平均としてはターゲットから離れる方向に進んでいても、瞬時的に正しい方向へと向かう状態が現れます。この挙動によって誤った周波数で適応が止まった状態から抜け出すことがあるかもしれません。

簡単に試した範囲では $\rho$ を 0.8 あたりまで下げれば、ターゲット周波数が高くても適応に成功します。

### 1 次オールパスによる改良
分析によって、ターゲット周波数がノッチのカットオフ周波数のより低いときは $(-\pi/2, \pi/2)$ 、高いときは $(\pi/2, 3\pi/2)$ の位相差を作ってやれば、 AM による $a$ の更新が上手くいきそうなことが分かりました。。

位相差を作るためには $[0, -/pi)$ の位相差を作ることができる 1 次オールパスフィルタが使えます。以下は 1 次オールパスフィルタの伝達関数です。

$$
H_{AP}(z) = \frac{k_{AP} + z^{-1}}{1 + k_{AP} z^{-1}}.
$$

以下は $k_{AP}$ を変えたときの位相特性です。

<figure>
<img src="img/allpass1_phase_response.svg" alt="Plot of phase response of first order allpass filter." style="padding-bottom: 12px;"/>
</figure>

$k_{AP}$ は、ノッチフィルタの $a$ から決めます。 $a$ からノッチのカットオフ周波数 $\omega_a$ を計算して、 $\omega_a$ をブレーク周波数 (break frequency) として $k_{AP}$ を計算します。

$$
\begin{aligned}
\omega_a &= \pi - \arccos(a / 2), \\
k_{AP} &= \frac{\tan(\omega_a / 2) - 1}{\tan(\omega_a / 2) + 1}. \\
\end{aligned}
$$

以下はブロック線図です。推定値に近くなったときに修正幅を狭めるためノッチの出力の絶対値を、入力とオールパス出力の AM 信号にさらに乗算しています。青の線で示された $\dot{a}$ はノッチの係数 $a$ の 1 サンプル当たりの変化量です。 $a$ を $k_{AP}$ に変換する処理をオレンジで示しています。

<figure>
<img src="img/adaptive_notch_improved_am_block_diagram.svg" alt="Block diagram of improved AM adaptive notch filter." style="padding-bottom: 12px;"/>
</figure>

以下は実装です。

```python
def adaptiveNotchAM2(x, rho=0.99, mu=1, initialGuess=0.5):
    out = np.zeros_like(x)
    a = -2 * np.cos(2 * np.pi * initialGuess)
    x1 = 0
    x2 = 0
    y1 = 0
    y2 = 0
    q1 = 0  # オールパスの 1 サンプル前の入力。
    r1 = 0  # オールパスの 1 サンプル前の出力。
    for i in range(len(x)):
        x0 = x[i]
        y0 = x0 + a * x1 + x2 - rho * a * y1 - rho * rho * y2
        out[i] = y0

        notchFreq = (np.pi - np.arccos(a / 2)) / (2 * np.pi)
        k_ap = breakFreqToA(notchFreq)

        r1 = k_ap * (x0 - r1) + q1
        q1 = x0
        a = np.clip(a - mu * np.abs(y0) * x0 * r1, -2, 2)

        x2 = x1
        x1 = x0
        y2 = y1
        y1 = y0
    return out
```

$\rho$ が 0.99 あたりで高い周波数への適応がいいです。出力が完全に 0 に収束しづらいという欠点があります。特に低い周波数では収束が遅くなります。また $\rho$ を 0 に近づけるほど収束し損ねた信号の振幅が大きくなります。

収束しづらいのはカットオフ周波数付近では位相差 $\phi$ が $-\pi/2$ に近づくため、直流成分 $\cos(-\phi)$ が小さくなるからだと考えられます。したがってターゲット周波数とカットオフ周波数の距離が遠いほど収束が速く、周波数の差が縮まるにつれて収束が遅くなります。

### 異なるノッチフィルタへの適用
AM による周波数の適応は $H_N$ とは異なるノッチフィルタとも組み合わせられます。計算量の面からみると利点は薄いですが、 [RBJ biquad](https://webaudio.github.io/Audio-EQ-Cookbook/audio-eq-cookbook.html) のノッチと組み合わせた実装を以下に掲載しています。

- [RBJ biquad と 1 次オールパスによる適応ノッチフィルタの実装 (github.com)](https://github.com/ryukau/filter_notes/blob/194d16d4568d05458e1df32adfc1c452227dc2fc/adaptive_notch/notch.py#L335-L400)

## 発散の抑制
$H_N$ を使った適応ノッチフィルタは、フィルタ係数 $a$ の挙動によって発散することがあります。

まずは $a$ の値が $[-2, 2]$ を超えないように制限します。以下のコード例の [`np.clip`](https://numpy.org/doc/stable/reference/generated/numpy.clip.html) は C++ などの [`clamp`](https://en.cppreference.com/w/cpp/algorithm/clamp) と同じ計算です。 `delta_a` は 1 サンプル当たりの $a$ の変動値です。

```python
a = np.clip(a - mu * delta_a, -2, 2)
```

$a$ が数サンプル以内で大きく動くと発散することがあります。 $a$ の動きを緩やかにする最も手軽な方法は $\mu$ の値を 0 に近づけることです。ただし副作用として収束が遅くなります。別の方法は $a$ の変動値にスルーレートリミッタあるいは EMA ローパスをかけることです。

スルーレートリミッタは `clamp` などで値の範囲を制限するだけで実装できます。ただし収束しきっていない波形が歪むことがあります。

```python
a = np.clip(a - mu * np.clip(delta_a, -1e-2, 1e-2), -2, 2)
```

EMA ローパスの実装には状態変数とフィルタ係数が必要です。状態変数はループの外で `w1 = 0` と定義しておきます。フィルタ係数は個別に設定できますが、 $\mu$ を流用すると楽です。

```python
w1 += mu * (delta_a - w1)
a = np.clip(a - mu * w1, -2, 2)
```

念を入れるのであれば、出力を `isfinite` に渡して、値が有限でなければフィルタの状態をリセットすると安全です。 C++ ではコンパイラの `-ffinite-math-only` や `/fp:fast` が有効だと [`std::isfinite`](https://en.cppreference.com/w/cpp/numeric/math/isfinite) は常に `true` を返すので注意してください。 `-ffinite-math-only` は `-ffast-math` によって有効になります。

- [Optimize Options (Using the GNU Compiler Collection (GCC))](https://gcc.gnu.org/onlinedocs/gcc/Optimize-Options.html#index-ffast-math)

## 参考文献
- Ishibashi, Satoru, Shunsuke Koshita, Masahide Abe, and Masayuki Kawamata. "[DSP implementation of adaptive notch filters with overflow avoidance in fixed-point arithmetic](http://www.apsipa.org/proceedings/2018/pdfs/0001355.pdf)." In 2018 Asia-Pacific Signal and Information Processing Association Annual Summit and Conference (APSIPA ASC), pp. 1355-1360. IEEE, 2018.
- [Constant Peak-Gain Resonator](https://ccrma.stanford.edu/~jos/filters/Constant_Peak_Gain_Resonator.html)
- [Complex Resonator](https://ccrma.stanford.edu/~jos/filters/Complex_Resonator.html)
