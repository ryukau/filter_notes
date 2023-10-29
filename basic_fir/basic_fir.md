# お手軽なFIRフィルタのレシピ
**注意**: 内容が怪しいです。

[FIR](https://en.wikipedia.org/wiki/Finite_impulse_response)のローパス、ハイパス、バンドパス、バンドリジェクトフィルタを作ります。

## 記号

- $f_{s}$ : [サンプリング周波数](https://en.wikipedia.org/wiki/Sampling_(signal_processing)#Sampling_rate)
- $f_l$ : 低いほうの[カットオフ周波数](https://en.wikipedia.org/wiki/Cutoff_frequency)
- $f_h$ : 高いほうのカットオフ周波数
- $\omega_l$ : $2 \pi f_l$
- $\omega_h$ : $2 \pi f_h$

[通過域](https://en.wikipedia.org/wiki/Passband)の大きさ（[Amplitude](https://en.wikipedia.org/wiki/Amplitude)）は1とします。

## ローパスフィルタ
<figure>
<img src="img/fir_filter_lowpass.png" alt="Image of FIR lowpass filter." style="width: 320px;padding-bottom: 12px;"/>
</figure>

$$
A(\omega) =
\begin{cases}
1, & (-\omega_l \leq \omega \leq \omega_l) \\
0, & \text{otherwise}
\end{cases}
$$

$A(\omega)$ を[逆離散時間フーリエ変換](https://en.wikipedia.org/wiki/Discrete-time_Fourier_transform)の式に入れて解きます。

$$
\begin{aligned}
\frac{1}{2\pi}\int^{\pi}_{-\pi} A(\omega) e^{j\omega n} d\omega
&= \frac{1}{2\pi}\int^{\omega_l}_{-\omega_l} e^{j\omega n} d\omega \\
&= \frac{\sin (\omega_l n)}{n \pi}
\end{aligned}
$$

コードに変えます。

```javascript
function makeLowpassCoefficient(length, cutoff) {
  // cutoff の範囲は [0, 1]
  var coefficient = new Array(length).fill(0)
  var half = (length % 2 === 0 ? length - 1 : length) / 2
  var omegaL = 2 * Math.PI * cutoff
  for (var i = 0; i < length; ++i) {
    var n = i - half
    coefficient[i] = (n === 0) ? 1 : Math.sin(omegaL * n) / (Math.PI * n)
  }
  return coefficient
}
```

## ハイパスフィルタ
<figure>
<img src="img/fir_filter_highpass.png" alt="Image of FIR highpass filter." style="width: 320px;padding-bottom: 12px;"/>
</figure>

$$
A(\omega) =
\begin{cases}
1, & (-\omega_s \leq \omega \leq -\omega_h) \\
1, & (\omega_h \leq \omega \leq \omega_s) \\
0, & \text{otherwise}
\end{cases}
$$

$A(\omega)$ を逆離散時間フーリエ変換の式に入れて解きます。

$$
\begin{aligned}
\frac{1}{2\pi}\int^{\pi}_{-\pi} A(\omega) e^{j\omega n} d\omega
&= \frac{1}{2\pi} \biggl ( \int^{-\omega_h}_{-\pi} e^{j\omega n} d\omega
 + \int^{\pi}_{\omega_h} e^{j\omega n} d\omega \biggr ) \\
&= {{\sin \left(n\,\pi\right)}\over{n\,\pi}}-{{\sin \left(\omega_h\,n\right)}\over{n\,\pi}}
\end{aligned}
$$

コードに変えます。

```javascript
function makeHighpassCoefficient(length, cutoff) {
  // cutoff の範囲は [0, 1]
  var coefficient = new Array(length).fill(0)
  var half = (length % 2 === 0 ? length - 1 : length) / 2
  var omegaH = 2 * Math.PI * cutoff
  for (var i = 0; i < length; ++i) {
    var n = i - half
    coefficient[i] = (n === 0)
      ? 1
      : (Math.sin(Math.PI * n) - Math.sin(omegaH * n)) / (Math.PI * n)
  return coefficient
}
```

## バンドパスフィルタ
<figure>
<img src="img/fir_filter_bandpass.png" alt="Image of FIR bandpass filter." style="width: 320px;padding-bottom: 12px;"/>
</figure>

$$
A(\omega) =
\begin{cases}
1, & (-\omega_h \leq \omega \leq -\omega_l) \\
1, & (\omega_l \leq \omega \leq \omega_h) \\
0, & \text{otherwise}
\end{cases}
$$

$A(\omega)$ を逆離散時間フーリエ変換の式に入れて解きます。

$$
\begin{aligned}
\frac{1}{2\pi}\int^{\pi}_{-\pi} A(\omega) e^{j\omega n} d\omega
&= \frac{1}{2\pi} \biggl ( \int^{-\omega_l}_{-\omega_h} e^{j\omega n} d\omega
 + \int^{\omega_h}_{\omega_l} e^{j\omega n} d\omega \biggr ) \\
&= {{\sin \left(\omega_h\,n\right)}\over{n\,\pi}}-{{\sin \left(\omega_l\,n\right)}\over{n\,\pi}}
\end{aligned}
$$

コードに変えます。

```javascript
function makeBandpassCoefficient(length, low, high) {
  // low, high の範囲は [0, 1]
  var coefficient = new Array(length).fill(0)
  var half = (length % 2 === 0 ? length - 1 : length) / 2
  var twoPi = 2 * Math.PI
  var omegaL = twoPi * low
  var omegaH = twoPi * high
  for (var i = 0; i < length; ++i) {
    var n = i - half
    coefficient[i] = (n === 0)
      ? 1
      : (Math.sin(omegaH * n) - Math.sin(omegaL * n)) / (Math.PI * n)
  return coefficient
}
```

## バンドリジェクトフィルタ
バンドストップフィルタとも呼ばれるようです。

<figure>
<img src="img/fir_filter_bandreject.png" alt="Image of FIR bandreject filter." style="width: 320px;padding-bottom: 12px;"/>
</figure>

$$
A(\omega) =
\begin{cases}
1, &  (-\omega_s \leq \omega \leq -\omega_h) \\
1, & (-\omega_l \leq \omega \leq \omega_l) \\
1, & (\omega_h \leq \omega \leq \omega_s) \\
0, & \text{otherwise}
\end{cases}
$$

$A(\omega)$ を逆離散時間フーリエ変換の式に入れて解きます。

$$
\begin{aligned}
\frac{1}{2\pi}\int^{\pi}_{-\pi} A(\omega) e^{j\omega n} d\omega
&= \frac{1}{2\pi} \biggl ( \int^{-\omega_h}_{-\pi} e^{j\omega n} d\omega
 + \int^{\omega_l}_{-\omega_l} e^{j\omega n} d\omega
 + \int^{\pi}_{\omega_h} e^{j\omega n} d\omega \biggr ) \\
&= {{\sin \left(n\,\pi\right)}\over{n\,\pi}}+{{\sin \left(\omega_l\,n\right)}\over{n\,\pi}}-{{\sin \left(\omega_h\,n\right)}\over{n\,\pi}}
\end{aligned}
$$

コードに変えます。

```javascript
function makeBandrejectCoefficient(length, low, high) {
  // low, high の範囲は [0, 1]
  var coefficient = new Array(length).fill(0)
  var half = (length % 2 === 0 ? length - 1 : length) / 2
  var twoPi = 2 * Math.PI
  var omegaL = twoPi * low
  var omegaH = twoPi * high
  for (var i = 0; i < length; ++i) {
    var n = i - half
    var piN = Math.PI * n
    coefficient[i] = (n === 0)
      ? 1
      : (Math.sin(piN) + Math.sin(omegaL * n) - Math.sin(omegaH * n)) / piN
  return coefficient
}
```

## FIRフィルタのかけ方
まずフィルタをかけるソースを用意します。

```javascript
var audioContext = new AudioContext()
var source = new Array(audioContext.sampleRate)
for (var i = 0; i < source.length; ++i) {
  source[i] = 2.0 * Math.random() - 1.0 // 適当なノイズ。
}
```

FIRフィルタは適切な[窓関数](https://en.wikipedia.org/wiki/Window_function)をかけることで特性が改善します。ここではお手軽でそれなりに特性がいい[Blackman–Harris窓](https://en.wikipedia.org/wiki/Window_function#Blackman%E2%80%93Harris_window)を使います。

```javascript
function blackmanHarrisWindow(length) {
  var window = new Array(length)
  var a0 = 0.35875
  var a1 = 0.48829
  var a2 = 0.14128
  var a3 = 0.01168
  var pi_N1 = Math.PI / (window.length - 1)
  var twopi_N1 = 2 * pi_N1
  var fourpi_N1 = 4 * pi_N1
  var sixpi_N1 = 6 * pi_N1
  for (var n = 0; n < window.length; ++n) {
    window[n] = a0
      - a1 * Math.cos(n * twopi_N1)
      + a2 * Math.cos(n * fourpi_N1)
      - a3 * Math.cos(n * sixpi_N1)
  }
  return window
}
```

フィルタを用意します。

```javascript
var filterLength = 1025
var cutoff = 1000 // Hz
var coefficient = makeLowpassCoefficient(filterLength, cutoff / audioContext.sampleRate)
var filter = blackmanHarrisWindow(filterLength).map((v, i) => v * coefficient[i])
```

[畳み込み（Convolution）](https://en.wikipedia.org/wiki/Convolution)を行ってソースにFIRフィルタをかけます。

```javascript
var destination = new Array(source.length).fill(0)
var buffer = new Array(filter.length).fill(0)
for (var i = 0; i < source.length; ++i) {
  buffer.push(source[i])
  buffer.shift()
  for (var j = 0; j < filter.length; ++j) {
    destination[i] += buffer[j] * filter[j]
  }
}
```

## Computer Algebra System の利用
手で式を解くと間違えることがあるので [Computer Algebra System (CAS)](https://en.wikipedia.org/wiki/Computer_algebra_system) を利用します。

今回のような簡単な式であれば[Wolfram Alpha](https://www.wolframalpha.com/)が便利です。Wolfram Alphaでは数字でない下付き文字が使えないようなので $\omega_l$ を $l$ に置き換えています。以降のCASのコードも同じ置き換えを使います。

```wolfram
(integral e^(i*omega*n) for omega from -l to l) / (2pi)
```

[Maxima](http://maxima.sourceforge.net/)は式の整理について指定する必要があります。コードの `demoivre` で `exp(%i*n)` を `%i * sin(n) + cos(n)` に置き換えています。

```maxima
expand(demoivre(integrate(exp(%i * omega * n) / (2 * pi), omega, -l, l)));
```

[SymPy](http://www.sympy.org/en/index.html)も使えますが少し長めです。 `rewrite(sin)` でオイラーの公式を適用しています。 `n = Symbol('n', positive=True)` が無いとコードの `integrate(...)` を解いてくれません。

```python
# SymPy 1.1.1
from sympy import *
n = Symbol('n', positive=True)
l = Symbol('l')
omega = Symbol('omega')
answer = simplify(integrate(exp(I * omega * n), (omega, -l, l)).rewrite(sin))
pprint(answer)
```

## フィルタ係数の計算について
フィルタ係数は周波数特性を逆離散時間フーリエ変換することで得られます。

$$
\frac{1}{2\pi}\int^{\pi}_{-\pi} A(\omega) e^{j\omega n} d\omega
$$

この逆離散時間フーリエ変換の式に、周波数$f$を対応させて使うときは$\pi{f}/{f_s}$と変換します。

逆離散時間フーリエ変換をするときは$[-f_s,\ f_s]$の範囲を考慮する必要がありますが、ここまでに出てきた周波数特性の図では$[-f_s,\ 0]$の範囲を省略していました。$[-f_s,\ 0]$の範囲での周波数特性は$[0,\ f_s]$の鏡像になっています。

ローパスフィルタを例に見ていきます。

<figure>
<img src="img/fir_filter_frequency_mirror.png" alt="Image of range [-f_s, f_s] of lowpass frequency responce." style="width: 600px;padding-bottom: 12px;"/>
</figure>

[偶関数の積分の性質](https://en.wikipedia.org/wiki/Even_and_odd_functions#Calculus_properties)を利用して式を変形できそうなので試してみます。

$$
\begin{aligned}
\frac{1}{2\pi}\int^{\pi}_{-\pi} A(\omega) e^{j\omega n} d\omega
&= \frac{1}{2\pi}\int^{\omega_l}_{-\omega_l} e^{j\omega n} d\omega \\
&= \frac{2}{2\pi}\int^{\omega_l}_{0} e^{j\omega n} d\omega \\
&= {{\sin \left(\omega_l\,n\right)}\over{n\,\pi}}-{{i\,\cos \left(\omega_l\,n\right)}\over{n\,\pi}}+{{i}\over{n\,\pi}}
\end{aligned}
$$

虚部が出てきました。実部は式を変形をしないときと同じになっています。

念のために定義どおり計算した方がよさそうです。

## その他
式を解かなくても周波数特性を逆離散フーリエ変換すればフィルタは作れます。

フィルタ係数が固定のときは[SciPy](https://www.scipy.org/)や[Octave](https://www.gnu.org/software/octave/)などを使って設計するほうが楽で確実です。

ここで作ったバンドパスフィルタを Banded Waveguides に使おうとしたのですが、切れ味が良すぎてWaveguide間の干渉がほとんど起こらず、面白い音になりませんでした。

## 参考サイト
- [The Ideal Lowpass Filter](https://ccrma.stanford.edu/~jos/sasp/Ideal_Lowpass_Filter.html)
- [Maxima: Expand e to cos and i sin? - Stack Overflow](https://stackoverflow.com/questions/42454464/maxima-expand-e-to-cos-and-i-sin)
