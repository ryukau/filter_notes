# n dB/octのスロープを持つフィルタ
n dB/oct のスロープを持つフィルタを設計します。

<figure>
<img src="img/fir_filter_n_dB_oct.png" alt="Image of the filter with n dB/oct slope." style="width: 300px; padding-bottom: 12px;"/>
</figure>

$$
\begin{aligned}
A(\omega) &=
\begin{cases}
(h/l)^{M}, &  (-\pi \leq \omega \leq -h) \\
(\omega/l)^{M}, &  (-h \lt \omega \lt -l) \\
1, & (-l \leq \omega \leq l) \\
(\omega/l)^{M}, & (l \lt \omega \lt h) \\
(h/l)^{M}, & (h \leq \omega \leq \pi)
\end{cases} \\\\
M &= {{\,m}\over{20}}\log_2 10
\end{aligned}
$$

## 逆離散フーリエ変換を使う方法
[逆離散フーリエ変換（IDFT）](https://en.wikipedia.org/wiki/Discrete_Fourier_transform)を使って楽に設計します。

コードを順番にPython3のインタープリタにコピペしていけば動作します。別ページでまとめて見ることもできます。

- [n dB/octのフィルタを作るコード (github.com)](https://github.com/ryukau/filter_notes/blob/master/ndboct/ndboct.py)

[SciPy](https://www.scipy.org/)、[NumPy](http://www.numpy.org/)、[Matplotlib](https://matplotlib.org/tutorials/introductory/pyplot.html)を import します。

```python
import matplotlib.pyplot as p
import numpy as np
from scipy import signal
from scipy.fftpack import *
```

プロットに使う関数を用意します。

```python
def showResponse(sample_rate, sig):
    # sig の周波数成分を dB に変換してプロット。
    res = fftshift(fft(sig, 2**16) / (len(sig) / 2.0))
    gain = np.array_split(
        20 * np.log10(abs(res / np.max(abs(res)))),
        2,
    )[0][::-1]
    freq = getFreq(sample_rate, gain)
    p.plot(freq, gain)
    p.grid(which='both')
    p.show()

def getFreq(sample_rate, spec):
    # 出力はプロットのx軸に使う。
    freq_step = sample_rate / 2 / (len(spec) - 1)
    return [i * freq_step for i in range(0, len(spec))]

def toDecibel(spec):
    return [20 * np.log10(v) for v in spec]
```

欲しい周波数特性を用意します。

```python
def makeFrequencySpecification(sample_rate, length, low, high, slope):
    # slope は負の値。単位は[dB/oct]。
    # low, high は周波数[Hz]。
    #
    # dB/oct を算出する式の解説。
    # https://femci.gsfc.nasa.gov/random/randomequations.html
    #
    # 6.020599913279624 = 20 * np.log10(2)
    #
    oct = (high / low)**(slope / 6.020599913279624)
    spec = [0] * length
    freq_step = sample_rate / 2 / (length - 1)
    for i in range(0, len(spec)):
        freq = i * freq_step
        if freq <= low:
            spec[i] = 1
        elif freq >= high:
            spec[i] = oct
        else:
            spec[i] = (freq / low)**(slope / 6.020599913279624)
    return spec


sample_rate = 44100
spec = makeFrequencySpecification(sample_rate, 1024, 100, 10000, -10)

# 設計した特性をプロット。
spec_db = toDecibel(spec)
freq = getFreq(sample_rate, spec)

p.title('Specification')
p.plot(freq, spec_db)
p.xscale('log')
p.grid(which='both')
p.show()
```

ここでは100Hzから10000Hz間において、-10dB/octの傾きで減衰する、窓長1024の特性を作りました。

<figure>
<img src="img/n_dB_oct_spec.png" alt="Image of filter specification." style="width: 480px;"/>
</figure>

フィルタ係数を出します。

```python
filt = fftshift(ifft(spec + spec[::-1]))

# フィルタ係数の見た目。
p.title('Filter Coefficients')
p.plot(filt)
p.show()

# フィルタの周波数応答。
p.title('Frequency Response of the filter')
p.xscale('log')
showResponse(sample_rate, filt)
```

フィルタ係数の見た目です。

<figure>
<img src="img/n_dB_oct_cofficients.png" alt="Image of the filter coefficients." style="width: 480px;"/>
</figure>

フィルタの周波数応答です。作った特性と一致します。

<figure>
<img src="img/n_dB_oct_response.png" alt="Image of frequency response of the filter." style="width: 480px;"/>
</figure>

不連続点を拡大すると[ギブス現象](https://en.wikipedia.org/wiki/Gibbs_phenomenon)が観察できます。

```python
showResponse(sample_rate, filt)
```

<figure>
<img src="img/n_dB_oct_response_gibbs.png" alt="Image of gibbs phenomenon." style="width: 593px;"/>
</figure>


図の左は周波数軸が対数でないときのフィルタの周波数応答です。右上は低域側の不連続点、右下は高域側の不連続点を拡大したものです。

フィルタができたので、ノイズを作ってかけてみます。

```python
# フィルタをノイズに適用。
noise = np.random.uniform(-1.0, 1.0, sample_rate)
p.title('Before')
p.xscale('log')
showResponse(sample_rate, noise)

filtered = signal.convolve(noise, filt, mode='same')
p.title('After')
p.xscale('log')
showResponse(sample_rate, filtered)
```

フィルタをかける前とかけた後のスペクトラムです。

<figure>
<img src="img/n_dB_oct_noise_and_filtered_response.png" alt="Image of frequency response of the filter." style="width: 612px;"/>
</figure>

## 式を途中まで解く
[解析解](https://math.stackexchange.com/questions/935405/what-s-the-difference-between-analytical-and-numerical-approaches-to-problems)について調べます。

問題を再掲します。

$$
\begin{aligned}
A(\omega) &=
\begin{cases}
(h/l)^{M}, &  (-\pi \leq \omega \leq -h) \\
(\omega/l)^{M}, &  (-h \lt \omega \lt -l) \\
1, & (-l \leq \omega \leq l) \\
(\omega/l)^{M}, & (l \lt \omega \lt h) \\
(h/l)^{M}, & (h \leq \omega \leq \pi)
\end{cases} \\\\
M &= {{\,m}\over{20}}\log_2 10
\end{aligned}
$$

$A(\omega)$ を[逆DTFT](https://en.wikipedia.org/wiki/Discrete-time_Fourier_transform#Inverse_transform)の式に入れてMaximaで解きます。

```maxima
/* Maxima */
expintrep: true;
expintexpand: true;
declare(n, integer, M, real, l, real, h, real);
assume(n > 0, l >= 0, l <= pi, h >= 0, h <= pi, h > l);

oct(m, l, h) := (h / l)**(M);

O: integrate(oct(m, l, h) * exp(%i*omega*n), omega, -pi, -h);
P: integrate(oct(m, l, omega) * exp(%i*omega*n), omega, -h, -l);
Q: integrate(exp(%i*omega*n), omega, -l, l);
R: integrate(oct(m, l, omega) * exp(%i*omega*n), omega, l, h);
S: integrate(oct(m, l, h) * exp(%i*omega*n), omega, h, pi);
answer: expand(demoivre((O + P + Q + R + S) / (2*pi)));

/* where */
M: log(10)*m/20/log(2);
```

$$
\begin{aligned}
\frac{1}{2\pi}\int^{\pi}_{-\pi} A(\omega) e^{j\omega n} d\omega
=
&{{h^{M}}\over{l^{M}}}\int_{-\pi}^{-h}{e^{i\,n\,\omega}\;d\omega}\\
&+{{1}\over{l^{M}}}\int_{-h}^{-l}{\omega^{M}\,e^{i\,n\,\omega}\;d\omega}\\
&+\int_{-l}^{l}{e^{i\,n\,\omega}\;d\omega}\\
&+{{1}\over{l^{M}}}\int_{l}^{h}{\omega^{M}\,e^{i\,n\,\omega}\;d\omega}\\
&+{{h^{M}}\over{l^{M}}}\int_{h}^{\pi}{e^{i\,n\,\omega}\;d\omega}\\
=
&{{\sin \left(l\,n\right)}\over{n\,\pi}}
+{{h^{M}\,\sin \left(n\,\pi\right)}\over{l^{M}\,n\,\pi}}
-{{h^{M}\,\sin \left(h\,n\right)}\over{l^{M}\,n\,\pi}}
\\
&+{{i\,\left(-1\right)^{{{M}\over{2}}} n^{-M-1}}\over{2\,\pi\,l^{M}}}
\Bigl(\Gamma\left(M+1 , i\,h\,n\right) - \Gamma\left(M+1 , i\,l\,n\right)\Bigr)
\\
&+{{i\,\left(-i\right)^{-M}\,n^{-M-1}}\over{2\,\pi\,l^{M}}}
\Bigl(\Gamma\left(M+1 , -i\,l\,n\right) - \Gamma\left(M+1 , -i\,h\,n\right)\Bigr)
\end{aligned}
$$

$\Gamma(s, x)$ は [upper incomplete gamma function](https://en.wikipedia.org/wiki/Incomplete_gamma_function#Definition) です。

$$
\Gamma(s, x) = \int_{x}^{\infty} t^{s - 1} e^{-t} dt
$$

$\Gamma(s, x)$ がどこから出てきたのかを確認するために、Maximaのコードで定義したP式とR式に注目します。

```maxima
oct(m, l, h) := (h / l)**(M);
assume(n > 0, l >= 0, l <= pi, h >= 0, h <= pi, h > l);
P: integrate(oct(m, l, omega) * exp(%i*omega*n), omega, -h, -l);
R: integrate(oct(m, l, omega) * exp(%i*omega*n), omega, l, h);
```

$$
\begin{aligned}
\int^{\pi}_{-\pi} A_P(\omega) e^{j\omega n} d\omega
=
&{{i\,\left(-1\right)^{{{M}\over{2}}} n^{-M-1}}\over{2\,\pi\,l^{M}}}
\Bigl(\Gamma\left(M+1 , i\,h\,n\right) - \Gamma\left(M+1 , i\,l\,n\right)\Bigr)
\\
\int^{\pi}_{-\pi} A_R(\omega) e^{j\omega n} d\omega
=
&{{i\,\left(-i\right)^{-M}\,n^{-M-1}}\over{2\,\pi\,l^{M}}}
\Bigl(\Gamma\left(M+1 , -i\,l\,n\right) - \Gamma\left(M+1 , -i\,h\,n\right)\Bigr)
\end{aligned}
$$

適当に $\Gamma\left(M+1 , i\,h\,n\right)$ を展開します。

$$
\Gamma\left(M+1 , i\,h\,n\right) = \int_{i\,h\,n}^{\infty} t^{M} e^{-t} dt
$$

[複素数の線積分](https://en.wikipedia.org/wiki/Line_integral#Complex_line_integral)に見えます。

## 参考サイト
- [Maxima 5.41.0 Manual: 11. Maximas Database](http://maxima.sourceforge.net/docs/manual/maxima_11.html#declare)
  - Maxima の declare。
- [Stop maxima from asking positive, negative or zero - Stack Overflow](https://stackoverflow.com/questions/14291119/stop-maxima-from-asking-positive-negative-or-zero)
  - Maxima の assume。
- [Maxima 5.41.0 Manual: 15. Special Functions](http://maxima.sourceforge.net/docs/manual/maxima_15.html)
- [Exponential integral - Wikipedia](https://en.wikipedia.org/wiki/Exponential_integral)
  - Incomplete gamma function と関連。
