<style>
body {
  max-width: 704px;
  margin: auto;
  padding: 32px 8px;
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

# DR-110風のシンバル
[DR-110](https://www.youtube.com/watch?v=-ap0ucC6QBU)のシンバルに似た音を作ります。

次のリンクから完成したコードを読むことができます。

- [コードを見る (github.com)](https://github.com/ryukau/filter_notes/tree/master/docs/demo/dr110_cymbal)

コードを実行するには[Python3](https://www.python.org/)と以下のライブラリが必要です。

- [SciPy](https://www.scipy.org/)
- [NumPy](http://www.numpy.org/)
- [Matplotlib](https://matplotlib.org/)
- [pyFFTW](https://hgomersall.github.io/pyFFTW/)
- [PySoundFile](http://pysoundfile.readthedocs.io/en/0.9.0/)

## 回路図の入手
[The Boss DR-110 Dr Rhythm Graphic Information Homepage](http://www.theninhotline.net/dr110/) の Schematics から回路図をダウンロードできます。

[richardc64のDR-110のページ](http://www.sdiy.org/richardc64/new_drums/dr110/dr110a1.html)に回路図の解説があります。

## 信号の流れ
DR-110のシンバルはオシレータから出た信号を High Metal と Low Metal の2つに分けて、それぞれでフィルタとエンベロープをかけてから足し合わせることで合成されます。High Metal と Low Metal という呼び方は[richardc64のDR-110のページ](http://www.sdiy.org/richardc64/new_drums/dr110/dr110a1.html)にならっています。

今回実装したプログラムの信号の流れです。

<figure>
<img src="img/dr110_cymbal/DR-110.svg" alt="Image of signal flow of DR-110 cymbal." style="width: 720px; padding-bottom: 12px;"/>
</figure>

## オシレータ
4つの矩形波と1つのノイズを足し合わせます。音量の比率は、矩形波がそれぞれ1、ノイズは0.3です。回路図では矩形波の周波数に317、465、820、1150[Hz]が使われています。

<figure>
<img src="img/dr110_cymbal/DR-110_osc.svg" alt="Image of signal flow of DR-110 cymbal oscillator." style="width: 360px; padding-bottom: 12px;"/>
</figure>

そのまま実装すると味気ない音になるので、オシレータの位相を進める時にノイズを加えて周波数を揺らしています。音量の比率も少しだけランダマイズしています。

```python
import numpy as np
import soundfile
import scipy.signal as signal

def normalize(data, peak=1.0):
    return peak * data / np.abs(np.max(data))


def pulse(time, frequency):
    noise = np.power(10, -(np.random.uniform(0, 5, len(time)) + 5))
    return signal.square(2 * np.pi * frequency * (time + noise))


def noise(size):
    return np.random.uniform(-1, 1, size)


def randomAmp(low, high):
    lowLog = np.log(low)
    diffLog = np.log(high) - lowLog
    return np.exp(lowLog + diffLog * np.random.random())


def cymbalSource(time, pulseFreq=[317, 465, 820, 1150]):
    sig = 0
    for freq in pulseFreq:
        sig += randomAmp(0.999, 1) * pulse(time, freq)
    return sig + noise(len(time)) / 3.3


samplerate = 44100
duration = 2.0

time = np.linspace(0, duration, int(samplerate * duration), endpoint=False)
sig = normalize(cymbalSource(time))

soundfile.write("cymbal_source.wav", sig, samplerate)
```

## フィルタ
[伝達関数](https://ccrma.stanford.edu/~jos/filters/Transfer_Function.html)を使ってDR-110の回路図から離散フィルタを設計します。

### バンドパスフィルタ
オシレータセクションで使われているバンドパスフィルタです。BP0とBP1という名前はこの文章で区別するためにつけました。

<figure>
<img src="img/dr110_cymbal/osc_crop.png" alt="Image of oscillator section of DR-110 schematic." style="width: 480px; padding-bottom: 12px;"/>
</figure>

単位の無い[キャパシタ](https://en.wikipedia.org/wiki/Capacitor)の値は、上位2桁が有効数字、下位1桁が10の指数で、単位は[[p](https://en.wikipedia.org/wiki/Pico-)[F](https://en.wikipedia.org/wiki/Farad)]です。

- 102 なら (10) * 10^(2) pF
- 333 なら (33) * 10^(3) pF

バンドパスフィルタの部分だけを抜き出した図です。このフィルタは2つのネットワークが組み合わさっています。

<figure>
<img src="img/dr110_cymbal/dr110_bp0.svg" alt="Image of DR-110 band-pass filter schematic." style="width: 480px; padding-bottom: 12px;"/>
</figure>

伝達関数 $H$ です。 $G$ は抵抗の値で単位は [[℧] (mho) あるいは [S] (siemens)](https://en.wikipedia.org/wiki/Siemens_(unit)) 、 $C$ はキャパシタの値で単位は [[F] (ファラド)](https://en.wikipedia.org/wiki/Farad) です。

$$
\begin{aligned}
H &= -\frac{y_a}{y_b}\\

y_a &= - \frac{Gs}{s + G / C}\\
y_b &= - k \frac{s^2 + b_0 s + c_0}{s + a_0}\\

a_0 &= \frac{G_2}{C_1 + C_2},\quad
b_0 = \frac{G_1 (C_1 + C_2)}{C_1 C_2},\quad
c_0 = \frac{G_1 G_2}{C_1 C_2},\quad
k =  \frac{C_1 C_2}{C_1 + C_2}\\
\end{aligned}
$$

伝達関数の式は以下を参考にしました。

- 電子回路、山本 外史、1994、第11刷、p.145
- ["Handbook of Operational Amplifier Active RC Networks" - sboa093a.pdf](http://www.ti.com/lit/an/sboa093a/sboa093a.pdf)、p.14

[Maxima](http://maxima.sourceforge.net/)で解きます。 例としてBP0の抵抗とキャパシタの値を使っています。

```maxima
pico: 10**-12;

G: 1 / 22000;
C: 10 * 10**2 * pico;

G_1: 1 / 82000;
G_2: 1 / 560;
C_1: 33 * 10**2 * pico;
C_2: 33 * 10**2 * pico;

a_0: G_2 / (C_1 + C_2);
b_0: G_1 * (C_1 + C_2) / (C_1 * C_2);
c_0: G_1 * G_2 / (C_1 * C_2);
k:  C_1 * C_2 /  (C_1 + C_2);
y_a: - G * s / (s + G / C);
y_b: - k * (s**2 + b_0 * s + c_0) / (s + a_0);

rat(- y_a / y_b);
```

出力です。

```maxima
/* bp0 */
-(94710000000*s^2+25625000000000000*s)/(3437973*s^3+181681500000*s^2+8030000000000000*s+312500000000000000000)

/* bp1 */
-(268345000000*s^2+35234375000000000*s)/(30108309*s^3+522707500000*s^2+15667187500000000*s+195312500000000000000)
```

得られた値を [`scipy.signal.cont2discrete`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.signal.cont2discrete.html) に渡すことで離散フィルタを設計できます。下のコードではMaximaで計算した値を `applyContinuousFilter` の `system` に渡しています。

```python
import numpy
import soundfile
import scipy.signal as signal


def applyContinuousFilter(samplerate, source, system):
    num, den, dt = signal.cont2discrete(system, 1.0 / samplerate)
    return signal.lfilter(num[0], den, source)


def highMetalFilter(samplerate, source):
    return applyContinuousFilter(samplerate, source, (
        [9.471e+10, 2.5625e+16, 0],
        [3437973.0, 1.816815e+11, 8.03e+15, 3.125e+20],
    ))


def lowMetalFilter(samplerate, source):
    return applyContinuousFilter(samplerate, source, (
        [2.68345e+11, 3.5234375e+16, 0],
        [30108309.0, 5.227075e+11, 1.56671875e+16, 1.953125e+20],
    ))


samplerate = 44100

source = numpy.random.uniform(-0.1, 0.1, samplerate * 2)
dest_high_metal = highMetalFilter(samplerate, source)
dest_low_metal = lowMetalFilter(samplerate, source)

soundfile.write("source.wav", source, samplerate)
soundfile.write("dest_high_metal.wav", dest_high_metal, samplerate)
soundfile.write("dest_low_metal.wav", dest_low_metal, samplerate)
```

BP0とBP1のボード線図です。

<figure>
<img src="img/dr110_cymbal/dr110_bandpass_bode.png" alt="Image of Bode plot of DR-110 band-pass filter." style="width: 480px; padding-bottom: 12px;"/>
</figure>

### ハイパスフィルタ
回路図のいたるところに出てくる1次のハイパスフィルタです。

下図はDR-110のシンバルのエンベロープに関する回路図です。網掛けの色ごとに別のハイパスフィルタを表しています。灰色の矢印は音の信号の流れを表しています。

<figure>
<img src="img/dr110_cymbal/env_crop.png" alt="Image of envelope section of DR-110 schematic." style="width: 600px; padding-bottom: 12px;"/>
</figure>

ハイパスフィルタの部分だけを抜き出した図です。

<figure>
<img src="img/dr110_cymbal/dr110_highpass.svg" alt="Image of DR-110 high-pass filter schematic." style="width: 320px; padding-bottom: 12px;"/>
</figure>

伝達関数 $H$ です。 $R$ は抵抗の値で単位は[Ω]、$C$はキャパシタの値で単位は[F]です。

$$
H = \frac{RC}{RCs + 1}
$$

`scipy.signal` を使った実装の例です。

```python
import numpy
import soundfile
import scipy.signal as signal

def rcHighpass(samplerate, source, r, c):
    rc = r * c
    num, den, dt = signal.cont2discrete(([rc, 0], [rc, 1]), 1.0 / samplerate)
    return signal.lfilter(num[0], den, source)

samplerate = 44100
resistance = 1e6
capacitance = 47e-9

source = numpy.random.uniform(-0.1, 0.1, samplerate * 2)
dest = rcHighpass(samplerate, source, resistance, capacitance)

soundfile.write("source.wav", source, samplerate)
soundfile.write("dest.wav", dest, samplerate)
```

### High Metal と Low Metal
High Metalはハイハットで使われるフィルタネットワークです。

<figure>
<img src="img/dr110_cymbal/DR-110_high_metal.svg" alt="Image of signal flow of DR-110 high metal section." style="width: 600px; padding-bottom: 12px;"/>
</figure>

Low Metalシンバルの音の低い部分の表現に使われるフィルタネットワークです。

<figure>
<img src="img/dr110_cymbal/DR-110_low_metal.svg" alt="Image of signal flow of DR-110 low metal section." style="width: 520px; padding-bottom: 12px;"/>
</figure>

エンベロープの手前での周波数応答です。

<figure>
<img src="img/dr110_cymbal/metal_filter.png" alt="Image of signal flow of DR-110 low metal section." style="width: 480px; padding-bottom: 12px;"/>
</figure>

## エンベロープ
DR-110の回路図から抜き出したシンバルのエンベロープの仕様です。

<figure>
<img src="img/dr110_cymbal/env_time_crop.png" alt="Image of specification of DR-110 cymbal envelopes." style="width: 320px; padding-bottom: 12px;"/>
</figure>

ここではエンベロープに指数的減衰 ([Exponential Decay](https://en.wikipedia.org/wiki/Exponential_decay)) を使います。

$$
y(t) = x(t)e^{at}
$$

$a < 0$ のときに係数 $e^{at}$ が指数的減衰を表します。 $x(t)$ は任意の入力信号です。

$t = 0$ のときの初期値が $E$ かつ、 $t = T$ のとき閾値 $\epsilon$ になるように $a$ を決めます。 $1000$ は $T$ の単位 [ms] のミリから来ています。

$$
\begin{aligned}
\epsilon &= E e^{aT / 1000}\\
\log \left( \frac{\epsilon}{E} \right) &= aT / 1000\\
\frac{1000}{T} \log \left( \frac{\epsilon}{E} \right) &= a
\end{aligned}
$$

NumPyでの実装は次のようになります。 `threshold` の値は適当に決めています。

```python
import numpy
import soundfile

def envelope(time, T, E, threshold=0.4):
    a = numpy.log(threshold / E) / T * 1000
    return numpy.exp(a * time)

samplerate = 44100
duration = 1

num_sample = int(duration * samplerate)

time = numpy.linspace(0, 1, num_sample)

source = numpy.random.uniform(-0.1, 0.1, num_sample)
dest = source * envelope(time, 700, 6)

soundfile.write("source.wav", source, samplerate)
soundfile.write("dest.wav", dest, samplerate)
```

## ミックス
DR-110ではクローズドハイハット (CH) 、オープンハイハット (OH) 、シンバル (CY) の3つの音が使えます。この3つの音の違いはエンベロープとミックスです。

ハイハットの合成には High Metal からの出力だけを使います。エンベロープの Check Point はクローズドが6、オープンは5です。

シンバルの合成には High Metal と Low Metal の出力を使います。High Metal は Check Point 7, 8 のエンベロープを 10:1 の比率で足し合わせます。 Low metal のエンベロープは Check Point 9です。High MetalとLow Metalは約2:1の比率で足し合わせます。

以下は完成したコードのミックスの部分へのリンクです。

- [ミックスを行う関数 metalMix (github.com)](https://github.com/ryukau/filter_notes/blob/de23e42693f5c95e0c2bb32a0796e3a270ef3ec4/docs/demo/dr110_cymbal/dr110_cymbal.py#L86)

## 結果
オシレータの音です。

<label>Oscillator output</label>
<audio controls>
    <source src="snd/dr110_osc.wav" type="audio/wav">
    Audio of simulated DR-110 cymbal oscillator.
</audio>

フィルタを通した音です。

<label>High Metal</label>
<audio controls>
    <source src="snd/dr110_high_metal.wav" type="audio/wav">
    Audio of simulated DR-110 high metal.
</audio>

<label>Low Metal</label>
<audio controls>
    <source src="snd/dr110_low_metal.wav" type="audio/wav">
    Audio of simulated DR-110 low metal.
</audio>

ミックスした音です。

<label>CY</label>
<audio controls>
    <source src="snd/dr110_cy_mixed.wav" type="audio/wav">
    Audio of simulated DR-110 cymbal sound.
</audio>

<label>CH</label>
<audio controls>
    <source src="snd/dr110_ch_mixed.wav" type="audio/wav">
    Audio of simulated DR-110 closed hihat sound.
</audio>

<label>OH</label>
<audio controls>
    <source src="snd/dr110_oh_mixed.wav" type="audio/wav">
    Audio of simulated DR-110 open hihat sound.
</audio>

オリジナルのDR-110の音です。

<label>CY</label>
<audio controls>
    <source src="snd/dr110_cy_original.wav" type="audio/wav">
    Audio of DR-110 cymbal sound.
</audio>

<label>CH</label>
<audio controls>
    <source src="snd/dr110_ch_original.wav" type="audio/wav">
    Audio of DR-110 closed hihat sound.
</audio>

<label>OH</label>
<audio controls>
    <source src="snd/dr110_oh_original.wav" type="audio/wav">
    Audio of DR-110 open hihat sound.
</audio>

作った音はノイズが足りていないように聞こえたので、オシレータのノイズの音量を0.3から6に変えてレンダリングした音です。

<label>CY</label>
<audio controls>
    <source src="snd/dr110_cy_noisy.wav" type="audio/wav">
    Audio of simulated DR-110 cymbal sound.
</audio>

<label>CH</label>
<audio controls>
    <source src="snd/dr110_ch_noisy.wav" type="audio/wav">
    Audio of simulated DR-110 closed hihat sound.
</audio>

<label>OH</label>
<audio controls>
    <source src="snd/dr110_oh_noisy.wav" type="audio/wav">
    Audio of simulated DR-110 open hihat sound.
</audio>

## その他
アタックが弱いです。エンベロープが単純な指数的減衰ではないのかもしれません。

シンバルのミックスの後に続く回路を無視しています。実機ではBalanceで大きく音が変わります。

トランジスタの回路がよく分からなかったのでハイパスフィルタを適当に試行錯誤して設計しました。

エンベロープを打ち切る電圧がわかっていません。

バンドパスフィルタを伝達関数から離散フィルタにするとき、ソルバが次のようなWarningを出します。

```
/__somewhere__/scipy/linalg/basic.py:40: RuntimeWarning: scipy.linalg.solve
Ill-conditioned matrix detected. Result is not guaranteed to be accurate.
Reciprocal condition number/precision: 1.5222285969682624e-17 / 1.1102230246251565e-16
  RuntimeWarning)
```

## 参考文献
- 電子回路、山本 外史、1994、第11刷、p.145、ISBN4-254-22111-8
- ["Handbook of Operational Amplifier Active RC Networks" - sboa093a.pdf](http://www.ti.com/lit/an/sboa093a/sboa093a.pdf)、p.14
- [What is the transfer function for a first order active high-pass filter - Electrical Engineering Stack Exchange](https://electronics.stackexchange.com/questions/125123/what-is-the-transfer-function-for-a-first-order-active-high-pass-filter)
- [Input Impedance of an Amplifier and How to Calculate it](https://www.electronics-tutorials.ws/amplifier/input-impedance-of-an-amplifier.html)
- ["Sine Wave Oscillator" - sloa060.pdf](http://www.ti.com/lit/an/sloa060/sloa060.pdf)
