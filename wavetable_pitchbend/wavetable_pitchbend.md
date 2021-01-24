# ウェーブテーブルのピッチベンド
[Yoshimi](http://yoshimi.sourceforge.net/) の ADDSynth ではユーザが生成した波形を FFT して得られた周波数成分を内部的に保持しています。そして、この周波数成分をノートオンごとに帯域制限して IFFT することでエイリアシングを低減しています。 ADDSynth の方法はノートの音程が一定ならエイリアシングノイズを消せますが、ピッチベンドがかかると高次の倍音の欠落あるいはエイリアシングノイズが出ます。ここではピッチベンドをかけても倍音の欠落やエイリアシングノイズが出ない方法を調べます。

ここで紹介するテクニックは Laurent de Soras さんによる [The Quest For The Perfect Resampler](http://ldesoras.free.fr/doc/articles/resampler-en.pdf) の内容を基にしています。

文中のコードは Python3 です。実行には [SciPy](https://www.scipy.org/), [NumPy](https://numpy.org/) が必要です。

## ウェーブテーブルの帯域制限
帯域制限の方法の要旨を再掲します。詳細は「[ウェーブテーブルのエイリアシングの低減](../wavetable/wavetable.html)」 にまとめています。

まず、波形の周波数成分を用意します。ユーザに直接入力させる方法と、ユーザが指定した波形を FFT する方法があります。

サンプリング周波数を $f_s$​ 、合成する音の周波数 (基本周波数) を $f_0$ とすると、ナイキスト周波数を超えない最大の倍音の次数は $\displaystyle N_h = \left\lfloor \frac{f_s}{2 f_0} \right\rfloor$ です。よって $N_h + 1$ より高い周波数成分の値を 0 にしてから IFFT によってウェーブテーブルを生成することで帯域制限できます。

逆に言えば、周波数成分を帯域制限せずに生成したウェーブテーブルの長さが $N$ のとき、再生時の周波数が $f_0 = \dfrac{f_s}{2N}$ までならエイリアシングが起こりません。 $N_h$ は 0 以上の整数なので計算に床関数が必要ですが、 $f_0$ は実数なので床関数がなくなっています。

## ピッチベンド
例を作ってピッチベンドの挙動を見ていきます。サンプリング周波数 $f_s = 48000\,\text{[Hz]}$ 、基本周波数 $f_0 = 880\,\text{[Hz]}$ として帯域制限したウェーブテーブルを生成します。

<figure>
<img src="img/table.svg" alt="Plot of wavetable and its spectrum." style="padding-bottom: 12px;"/>
</figure>

生成したウェーブテーブルを指定したピッチで再生します。

以下はウェーブテーブルを一定の周波数 $f_0$ で再生したときのスペクトログラムです。明るい線がウェーブテーブルの倍音です。再生時の周波数が一定なので横一直線に伸びています。

<figure>
<img src="img/ConstantPitch.png" alt="Plot of a spectrogram of wavetable oscillator with constant pitch." style="padding-bottom: 12px;"/>
</figure>

以下は基本周波数より高くなるようにウェーブテーブルオシレータを $f_0$ から $2 f_0$ までピッチベンドしたときのスペクトログラムです。明るい線が上端で跳ね返っているように見えるのがエイリアシングノイズです。

<figure>
<img src="img/RisingPitch.png" alt="Plot of a spectrogram of wavetable oscillator with rising pitch." style="padding-bottom: 12px;"/>
</figure>

以下は基本周波数より低くなるようにウェーブテーブルオシレータを $f_0$ から $0.5f_0$ までピッチベンドしたときのスペクトログラムです。高次の倍音が欠けています。

<figure>
<img src="img/FallingPitch.png" alt="Plot of a spectrogram of wavetable oscillator with falling pitch." style="padding-bottom: 12px;"/>
</figure>

倍音の欠落とエイリアシングノイズの問題を解決するためにはオーバーサンプリングが使えます。オーバーサンプリングの倍率を $L$ とするとナイキスト周波数を超えない最大の倍音の次数は $\displaystyle N_h = \left\lfloor \frac{f_s}{2 L f_0} \right\rfloor$ です。例として 2 倍のオーバーサンプリングを行います。

<figure>
<img src="img/table2x.svg" alt="Plot of 2 times oversampled wavetable and its spectrum." style="padding-bottom: 12px;"/>
</figure>

以下は基本周波数より高くなるようにウェーブテーブルオシレータを $f_0$ から $2 f_0$ までピッチベンドしたときのスペクトログラムです。オーバーサンプリングによってできたマージンがあるのでエイリアシングが出ていません。

<figure>
<img src="img/RisingPitch2x.png" alt="Plot of a spectrogram of 2 times oversampled wavetable oscillator with rising pitch." style="padding-bottom: 12px;"/>
</figure>

マージンに現れた成分はダウンサンプリング時にかけるローパスフィルタで低減できます。以下はローパスフィルタをかけた後の信号のスペクトログラムです。

<figure>
<img src="img/FilteredRisingPitch2x.png" alt="Plot of a spectrogram of wavetable oscillator after passing down-sampling lowpass filter." style="padding-bottom: 12px;"/>
</figure>

以下はフィルタをかけた信号をデシメーションした後のスペクトログラムです。オーバーサンプリングなしで高いほうに向かってピッチベンドしたときと比べると、エイリアシングノイズが消えていることが確認できます。また、ダウンサンプリングのフィルタによって 20000 Hz 以上の成分が低減されているので、図の上部が真っ黒になっています。

<figure>
<img src="img/FilteredRisingPitch2xDecimated.png" alt="Plot of a spectrogram of down-sampled wavetable oscillator output." style="padding-bottom: 12px;"/>
</figure>

アップサンプリングの倍率を $L$ とするとエイリアシングノイズが出ないピッチベンドの範囲は $[f_0, L f_0)$ です。ダウンサンプリングに使うローパスフィルタが理想的な特性ならピッチベンドの範囲の上限を $\dfrac{3}{2} L f_0$ まで広くできます。あとは必要な範囲に応じて複数のウェーブテーブルを用意して、入力されたピッチによって再生するウェーブテーブルを切り替えることで、ピッチベンドしてもエイリアシングノイズが出ないオシレータが作れます。

基本周波数より低くなるようにピッチベンドするとオーバーサンプリングの有無にかかわらず倍音の欠落が防げません。よって切り替えを行うときには、基本周波数が入力された周波数 $f$ 以下のウェーブテーブルの中から、基本周波数が $f$ に最も近いものを選びます。

### コード
ここまで掲載した図の作成に使ったウェーブテーブルの生成とピッチベンドのテストコードです。 `table` がウェーブテーブルです。プロットの部分は省略しています。

```python
import numpy as np
import scipy.signal as signal

def saw_spectrum(size):
    spec = [0] + [(-1)**k / k for k in range(1, int(size / 2 + 1))]
    spec = -1j * size / np.pi * np.array(spec)
    return spec

def renderTable(table, freqStart, freqEnd, fs):
    phase = np.linspace(freqStart / fs, freqEnd / fs, fs).cumsum() % 1.0
    x = np.linspace(0, 1, len(table))
    return np.interp(phase, x, table, period=1)

def testNaiveTable(fs=48000, f0=880, size=8192):
    """オーバーサンプリングなし。"""
    source = saw_spectrum(size)  # のこぎり波の周波数成分を生成。

    spectrum = source.copy()
    cutoff = int(fs / (2 * f0))  # 数式の N_h に対応。
    spectrum[cutoff:] = 0
    table = np.fft.irfft(spectrum)

    sigF0 = renderTable(table, f0, f0, fs)          # ピッチ一定
    sigUp = renderTable(table, f0, 2 * f0, fs)      # ピッチ上昇
    sigDown = renderTable(table, f0, 0.5 * f0, fs)  # ピッチ下降

    # プロットは省略

def test2xTable(fs=48000, f0=880, size=8192):
    """2倍のオーバーサンプリング。"""
    oversample = 2  # 数式の L に対応。
    upRate = oversample * fs

    source = saw_spectrum(size)

    spectrum = source.copy()
    cutoff = int(upRate / (2 * oversample * f0))  # 数式の N_h に対応。
    spectrum[cutoff:] = 0
    table = np.fft.irfft(spectrum)

    lowpass = signal.ellip(12, 0.01, 100, 0.4, "low", output="sos", fs=2)

    sigUp = renderTable(table, f0, 2 * f0, upRate)
    sigFiltered = signal.sosfilt(lowpass, sigUp)
    sigDecimated = sigFiltered[::2]

    # プロットは省略

testNaiveTable()
test2xTable()
```

## 設計のポイント
以下はウェーブテーブルオシレータの仕様を決める変数の一覧です。

- オーバーサンプリングの倍率
- 1 つのウェーブテーブルのピッチベンドの範囲
- ウェーブテーブルの長さ
- 倍音が欠落しない最低周波数

まずオーバーサンプリングの倍率を決めます。すると 1 つのウェーブテーブルのピッチベンドの範囲も決まります。 FM などのモジュレーションのかかり方を細かくチューニングしたいときはオーバーサンプリングの倍率と 1 つのウェーブテーブルのピッチベンドの範囲を分けて決めることも考えられます。

次にウェーブテーブルの長さ $N$ を決めます。するとオシレータ全体で倍音の欠落が起こらない最低周波数は $f_s / N$ となります。逆に最低周波数 $f_L$ から $N$ を決めたいときは以下の式が使えます。

$$
N = \frac{f_s}{f_L}
$$

最低周波数 $f_L$ は、ナイキスト周波数 $f_s / 2$ から 1 つのウェーブテーブルのピッチベンドの範囲を減算していくように決めると帯域制限の計算が簡単になります。例えば 1 つのウェーブテーブルのピッチベンドの範囲が 1 オクターブのときは、 $N$ を 2 のべき乗にして、波形の周波数成分の高いほうから全体の $\dfrac{N}{2}, \dfrac{N}{4}, \dfrac{N}{8}, \dots$ 個の要素を 0 埋めすることで帯域制限できます。

以下は 1 オクターブ間隔で生成したウェーブテーブルと基本周波数、インデックスの対応表です。インデックスがそのまま基本周波数のオクターブ差になっています。表の半音 (semitone) は $f_L$ を 0 とした相対的な値です。

<figure>
<img src="img/multi_wavetable.svg" alt="Image of ." style="padding-bottom: 12px;"/>
</figure>

## オーバーサンプリングを行う実装
オーバーサンプリングを行うウェーブテーブルオシレータのブロック線図です。

<figure>
<img src="img/simple_wavetable_block_diagram.svg" alt="Block diagram of simple wavetable oscillator." style="padding-bottom: 12px;"/>
</figure>

ピッチベンドの範囲は $[f_0, 2 f_0)$ とします。こうすると事前に用意するウェーブテーブルが 1 オクターブ間隔になります。

### コード
以下は 2 倍のオーバーサンプリングを行うウェーブテーブルオシレータのコードです。波形はのこぎり波で固定しています。

```python
import numpy as np
import scipy.signal as signal

def frequencyToMidinote(freq):
    return 12 * np.log2(freq / 440) + 69

def midinoteToFrequency(note):
    return 440 * np.exp2((note - 69) / 12)

def saw_spectrum(size):
    spec = [0] + [(-1)**k / k for k in range(1, int(size / 2 + 1))]
    spec = -1j * size / np.pi * np.array(spec)
    return spec

class TableOsc:
    def __init__(self, samplerate, oversample=2):
        self.fs = oversample * samplerate
        self.phase = 0

        exponent = int(np.log2(self.fs / 10))  # 倍音が欠けない最低の周波数は 10 Hz.
        size = 2**exponent
        spectrum = saw_spectrum(size)

        self.basefreq = self.fs / size
        self.basenote = frequencyToMidinote(self.basefreq)

        self.table = []
        for idx in range(2, exponent + 1):
            spec = np.zeros_like(spectrum)
            cutoff = int(size / 2**(idx))
            spec[:cutoff] = spectrum[:cutoff]
            tbl = np.fft.irfft(spec)
            tbl = np.hstack((tbl, tbl[0]))
            self.table.append(tbl)

    def process(self, note):
        """
        note は MIDI ノート番号。
        """
        self.phase += midinoteToFrequency(note) / self.fs
        self.phase -= np.floor(self.phase)

        octave = int(np.clip((note - self.basenote) / 12, 0, len(self.table) - 1))
        pos = (len(self.table[0]) - 1) * self.phase
        idx = int(pos)
        frac = pos - idx

        x0 = self.table[octave][idx]
        x1 = self.table[octave][idx + 1]
        return x0 + frac * (x1 - x0)

gain = 0.25
samplerate = 48000
lowpass = signal.ellip(12, 0.01, 100, 0.4, "low", output="sos", fs=2)

osc = TableOsc(samplerate)

sig = np.empty(16 * samplerate)
note = np.linspace(0, 128, len(sig) // 2)
note = np.hstack((note, np.flip(note)))
for i in range(len(sig)):
    sig[i] = osc.process(note[i])
sig = gain * signal.sosfilt(lowpass, sig)[::2]
```

`sig` がレンダリング結果です。

`frequencyToMidinote` と `midinoteToFrequency` は MIDI ノート番号と周波数を変換する関数です。以下のページに計算式が載っています。

- [Note names, MIDI numbers and frequencies](http://newt.phys.unsw.edu.au/jw/notes.html)

`saw_spectrum` はのこぎり波のスペクトラムを生成する関数です。以下のリンク先の式を参考に、位相を 180° ずらしています。

- [Fourier Series--Sawtooth Wave -- from Wolfram MathWorld](https://mathworld.wolfram.com/FourierSeriesSawtoothWave.html)

`__init__` ではウェーブテーブルの長さが $2^k$ になるようにしています。 $k$ はコードの `exponent` と対応しています。

ウェーブテーブルの補間には線形補間を使っています。インデックスの余り演算を避けるために `self.table` の先頭の値を最後尾に追加して、長さを $2^k + 1$ にしています。例えばウェーブテーブルの値が `[0, 1, 2]` のときは `[0, 1, 2, 0]` に変換されます。

C++ による実装は別ページに掲載しています。

- [TableOsc の C++ 実装を読む (github.com)](https://github.com/ryukau/filter_notes/blob/f466fff89bb16ce29d9a2e0cbf7d538e0a8be7f8/wavetable_pitchbend/cpp/bench/bench.cpp#L232)

### 出力の癖
以下の図は `TableOsc` の出力のスペクトログラムです。MIDI ノート番号で 0, 128, 0 とピッチベンドしています。

<figure>
<img src="img/TableOsc.png" alt="Spectrogram of an output of 2 times oversampled wavetable oscillator." style="padding-bottom: 12px;"/>
</figure>

エイリアシングは抑えられていますが、垂直な縦線が見えるのでポップノイズが出ています。このポップノイズはウェーブテーブルの切り替えによって発生しています。

ポップノイズの低減方法としては以下の 2 つのアプローチがあります。

- オーバーサンプリングの倍率を上げる。
- ウェーブテーブルを周波数軸に沿った方向でも補間する。

オーバーサンプリングを 8 倍まで上げるとポップノイズが聞こえなくなりました。

ウェーブテーブルを周波数軸に沿った方向でも線形補間すると、 2 倍のオーバーサンプリングでもポップノイズが聞こえなくなります。この補間は位相と周波数の軸に沿った[バイリニア補間](https://en.wikipedia.org/wiki/Bilinear_interpolation)です。ただし、周波数方向の補間によって、高次の倍音が弱まる区間が現れる欠点があります。以下は 2 倍のオーバーサンプリングで周波数方向の補間を追加したときのスペクトログラムです。図の上側の周期的に暗くなっているところで倍音が弱まっています。

<figure>
<img src="img/TableOscBilinear.png" alt="Spectrogram of an output of bilinear interpolated wavetable oscillator." style="padding-bottom: 12px;"/>
</figure>

倍音の弱まりを抑えるには、結局のところオーバーサンプリングの倍率を上げて対処することになります。ただし、周波数方向の補間なしのときと比べると、低めの倍率が使えるかもしれません。

- [バイリニア補間を使った Python3 の実装を読む (github.com)](https://github.com/ryukau/filter_notes/blob/f466fff89bb16ce29d9a2e0cbf7d538e0a8be7f8/wavetable_pitchbend/pitchbend.py#L78)

## オーバーサンプリングを行わない実装
ウェーブテーブルの数を増やして周波数方向の間隔を狭くすることで、オーバーサンプリング無しでもエイリアシングを低減できます。この方法は [CubicPadSynth](https://ryukau.github.io/VSTPlugins/manual/CubicPadSynth/CubicPadSynth_ja.html) で使いました。

### コード
以下はオーバーサンプリングを行わないウェーブテーブルオシレータのコードです。 1 半音ごとに、合計 137 のウェーブテーブルを用意しています。波形はのこぎり波で固定しています。

```python
class LpsOsc:
    def __init__(self, samplerate):
        self.fs = samplerate
        self.phase = 0

        minFreq = midinoteToFrequency(0)
        expFloat = np.log2(self.fs / np.floor(minFreq))
        exponent = int(np.ceil(expFloat))
        self.size = 2**exponent
        spectrum = saw_spectrum(self.size)

        self.table = []
        self.maxNote = np.ceil(frequencyToMidinote(20000))
        spec = np.zeros_like(spectrum)
        for idx in range(int(self.maxNote)):
            freq = midinoteToFrequency(idx)
            cutoff = int(len(spectrum) * minFreq / freq)

            spec[:cutoff] = spectrum[:cutoff]
            spec[cutoff:] = 0

            tbl = np.fft.irfft(spec)
            tbl = np.hstack((tbl, tbl[0]))
            self.table.append(tbl)
        self.table.append(np.zeros_like(self.table[0]))

    def processPhase(self, note):
        tick = midinoteToFrequency(note) * self.size / self.fs
        if tick >= self.size or tick < 0:
            tick = 0

        self.phase += tick
        if self.phase >= self.size:
            self.phase -= self.size

    def processLow(self, note):
        """
        note が 0 以下のときは self.table[0] だけでレンダリング。
        """
        idx = int(self.phase)
        a0 = self.table[0][idx]
        a1 = self.table[0][idx + 1]
        fracX = self.phase - np.floor(self.phase)
        return a0 + fracX * (a1 - a0)

    def process(self, note):
        """
        note は MIDI ノート番号。
        """
        self.processPhase(note)

        if note < 0:
            return processLow(note)

        nn = int(note)
        if nn >= self.maxNote:
            nn = self.maxNote - 1

        idx = int(self.phase)
        a0 = self.table[nn][idx]
        a1 = self.table[nn][idx + 1]
        b0 = self.table[nn + 1][idx]
        b1 = self.table[nn + 1][idx + 1]

        fracX = self.phase - np.floor(self.phase)
        x0 = a0 + fracX * (a1 - a0)
        x1 = b0 + fracX * (b1 - b0)

        fracY = note - nn
        return x0 + fracY * (x1 - x0)

gain = 0.25
samplerate = 48000

osc = LpsOsc(samplerate)

sig = np.empty(8 * samplerate)
note = np.linspace(0, 128, len(sig) // 2)
note = np.hstack((note, np.flip(note)))
for i in range(len(sig)):
    sig[i] = osc.process(note[i])
sig = gain * sig
```

補間にはバイリニア補間を使っています。

`self.table` のインデックスは MIDI ノート番号と対応しています。 `process` に入力された `note` が 0 より小さいときは、 `processLow` に飛んでインデックス 0 のウェーブテーブルだけを使って出力を計算しています。

C++ による実装は別ページに掲載しています。

- [LpsOsc の C++ 実装を読む (github.com)](https://github.com/ryukau/filter_notes/blob/f466fff89bb16ce29d9a2e0cbf7d538e0a8be7f8/wavetable_pitchbend/cpp/bench/bench.cpp#L460)

### 出力の癖
以下の図は `LpsOsc` の出力のスペクトログラムです。MIDI ノート番号で 0, 128, 0 とピッチベンドしています。

<figure>
<img src="img/LpsOsc.png" alt="Spectrogram of an output of wavetable oscillator without oversampling." style="padding-bottom: 12px;"/>
</figure>

図の上側を見ると少しだけエイリアシングノイズが出ていることが確認できます。原因はわからなかったのですが、以下の変更では解決しませんでした。

- `process` 内の `self.table` の一つ目のインデックスの変更。
- `cutoff` の変更。
- バイリニア補間から位相方向のみの線形補間に変更。

全体としてはエイリアシングノイズが低減できているので使える場面はありそうです。速度はオーバーサンプリングがない分だけ速いです。ただし、ウェーブテーブルが長くなると遅くなる可能性があります。

### ミップマップを使う実装
The Quest For The Perfect Resampler で紹介されていた方法です。

[ミップマップ](https://en.wikipedia.org/wiki/Mipmap)は3DCG のテクスチャで使われる手法です。事前にテクスチャの画素数を $1/4^n$ に減らしたデータを用意しておいて、カメラからの距離に応じて十分に解像度の低いデータを選ぶことで計算量を減らすことができます。

音のデータは 1 次元なので 1 オクターブごとにウェーブテーブルの長さを半分にしていきます。言い換えると 1 オクターブ上がるごとにサンプリング周波数を半分にしていきます。

ミップマップを使ったウェーブテーブルの利点は、メモリが節約できる点です。ウェーブテーブルの数を $M$ 、最も倍音を多く含むウェーブテーブルの長さを $N$ とすると、必要なメモリの量はミップマップありのときで $2N$ 、 ミップマップなしのときで $NM$ です。

ミップマップの欠点は補間の計算量です。入力されたピッチに応じてフィルタ係数を計算しなおしつつアップサンプリングする必要があるので、線形補間と比べるとかなり重たそうです。補間以外の部分を実装して速度を比較してみたのですが、オーバーサンプリングありの実装とほぼ同じだったのでテスト対象から外しました。

以下は線形補間を使ったミップマップオシレータのスペクトログラムです。補間の質が悪いので、全体に薄くノイズが乗っていることが見て取れます。

<figure>
<img src="img/MipmapOsc.png" alt="Spectrogram of an output of mipmap wavetable oscillator with linear interpolation." style="padding-bottom: 12px;"/>
</figure>

以下は C++ と Python3 による実装へのリンクです。

- [MipmapOsc の C++ 実装を読む (github.com)](https://github.com/ryukau/filter_notes/blob/f466fff89bb16ce29d9a2e0cbf7d538e0a8be7f8/wavetable_pitchbend/cpp/bench/bench.cpp#L381)
- [MipmapOsc の Python3 実装を読む (github.com)](https://github.com/ryukau/filter_notes/blob/f466fff89bb16ce29d9a2e0cbf7d538e0a8be7f8/wavetable_pitchbend/pitchbend.py#L124)

## 音のサンプル
各実装でのこぎり波をレンダリングした結果です。

<figure>
  <figcaption>オーバーサンプリングを行う実装、2 倍のオーバーサンプリング</figcaption>
  <audio controls>
    <source src="snd/simple.wav" type="audio/wav">
  </audio>
</figure>

<figure>
  <figcaption>オーバーサンプリングを行わない実装</figcaption>
  <audio controls>
    <source src="snd/lpsosc.wav" type="audio/wav">
  </audio>
</figure>

<figure>
  <figcaption>ミップマップを使う実装、線形補間</figcaption>
  <audio controls>
    <source src="snd/mipmap.wav" type="audio/wav">
  </audio>
</figure>

オーバーサンプリングを行う実装は、音が高くなったときにウェーブテーブルの切り替えによるポップノイズが聞こえます。ミップマップを使う実装では音が高くなったときに線形補間によるノイズが聞こえます。

### オーバーサンプリングを行う実装のポップノイズ
以下はオーバーサンプリングを行う実装で、オーバーサンプリングの倍率と、周波数方向の補間の有無を変えて、のこぎり波をレンダリングした結果です。

<figure>
  <figcaption>2 倍のオーバーサンプリング、位相方向のみ線形補間</figcaption>
  <audio controls>
    <source src="snd/TableOsc_x2_switch.wav" type="audio/wav">
  </audio>
</figure>

<figure>
  <figcaption>4 倍のオーバーサンプリング、位相方向のみ線形補間</figcaption>
  <audio controls>
    <source src="snd/TableOsc_x4_switch.wav" type="audio/wav">
  </audio>
</figure>

<figure>
  <figcaption>8 倍のオーバーサンプリング、位相方向のみ線形補間</figcaption>
  <audio controls>
    <source src="snd/TableOsc_x8_switch.wav" type="audio/wav">
  </audio>
</figure>

<figure>
  <figcaption>2 倍のオーバーサンプリング、位相方向と周波数方向でバイリニア補間</figcaption>
  <audio controls>
    <source src="snd/TableOscBilinear_x2_switch.wav" type="audio/wav">
  </audio>
</figure>

<figure>
  <figcaption>4 倍のオーバーサンプリング、位相方向と周波数方向でバイリニア補間</figcaption>
  <audio controls>
    <source src="snd/TableOscBilinear_x4_switch.wav" type="audio/wav">
  </audio>
</figure>

<figure>
  <figcaption>8 倍のオーバーサンプリング、位相方向と周波数方向でバイリニア補間</figcaption>
  <audio controls>
    <source src="snd/TableOscBilinear_x8_switch.wav" type="audio/wav">
  </audio>
</figure>

4 倍のオーバーサンプリング、位相方向のみ線形補間では、音量を上げるとポップノイズを聞き取ることができます。

2 倍のオーバーサンプリング、位相方向と周波数方向でバイリニア補間での、高次倍音の周期的な弱まりは聞き取れなかったです。波形によって変わるかもしれません。

## 参考文献
- [The Quest For The Perfect Resampler - resampler-en.pdf](http://ldesoras.free.fr/doc/articles/resampler-en.pdf)
- [Note names, MIDI numbers and frequencies](http://newt.phys.unsw.edu.au/jw/notes.html)
- [Fourier Series--Sawtooth Wave -- from Wolfram MathWorld](https://mathworld.wolfram.com/FourierSeriesSawtoothWave.html)
