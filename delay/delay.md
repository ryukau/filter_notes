# ディレイの実装
ディレイは信号を遅らせる処理で、1つのリングバッファと書き込みと読み込みの2つポインタで表現できます。

<figure>
<img src="img/ringbuffer.svg" alt="Image of ringbuffer." style="width: 640px; padding-bottom: 12px;"/>
</figure>

時間が1サンプル進むたびにリングバッファの書き込みポインタと読み取りポインタのインデックスを1つ進めます。バッファの終端から1サンプル進むとインデックス0に戻ります。

ディレイ時間は実装の内部ではサンプル数で表現されます。
秒数 $t$ とサンプリング周波数 $f_s$ を掛け合わせることでサンプル数で表された時間 $d$ が計算できます。

$$
d = tf_s
$$

例として $t = 0.23,\,f_s = 10$ とすると $d = 2.3$ となります。 $d$ はリングバッファのインデックスを表していますが、実数なので `buffer[d]` などと書くとエラーになります。そこで補間が必要になってきます。

補間を行うディレイは、整数サンプルのディレイと、その後に続く補間に処理を分けることができます。補間は分数サンプルのディレイを行うフィルタとして扱えるので、分数ディレイフィルタ (fractional delay filter) と呼ばれることがあります。

<figure>
<img src="img/delay_outline.svg" alt="Image of signal flow of a delay." style="width: 720px; padding-bottom: 12px;"/>
</figure>

## 整数ディレイ
まずは補間のないディレイを実装します。

```javascript
class IntDelay {
  constructor(sampleRate, time) {
    this.sampleRate = sampleRate
    this.buf = new Array(Math.ceil(sampleRate * time) + 1).fill(0)

    this.wptr = 0
    this.time = time
  }

  set time(time) {
    var delay = time * this.sampleRate
    this.rptr = Math.floor(this.wptr - delay)
    while (this.rptr < 0) this.rptr += this.buf.length
  }

  process(input) {
    this.buf[this.wptr] = input
    this.wptr = (this.wptr + 1) % this.buf.length

    var output = this.buf[this.rptr]
    this.rptr = (this.rptr + 1) % this.buf.length

    return output
  }
}
```

`time` の単位は秒です。 `this.buf` の長さはコンストラクタの `time` に応じて設定されるので、後でディレイタイムを変更するときは予想される最大値をコンストラクタの `time` に渡しておきます。

`while` の行で `this.rptr` が [0, this.buf.length] の範囲に収まるようにしています。 `delay >= this.buf.length || delay < 0` のときは正しい値が設定されないので注意してください。

テストします。

```javascript
var sampleRate = 44100

var source = new Array(8).fill(0)
source[0] = 1

var result = [source]
for (var time = 1; time <= 4; ++time) {
  var delay = new IntDelay(sampleRate, time / sampleRate)

  var dest = new Array(source.length).fill(0)
  for (var i = 0; i < source.length; ++i) {
    dest[i] = delay.process(source[i])
  }
  result.push(dest)
}

console.log(result)
```

出力です。

```javascript
[
  [1, 0, 0, 0, 0, 0, 0, 0], // source
  [0, 1, 0, 0, 0, 0, 0, 0], // 1 sample delay
  [0, 0, 1, 0, 0, 0, 0, 0], // 2 sample delay
  [0, 0, 0, 1, 0, 0, 0, 0], // 3 sample delay
  [0, 0, 0, 0, 1, 0, 0, 0]  // 4 sample delay
]
```

## 分数ディレイフィルタ
ここでは分数ディレイフィルタの実装に、テイラー展開に基づくラグランジュ補間を使います。次の図はテイラー展開に基づくラグランジュ補間のブロック線図です。

<figure>
<img src="img/lagrange_taylor.svg" alt="Image of a implementation of lagrange interpolation based on taylor expansion." style="width: 640px;padding-bottom: 12px;"/>
</figure>

$z^{-n}$ は $n$ サンプルディレイを表しています。

$\Delta$ はフィルタを通したときに追加されるディレイで、単位はサンプル数、範囲は $[0, N]$ の実数です。 $N$ 次のときの $\Delta$ の値は $[(N-1)/2, (N+1)/2]$ の範囲で設定するといいようなので、次の式のように書けます。

$$
\Delta = \frac{N - 1}{2} + \delta, \quad 0 \leq \delta < 1
$$

例えば3次のラグランジュ補間では $\Delta$ の値が1.0 から 2.0 の範囲に収まるように使えば良いことがわかります。

上のブロック線図をそのまま実装すると $\delta=1$ のときに時間が1サンプル戻るようになります。

<figure>
<img src="img/fd_filter_direction.svg" alt="Image of ." style="width: 480px;padding-bottom: 12px;"/>
</figure>

3次のラグランジュ補間を実装します。

```javascript
class Lagrange3Interp {
  constructor() {
    this.xd = new Array(3).fill(0)
    this.diff = new Array(3).fill(0)
  }

  process(input, fraction) {
    var delta = 1 + fraction

    this.diff[0] = input - this.xd[0]
    this.diff[1] = this.diff[0] - this.xd[1]
    this.diff[2] = this.diff[1] - this.xd[2]

    this.xd[0] = input
    this.xd[1] = this.diff[0]
    this.xd[2] = this.diff[1]

    var sig = 0
    sig = (2 - delta) / 3 * (this.diff[2] + sig)
    sig = (1 - delta) / 2 * (this.diff[1] + sig)
    sig = (0 - delta) / 1 * (this.diff[0] + sig)
    sig = input + sig

    return sig
  }
}
```

`fraction` は補間の位置 $\delta$ で、範囲 `[0.0, 1.0)` の値です。 `sig` の計算はわかりやすさのために `(0 - delta) / 1` のような意味のない計算も明示的に書いています。

$N$ 次のラグランジュ補間が計算できるように一般化します。

```javascript
class LagrangeNInterp {
  constructor(order) {
    this.N = order
    this.fix = (this.N - 1) / 2
    this.xd = new Array(this.N).fill(0)
    this.diff = new Array(this.N).fill(0)
  }

  reset() {
    this.xd.fill(0)
    this.diff.fill(0)
  }

  push(input) {
    this.diff[0] = input - this.xd[0]
    this.xd[0] = input
    for (var i = 1; i < this.N; ++i) {
      this.diff[i] = this.diff[i - 1] - this.xd[i]
      this.xd[i] = this.diff[i - 1]
    }
  }

  at(fraction) {
    var delta = fraction + this.fix
    var sig = 0

    var i = this.N
    while (i > 0) {
      var next = i - 1
      sig = (next - delta) / i * (this.diff[next] + sig)
      i = next
    }

    return sig + this.xd[0]
  }

  process(input, fraction) {
    this.push(input)
    return this.at(fraction)
  }
}
```

### 奇数の $N$ と偶数の $N$
サンプル数で表したディレイ時間 $d$ を整数部 $n$ と分数部 $\delta$ に分けます。

$$
d = n + \delta
\quad \text{where} \quad n = \lfloor d \rfloor
, \quad 0 \leq \delta < 1.
$$

ラグランジュ補間では、計算のために補間したい位置の周りにある $N$ サンプルを参照するのですが、 $N$ が偶数のときは $0 \leq \delta \leq 0.5$ の場合と $0.5 < \delta \leq 1$ の場合で参照するサンプルの組が変わってしまいます。

<figure>
<img src="img/lagrange_fd_filter_odd_even_consideration.svg" alt="Image of several pair of samples that is needed to compute lagrange interpolation." style="width: 640px;padding-bottom: 12px;"/>
</figure>

この文章では $N$ が奇数の実装だけを掲載しています。$N$ が偶数のときは入力を1サンプルずらした2つの分数ディレイフィルタを用意して、 $\Delta$ の値によって通過させるフィルタを入れ替える実装が考えられます。

## オーバーサンプリング
整数ディレイと分数ディレイフィルタをそのままつないだだけだと、ディレイ時間を変更したときにノイズが乗ります。ノイズはサンプル数で表されたディレイ時間の整数部が変わるときに乗るようです。

いろいろ試したところ、オーバーサンプリングによってノイズを減らすことができました。オーバーサンプリングは、信号のサンプリング周波数を $K$ 倍に増やして処理を行うことでノイズを低減するテクニックです。ここではサンプルとサンプルの間の値は分かっていないので、補間を使って近似します。

実装ではオーバーサンプリングしたい区間をアップサンプリングとダウンサンプリングの処理で挟みます。アップサンプリングはサンプリング周波数を $K$ 倍に増やす処理、ダウンサンプリングはサンプリング周波数を $1 / K$ 倍に減らす処理です。ここではアップサンプリングの補間に分数ディレイフィルタを使います。

オーバーサンプリングの処理は次のように書けます。

```javascript
var interp = new LagrangeNInterp(3)

function overSampling(section, input, K) {
  interp.push(input)

  for (var i = 0; i < K; ++i) {
    var fraction = (K - i) / K
    var output = section.process(interp.at(fraction))
  }

  return output
}
```

## ディレイ時間の補正
整数ディレイの前後に取り付けられた分数ディレイフィルタの遅延を補正します。次の図はオーバーサンプリングを追加したときのディレイの処理の流れです。

<figure>
<img src="img/oversampling_outline.svg" alt="Image of ." style="width: 720px;padding-bottom: 12px;"/>
</figure>

FD Filter は Fractional Delay Filter の略です。

例として、オーバーサンプリングを4倍、書き込みと読み取りのラグランジュ補間を3次に設定したディレイを使って、入力信号 を3.375 (3 + 3/8) サンプル遅らせることを考えます。

入力信号です。

<figure>
<img src="img/oversampling_input.svg" alt="Image of ." style="width: 640px;padding-bottom: 12px;"/>
</figure>

今、図の一番右端のサンプルがディレイに入力されました。右端のサンプルを入力を補間する分数ディレイフィルタに入力して信号を補間します。

<figure>
<img src="img/oversampling_oversampled.svg" alt="Image of ." style="width: 640px;padding-bottom: 12px;"/>
</figure>

入力信号を補間した状態です。上段の time は入力信号のサンプリング周波数に基づいてサンプル数で表された時間です。下段の time in buffer はディレイのリングバッファ上でのサンプル数で表された時間で、オーバーサンプリングされた信号のサンプリング周波数に基づいています。

入力信号のサンプリング周波数で 3.375 サンプルのディレイは 4 倍にオーバーサンプリングされると 4 * 3.375 = 13.5 サンプルのディレイになります。

`wptr` はリングバッファへの書き込みポインタです。上の図では書き込みを終えたときの位置を指しています。ここでの目的は `wptr` の位置から、分数ディレイフィルタで追加される遅れを考慮した整数ディレイの読み取りポインタ `rptr` の位置を計算することです。

`input` と `wptr` のずれは書き込み時の補間によるものです。読み取りでも補間を行うのでディレイが追加されます。書き込みの補間によるディレイを `wDelay` 、読み取りの補間によるディレイを `rDelay` とすると指定したディレイ時間 `d_buf` との関係は次の図のようになります。

<figure>
<img src="img/oversampling_pointer_offsets.svg" alt="Image of ." style="width: 640px;padding-bottom: 12px;"/>
</figure>

各変数の計算式をまとめます。 `wOrder` は入力の補間の次数、 `rOrder` は出力の補間の次数です。

```
wDelay = oversampling * (wOrder - 1) / 2
rDelay = (rOrder - 1) / 2
d_buf = oversampling * input_samplerate * time
rptr = ceil(wptr - d_buf) + wDelay + rDelay
```

## まとめて実装
ここまでのテクニックを一つにまとめて実装します。

```javascript
class Delay {
  // readOrder と writeOrder は奇数。
  constructor(sampleRate, time, overSample = 2, readOrder = 1, writeOrder = 1) {
    this.overSample = overSample
    this.sampleRate = this.overSample * sampleRate

    // 補間で追加されるディレイの補正値。
    this.fix = (readOrder - 1 + this.overSample * (writeOrder - 1)) / 2

    this.buf = new Array(Math.ceil(time * this.sampleRate) + 1).fill(0)

    this.rInterp = new LagrangeNInterp(readOrder)
    this.wInterp = new LagrangeNInterp(writeOrder)

    this.wptr = 0
    this.time = time
  }

  set time(value) {
    var delayTime = Math.max(
      this.fix, Math.min(this.sampleRate * value, this.buf.length))
    this.rFraction = delayTime % 1
    this.rptr = Math.ceil(this.wptr - delayTime) + this.fix

    // mod(rptr, buf.length). 正の余りを計算。
    while (this.rptr < 0) this.rptr += this.buf.length
    this.rptr %= this.buf.length
  }

  reset() {
    this.buf.fill(0)
    this.rInterp.reset()
    this.wInterp.reset()
  }

  process(input) {
    this.wInterp.push(input)

    for (var i = 1; i <= this.overSample; ++i) {
      var wFraction = (this.overSample - i) / this.overSample

      this.buf[this.wptr] = this.wInterp.at(wFraction)
      this.wptr = (this.wptr + 1) % this.buf.length

      this.rInterp.push(this.buf[this.rptr])
      this.rptr = (this.rptr + 1) % this.buf.length
    }

    return this.rInterp.at(this.rFraction)
  }
}
```

音のサンプルです。220Hzのサイン波を入力して、ディレイ時間を2Hzのサイン波で変調しています。入力 と出力の補間の次数は同じです。

入力したサイン波です。

<label>220Hzのサイン波</label>
<audio controls>
    <source src="snd/sin_source.wav" type="audio/wav">
    Audio of 220Hz sine wave.
</audio>

オーバーサンプリングなしの音です。ノイズがのっています。

<label>1倍、1次補間</label>
<audio controls>
    <source src="snd/sin_01x_r1_w1.wav" type="audio/wav">
    Audio of delay. 1x oversampling, 1st order lagrange interpolation.
</audio>

<label>1倍、3次補間</label>
<audio controls>
    <source src="snd/sin_01x_r3_w3.wav" type="audio/wav">
    Audio of delay. 1x oversampling, 3rd order lagrange interpolation.
</audio>

<label>1倍、5次補間</label>
<audio controls>
    <source src="snd/sin_01x_r5_w5.wav" type="audio/wav">
    Audio of delay. 1x oversampling, 5th order lagrange interpolation.
</audio>

<label>1倍、7次補間</label>
<audio controls>
    <source src="snd/sin_01x_r7_w7.wav" type="audio/wav">
    Audio of delay. 1x oversampling, 7th order lagrange interpolation.
</audio>

<label>1倍、9次補間</label>
<audio controls>
    <source src="snd/sin_01x_r9_w9.wav" type="audio/wav">
    Audio of delay. 1x oversampling, 9th order lagrange interpolation.
</audio>

オーバーサンプリング2倍の音です。1次補間ではノイズが大きく減ったように聞こえます。

<label>2倍、1次補間</label>
<audio controls>
    <source src="snd/sin_02x_r1_w1.wav" type="audio/wav">
    Audio of delay. 2x oversampling, 1st order lagrange interpolation.
</audio>

<label>2倍、3次補間</label>
<audio controls>
    <source src="snd/sin_02x_r3_w3.wav" type="audio/wav">
    Audio of delay. 2x oversampling, 3rd order lagrange interpolation.
</audio>

<label>2倍、5次補間</label>
<audio controls>
    <source src="snd/sin_02x_r5_w5.wav" type="audio/wav">
    Audio of delay. 2x oversampling, 5th order lagrange interpolation.
</audio>

<label>2倍、7次補間</label>
<audio controls>
    <source src="snd/sin_02x_r7_w7.wav" type="audio/wav">
    Audio of delay. 2x oversampling, 7th order lagrange interpolation.
</audio>

<label>2倍、9次補間</label>
<audio controls>
    <source src="snd/sin_02x_r9_w9.wav" type="audio/wav">
    Audio of delay. 2x oversampling, 9th order lagrange interpolation.
</audio>

オーバーサンプリング4倍の音です。5次補間以降は音量を上げるとノイズが聞き取れます。

<label>4倍、1次補間</label>
<audio controls>
    <source src="snd/sin_04x_r1_w1.wav" type="audio/wav">
    Audio of delay. 4x oversampling, 1st order lagrange interpolation.
</audio>

<label>4倍、3次補間</label>
<audio controls>
    <source src="snd/sin_04x_r3_w3.wav" type="audio/wav">
    Audio of delay. 4x oversampling, 3rd order lagrange interpolation.
</audio>

<label>4倍、5次補間</label>
<audio controls>
    <source src="snd/sin_04x_r5_w5.wav" type="audio/wav">
    Audio of delay. 4x oversampling, 5th order lagrange interpolation.
</audio>

<label>4倍、7次補間</label>
<audio controls>
    <source src="snd/sin_04x_r7_w7.wav" type="audio/wav">
    Audio of delay. 4x oversampling, 7th order lagrange interpolation.
</audio>

<label>4倍、9次補間</label>
<audio controls>
    <source src="snd/sin_04x_r9_w9.wav" type="audio/wav">
    Audio of delay. 4x oversampling, 9th order lagrange interpolation.
</audio>

オーバーサンプリング8倍の音です。ノイズは聞き取れません。

<label>8倍、1次補間</label>
<audio controls>
    <source src="snd/sin_08x_r1_w1.wav" type="audio/wav">
    Audio of delay. 8x oversampling, 1st order lagrange interpolation.
</audio>

<label>8倍、3次補間</label>
<audio controls>
    <source src="snd/sin_08x_r3_w3.wav" type="audio/wav">
    Audio of delay. 8x oversampling, 3rd order lagrange interpolation.
</audio>

<label>8倍、5次補間</label>
<audio controls>
    <source src="snd/sin_08x_r5_w5.wav" type="audio/wav">
    Audio of delay. 8x oversampling, 5th order lagrange interpolation.
</audio>

<label>8倍、7次補間</label>
<audio controls>
    <source src="snd/sin_08x_r7_w7.wav" type="audio/wav">
    Audio of delay. 8x oversampling, 7th order lagrange interpolation.
</audio>

<label>8倍、9次補間</label>
<audio controls>
    <source src="snd/sin_08x_r9_w9.wav" type="audio/wav">
    Audio of delay. 8x oversampling, 9th order lagrange interpolation.
</audio>

ディレイ時間の変更でリングバッファの内容が低速再生される状態ではどうなるか試します。

220Hzのサイン波を入力して低速再生したときの音です。音量を上げると1次補間でノイズが聞こえます。

<label>低速再生サイン波、8倍、1次補間</label>
<audio controls>
    <source src="snd/sin_slow_08x_r1_w1.wav" type="audio/wav">
    Audio of delay with low speed playback. 8x oversampling, 1st order lagrange interpolation.
</audio>

<label>低速再生サイン波、8倍、3次補間</label>
<audio controls>
    <source src="snd/sin_slow_08x_r3_w3.wav" type="audio/wav">
    Audio of delay with low speed playback. 8x oversampling, 3rd order lagrange interpolation.
</audio>

<label>低速再生サイン波、8倍、5次補間</label>
<audio controls>
    <source src="snd/sin_slow_08x_r5_w5.wav" type="audio/wav">
    Audio of delay with low speed playback. 8x oversampling, 5th order lagrange interpolation.
</audio>

<label>低速再生サイン波、8倍、7次補間</label>
<audio controls>
    <source src="snd/sin_slow_08x_r7_w7.wav" type="audio/wav">
    Audio of delay with low speed playback. 8x oversampling, 7th order lagrange interpolation.
</audio>

<label>低速再生サイン波、8倍、9次補間</label>
<audio controls>
    <source src="snd/sin_slow_08x_r9_w9.wav" type="audio/wav">
    Audio of delay with low speed playback. 8x oversampling, 9th order lagrange interpolation.
</audio>

`Math.random()` で生成したノイズを入力します。

<label>ノイズ</label>
<audio controls>
    <source src="snd/noise_source.wav" type="audio/wav">
    Audio of 220Hz sine wave.
</audio>

ノイズを入力して低速再生したときの音です。補間の次数が高くなるほどローパスフィルタがかかったような音に近くなります。

<label>低速再生ノイズ、16倍、1次補間</label>
<audio controls>
    <source src="snd/noise_16x_r1_w1.wav" type="audio/wav">
    Audio of delay with low speed playback. 16x oversampling, 1st order lagrange interpolation.
</audio>

<label>低速再生ノイズ、16倍、3次補間</label>
<audio controls>
    <source src="snd/noise_16x_r3_w3.wav" type="audio/wav">
    Audio of delay with low speed playback. 16x oversampling, 3rd order lagrange interpolation.
</audio>

<label>低速再生ノイズ、16倍、5次補間</label>
<audio controls>
    <source src="snd/noise_16x_r5_w5.wav" type="audio/wav">
    Audio of delay with low speed playback. 16x oversampling, 5th order lagrange interpolation.
</audio>

<label>低速再生ノイズ、16倍、7次補間</label>
<audio controls>
    <source src="snd/noise_16x_r7_w7.wav" type="audio/wav">
    Audio of delay with low speed playback. 16x oversampling, 7th order lagrange interpolation.
</audio>

<label>低速再生ノイズ、16倍、9次補間</label>
<audio controls>
    <source src="snd/noise_16x_r9_w9.wav" type="audio/wav">
    Audio of delay with low speed playback. 16x oversampling, 9th order lagrange interpolation.
</audio>

レンダリングした音を聞く限りでは、高次の補間を使うときはオーバーサンプリングの倍率を補間の次数以上にするとノイズが大きく減るようです。

ディレイ時間があまり変更されないときは、1次補間で2倍のオーバーサンプリングにすれば十分な気がします。テープエコーのシミュレーションなどでリングバッファ内の音を低速再生する状態を作るときは高次の補間のほうが良さそうです。

## その他
1次のラグランジュ補間は線形補間です。 $k$ 次のラグランジュ補間の定義です。

$$
L_k(x) = \sum_{j=0}^{k} y_j \prod_{m=0,\,m \neq j}^{k} \frac{x - x_m}{x_j - x_m}
$$

$k = 1$ を代入します。

$$
\begin{aligned}
L_1(x)
&= \sum_{j=0}^{1} y_j \prod_{m=0,\,m \neq j}^{1} \frac{x - x_m}{x_j - x_m}\\
&= y_0 \left(
    \prod_{m=0,\,m \neq 0}^{1} \frac{x - x_m}{x_0 - x_m}
  \right)
  + y_1 \left(
    \prod_{m=0,\,m \neq 1}^{1} \frac{x - x_m}{x_1 - x_m}
  \right)\\
&= y_0 \frac{x - x_1}{x_0 - x_1}
  + y_1 \frac{x - x_0}{x_1 - x_0}
\end{aligned}
$$

ここで $x$ が等間隔にサンプリングされているとき $x_0 = 0,\,x_1 = 1$ と代入できます。

$$
y_0 \frac{x - x_1}{x_0 - x_1} + y_1 \frac{x - x_0}{x_1 - x_0}
= - y_0 (x - 1) + y_1 x
$$

$0 \leq x \leq 1$ で $y_0$ と $y_1$ の間を線形補間する式になっています。

## 参考文献
- [MUS420 Lecture 4A Interpolated Delay Lines, Ideal Bandlimited Interpolation, and Fractional Delay Filter Design](https://ccrma.stanford.edu/~jos/Interpolation/)
