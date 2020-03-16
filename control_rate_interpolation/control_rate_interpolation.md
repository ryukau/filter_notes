# オーディオプラグインの UI から入力された値の補間
UI からプラグインへと投げられてきたパラメータの補間について見ていきます。

## 素朴なパラメータの受け取り
VST 3 や LV2 では、オーディオデバイスのバッファの長さ $L$ サンプルごとにオーディオ処理のメソッドが呼び出されます。

オーディオ処理のメソッドが呼び出される頻度のことをコントロールレートと呼ぶことがあります。サンプリング周波数を $f_s$ とするとコントロールレートは $f_s / L$ です。

次のコードはバッファ単位で `process` を呼び出す例です。バッファに書き込み直す処理を `main` 内の `for` ループで回していますが、実際にはオーディオデバイスによる割り込みに応じて DAW のオーディオスレッドから呼び出されます。

```cpp
// naive.cpp
#include <cstring>
#include <random>
#include <vector>

constexpr float sampleRate = 48000.0f;
constexpr size_t nFrame = 512; // 式中の L 。
constexpr size_t nBuffer = 32;

struct DSP {
  float gain = 1.0f; // 適当なパラメータ。

  void process(const float frame, float *out) // オーディオ処理のメソッド
  {
    for (size_t i = 0; i < frame; ++i) out[i] = gain;
  }
};

int main()
{
  std::vector<float> wav(nFrame * nBuffer);

  float out[nFrame]; // VST 3 や LV2 にならってバッファに C の配列を使う。

  std::minstd_rand rng{0};
  std::uniform_real_distribution<float> gainDist(0.0f, 1.0f);
  DSP dsp;

  float gainValue = 1.0f;
  for (size_t idx = 0; idx < nBuffer; ++idx) {
    if (idx % 4 == 0) gainValue = gainDist(rng); // パラメータの更新頻度を下げる。
    dsp.gain = gainValue;
    dsp.process(nFrame, out);
    std::memcpy(&wav[idx * nFrame], out, sizeof(float) * nFrame);
  }

  return writeWave("snd/naive.wav", wav, sampleRate);
}
```

`writeWave` の実装を含む完成したコードは次のリンクに掲載しています。コンパイルには [libsndfile]([libsndfile](http://www.mega-nerd.com/libsndfile/)) が必要です。

- [filter_notes/naive.cpp at master · ryukau/filter_notes · GitHub](https://github.com/ryukau/filter_notes/blob/master/control_rate_interpolation/demo/naive.cpp)

コンパイルして実行します。

```bash
g++ -lsndfile -O3 naive.cpp
./a.out
```

実行結果の波形です。

<figure>
<img src="img/naive.png" alt="Image of naively received parameter." style="padding-bottom: 12px;"/>
</figure>

図の波形に 100 Hz の正弦波を掛け合わせた音です。

<figure>
  <figcaption>Naive</figcaption>
  <audio controls>
    <source src="snd/naive.wav" type="audio/wav">
  </audio>
</figure>

バッファの境界で値が不連続に変化しているので、そのまま使うとノイズが乗ります。

不連続点が目立ってしまうパラメータの例としては、音量、フィルタのカットオフ周波数、ディレイ時間、オシレータの周波数などが挙げられます。特に IIR フィルタのカットオフ周波数を急激に変更すると発散することがあるので、パラメータを補間せずに使うことは避けたいです。

## バッファ内で線形補間
サンプリング周波数とバッファサイズが固定なら、パラメータが変化したときだけ次のバッファまで線形補間することでノイズが大きく減ります。

デスクトップ環境では設定によってサンプリング周波数とバッファサイズが変わります。そこで、古い値から新しい値へと更新される時間を指定できるようにします。

まずはサンプリング周波数 $f_s$ とユーザが指定した新しい値へと更新される秒数 $t$ から、新しい値へと更新されるサンプル数 $n$ を求めます。

$$
n = t f_s
$$

次に、新しく受け取ったパラメータの値 $p$ 、 1 つ前のバッファの終端でのパラメータの値 $p_1$ 、 バッファの長さ $L$ を使って、現在のバッファの終端でのパラメータの値 $p_0$ を計算します。

$$
p_0 = \frac{L}{n} (p - p_1) + p_1
$$

そして、バッファ内でのパラメータの値 $p_i$ を、バッファの先頭から経過したサンプル数 $i$ を使って計算します。この式は $p_1$ と $p_0$ の線形補間です。

$$
p_i = p_1 + \frac{i}{L} (p_0 - p_1)
$$

この計算は指定した時間で目的の値に到達しませんが、シンセサイザで使うには十分です。もう少し適切に指定したいときはパラメータチューニングの項を参照してください。

実装します。

```cpp
constexpr float sampleRate = 48000.0f;
constexpr size_t nFrame = 512;
constexpr size_t nBuffer = 32;

template<typename Sample> class Smoother {
public:
  void setSampleRate(Sample sampleRate, Sample time = 0.04)
  {
    this->sampleRate = sampleRate;
    setTime(time);
  }

  void setTime(Sample seconds) { timeInSamples = seconds * sampleRate; }
  void setBufferSize(Sample bufferSize) { this->bufferSize = bufferSize; }
  inline Sample getValue() { return value; }

  void push(Sample newTarget) // newTarget は式中の p 。
  {
    p1 = p0;
    p0 = (timeInSamples >= bufferSize) && (fabs(p0 - newTarget) >= Sample(1e-5))
      ? (newTarget - p0) * bufferSize / timeInSamples + p0
      : newTarget;
  }

  Sample process(float index) { return value = p1 + index / bufferSize * (p0 - p1); }

protected:
  Sample sampleRate = 44100;
  Sample timeInSamples = -1; // 式中の n 。
  Sample bufferSize = 0;     // 式中の L 。
  Sample p0 = 1;
  Sample p1 = 1;
  Sample value = 0;          // 1 フレーム内で計算結果を 2 回以上使うときに必要。
};

struct DSP {
  Smoother<float> gain;

  void process(const size_t frame, float *out)
  {
    for (size_t i = 0; i < frame; ++i) out[i] = gain.process(i);
  }
};

int main()
{
  std::vector<float> wav(nFrame * nBuffer);

  float out[nFrame];

  std::minstd_rand rng{0};
  std::uniform_real_distribution<float> gainDist(0.0f, 1.0f);
  DSP dsp;
  dsp.gain.setSampleRate(sampleRate, 0.02f);

  float gainValue = 1.0f;
  for (size_t idx = 0; idx < nBuffer; ++idx) {
    if (idx % 4 == 0) gainValue = gainDist(rng); // パラメータの更新頻度を下げる。
    dsp.gain.setBufferSize(nFrame);
    dsp.gain.push(gainValue);
    dsp.process(nFrame, out);
    std::memcpy(&wav[idx * nFrame], out, sizeof(float) * nFrame);
  }

  return writeWave("snd/smoother.wav", wav, sampleRate);
}
```

完成したコードへのリンクです。

- [filter_notes/smoother.cpp at master · ryukau/filter_notes · GitHub](https://github.com/ryukau/filter_notes/blob/master/control_rate_interpolation/demo/smoother.cpp)

実行結果です。

<figure>
<img src="img/smoother.png" alt="Image of smoothed parameter." style="padding-bottom: 12px;"/>
</figure>

図の波形に 100 Hz のサイン波を掛け合わせた音です。素朴な実装と比べるとプチノイズが減っています。

<figure>
  <figcaption>Smoother</figcaption>
  <audio controls>
    <source src="snd/smoother.wav" type="audio/wav">
  </audio>
</figure>

比較のために素朴な実装の結果を再掲します。

<figure>
<img src="img/naive.png" alt="Image of naively received parameter." style="padding-bottom: 12px;"/>
</figure>

<figure>
  <figcaption>Naive</figcaption>
  <audio controls>
    <source src="snd/naive.wav" type="audio/wav">
  </audio>
</figure>

パラメータが変更されたバッファのみ線形補間したときの結果も掲載します。手間を省きたいときは、これで十分です。

<figure>
<img src="img/linterp.png" alt="Image of linear interpolated parameter." style="width: 480px;padding-bottom: 12px;"/>
</figure>

<figure>
  <figcaption>Linterp</figcaption>
  <audio controls>
    <source src="snd/linterp.wav" type="audio/wav">
  </audio>
</figure>

コードへのリンク: [filter_notes/linterp.cpp at master · ryukau/filter_notes · GitHub](https://github.com/ryukau/filter_notes/blob/master/control_rate_interpolation/demo/linterp.cpp)


### パラメータチューニング
$p_0$ の振る舞いは指数関数的減衰 (exponential decay) になります。経過したバッファ数を $j$ として、 $L, n, p$ を固定します。さらに $\alpha = L / n,\ p_1 = 0$ とすれば次のように書き直せます。

$$
p_0 = \alpha^j p
$$

$\alpha^j$ が適当なしきい値 $\epsilon$ に到達する時間を求めます。

$$
\begin{aligned}
\epsilon &= \alpha^j \\
j &= \frac{\log \epsilon}{\log \alpha}
\end{aligned}
$$

また $\epsilon$ と $j$ が分かっているとき、 $\alpha = \epsilon^{\frac{1}{j}}$ です。

## 異なる計算方法
Uhhyou Plugins では以前、異なる計算方法を使っていたので記録しておきます。

`Smoother::process()` での計算式を再掲します。

```cpp
value = p1 + index / bufferSize * (p0 - p1);
```

`index` はバッファ内のインデックスです。 1 サンプルあたりに `value` が増える量を事前に計算しておけば、バッファ内のインデックスを使わずに計算できます。

```cpp
template<typename Sample> class SmootherRamp {
public:
  void setSampleRate(Sample sampleRate, Sample time = 0.04)
  {
    this->sampleRate = sampleRate;
    timeInSamples = seconds * sampleRate;
  }

  void setBufferSize(Sample bufferSize) { this->bufferSize = bufferSize; }

  void push(Sample newTarget)
  {
    target = newTarget;
    if (timeInSamples < bufferSize)
      value = target;
    else
      ramp = (target - value) / timeInSamples;
  }

  Sample process()
  {
    if (value == target) return value;
    value += ramp;
    if (fabs(value - target) < 1e-5) value = target;
    return value;
  }

protected:
  Sample sampleRate = 44100;
  Sample timeInSamples = -1;
  Sample bufferSize = 0;
  Sample value = 1.0;
  Sample target = 1.0;
  Sample ramp = 0.0; // 1 サンプルあたりに value が増える量。
};
```

指定されたスムーシング時間がバッファ長よりも短いときは `push()` の条件分岐によってスムーシングを無効にしています。

`process()` の条件分岐 `fabs(value - target) < 1e-5` でオーバーシュートを防いでいます。 `1e-5` は適当に決めた値です。パラメータの値の範囲が広いときは `1e-5` をより大きな値に変えないと上手く動かない可能性があります。

1 サンプルあたりの増加量を使う方法は、バッファ内のインデックスを使う方法よりも遅いです。スムーシング時間が固定で、オーバーシュートがあってもいいなら `process()` 内の条件分岐を消せるので使えるかもしれません。

ベンチマークに使ったコードを次のリンクに掲載しています。

- [filter_notes/bench_smoother.cpp at master · ryukau/filter_notes · GitHub](https://github.com/ryukau/filter_notes/blob/master/control_rate_interpolation/demo/bench_smoother.cpp)
