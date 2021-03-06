# 直線の ADSR エンベロープ
シンセサイザに使う直線を組み合わせた ADSR エンベロープを実装します。計算はシンプルですが、場合分けの塊なのでテストのポイントをまとめています。

## ADSR エンベロープ
シンセサイザでよく使われるエンベロープは、トリガされると low から high に向かって増加したあと high から low に向かって減少するような信号を出力します。例えばエンベロープの出力をオシレータの出力やサンプルデータと掛け合わせることで音量の変化を作ることができます。

よくあるシンセサイザのエンベロープはアタック、ディケイ、サステイン、リリース (Attack, Decay, Sustain, Release) の4つの区間が組み合わさっています。4つの区間をそれぞれ1文字に略して ADSR と言うことがあります。次の図は ADSR エンベロープを表しています。

<figure>
<img src="img/adsr_envelope.svg" alt="Image of ADOsc audio graph." style="padding-bottom: 12px;"/>
</figure>

アタックは low から high までの立ち上がりにかかる時間、 ディケイは high からサステインの値までの減衰にかかる時間です。

サステインは、ディケイの終わりから指が鍵盤が離される(ノートオフ)までの区間ですが、ユーザが操作するパラメータとしてはその区間でエンベロープが出力する値の大きさになっています。

リリースはノートオフの後、エンベロープが low まで減衰するのにかかる時間です。アタックやディケイの途中でノートオフされたときもリリースに移行します。

## 直線の計算
一つの区間での直線の計算の実装例です。出力の範囲は [0, 1] です。移動距離 1 と、区間の長さ `seconds` から、速度 `delta` [gain/sample] を計算しています。距離の単位はメートルなどではない抽象的な値なので適当に gain としています。

```cpp
struct LinearAttack {
  float delta;
  float value = 0;

  LinearAttack(float sampleRate, float seconds) : delta(1.0f / (seconds * sampleRate)) {}

  float process()
  {
    value += delta;
    if (value > 1.0f) return 1.0f;
    return value;
  }
};
```

`process` の出力を表した図です。

<figure>
<img src="img/line.svg" alt="Image of a section of line." style="padding-bottom: 12px;"/>
</figure>

## 挙動
### ノートオン中の ADSR の変更
ノートオン中に ADSR のパラメータを変更可能にするかどうかの選択があります。パラメータを変更可能にしたほうが便利ですが、傾きを再計算するコストがかかります。全てのサンプルで傾きを再計算すると無視できないほど重くなるので、この文章ではバッファの先頭でのみ傾きの再計算を行うようにしています。

リリースについてはサステインが変更されている途中でノートオフされると、ノートオフの時点でさらに傾きを計算し直すかどうかの選択があります。今回の実装では傾きを再計算しています。

### アタック中のリリース
ノートオフの時点の高さから減衰を始めます。

<figure>
<img src="img/release_while_attack.svg" alt="Image of envelope receiving note-off while attack." style="padding-bottom: 12px;"/>
</figure>

### ディケイ中のリリース
ノートオフの時点の高さから減衰を始めます。

<figure>
<img src="img/release_while_decay.svg" alt="Image of envelope receiving note-off while decay." style="padding-bottom: 12px;"/>
</figure>


### リリース中のノートオン
左の図のようにノートオンの時点での高さを維持してアタックに移る挙動と、右の図のようにノートオンごとに高さを 0 にリセットする挙動の 2 種類が考えられます。

<figure>
<img src="img/trigger_while_release.svg" alt="Image of 2 variation of envelope when receiving note-on while release." style="padding-bottom: 12px;"/>
</figure>

特殊な場合を除けば、高さを維持する挙動を使うほうがポップノイズが減ります。 [CubicPadSynth](https://ryukau.github.io/VSTPlugins/manual/CubicPadSynth/CubicPadSynth_ja.html) では、ノートオンを受け取ったボイスから、ポップノイズを減らすために用意された専用のボイスにパラメータをコピーして短いリリースをレンダリングすることでポップノイズを低減しているので、音量エンベロープに高さをリセットする挙動を使っています。

## 実装
C++ です。 `std::clamp` を使っているので C++17 以降で動作します。

```cpp
// linear_envelope.hpp
#include <algorithm>
#include <cstdint>

template<typename Sample> class LinearADSREnvelope {
public:
  void setup(Sample sampleRate) { this->sampleRate = sampleRate; }

  Sample adaptTime(Sample seconds, Sample noteFreq)
  {
    const Sample cycle = Sample(1) / noteFreq;
    return seconds >= cycle ? seconds : cycle > Sample(0.1) ? Sample(0.1) : cycle;
  }

  Sample secondToDelta(Sample seconds) { return Sample(1) / (sampleRate * seconds); }

  void trigger(
    Sample attackTime,
    Sample decayTime,
    Sample sustainLevel,
    Sample releaseTime,
    Sample noteFreq)
  {
    state = stateAttack;
    value = Sample(1);
    atkOffset = state == stateTerminated ? 0 : out;
    atkRange = Sample(1) - atkOffset;
    set(attackTime, decayTime, sustainLevel, releaseTime, noteFreq);
  }

  void set(
    Sample attackTime,
    Sample decayTime,
    Sample sustainLevel,
    Sample releaseTime,
    Sample noteFreq)
  {
    sus = std::clamp<Sample>(sustainLevel, Sample(0), Sample(1));
    atk = secondToDelta(adaptTime(attackTime, noteFreq));
    dec = secondToDelta(adaptTime(decayTime, noteFreq));
    rel = secondToDelta(adaptTime(releaseTime, noteFreq));
  }

  void release()
  {
    state = stateRelease;
    value = Sample(1);
    relRange = out;
  }

  Sample process()
  {
    if (value <= Sample(0)) {
      state = state + 1;
      value = Sample(1);
    }

    switch (state) {
      case stateAttack:
        value -= atk;
        out = atkOffset + atkRange * (Sample(1) - value);
        break;

      case stateDecay:
        value -= dec;
        out = (Sample(1) - sus) * value + sus;
        break;

      case stateSustain:
        out = sus;
        break;

      case stateRelease:
        value -= rel;
        out = relRange * value;
        break;

      default:
        return Sample(0);
    }
    return std::clamp<Sample>(out, Sample(0), Sample(1));
  }

protected:
  enum State : int32_t {
    stateAttack,
    stateDecay,
    stateSustain,
    stateRelease,
    stateTerminated
  };

  Sample sampleRate = 44100;
  Sample sus = 0.5;
  Sample atk = 0.01;
  Sample atkOffset = 0;
  Sample atkRange = 1;
  Sample dec = 0.01;
  Sample rel = 0.01;
  Sample relRange = 0.5;
  Sample value = 0;
  Sample out = 0;
  int32_t state = stateTerminated;
};
```

`adaptTime` は信号の周波数に応じて区間の時間を少し長くすることで、ポップノイズを低減するメソッドです。

`secondToDelta` は秒数から 1 サンプルあたりのエンベロープの値の変化量を計算します。各区間の内部的な出力は [0, 1] の範囲に正規化されています。 `process` から出力されるときに適切な値にマッピングしています。

`trigger` の引数の `*Time` の単位は秒です。リリース中のノートオンは高さを維持する挙動にしています。

`set` は UI で変更されたパラメータを更新するために呼び出します。今回の実装ではサステイン中に `sustainLevel` を変更するとノイズが乗ります。ノイズを減らすときは次のページで紹介している `Smoother` クラスが使えます。

- [オーディオプラグインの UI から入力された値の補間](https://ryukau.github.io/filter_notes/control_rate_interpolation/control_rate_interpolation.html)

`process` はエンベロープの出力を計算します。出力が [0, 1] の範囲から少しだけはみ出してもいいなら `clamp` を消すことができます。このとき、はみ出す量の絶対値の上限は `1 / (sampleRate * seconds)` です。

## テスト
wav ファイルへの書き出しに [libsndfile](http://www.mega-nerd.com/libsndfile/) を使っています。

```cpp
// test.cpp
#include "linear_envelope.hpp"

#include <cmath>
#include <iostream>
#include <sndfile.h>
#include <stdlib.h>
#include <string.h>
#include <vector>

using Sample = float;

constexpr Sample sampleRate = 44100.0;

int32_t
writeWave(const char *filename, std::vector<float> &buffer, const size_t &samplerate)
{
  SF_INFO sfinfo;
  memset(&sfinfo, 0, sizeof(sfinfo));
  sfinfo.samplerate = samplerate;
  sfinfo.frames = buffer.size();
  sfinfo.channels = 1;
  sfinfo.format = (SF_FORMAT_WAV | SF_FORMAT_FLOAT);

  SNDFILE *file = sf_open(filename, SFM_WRITE, &sfinfo);
  if (!file) {
    std::cout << "Error: sf_open failed." << std::endl;
    return 1;
  }

  size_t length = sfinfo.channels * buffer.size();
  if (sf_write_float(file, &buffer[0], length) != length)
    std::cout << sf_strerror(file) << std::endl;

  sf_close(file);

  return 0;
}

void testADSR(
  float attackTime = 1.0f,
  float decayTime = 1.0f,
  float sustainLevel = 0.5f,
  float releaseTime = 2.0f,
  float releaseAt = 3.0f,
  float noteFreq = 100.0f)
{
  LinearADSREnvelope<float> envelope;
  envelope.setup(sampleRate);
  envelope.trigger(attackTime, decayTime, sustainLevel, releaseTime, noteFreq);

  std::vector<float> buffer(size_t(8 * sampleRate));
  for (size_t idx = 0; idx < buffer.size(); ++idx) {
    if (idx == size_t(releaseAt * sampleRate)) envelope.release();
    buffer[idx] = envelope.process();
  }
  writeWave("ADSR.wav", buffer, sampleRate);
}

void testReleaseWhileAttack()
{
  LinearADSREnvelope<float> envelope;
  envelope.setup(sampleRate);
  envelope.trigger(4.0f, 1.0f, 0.5f, 2.0f, 100.0f);

  std::vector<float> buffer(size_t(8 * sampleRate));
  for (size_t idx = 0; idx < buffer.size(); ++idx) {
    if (idx == size_t(3.0f * sampleRate)) envelope.release();
    buffer[idx] = envelope.process();
  }
  writeWave("ReleaseWhileAttack.wav", buffer, sampleRate);
}

void testReleaseWhileDecay()
{
  LinearADSREnvelope<float> envelope;
  envelope.setup(sampleRate);
  envelope.trigger(1.0f, 1.0f, 0.5f, 0.5f, 100.0f);

  std::vector<float> buffer(size_t(8 * sampleRate));
  for (size_t idx = 0; idx < buffer.size(); ++idx) {
    if (idx == size_t(1.5f * sampleRate)) envelope.release();
    buffer[idx] = envelope.process();
  }
  writeWave("ReleaseWhileDecay.wav", buffer, sampleRate);
}

void testTriggerWhileRelease()
{
  LinearADSREnvelope<float> envelope;
  envelope.setup(sampleRate);
  envelope.trigger(1.0f, 1.0f, 0.5f, 2.0f, 100.0f);

  std::vector<float> buffer(size_t(8 * sampleRate));
  for (size_t idx = 0; idx < buffer.size(); ++idx) {
    if (idx == size_t(3.0f * sampleRate)) envelope.release();
    if (idx == size_t(4.0f * sampleRate))
      envelope.trigger(0.5f, 0.5f, 0.3f, 2.0f, 100.0f);
    if (idx == size_t(6.0f * sampleRate)) envelope.release();
    buffer[idx] = envelope.process();
  }
  writeWave("TriggerWhileRelease.wav", buffer, sampleRate);
}

int main()
{
  testADSR();
  testReleaseWhileAttack();
  testReleaseWhileDecay();
  testTriggerWhileRelease();
  return 0;
}
```

ファイルの配置です。

```bash
$ tree
.
├── linear_envelope.hpp
└── test.cpp
```

コンパイルして実行します。

```bash
$ g++ -std=c++17 -lsndfile -O3 test.cpp
$ ./a.out
```

結果です。この画像は matplotlib でプロットしています。縦軸は振幅、横軸は時間で単位は秒数です。大まかなチェックなら [Audacity](https://www.audacityteam.org/) や [Sonic Visualiser](https://sonicvisualiser.org/) も使えますが、どちらも拡大縮小がしづらいです。

<figure>
<img src="img/testADSR.png" alt="Image of test result of linear ADSR envelope." style="padding-bottom: 12px;"/>
</figure>

<figure>
<img src="img/testReleaseWhileAttack.png" alt="Image of test result of release while attack." style="padding-bottom: 12px;"/>
</figure>

<figure>
<img src="img/testReleaseWhileDecay.png" alt="Image of test result of release while decay." style="padding-bottom: 12px;"/>
</figure>

<figure>
<img src="img/testTriggerWhileRelease.png" alt="Image of test result of trigger while release." style="padding-bottom: 12px;"/>
</figure>

## SIMD による並列化
CubicPadSynth に [vector class library](https://github.com/vectorclass/version2) を使って並列化した実装があります。この実装ではリリース時に傾きの再計算をしていません。

- [VSTPlugins/envelope.hpp at master · ryukau/VSTPlugins · GitHub](https://github.com/ryukau/VSTPlugins/blob/master/CubicPadSynth/source/dsp/envelope.hpp#L155)

## 記事の変更点
- 2020-03-23:
  - `P Controller` の節を `オーディオプラグインの UI から入力された値の補間` に移動。
  - `カットオフ周波数から $k_p$ を求める` で式の求め方が間違っていたので修正。
  - `複素数を使わずにカットオフ周波数から $k_p$ 求める式` を追加。
