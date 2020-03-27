# 指数曲線のエンベロープ
指数曲線 (Exponential Curve) を使ったエンベロープを作ります。

## 減衰する指数曲線
値が 1 から 0 に向かって減衰する指数曲線 $E_d(t)$ の式です。

$$
E_d(t) = \alpha^t, \quad 0 \leq \alpha \leq 1
$$

$\alpha$ は減衰の速さを決める任意の値、 $t$ は単位が秒数の時間です。

ユーザから指定された減衰時間 $\tau$ から $\alpha$ を決めます。 $E_d$ は $0 < \alpha$ のとき $t = +\infty$ でようやく 0 になります。言い換えると $E_d$ はいつまで経っても 0 になりません。そこで 0 の代わりに十分に小さな値 $\epsilon$ に到達する時間を求めます。

$$
E_d(\tau) = \alpha^\tau = \epsilon
$$

時間の単位を秒数 $\tau$ からサンプル数 $n_\tau$ に置き換えます。

$$
n_\tau = \tau f_s
$$

$f_s$ はサンプリング周波数です。

$\tau$ と $\epsilon$ が与えられたとき $\alpha^\tau = \epsilon$ の関係と $n_\tau$ の式より $\alpha$ が求められます。

$$
\alpha = \epsilon^\Tau, \quad \Tau = \frac{1}{\tau f_s}.
$$

$\epsilon$ の値は任意です。この文章では `1e-5` を使っています。値を決めるときの参考までに、 10 進数では `float` は約 7 桁、 `double` は約 16 桁の精度があります。

- [floating point - How many significant digits do floats and doubles have in java? - Stack Overflow](https://stackoverflow.com/questions/13542944/how-many-significant-digits-do-floats-and-doubles-have-in-java)

実装します。

```cpp
#include <cmath>

template<typename Sample> class ExpDecayCurve {
public:
  void reset(Sample sampleRate, Sample seconds)
  {
    value = Sample(1);
    set(sampleRate, seconds);
  }

  void set(Sample sampleRate, Sample seconds)
  {
    alpha = pow(threshold, Sample(1) / (seconds * sampleRate));
  }

  bool isTerminated() { return value <= threshold; }

  Sample process()
  {
    if (value <= threshold) return Sample(0);
    value *= alpha;
    return value - threshold;
  }

protected:
  const Sample threshold = 1e-5;
  Sample value = 0;
  Sample alpha = 0;
};
```

`reset()` は再トリガ時、 `set()` はコントロールレートで呼び出されます。

`isTerminated()` は複数の曲線を組み合わせて ADSR エンベロープなどを作るとき、状態遷移を行うために使います。

`process()` の最後で `threshold` を引いているのは、エンベロープの終端の不連続点を小さくするためです。この処理によって出力の範囲は `[0, 1 - threshold]` になります。

<figure>
<img src="img/discontinuity_at_termination.svg" alt="Image of ." style="padding-bottom: 12px;"/>
</figure>

値の範囲が `[0, 1]` かつ終端に崖を作りたくないときはルックアップテーブルが使えます。ただし、滑らかさは劣ります。

## 増加する指数曲線
増加する指数曲線 $E_a(t)$ の式です。

$$
E_d(\tau) = \epsilon \alpha^\tau = 1
$$

$\alpha$ について解きます。

$$
\alpha = \left( \frac{1}{\epsilon} \right)^\Tau, \quad \Tau = \frac{1}{\tau f_s}.
$$

実装します。

```cpp
template<typename Sample> class ExpAttackCurve {
public:
  void reset(Sample sampleRate, Sample seconds)
  {
    value = threshold;
    set(sampleRate, seconds);
  }

  void set(Sample sampleRate, Sample seconds)
  {
    alpha
      = pow(Sample(1) / threshold, Sample(1) / (seconds * sampleRate));
  }

  bool isTerminated() { return value >= Sample(1); }

  Sample process()
  {
    value *= alpha;
    if (value >= Sample(1)) return Sample(1 - threshold);
    return value - threshold;
  }

protected:
  const Sample threshold = 1e-5;
  Sample value = 0;
  Sample alpha = 0;
};
```

## 減衰する指数曲線の反転
$1 - E_d(t)$ として減衰する指数曲線を上下反転して使う方法があります。

実装です。

```cpp
// 1 - ExpDecayCurve.process();
template<typename Sample> class NegativeExpAttackCurve {
public:
  void reset(Sample sampleRate, Sample seconds)
  {
    value = Sample(1);
    set(sampleRate, seconds);
  }

  void set(Sample sampleRate, Sample seconds)
  {
    alpha = pow(threshold, Sample(1) / (seconds * sampleRate));
  }

  bool isTerminated() { return value <= threshold; }

  Sample process()
  {
    if (value <= threshold) return Sample(1 - threshold);
    value *= alpha;
    return Sample(1 - threshold) - value;
  }

protected:
  const Sample threshold = 1e-5;
  Sample value = 0;
  Sample alpha = 0;
};
```

## 指数曲線を使った ADSR エンベロープ
C++ での実装例です。

このエンベロープは [IterativeSinCluster](https://ryukau.github.io/VSTPlugins/manual/IterativeSinCluster/IterativeSinCluster_ja.html) で使われています。この文章に掲載している実装では `sustain` をバッファ内で補間していないので、サステイン中にサステイン音量を変更するとノイズが乗ります。また、エンベロープが終了していない状態で再トリガすると、デクリックの副エンベロープが掛け合わされることによってノイズが乗ります。

`reset()` の `curve` の値で `ExpAttackCurve` と `NegativeExpAttackCurve` を入れ替えられるようにしています。 `curve` の値はトリガ時に渡された値をエンベロープの終了までに使い続けることを想定しています。

```cpp
#include <algorithm>
#include <cmath>

// t in [0, 1].
template<typename Sample> inline Sample cosinterp(Sample t)
{
  return 0.5 * (1.0 - cos(pi * t));
}

template<typename Sample> class ExpADSREnvelope {
public:
  void setup(Sample sampleRate)
  {
    this->sampleRate = sampleRate;
    declickLength = int32_t(0.001 * sampleRate);
  }

  Sample adaptTime(Sample seconds, Sample noteFreq)
  {
    const Sample cycle = Sample(1) / noteFreq;
    return seconds < cycle ? cycle : seconds;
  }

  void reset(
    Sample attackTime,
    Sample decayTime,
    Sample sustainLevel,
    Sample releaseTime,
    Sample noteFreq,
    Sample curve)
  {
    if (declickCounter >= declickLength || state == State::terminated) declickCounter = 0;
    state = State::attack;

    sustain = std::clamp<Sample>(sustainLevel, Sample(0), Sample(1));

    offset = value;
    range = Sample(1) - value;

    this->curve = std::clamp<Sample>(curve, Sample(0), Sample(1));

    attackTime = adaptTime(attackTime, noteFreq);
    atk.reset(sampleRate, attackTime);
    atkNeg.reset(sampleRate, attackTime);
    dec.reset(sampleRate, decayTime);
    rel.reset(sampleRate, adaptTime(releaseTime, noteFreq));
  }

  void set(
    Sample attackTime,
    Sample decayTime,
    Sample sustainLevel,
    Sample releaseTime,
    Sample noteFreq)
  {
    switch (state) {
      case State::attack:
        attackTime = adaptTime(attackTime, noteFreq);
        atk.set(sampleRate, attackTime);
        atkNeg.set(sampleRate, attackTime);
        // Fall through.

      case State::decay:
        dec.set(sampleRate, decayTime);
        // Fall through.

      case State::sustain:
        sustain = std::clamp<Sample>(sustainLevel, Sample(0), Sample(1));
        // Fall through.

      case State::release:
        rel.set(sampleRate, adaptTime(releaseTime, noteFreq));
        // Fall through.

      default:
        break;
    }
  }

  void release()
  {
    range = value;
    state = State::release;
  }

  bool isAttacking() { return state == State::attack; }
  bool isReleasing() { return state == State::release; }
  bool isTerminated() { return state == State::terminated; }

  inline Sample declickIn(Sample input)
  {
    if (declickCounter >= declickLength) return input;
    declickCounter += 1;
    return input * cosinterp<Sample>(declickCounter / Sample(declickLength));
  }

  Sample process()
  {
    switch (state) {
      case State::attack: {
        const auto atkPos = atk.process();
        const auto atkMix = atkPos + curve * (atkNeg.process() - atkPos);
        value = range * declickIn(atkMix) + offset;
        if (atk.isTerminated()) {
          state = State::decay;
          range = Sample(1) - sustain;
        }
      } break;

      case State::decay:
        value = range * declickIn(dec.process()) + sustain;
        if (value <= sustain) state = State::sustain;
        break;

      case State::sustain:
        value = declickIn(sustain);
        break;

      case State::release:
        value = range * declickIn(rel.process());
        if (rel.isTerminated()) state = State::terminated;
        break;

      default:
        return 0;
    }
    return value;
  }

protected:
  enum class State : int32_t { attack, decay, sustain, release, terminated };

  int32_t declickLength;
  int32_t declickCounter = 0;

  ExpAttackCurve<Sample> atk{};
  NegativeExpAttackCurve<Sample> atkNeg{};
  ExpDecayCurve<Sample> dec{};
  ExpDecayCurve<Sample> rel{};

  State state = State::terminated;
  Sample value = 0;
  Sample curve = 0;
  Sample sampleRate = 44100;
  Sample offset = 0;
  Sample range = 1;
  Sample sustain = 1;
};
```

テストに使ったコードへのリンクです。

- [filter_notes/test.cpp at master · ryukau/filter_notes · GitHub](https://github.com/ryukau/filter_notes/blob/master/exponential_envelope/demo/test.cpp)

テスト結果です。図の縦軸は振幅、横軸は秒数です。音のサンプルはエンベロープを 100 Hz のサイン波の音量に適用しています。

<figure>
<img src="img/testADSR.png" alt="Image of test result of exponential ADSR envelope." style="padding-bottom: 12px;"/>
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

<figure>
<img src="img/testChangeSustain.png" alt="Image of test result of changing sustain." style="padding-bottom: 12px;"/>
</figure>

<figure>
  <figcaption>ADSR</figcaption>
  <audio controls>
    <source src="snd/tone_ADSR.wav" type="audio/wav">
  </audio>
</figure>

<figure>
  <figcaption>ReleaseWhileAttack</figcaption>
  <audio controls>
    <source src="snd/tone_ReleaseWhileAttack.wav" type="audio/wav">
  </audio>
</figure>

<figure>
  <figcaption>testReleaseWhileDecay</figcaption>
  <audio controls>
    <source src="snd/tone_ReleaseWhileDecay.wav" type="audio/wav">
  </audio>
</figure>

<figure>
  <figcaption>TriggerWhileRelease</figcaption>
  <audio controls>
    <source src="snd/tone_TriggerWhileRelease.wav" type="audio/wav">
  </audio>
</figure>

<figure>
  <figcaption>ChangeSustain</figcaption>
  <audio controls>
    <source src="snd/tone_ChangeSustain.wav" type="audio/wav">
  </audio>
</figure>

### デクリック (declick)
エンベロープのアタック時間やディケイ時間などがとても短いときにでるプチノイズを低減することをデクリック (declick) と呼びます。 Image-Line の [Sytrus](https://www.image-line.com/plugins/Synths/Sytrus/) というシンセサイザでの用例にならっています。

ここでは 1 ミリ秒の短いアタックを持つ副エンベロープを用意して指数曲線を使った主エンベロープに掛け合わせています。また、[直線の ADSR エンベロープ](https://ryukau.github.io/filter_notes/linear_envelope/linear_envelope.html)で紹介した `adaptTime()` も組み合わせて使っています。

モノフォニックのときに副エンベロープを使うと再トリガの処理が複雑になるのでお勧めしません。出力を slew limiter に通すほうが楽に実装できます。[オーディオプラグインの UI から入力された値の補間](../control_rate_interpolation/control_rate_interpolation.html)も参考にしてみてください。

図は副エンベロープによるデクリックを行ったエンベロープ出力の例です。副エンベロープのアタックカーブは $0.5 + 0.5 \cos(n)$ で、 $n$ は $-\pi$ から $0$ に向かって増加しています。

<figure>
<img src="img/declick.png" alt="Image of comparison of raw attack and declicked attack." style="padding-bottom: 12px;"/>
</figure>

1000 Hz のサイン波の音量に適用した音のサンプルです。プチノイズが聞き取りやすいようにサイン波の初期位相を $\dfrac{\pi}{2}$ にしています。

<figure>
  <figcaption>declick_on</figcaption>
  <audio controls>
    <source src="snd/declick_on.wav" type="audio/wav">
  </audio>
</figure>

<figure>
  <figcaption>declick_off</figcaption>
  <audio controls>
    <source src="snd/declick_off.wav" type="audio/wav">
  </audio>
</figure>

## P Controller による実装
P controller については[オーディオプラグインの UI から入力された値の補間](../control_rate_interpolation/control_rate_interpolation.html)を参照してください。

ステップ状の不連続点を始点として指数曲線を描く出力が得られる P controller の特性を利用してエンベロープを実装します。

ユーザから指定された時間 $T$ の逆数をカットオフ周波数 $f_c$ に使っています。

$$
f_c = \frac{1}{T}
$$

アタックはカウンタを使って時間経過で次の状態に進みます。

リリースは出力がしきい値 `threshold` 以下になったときに次の状態に進みます。

状態 `State::tail` では出力が `threshold` から 0 に到達するように直線を描きます。

```cpp
#include <algorithm>
#include <cmath>

constexpr double twopi = 6.283185307179586;

template<typename Sample> class PController {
public:
  // float 型での cutoffHz の下限は 3~4 Hz 程度。
  static Sample cutoffToP(Sample sampleRate, Sample cutoffHz)
  {
    auto omega_c = Sample(twopi) * cutoffHz / sampleRate;
    auto y = Sample(1) - cos(omega_c);
    return -y + sqrt((y + Sample(2)) * y);
  }

  void setP(Sample p) { kp = std::clamp<Sample>(p, Sample(0), Sample(1)); }
  void reset(Sample value = 0) { this->value = value; }
  Sample process(Sample input) { return value += kp * (input - value); }

  Sample kp = 1; // Range in [0, 1].
  Sample value = 0;
};

template<typename Sample> class ExpADSREnvelopeP {
public:
  void setup(Sample sampleRate)
  {
    this->sampleRate = sampleRate;
    tailLength = uint32_t(0.01 * sampleRate);
  }

  void reset(Sample attackTime, Sample decayTime, Sample sustainLevel, Sample releaseTime)
  {
    state = State::attack;
    sustain = std::clamp<Sample>(sustainLevel, Sample(0), Sample(1));
    atk = int32_t(sampleRate * attackTime);
    decTime = decayTime;
    relTime = releaseTime;
    pController.setP(PController<Sample>::cutoffToP(sampleRate, Sample(1) / attackTime));
  }

  void set(Sample attackTime, Sample decayTime, Sample sustainLevel, Sample releaseTime)
  {
    switch (state) {
      case State::attack:
        atk = int32_t(sampleRate * attackTime);
        // Fall through.

      case State::decay:
        decTime = decayTime;
        sustain = std::clamp<Sample>(sustainLevel, Sample(0), Sample(1));
        // Fall through.

      case State::release:
        relTime = releaseTime;

      default:
        break;
    }

    if (state == State::attack)
      pController.setP(
        PController<Sample>::cutoffToP(sampleRate, Sample(1) / attackTime));
    else if (state == State::decay)
      pController.setP(PController<Sample>::cutoffToP(sampleRate, Sample(1) / decayTime));
    else if (state == State::release)
      pController.setP(
        PController<Sample>::cutoffToP(sampleRate, Sample(1) / releaseTime));
  }

  void release()
  {
    state = State::release;
    pController.setP(PController<Sample>::cutoffToP(sampleRate, Sample(1) / relTime));
  }

  bool isAttacking() { return state == State::attack; }
  bool isReleasing() { return state == State::release; }
  bool isTerminated() { return state == State::terminated; }

  Sample process()
  {
    switch (state) {
      case State::attack: {
        value = pController.process(Sample(1));
        --atk;
        if (atk == 0) {
          state = State::decay;
          pController.setP(
            PController<Sample>::cutoffToP(sampleRate, Sample(1) / decTime));
        }
      } break;

      case State::decay:
        value = pController.process(sustain);
        break;

      case State::release:
        value = pController.process(0);
        if (value < threshold) {
          value = threshold;
          state = State::tail;
          tailCounter = tailLength;
        }
        break;

      case State::tail:
        --tailCounter;
        value = threshold * tailCounter / float(tailLength);
        if (tailCounter == 0) {
          state = State::terminated;
          pController.reset(0);
        } else {
          pController.reset(value);
        }
        break;

      default:
        return 0;
    }
    return value;
  }

private:
  enum class State : int32_t { attack, decay, release, tail, terminated };
  const Sample threshold = 1e-5;

  uint32_t tailLength = 32;
  uint32_t tailCounter = tailLength;

  PController<Sample> pController;
  State state = State::terminated;
  uint32_t atk = 0;
  Sample decTime = 0;
  Sample relTime = 0;
  Sample sampleRate = 44100;
  Sample sustain = 1;
  Sample value = 0;
};
```

テストに使ったコードへのリンクです。

- TODO リンク

テスト結果です。

<figure>
<img src="img/testP_ADSR.png" alt="Image of test result of exponential ADSR envelope with P controller." style="padding-bottom: 12px;"/>
</figure>

<figure>
<img src="img/testP_ReleaseWhileAttack.png" alt="Image of test result of release while attack." style="padding-bottom: 12px;"/>
</figure>

<figure>
<img src="img/testP_ReleaseWhileDecay.png" alt="Image of test result of release while decay." style="padding-bottom: 12px;"/>
</figure>

<figure>
<img src="img/testP_TriggerWhileRelease.png" alt="Image of test result of trigger while release." style="padding-bottom: 12px;"/>
</figure>

<figure>
<img src="img/testP_ChangeSustain.png" alt="Image of test result of changing sustain." style="padding-bottom: 12px;"/>
</figure>

`P_ADSR` の `State::tail` の部分を拡大した図です。

<figure>
<img src="img/P_tail.png" alt="Image of tail state of envelope." style="padding-bottom: 12px;"/>
</figure>


## 減衰する指数曲線とその反転の乗算
減衰する指数曲線とその反転を掛け合わせて合成したエンベロープ $E_{\mathtt{AD}}(t)$ を作ります。

$$
E_{\mathtt{AD}}(t) = (1 - a^{t}) d^{t}
$$

ピークの位置 $t_p$ とピークの大きさ $E_{\mathtt{AD}}(t_p)$ を求めます。 $t_p$ はアタック時間です。 $E_{\mathtt{AD}}(t_p)$ は出力範囲を $[0, 1]$ に正規化するために使えます。

$t_p$ の時点で $E_{\mathtt{AD}}(t)$ を $t$ について微分した関数の値が 0 になるはずです。 Maxima で解きます。

```maxima
expr: diff((1-a^t) * d^t, t);
solve(0 = expr, t);
```

出力です。

$$
t_p = \frac{\log{\left( \dfrac{\log{(d)}}{\log{(d)}+\log{(a)}}\right) }}{\log{(a)}}
$$

$a$ と $d$ を求めます。ユーザが指定したアタック時間を $A$ 、 ユーザが指定したディケイ時間を $D$ とします。 $[0, 1)$ の範囲の適当なしきい値 $\epsilon$ を用意して、 $a^{A} = d^{D} = \epsilon$ とすると $a, d$ は次の式で計算できます。

$$
a = \epsilon^{1/A}, \quad b = \epsilon^{1/B}
$$

C++ での実装です。

```cpp
#include <cmath>

class ExpAD {
public:
  void setup(float sampleRate) { this->sampleRate = sampleRate; }
  bool isTerminated() { return valueD <= threshold; }

  // attack and decay in seconds.
  void reset(float attack, float decay)
  {
    valueA = 1.0f;
    if (attack < 1e-5) attack = 1e-5;
    alphaA = powf(threshold, 1.0f / (attack * sampleRate));

    valueD = 1.0f;
    if (decay < 1e-5) decay = 1e-5;
    alphaD = powf(threshold, 1.0f / (decay * sampleRate));

    if (attack <= 0.0f) {
      gain = 1.0f;
    } else if (decay <= 0.0f) {
      gain = 0.0f;
    } else {
      auto log_a = logf(alphaA);
      auto log_d = logf(alphaD);
      auto t_p = logf(log_d / (log_a + log_d)) / log_a;
      gain = 1.0f / ((1.0f - powf(alphaA, t_p)) * powf(alphaD, t_p));
    }
  }

  float process()
  {
    valueA *= alphaA;
    valueD *= alphaD;
    return gain * (1.0f - threshold - valueA) * (valueD - threshold);
  }

protected:
  const float threshold = 1e-5;
  float sampleRate = 44100;
  float gain = 0;
  float valueA = 0;
  float alphaA = 0;
  float valueD = 0;
  float alphaD = 0;
};
```

テストコードへのリンクです。

- [filter_notes/expAD.cpp at master · ryukau/filter_notes · GitHub](https://github.com/ryukau/filter_notes/blob/master/exponential_envelope/demo/expAD.cpp)

テスト結果です。

<figure>
<img src="img/ExpAD_cpp.png" alt="Image of ExpAD envelope. C++ implementation." style="padding-bottom: 12px;"/>
</figure>

### 別解
この計算方法は [EnvelopedSine](https://ryukau.github.io/VSTPlugins/manual/EnvelopedSine/EnvelopedSine_ja.html) で使っています。エンベロープ $\tilde{E}_{\mathtt{AD}}$ の式です。

$$
\tilde{E}_{\mathtt{AD}}(t) = (1 - e^{-at}) e^{-bt}
$$

$E_{\mathtt{AD}}(t)$ の式について $a \to e^{-a},\ d \to e^{-b}$ と置き換えています。

ピークの位置 $t_p$ とピークの大きさ $\tilde{E}_{\mathtt{AD}}(t_p)$ を求めます。

```maxima
expr: diff((1-exp(-a*t)) * exp(-b*t), t);
solve(0 = expr, t);
```

$$
t_p = \frac{\log{\left( \dfrac{a}{b}+1\right) }}{a}
$$

$a$ と $b$ を求めます。アタック時間を $A$ 、 ディケイ時間を $B$ 、 適当なしきい値を $\epsilon \in [0, 1)$ とします。 $e^{-a A},\ e^{-b B},\ \epsilon$ が等しくなるような $a, b$ は次の式で計算できます。

$e^{-a A} = \epsilon$ より、

$$
a = - \dfrac{\log(\epsilon)}{A}.
$$

$e^{-b B} = \epsilon$ より、

$$
\quad b = - \dfrac{\log(\epsilon)}{B}.
$$

コード例です。

```python
import numpy as np

def envelopeExpr(a, b, time):
    return (1 - np.exp(-a * time)) * np.exp(-b * time)

def envelope(attack, decay, eps=1e-5):
    _a = -np.log(eps) / attack
    _b = -np.log(eps) / decay
    peakTime = np.log(_a / _b + 1) / _a
    gain = 1 / envelopeExpr(_a, _b, peakTime)

    samplerate = 48000
    duration = 1
    time = np.linspace(0, duration, duration * samplerate)

    return envelopeExpr(_a, _b, time)

output = envelope(1.0, 2.0)
```

<figure>
<img src="img/ExpAD.png" alt="Image of ExpAD envelope. Alternative implementation." style="padding-bottom: 12px;"/>
</figure>

## 変更点
- 2020-03-27
  - `P Controller による実装` を追加。
