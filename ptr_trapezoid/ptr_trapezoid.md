# PTR 台形オシレータ
[TrapezoidSynth](https://ryukau.github.io/VSTPlugins/manual/TrapezoidSynth/TrapezoidSynth_ja.html) の[オシレータのコード](https://github.com/ryukau/VSTPlugins/blob/91d663e7ed916bc7c5e372990ad40eb9736ad84e/TrapezoidSynth/source/dsp/oscillator.hpp#L106)があまりにも場当たり的で自分でも読めなくなったので整理しました。

## PTR 台形の特性
TrapezoidSynth の PTR 台形オシレータは 4 つの PTR [ランプ関数](https://en.wikipedia.org/wiki/Ramp_function)をつぎはぎすることで台形を合成しています。このつぎはぎが複雑さの原因です。

ランプ関数の定義です。

$$
R(x) = \begin{cases}
0, & \text{if}\ x < 0,\\
x, & \text{otherwise.}
\end{cases}
$$

PTR ランプ関数の導出は次のリンクにまとめています。

- [PTR オシレータ](../ptr_oscillator/ptr_oscillator.html)

次の図はランプ関数のつぎはぎの様子を表しています。下側の矢印はランプ関数の $x$ が増える方向です。

<figure>
<img src="img/ramp_patch.svg" alt="Image of trapezoid made from 4 ramp functions." style="width: 480px;padding-bottom: 12px;"/>
</figure>

- [Ramp function - Wikipedia](https://en.wikipedia.org/wiki/Ramp_function)

まずは台形の各部に名前をつけます。

<figure>
<img src="img/trapezoid_dc1.svg" alt="Image of trapezoid with name of parts." style="width: 480px;padding-bottom: 12px;"/>
</figure>

### $y$ の値
PTR の次数を $N$ 、 1 サンプルあたりに進む位相を $T$ とすると $y = 1 - 2 K N T$ です。 $- 2 K N T$ という値は試行錯誤して見つけた適当な値です。

### $A_1$ と $A_2$ の関係
$A_1$ と $K$ はユーザが決定する定数です。 $A_1$ は1つ目の台形の上辺の長さ、 $K$ は台形の斜線の傾きです。

$B_1 + B_2$ は 1 周期の波長です。ここでは正規化して $B_1 + B_2 = 1$ と定義します。この定義を使って、2つ目の台形の上辺の長さ $A_2$ を計算できます。

$$
\begin{aligned}
A_2 &= B_1 + B_2 - A_1 - \frac{1}{K}\\
  &= 1 - A_1 - \frac{1}{K} & \text{(1)}
\end{aligned}
$$

### $A_1$ の範囲
$A_1$ の最小値は $0$ です。

$A_1$ が最大となるとき、 $A_2 = 0$ です。 $A_2 = 0$ を式 (1) に代入すると $A_1$ の最大値 $1 - \dfrac{1}{K}$ が得られます。

### $K$ の範囲
$A_1 = 0,\ A_2 = 0$ のとき、式 (1) より $K = 1$ です。このとき三角波に近い波形になります。 $K < 1$ のときは三角波が途中でハードシンクされたような滑らかでない波形になります。

$A_1 = 0,\ A_2 = 1$ のとき、式 (1) より $K = \infty$ となります。ただし、実際に $K = \infty$ とすると PTR の遷移領域を飛び越えてしまうので波形が不連続になってしまいます。そこで、実装では $A_1 = 0,\ A_2 = 1 - 4 N T$ を式 (1) に代入して得られる $K = \dfrac{1}{4 N T}$ を $K$ の最大値としています。 $4 N T$ は 4 つの PTR ランプ関数の遷移領域の長さです。

### 直流補正
図では $h$ を引いて直流成分を補正しています。 $A_1, B_1$ からなる台形 1 と、 $A_2, B_2$ からなる台形 2 の面積が等しいときに直流が 0 になります。 $h$ を求める方程式を立てます。

$$
\begin{aligned}
B_1 &= A_1 + \frac{2(y - h)}{K}\\
B_2 &= A_2 + \frac{2h}{K}\\
\frac{(A_1 + B_1) (y - h)}{2} &= \frac{(A_2 + B_2) h}{2}\\
\end{aligned}
$$

Maxima で解きます。

```maxima
A2: 1 - A1 - 1 / K;
B1: A1 + 2 * (y - h) / K;
B2: A2 + 2 * h / K;
eq: (A1 + B1) * (y - h) / 2 = (A2 + B2) * h / 2;
expand(solve(eq, h));
```

出力です。

$$
h=\frac{y^2 + A_1 K y}{2 y+K-1}
$$

これで直流成分の大きさ $h$ が分かりました。

## 音程が高いときのノイズの低減
記号の定義です。

- $f$ : オシレータの周波数
- $f_s$ : サンプリング周波数
- $N$ : PTR の次数

PTR 台形オシレータの次数 $N$ を大きくするほどエイリアシングノイズが低減されますが、周波数が $\dfrac{f_s}{4 N}$ を超えるとノイズが増えてしまうというトレードオフがあります。ノイズが乗るのは 4 つのランプ関数のすべての遷移領域を通過する前に波形の 1 周期が終わってしまうからです。

$f_s = 44100$ として $N$ を変えたときのノイズが乗らない周波数の上限は次のようになります。

- $N = 2 \to 5512.5\,\mathrm{Hz}$
- $N = 3 \to 3675\,\mathrm{Hz}$
- $N = 4 \to 2756.25\,\mathrm{Hz}$
- $N = 5 \to 2205\,\mathrm{Hz}$

$f_s \geq 44100$ のときは $N = 4$ にして 8 倍のオーバーサンプリングをかけると、ノイズが乗らない周波数の上限が 22050 Hz になるので十分な音域が使えるようになります。ただし、オーバーサンプリングを使わずにエイリアシングノイズを効率よく低減できるという PTR の利点を無視することになります。

オーバーサンプリングを極力避けるためには次の式のように PTR の次数を切り替えることが考えられます。例として $N$ が 2 から 5 の間で切り替わっています。

$$
\begin{aligned}
N_{\mathtt{TrapezoidSynth}}(f) &= \begin{cases}
5, & \text{if} \enspace f < f_u(5)\\
4, & \text{if} \enspace f_u(5) \leq f < f_u(4)\\
3, & \text{if} \enspace f_u(4) \leq f < f_u(3)\\
2, & \text{if} \enspace f_u(3) \leq f\\
\end{cases}\\
\text{where}\enspace
f_u(N) &= \frac{f_s}{4 N}
\end{aligned}
$$

分岐処理がオーバーサンプリングよりも速いかどうかはベンチマークを取ってみないとわかりません。また $N$ の切り替えのタイミングで不連続になるので、ピッチベンドなどで周波数が動くとノイズが乗ります。

TrapezoidSynth では何を思ったのか $N$ の切り替えと 8 倍のオーバーサンプリングの両方が実装されています。

## 実装
C++ です。

定数をキャストしているのはテンプレートで `float` と `double` を切り替えることを想定していたからです。結局オーバーサンプリングしてごまかしたのでテンプレート化は見送りました。

```cpp
#include <algorithm>

class PTRTrapezoidOsc {
public:
  PTRTrapezoidOsc(float sampleRate, float frequencyHz) : sampleRate(sampleRate)
  {
    setFreq(frequencyHz);
  }

  void setFreq(float hz)
  {
    if (hz >= 0) tick = hz / sampleRate;
  }

  void setPhase(float phase) { this->phase = phase; }
  void addPhase(float phase) { this->phase += phase; }
  void setSlope(float slope) { this->slope = slope; }
  void setPulseWidth(float pw) { this->pw = pw; }
  void reset() { phase = 0; }

  float process()
  {
    if (tick <= 0) return 0;
    phase += tick;
    phase -= somefloor<float>(phase);
    return ptrTpz5(phase, tick, slope, pw);
  }

  float sampleRate;
  float phase = 0;
  float tick = 0;
  float slope = 8;       // 台形の斜辺の傾き。この文章の K に相当。
  float pw = float(0.5); // Pulse width. A1 に相当。

protected:
  // tick must be greater than 0.
  static float ptrTpz5(float phase, float tick, float slope, float pw)
  {
    uint32_t order = 5; // PTR の次数。 N に相当。
    if (order > float(0.25) / tick) order = int(float(0.25) / tick);

    const float ptrLen = order * tick; // NT に相当。

    // A1 の範囲を制限。
    const float maxSlope = float(0.25) / ptrLen;
    if (slope > maxSlope) {
      slope = maxSlope;
    } else {
      const float minSlope = float(1);
      if (slope < minSlope) slope = minSlope;
    }

    // K の範囲を制限。
    const float maxPw = float(1) - float(1) / slope;
    if (pw > maxPw) pw = std::max<float>(float(0), maxPw);

    // 直流の計算。 y は適当に見つけた上手く動く値。
    const float y = float(1) - float(2) * slope * ptrLen;
    const float dc = (y * y + pw * slope * y) / (float(2) * y + slope - float(1));

    // N による分岐。
    if (order == 5)
      return branch<PTRTrapezoidOsc::ptrRamp5>(slope, pw, y, phase, tick) - dc;
    else if (order == 4)
      return branch<PTRTrapezoidOsc::ptrRamp4>(slope, pw, y, phase, tick) - dc;
    else if (order == 3)
      return branch<PTRTrapezoidOsc::ptrRamp3>(slope, pw, y, phase, tick) - dc;
    return branch<PTRTrapezoidOsc::ptrRamp2>(slope, pw, y, phase, tick) - dc;
  }

  template<float (*ptrfunc)(float, float)>
  static float branch(float slope, float pw, float y, float phase, float tick)
  {
    // ランプ関数のつぎはぎを行う分岐。
    if (phase <= float(0.25) / slope) // Ramp 1
      return slope * ptrfunc(phase, tick);
    else if (phase <= float(0.5) / slope) // Ramp 2
      return y - slope * ptrfunc(float(0.5) / slope - phase, tick);
    else if (phase <= float(0.5) / slope + pw) // 台形の上辺。
      return y;
    else if (phase <= float(0.75) / slope + pw) // Ramp 3
      return y - slope * ptrfunc(phase - float(0.5) / slope - pw, tick);
    else if (phase <= float(1) / slope + pw) // Ramp 4
      return slope * ptrfunc(float(1) / slope + pw + -phase, tick);
    return float(0); // 残りを 0 埋め。台形の各部の名前の図の A2 にあたる領域。
  }

  // 以降は PTR ランプ関数。
  static float ptrRamp2(float phi, float T)
  {
    float n = phi / T;
    if (n >= float(1)) return float(2) * T * n - float(2) * T;
    if (n < float(1)) return (T * n * n * n) / float(3);
    if (n < float(2))
      return -(T * n * n * n) / float(3) + float(2) * T * n * n - float(2) * T * n
        + (float(2) * T) / float(3);
    return 0.0; // Just in case.
  }

  static float ptrRamp3(float phi, float T)
  {
    float n = phi / T;
    if (n >= float(2)) return float(2) * T * n - float(3) * T;
    if (n < float(1)) return (T * n * n * n * n) / float(12);
    if (n < float(2))
      return -(T * n * n * n * n) / float(6) + T * n * n * n
        - (float(3) * T * n * n) / float(2) + T * n - T / float(4);
    if (n < float(3))
      return (T * n * n * n * n) / float(12) - T * n * n * n
        + (float(9) * T * n * n) / float(2) - float(7) * T * n
        + (float(15) * T) / float(4);
    return 0.0; // Just in case.
  }

  static float ptrRamp4(float phi, float T)
  {
    float n = phi / T;
    if (n >= float(3)) return float(2) * T * n - float(4) * T;
    if (n < float(1)) return (T * n * n * n * n * n) / float(60);
    if (n < float(2))
      return -(T * n * n * n * n * n) / float(20) + (T * n * n * n * n) / float(3)
        - (float(2) * T * n * n * n) / float(3) + (float(2) * T * n * n) / float(3)
        - (T * n) / float(3) + T / float(15);
    if (n < float(3))
      return (T * n * n * n * n * n) / float(20)
        - (float(2) * T * n * n * n * n) / float(3)
        + (float(10) * T * n * n * n) / float(3) - (float(22) * T * n * n) / float(3)
        + (float(23) * T * n) / float(3) - (float(47) * T) / float(15);
    if (n < float(4))
      return -(T * n * n * n * n * n) / float(60) + (T * n * n * n * n) / float(3)
        - (float(8) * T * n * n * n) / float(3) + (float(32) * T * n * n) / float(3)
        - (float(58) * T * n) / float(3) + (float(196) * T) / float(15);
    return 0.0; // Just in case.
  }

  static float ptrRamp5(float phi, float T)
  {
    float n = phi / T;
    if (n >= float(4)) return float(2) * T * n - float(5) * T;
    if (n < float(1)) return (T * n * n * n * n * n * n) / float(360);
    if (n < float(2))
      return -(T * n * n * n * n * n * n) / float(90)
        + (T * n * n * n * n * n) / float(12) - (float(5) * T * n * n * n * n) / float(24)
        + (float(5) * T * n * n * n) / float(18) - (float(5) * T * n * n) / float(24)
        + (T * n) / float(12) - T / float(72);
    if (n < float(3))
      return (T * n * n * n * n * n * n) / float(60) - (T * n * n * n * n * n) / float(4)
        + (float(35) * T * n * n * n * n) / float(24)
        - (float(25) * T * n * n * n) / float(6) + (float(155) * T * n * n) / float(24)
        - (float(21) * T * n) / float(4) + (float(127) * T) / float(72);
    if (n < float(4))
      return -(T * n * n * n * n * n * n) / float(90) + (T * n * n * n * n * n) / float(4)
        - (float(55) * T * n * n * n * n) / float(24)
        + (float(65) * T * n * n * n) / float(6) - (float(655) * T * n * n) / float(24)
        + (float(141) * T * n) / float(4) - (float(1331) * T) / float(72);
    if (n < float(5))
      return (T * n * n * n * n * n * n) / float(360)
        - (T * n * n * n * n * n) / float(12)
        + (float(25) * T * n * n * n * n) / float(24)
        - (float(125) * T * n * n * n) / float(18) + (float(625) * T * n * n) / float(24)
        - (float(601) * T * n) / float(12) + (float(2765) * T) / float(72);
    return 0.0; // Just in case.
  }
};
```
