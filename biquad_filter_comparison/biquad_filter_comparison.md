# Biquad フィルタの比較
以下の Biquad フィルタの実装を比較します。

1. Robert Bristow-Johnson (RBJ) さんによる[レシピ](https://www.w3.org/TR/audio-eq-cookbook/)。
2. Vadim Zavalishin さんによる [The Art of VA Filter Design](https://archive.org/details/the-art-of-va-filter-design-rev.-2.1.2) の SVF。 Chapter 4 に詳細。

ここでは 1 を RBJ フィルタ、 2 を TPT フィルタと呼ぶことにします。 TPT は topology-preserving transform の略で、 The Art of VA Filter Design で紹介されているアナログフィルタの離散化の方法です。

## 実装
C++ 20 で実装します。標準ライブラリの `<array>`, `<cmath>`, `<numbers>` をインクルードしています。

三角関数の計算は重たいのでフィルタ係数はコントロールレートでのみ計算しています。オーディオレートでの補間には以下の `ExpSmoother` を使っています。 `ExpSmoother` は一次ローパスあるいは exponential moving average (EMA) と呼ばれるフィルタです。

```c++
template<typename Sample> struct ExpSmoother {
  Sample value = 0;
  Sample target = 0;
  Sample kp = 1; // Exponential moving average coefficient.

  void setCutoff(Sample sampleRate, Sample cutoffHz)
  {
    // `double` is used for accuracy.
    double y = double(1)
      - std::cos(double(2) * std::numbers::pi_v<double> * cutoffHz / sampleRate);
    kp = Sample(std::sqrt((y + double(2)) * y) - y);
  }

  void reset(Sample newTarget = 0)
  {
    value = newTarget;
    target = newTarget;
  }

  void push(Sample newTarget = 0) { target = newTarget; }
  inline Sample process() { return value += kp * (target - value); }
};
```

### RBJ フィルタ
[Audio EQ Cookbook](https://www.w3.org/TR/audio-eq-cookbook/) に基づいたローパスのみの実装例です。

```c++
#include <array>
#include <cmath>
#include <numbers>

template<typename Sample> class RbjBiquad {
  Sample x1 = 0;
  Sample x2 = 0;
  Sample y1 = 0;
  Sample y2 = 0;
  std::array<ExpSmoother<Sample>, 5> co; // Coefficients {b0, b1, b2, -a1, -a2}.

public:
  const std::string name{"RbjBiquad"};

  RbjBiquad()
  {
    for (auto &x : co) x.setCutoff(Sample(1), Sample(0.02));
  }

#define ASSIGN_COEFFICINETS(method)                                                      \
  /* normalizedFreq = cutoffHz / sampleRate. */                                          \
  void method##Lowpass(Sample normalizedFreq, Sample Q)                                  \
  {                                                                                      \
    auto omega = std::numbers::pi_v<Sample> * Sample(2) * normalizedFreq;                \
    auto cos = std::cos(omega);                                                          \
    auto sin = std::sin(omega);                                                          \
    auto alpha = sin / (Sample(2) * Q);                                                  \
    auto a0_inv = Sample(1) / (Sample(1) + alpha);                                       \
    co[0].method(a0_inv *((Sample(1) - cos) / Sample(2))); /* b0 */                      \
    co[1].method(a0_inv *(Sample(1) - cos));               /* b1 */                      \
    co[2].method(a0_inv *((Sample(1) - cos) / Sample(2))); /* b2 */                      \
    co[3].method(-a0_inv *(Sample(-2) * cos));             /* a1 */                      \
    co[4].method(-a0_inv *(Sample(1) - alpha));            /* a2 */                      \
  }

  ASSIGN_COEFFICINETS(reset);
  ASSIGN_COEFFICINETS(push);
#undef ASSIGN_COEFFICINETS

  void reset(Sample normalizedFreq, Sample Q, Sample value = 0)
  {
    x1 = value;
    x2 = value;
    y1 = value;
    y2 = value;

    resetLowpass(normalizedFreq, Q);
  }

  void prepare(Sample normalizedFreq, Sample Q) { pushLowpass(normalizedFreq, Q); }

  Sample process(Sample x0)
  {
    for (auto &x : co) x.process();

    auto y0 = co[0].value * x0 + co[1].value * x1 + co[2].value * x2 + co[3].value * y1
      + co[4].value * y2;

    x2 = x1;
    x1 = x0;
    y2 = y1;
    y1 = y0;

    return y0;
  }
};
```

`process()` の `y0` の計算で加算と減算を混ぜないようにするため、フィルタ係数 `co` の $a_1, a_2$ の符号を反転しています。

マクロで定義しているのでややこしいですが、 `resetLowpass()` と `pushLowpass()` でフィルタ係数を設定します。 `reset*` は音の再生開始前などに内部状態をリセットするときに呼び出します。 `push*` と `prepare` は音の再生中にコントロールレートでフィルタ係数を設定するときに使います。 `prepare` は個人的に使っているコントロールレートでのパラメータ設定を行うメソッド名で、消しても問題ありません。

テンプレートで並列処理する信号の数を指定できるようにした実装も行いました。これは `cl.exe` のオプションに `/Qvec-report:1` を指定して SIMD がどれくらい入るかを検証するためです。

- [RbjBiquadParallel](https://github.com/ryukau/filter_notes/blob/3b830660c7a13fbb730ffd14e0a251fc53acdc54/biquad_filter_comparison/code/cpp/benchmark.cpp#L283-L377)

### TPT フィルタ
The Art of VA Filter Design で紹介されている手法に基づいたローパスのみの実装例です。 [Faust Libraries の svf](https://github.com/grame-cncm/faustlibraries/blob/ea4fae33bb57357c7e2880c094079018bb4bd769/filters.lib#L2628-L2707) を参考にしています。

```c++
template<typename Sample> class TptBiquad {
private:
  ExpSmoother<Sample> g;
  ExpSmoother<Sample> R;
  ExpSmoother<Sample> d;

  Sample s1 = 0;
  Sample s2 = 0;

public:
  const std::string name{"TptBiquad"};

  TptBiquad()
  {
    constexpr Sample cutoff = Sample(0.02);

    g.setCutoff(Sample(1), cutoff);
    R.setCutoff(Sample(1), cutoff);
    d.setCutoff(Sample(1), cutoff);
  }

#define ASSIGN_COEFFICINETS(method)                                                      \
  /* normalizedFreq = cutoffHz / sampleRate. */                                          \
  void method##Coefficient(Sample normalizedFreq, Sample Q)                              \
  {                                                                                      \
    g.method(std::tan(std::numbers::pi_v<Sample> *normalizedFreq));                      \
    R.method(Sample(1) / Q);                                                             \
    d.method(Sample(1) / (Sample(1) + g.target * g.target + g.target * R.target));       \
  }

  ASSIGN_COEFFICINETS(push);
  ASSIGN_COEFFICINETS(reset);
#undef ASSIGN_COEFFICINETS

  void reset(Sample normalizedFreq, Sample Q, Sample value = 0)
  {
    s1 = value;
    s2 = value;
    resetCoefficient(normalizedFreq, Q);
  }

  void prepare(Sample normalizedFreq, Sample Q) { pushCoefficient(normalizedFreq, Q); }

  // Lowpass.
  Sample process(Sample x0)
  {
    g.process();
    d.process();

    auto v1 = (s1 + g.value * (x0 - s2)) * d.value;
    auto v2 = s2 + g.value * v1;
    s1 = Sample(2) * v1 - s1;
    s2 = Sample(2) * v2 - s2;

    return v2;
  }
};
```

以下はテンプレートで並列処理する信号の数を指定できるようにした実装です。

- [TptBiquadParallel](https://github.com/ryukau/filter_notes/blob/3b830660c7a13fbb730ffd14e0a251fc53acdc54/biquad_filter_comparison/code/cpp/benchmark.cpp#L379-L450)

以下はオーディオレートで計算を行う素朴な実装です。

- [filter_notes/biquad_filter_comparison/code/svf.cpp at master · ryukau/filter_notes · GitHub](https://github.com/ryukau/filter_notes/blob/master/biquad_filter_comparison/code/svf.cpp)

以下は動作確認に使った Python 3 での実装です。

- [filter_notes/biquad_filter_comparison/code/svf.py at master · ryukau/filter_notes · GitHub](https://github.com/ryukau/filter_notes/blob/master/biquad_filter_comparison/code/svf.py)

## 比較
使い勝手については以下の違いがあります。

- RBJ フィルタはフィルタの種類によって係数の設定が異なるが、 `process` の計算は共通。
- TPT フィルタはフィルタの種類によって `process` の計算が異なるが、係数の設定はおおよそ共通。

コードの長さについては TPT フィルタのほうが短くなる傾向があります。ただし TPT フィルタは bell 、 low-shelf 、 high-shelf については係数の計算が変わります。

RBJ フィルタはフィルタ係数が 5 つで、補間を行うときは安定性に疑問が残ります。 TPT フィルタの係数は 3 つと比較的少ないですが、こちらも `d` については補間による安定性がはっきりしません。 `d` については補間を行わずに `g` と `R` からオーディオレートで計算することで安定しますが、除算が出てくるので速度は落ちます。

今回の実装の計算速度について簡単なベンチマークをとってみたところ、以下のような結果が得られました。`Sample Count` は計算したサンプル数です。 `Last Output` は最後の出力の総和で、実装間の出力に差がないことを確認するためにプリントしています。 `TptSvf` はローパス、バンドパス、ハイパスを同時に計算する `TptBiquad` の実装の派生です。 `RbjBiquadParallel` と `TptBiquadParallel` は 16 出力を同時に計算しています。

以下は `float` でのベンチマーク結果です。

| Filter Type       |  Time Elapsed [ms] | Sample Count |  Last Output |
|-------------------|-------------------:|-------------:|-------------:|
| RbjBiquad         |  17.73569999997201 |       960000 |  -0.03840086 |
| TptBiquad         | 17.528799999972524 |       960000 | -0.038400784 |
| TptSvf            | 17.455699999972666 |       960000 | -0.038400736 |
| RbjBiquadParallel |  290.3656000015933 |       960000 |   0.15570019 |
| TptBiquadParallel | 148.51210000301913 |       960000 |   0.15563208 |

以下は `double` でのベンチマーク結果です。

| Filter Type       |  Time Elapsed [ms] | Sample Count |         Last Output |
|-------------------|-------------------:|-------------:|--------------------:|
| RbjBiquad         |  19.11819999996879 |       960000 | 0.06474569475653402 |
| TptBiquad         |  22.07949999996203 |       960000 | 0.06474569475653402 |
| TptSvf            |  21.63369999996292 |       960000 | 0.06474569475653393 |
| RbjBiquadParallel | 293.51290000144576 |       960000 |    0.81283637843879 |
| TptBiquadParallel |  150.7567000030868 |       960000 |  0.8128361199584616 |

以下はベンチマークに使ったコードへのリンクです。

- [filter_notes/biquad_filter_comparison/code/cpp/benchmark.cpp at master · ryukau/filter_notes · GitHub](https://github.com/ryukau/filter_notes/blob/master/biquad_filter_comparison/code/cpp/benchmark.cpp)

結論としては TPT フィルタだけを使えばよさそうです。特に多出力では TPT フィルタが圧倒的に速いです。 `double` かつ 1 出力であれば RBJ フィルタがわずかに速いので、普段は TPT フィルタを使い、どうしても速度が要るときに RBJ フィルタへの切り替えも検討するという方針がよさそうです。

また `cl.exe` の `/Qvec-report:2` で検証したところ、今回の実装のすべてのフィルタのコードについて[自動ベクトル化](https://learn.microsoft.com/en-us/cpp/parallel/auto-parallelization-and-auto-vectorization?view=msvc-170#auto-vectorizer) (auto vectorization) は行われませんでした。より速くするには SIMD を直接、あるいは [simde](https://github.com/simd-everywhere/simde) や [VCL](https://github.com/vectorclass/version2) などのライブラリを通して間接的に書くしかなさそうです。

## その他
### Butterworth フィルタ
TPT フィルタのレゾナンスが $R = \dfrac{1}{\sqrt{2}}$ のとき Butterworth フィルタになります。

Butterworth フィルタは maximally flat なフィルタです。つまり、ローパスであればカットオフ周波数から 0 Hz までの範囲の振幅特性が最も平坦となります。また、振幅特性が 1 倍 (あるいは 0 dB) を超えない最大のレゾナンスが設定されたフィルタとも言えます。

### レゾナンスのゲイン
The Art of VA Filter Design (rev. 2.1.2) の 4.2 Resonance にローパスのレゾナンスが $R$ のときのゲイン $A_{\mathrm{LP}}$ の計算方法が載っています。

$$
\begin{align*}
A_{\mathrm{LP}} &= \frac{1}{2 R \sqrt{1 - R^2}}, \quad \text{where} \quad  R < \dfrac{1}{\sqrt{2}}. \\
\end{align*}
$$

$R \geq \dfrac{1}{\sqrt{2}}$ のときはレゾナンスによるゲインのピークが現れません。未定義としないときは周波数 0 でのゲインを使って $A_{\mathrm{LP}} = 1$ とすることが考えられます。

ハイパスのゲイン $A_{\mathrm{HP}}$ も導出します。以下はハイパスのアナログプロトタイプの伝達関数です。

$$
\begin{align*}
H_{\mathrm{HP}}(s)
= \frac{s^2}{s^2 + 2 R s + 1}
= \frac{1}{1 + 2 R s^{-1} + s^{-2}}
\end{align*}
$$

$s = j\omega$ を代入して絶対値を取ります。伝達関数の絶対値はゲインを表します。

$$
\begin{align*}
|H_{\mathrm{HP}}(j\omega)| = \frac{1}{|1 + 2 R (j\omega)^{-1} + (j\omega)^{-2}|}
\end{align*}
$$

ここで分母が最小となるとき、ゲインが最大となります。分母の最小値を求めます。

$$
\begin{align*}
0
&= |1 + 2 R (j\omega)^{-1} + (j\omega)^{-2}| \\
&= \sqrt{(1 - \omega^{-2})^2 + (2 R \omega^{-1})^2} \\
&= \sqrt{1 + (4 R^2 - 2) \omega^{-2} + \omega^{-4}} \\
\\
0
&= \omega^{4} + (4 R^2 - 2) \omega^{2} + 1
\end{align*}
$$

この式は The Art of VA Filter Design (rev. 2.1.2), p.102 の $|H_{\mathrm{LP}}(j\omega)|^{-2}$ の式と同じです。左辺が 0 なので、平方根に入っている $\omega^{-2}$ の多項式を二次方程式の解の公式で解くとゲインが最大となる周波数 $\omega_{\mathrm{peak}}$ が得られます。

$$
\omega_{\mathrm{peak}} = \sqrt{1 - 2 R^2}
$$

$\omega_{\mathrm{peak}}$ を伝達関数に代入するとハイパスのレゾナンスのゲイン $A_{\mathrm{HP}}$ が得られます。

$$
\begin{align*}
|H_{\mathrm{HP}}(\omega_{\mathrm{peak}})|^2
&= \frac{1}{1 + (4 R^2 - 2) \omega^{-2} + \omega^{-4}} \\
&= \frac{1}{1 + (4 R^2 - 2) (1 - 2 R^2)^{-1} + (1 - 2 R^2)^{-2}} \\
&= \frac{(1 - 2 R^2)^2}{4 R^2 - 4 R^4} \\
\\
A_{\mathrm{HP}} = |H_{\mathrm{HP}}(\omega_{\mathrm{peak}})|
&=\frac{1 - 2 R^2}{2 R \sqrt{1 - R^2}}, \quad \text{where} \quad R < \dfrac{1}{\sqrt{2}}.\\
\end{align*}
$$

## 参考文献
- [Audio EQ Cookbook](https://www.w3.org/TR/audio-eq-cookbook/)
- [The Art of VA Filter Design (rev. 2.1.2) (February 14, 2020) : Vadim Zavalishin : Free Download, Borrow, and Streaming : Internet Archive](https://archive.org/details/the-art-of-va-filter-design-rev.-2.1.2)
- [faustlibraries/filters.lib at master · grame-cncm/faustlibraries · GitHub](https://github.com/grame-cncm/faustlibraries/blob/ea4fae33bb57357c7e2880c094079018bb4bd769/filters.lib#L2628-L2707)
- [Auto-Parallelization and Auto-Vectorization | Microsoft Learn](https://learn.microsoft.com/en-us/cpp/parallel/auto-parallelization-and-auto-vectorization?view=msvc-170#auto-vectorizer)
- [/Qvec-report (Auto-Vectorizer Reporting Level) | Microsoft Learn](https://learn.microsoft.com/en-us/cpp/build/reference/qvec-report-auto-vectorizer-reporting-level?view=msvc-170)
- [Vectorizer and parallelizer messages | Microsoft Learn](https://learn.microsoft.com/en-us/cpp/error-messages/tool-errors/vectorizer-and-parallelizer-messages?view=msvc-170)
