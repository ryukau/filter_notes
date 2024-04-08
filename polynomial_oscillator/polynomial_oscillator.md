# 多項式オシレータ
[GlitchSprinkler](https://ryukau.github.io/UhhyouWebSynthesizers/GlitchSprinkler/synth.html) と [IntegerArpeggio](https://ryukau.github.io/UhhyouWebSynthesizers/IntegerArpeggio/synth.html) で使った多項式オシレータについて紹介します。

多項式とは以下の形の関数のことです。

$$
\begin{aligned}
P(x)
&= a_0 + a_1 x + a_2 x^2 + a_3 x^3 \dots \\
&= \sum_{n=0} a_n x^n
\end{aligned}
$$

- $x$ : 変数。
- $a_n$ : $n$ 次の項に乗算される定数。

## オシレータの仕様
楽に速く計算したいので以下のように仕様を定めます。

多項式から波形を取り出す範囲を $[0, 1)$ とします。これで位相のスケーリングを行わずに済みます。詳細は「[その他 -> 位相の計算](#位相の計算)」を参照してください。

波形の始点は 0 とします。これで多項式の `a_0` が常に 0 となり、計算量が減ります。

多項式はユーザが指定した $k$ 個の制御点を通過します。言い換えると、制御点によってユーザが波形を変更できるようにします。

## 指定された制御点を通る多項式の近似
以下は $k$ 個の制御点があるときに使える多項式です。

$$
\hat{P}(x) = \sum_{n=1}^{k} a_n x^n
$$

この節では、指定した制御点を通るように $a_n$ を求めます。

下の図のように制御点を定義します。

<figure>
<img src="img/control_point.svg" alt="" style="padding-bottom: 12px;"/>
</figure>

$x$ が小さい制御点から順にインデックスを与えて $x_0 < x_1 < \dots < x_{k+1}$ となるようにします。

インデックス $0$ と インデックス $k + 1$ の制御点は波形の両端です。インデックス $0$ の制御点の値を $(x_0, y_0) = (0, 0)$ 、インデックス $k + 1$ の制御点の値を $(x_{k+1}, y_{k+1}) = (1, 0)$ と固定します。 $x_0$ と $x_{k+1}$ の値は多項式から波形を取り出す範囲についての仕様に基づきます。 $y_0$ は波形の始点を 0 とする仕様から 0 になります。 $y_{k+1}$ はユーザが設定できるようにもできますが、簡略化のために 0 に固定しています。

ここまでの定義を使えば $\mathbf{A} \mathbf{x} = \mathbf{b}$ の形にして $a_n$ を求めることができます。

多項式 $\hat{P}$ はベクトルを使って以下のように書き直すことができます。

$$
\begin{aligned}
\hat{P}(x)
&= \sum_{n=1}^{k} a_n x^n
&=
\begin{bmatrix} x   &  x^2 &  x^3 &  \cdots &  x^k\end{bmatrix}
\begin{bmatrix} a_1 \\ a_2 \\ a_3 \\ \vdots \\ a_k\end{bmatrix}
\end{aligned}
$$

ここで $k$ 個の制御点の組 $(x_n, y_n)$ が分かっているので、以下の行列式を組み立てられます。

$$
\begin{bmatrix} y_1 \\ y_2 \\ y_3 \\ \vdots \\ y_k\end{bmatrix}
=
\begin{bmatrix}
  x_1    &  x_1^2 &  x_1^3 & \cdots &  x_1^k  \\
  x_2    &  x_2^2 &  x_2^3 &        &  x_2^k  \\
  x_3    &  x_3^2 &  x_3^3 &        &  x_3^k  \\
  \vdots &        &        & \ddots &  \vdots \\
  x_k    &  x_k^2 &  x_k^3 & \cdots &  x_k^k  \\
\end{bmatrix}
\begin{bmatrix} a_1 \\ a_2 \\ a_3 \\ \vdots \\ a_k\end{bmatrix}
$$

上の式は $\mathbf{A} \mathbf{x} = \mathbf{b}$ の形になっているので [`scipy.linalg.solve`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.linalg.solve.html) のようなソルバで解けます。

この式の組み立て方は [Fast fixed-point sine and cosine approximation with Julia – Nextjournal](https://nextjournal.com/zorn/fast-fixed-point-sine-and-cosine-approximation-with-julia) で紹介されていました。 sin と cos の多項式近似にも使えますが、ここで紹介している式の形だと精度はいまいちです。詳細はリンク先の記事や [Faster Math Functions](https://basesandframes.wordpress.com/2016/05/17/faster-math-functions/) などを参照してください。

今回の応用では $\mathbf{A}$ がそれほど大きくならないので、ソルバを自前で実装するときは[ガウスの消去法](https://ocw.mit.edu/courses/1-204-computer-algorithms-in-systems-engineering-spring-2010/ba1e1301f80d1f0351fa1746aff486e1_MIT1_204S10_lec20.pdf)を使えば間に合います。また [Wikipedia の LU decomposition の記事](https://en.wikipedia.org/wiki/LU_decomposition)にもソルバのコードが掲載されていますが、対角成分が 0 となるときのピボッティングなどが省略されているので注意してください。以下は IntegerArpeggio で使用した、JavaScript で書かれたソルバのコードへのリンクです。

- https://github.com/ryukau/UhhyouWebSynthesizers/blob/614b81f85c0a99f4a44a169003f890d29825d0ba/IntegerArpeggio/waveformxypad.js#L10-L59

## 振幅の正規化
以下の条件を満たすゲイン $g$ を求めて最大振幅を 1 に正規化します。

$$
\begin{aligned}
& \text{Find}      && g, \\
& \text{such that} && \max(|g P(x)|) = 1, \quad P(x) = \sum_{n=1}^{k} a_n x^n, \\
& \text{where}     && x \in [0, 1).
\end{aligned}
$$

この問題は $x \in [0, 1)$ の範囲における $P(x)$ のピークを探して、最大のピークの絶対値の逆数を $g$ とすれば解けます。

$P(x)$ のピークは、以下のように $P(x)$ を微分した式 $P'(x)$ が 0 になる点と定義します。

$$
\mathtt{peaks} = \left\{ P(x) \mid \forall x \enspace \text{such that} \enspace P'(x) = 0 \right\}, \quad
P'(x) = \frac{d P(x)}{dx}.
$$

一般的には root finding アルゴリズムに $P'(x)$ を渡せばピークが求められます。ただし root finding アルゴリズムは多数あり、選定も実装も手間です。そこで何とかならないかと試行錯誤していたところ、今回の応用には制御点があり、隣り合う制御点の間には最大でも 1 つのピークしか現れないことに気が付きました。以降ではこの観察に基づいて話を進めていきます。ただし、制御点の間のピークの数はあくまでも観察に基づく仮説であって、証明された事実ではないので例外があるかもしれません。

制御点の間のピークが 1 つ以下であれば、二分探索で安定して $P'(x)$ の解を探すことができます。以下に $g$ を求める手続きを書き下します。途中で出てくる $\mathrm{sgn}$ は[符号関数](https://en.wikipedia.org/wiki/Sign_function)です。

0. $P'(x)$ の定数 $a'_n$ をすべて求めて、任意の $x$ について計算できるようにする。
1. インデックスを $i = 0$ と定義。
2. 二分探索の初期状態として、ピークを求める区間の左端を $x_L = x_i$ 、右端を $x_R = x_{i+1}$ とする。
3. $P'(x_L)$ と $P'(x_R)$ の符号を比較。
   - 符号が等しければ $\max \big( |P(x_L)|, \, |P(x_R)| \big)$ をピークとして次の区間へ移動。 7 に進む。
   - 符号が等しくなければ 4 に進む。
4. 区間の中間点 $x_m = \dfrac{x_L + x_R}{2}$ を計算。
5. $P'(x_m)$ を計算する。
6. $P'(x_m)$ の符号を $P'(x_L)$ 、 $P'(x_R$ の符号と比較する。
   - $\mathrm{sgn}(P'(x_m)) = 0$ のときは $P(x_m)$ をピークとして、次の区間に移動。 7 に進む。
   - $\mathrm{sgn}(P'(x_m)) = \mathrm{sgn}(P'(x_L))$ のときは左端を $x_m$ として 2 に戻る。
   - $\mathrm{sgn}(P'(x_m)) = \mathrm{sgn}(P'(x_R))$ のときは右端を $x_m$ として 2 に戻る。
7. $i = i + 1$ としてインデックスを一つ進める。
   - $i < k$ なら 2 に戻る。
   - $i = k$ なら 8 に進む。
8. 得られたピークの中から最大のものを取り出して逆数を取ることで $g$ が得られる。

以下は上の手続きを JavaScript で実装したコードへのリンクです。

- https://github.com/ryukau/UhhyouWebSynthesizers/blob/58ecdb1bb85ee5747dfd9e6ea2aa232498e4dcab/IntegerArpeggio/waveformxypad.js#L168-L219

関数 $P'(x)$ はごく普通の多項式であり、特に精度が必要な応用でもないので、二分探索よりもニュートン法のほうが適しているかもしれません。

## その他
### 位相の計算
オシレータの位相の計算は様々なバリエーションがあるのでまとめておきます。コードは C++ です。

#### `floor`
おすすめの方法です。 `floor` を使う利点は `phase` が負の値になっても `[0, 1)` の範囲に丸められることです。

```c++
phase += frequencyHz / sampleRateHz;
phase -= std::floor(phase);
```

#### `fmod`
[`std::fmod`](https://en.cppreference.com/w/cpp/numeric/math/fmod) の利点は位相の範囲の上限を任意に設定できることです。欠点は割られる数 (1 つ目の引数) が負の値のときに出力も負の値となることです。したがって FM などによって負の周波数が出てくるときは対応が必要です。

```c++
phase = std::fmod(phase + frequencyHz / sampleRateHz, upperLimit);
if (frequencyHz < 0) phase = -phase;

// あるいは

phase = std::fmod(phase + frequencyHz / sampleRateHz, upperLimit);
phase *= std::copysign(float(1), frequencyHz);
```

また、 `sin` や `cos` の位相を回すときは `fmod` を使って `upperLimit = 2 * pi` と設定するとよさそうに思えますが、私が試した範囲では `floor` を使ったほうが速いです。

```c++
constexpr auto twoPi = float(2) * std::numbers::pi_v<float>;

// fmod
phase = std::fmod(phase + twoPi * frequencyHz / sampleRateHz, twoPi);
output = std::sin(phase);

// floor
phase += frequencyHz / sampleRateHz;
phase -= std::floor(phase);
output = std::sin(twoPi * phase);
```

#### `if` または `while`
位相の値が上限を超えたときに分岐を設ける方法があります。

```c++
phase += std::clamp(frequencyHz / sampleRateHz, float(0), float(0.5));
if (phase >= float(1)) phase -= float(1);

// あるいは

phase += std::max(frequencyHz / sampleRateHz, float(0));
while (phase >= float(1)) phase -= float(1);
```

`while` を使う方法は `frequencyHz` の大きさに比例して計算量が増えるので、使わないでください。似たようなコードを見たことがあるので一応紹介しています。

#### `remainder`
[`std::remainder`](https://en.cppreference.com/w/cpp/numeric/math/remainder) は引数が両方とも正の値のとき、出力が負の値となるので位相の回転には不適です。

### 符号なし整数
符号なし整数を使えば丸めの処理はコンパイラが何とかしてくれます。 DSP のコードがすべて整数演算なら高速です。浮動小数点数を使っているときはキャストやスケーリングが入るので速度面での利点は微妙です。

```c++
constexpr auto scaler = float(std::numeric_limits<unsigned>::max());
unsigned frequencyInt = unsigned(scaler * frequencyHz / sampleRateHz);
phase += frequencyInt;

// 使用例。
output = std::sin(twoPi * float(phase) / scaler);
```

標準ライブラリの [`<cinttypes>`](https://en.cppreference.com/w/cpp/header/cinttypes) に様々な大きさの整数が定義されています。

### 符号なし整数 + ビットマスク
波形をルックアップテーブルから読み出すときは符号なし整数とビットマスクを組み合わせる実装があります。ただし、ルックアップテーブルの長さは 2^n に制限されます。

```c++
constexpr auto scaler = unsigned(0xfffff); // 任意の 2^n - 1 の定数。
unsigned frequencyInt = unsigned(float(scaler) * frequencyHz / sampleRateHz);
phase += frequencyInt;
phase &= scaler;

// 使用例。
output = wavetable[phase];
```

### 固定小数点数
以下は前述の `floor` による位相の回転を固定小数点数で実装した例です。

```c++
// 定義。
constexpr auto fractionBits = int32_t(16);
constexpr auto fractionMult = int32_t(1) << fractionBits;
constexpr auto fractionMask = fractionMult - int32_t(1);
int32_t phase = 0;
int32_t sampleRateFixed = int32_t(sampleRateHz) << fractionBits;
int32_t frequencyFixed = int32_t(frequencyHz) << fractionBits;

auto div = [](int32_t x, int32_t y) -> int32_t {
  return int32_t((int64_t(x) * fractionMult) / y);
};

// 位相の回転。
phase += div(frequencyFixed, sampleRateFixed);
phase &= fractionMask; // phase -= floor(phase) と同じ意味。
```

上記のコードの固定小数点数のフォーマットは符号 1 ビット、整数部 15 ビット、小数部 16 ビットで、 Q15.16 あるいは Q16.16 と表記されます。 Qm.n という表記は Q フォーマットと呼ばれ、小数点より上の桁数が m ビット 、小数点より下の桁数が n ビットの固定小数点数という意味です。

Q フォーマットは [Texas Instruments の資料](https://www.ti.com/lit/ug/spru565b/spru565b.pdf) (Appendix A.2, p.A-3) と [ARM の資料](https://web.archive.org/web/20171104105632/http://infocenter.arm.com/help/topic/com.arm.doc.dui0066g/DUI0066.pdf) (4.7.9 Q-format, p.4-24) に掲載されていますが、 m に符号ビットを含むかどうかで違いがあります。 Texas Instruments では符号ビットを含まない、 ARM では符号ビットを含む、としています。

固定小数点数を自前で実装すると数学関数の用意が大変なので、特に理由がなければ fpm のようなライブラリを使うことを推奨します。

- [GitHub - MikeLankamp/fpm: C++ header-only fixed-point math library](https://github.com/MikeLankamp/fpm)

## 参考文献
- [Fast fixed-point sine and cosine approximation with Julia – Nextjournal](https://nextjournal.com/zorn/fast-fixed-point-sine-and-cosine-approximation-with-julia)
- [OlliW's Bastelseiten » Fast Computation of Functions on Microcontrollers](http://www.olliw.eu/2014/fast-functions/)
- [1.204 Lecture 20, Linear systems - ba1e1301f80d1f0351fa1746aff486e1_MIT1_204S10_lec20.pdf](https://ocw.mit.edu/courses/1-204-computer-algorithms-in-systems-engineering-spring-2010/ba1e1301f80d1f0351fa1746aff486e1_MIT1_204S10_lec20.pdf)
- [LU decomposition - Wikipedia](https://en.wikipedia.org/wiki/LU_decomposition)
- [TMS320C64x DSP Library Programmer's Reference (Rev. B) - spru565b.pdf](https://www.ti.com/lit/ug/spru565b/spru565b.pdf)
- [Wayback Machine](https://web.archive.org/web/20171104105632/http://infocenter.arm.com/help/topic/com.arm.doc.dui0066g/DUI0066.pdf)
