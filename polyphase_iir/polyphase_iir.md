# Polyphase IIR フィルタ
ダウンサンプリングにつかう IIR (infinite impulse response) フィルタをポリフェイズ (polyphase) 化して direct form による実装と比べます。

## 理論
IIR フィルタの伝達関数から始めます。

$$
H(z) = \frac{B(z)}{A(z)} = \frac{\sum_{i=0}^{L} b_i z^{-i}}{1 + \sum_{i=1}^{N} a_i z^{-i}}.
$$

分母の多項式 $A$ を因数分解します。 $p_i$ は $H$ の極です。 $H$ は有理関数なので、その極とは $A(z)=0$ の解のことです。

$$
A(z) = \prod_{i=1}^{N} (1 - p_i z^{-1}).
$$

$z^M$ にあわせて極をスケーリングした伝達関数 $A_M$ を定義します。 $M$ はダウンサンプリングの倍率です。

$$
A_M(z^M) = \prod_{i=1}^{N} (1 - p_i^M z^{-M}).
$$

伝達関数に $A_M$ をねじこみます。

$$
H(z)
= \frac{B(z)}{A(z)} \cdot \frac{A_M(z^{-M})}{A_M(z^{-M})}
= \frac{Q(z)}{A_M(z^{-M})}, \quad \text{where} \enspace Q(z) = B(z) \dfrac{A_M(z^{-M})}{A(z)}.
$$

ここで $Q$ について以下の変形ができます。

$$
\begin{aligned}
B(z) \frac{A_M(z^M)}{A(z)}
&= B(z) \prod_{i=1}^{N} \frac{1 - (p_i z^{-1})^M}{1 - p_i z^{-1}} \\
&= B(z) \prod_{i=1}^{N} \frac{
  \cancel{(1 - (p_i z^{-1}))}(1 + (p_i z^{-1}) + (p_i z^{-1})^2 + (p_i z^{-1})^3 + \dots + (p_i z^{-1})^{M-1})
}{
  \cancel{1 - p_i z^{-1}}
} \\
&= B(z) \prod_{i=1}^{N} \left( \sum_{k=0}^{M-1} p_i^k z^{-k} \right).
\end{aligned}
$$

分母が 1 となって消えたので FIR です。 FIR である $Q$ はポリフェイズ分解できます。

$$
Q(z) = \sum_{k=0}^{M-1} z^{-k} Q_k(z^M).
$$

ここでは FIR のポリフェイズ分解の詳細は省略します。ざっくり言えばインデックス $n$ について $z^{-(nM+k)}$ の項を集めたものがポリフェイズ分解後の FIR となる $Q_k$ です。例えば $M=2$ で係数を $q_n$ とすると、 $Q_0 = q_0 z + q_2 z^{-2} + q_4 z^{-4}, \dots$ 、 $Q_1 = q_1 z^{-1} + q_3 z^{-3} + q_5 z^{-5}, \dots$ となります。

$H$ に戻ります。

$$
H(z) = \sum_{k=0}^{M-1} z^{-k} \frac{Q_k(z^M)}{A_M(z^M)}.
$$

分母を含めてポリフェイズ分解できています。

## 設計
プロトタイプとなる [zpk 形式の伝達関数](https://www.mathworks.com/help/control/ref/zpk.html)から $A_M$ と $Q_k$ を求めれば計算できる形になります。 zpk は MATLAB や SciPy で使われている伝達関数のゼロ、極、ゲイン (zero, pole, gain) の組です。 SciPy では [`scipy.signal.butter`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.signal.butter.html) や [`scipy.signal.ellip`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.signal.ellip.html) の引数に `output="zpk"` を指定すれば zpk の組が得られます。

「[実装](#実装)」の節にコード例と、完全な設計コードへのリンクを掲載しています。

## 計算方法
高速で正確な計算方法を検討します。計算方法のバリエーションには以下の軸があります。

- 軸 1 : フィルタ係数の表現をすべて 2 次セクション ([sos](https://docs.scipy.org/doc/scipy/tutorial/signal.html#second-order-sections-representation)) にするか、分子と分母を分けて計算 (ハイブリッド) するか。
- 軸 2 : 単純加算するか、 [Kahan summation](https://en.wikipedia.org/wiki/Kahan_summation_algorithm) を使うか。

以降では、軸 1 の分子と分母を分けて計算する方法のことをハイブリッド (hybrid) と呼ぶことにします。ハイブリッド形式では分子を FIR 、分母を sos として計算します。分子を別に計算するため、分母の sos は 1 セクションあたりの係数を $a_1, a_2$ の 2 つに減らせます。

軸 2 の Kahan summation はポリフェイズ出力と、ハイブリッドの FIR の計算に適用できます。単純加算に比べると遅くなりますが、その遅さが許容範囲内か、また意味のある正確さの向上が得られるかについて調べます。

以下は C++20 での実装へのリンクです。リンク先のコードには、いくつか FMA の無しの実装が含まれています。この文章では FMA 有りの実装のみについて検討します。

- [filter_notes/polyphase_iir/polyphaseiir.hpp at master · ryukau/filter_notes · GitHub](https://github.com/ryukau/filter_notes/blob/master/polyphase_iir/polyphaseiir.hpp)

### テスト結果
数が多いので折りたたんでいます。

##### Butterworth, order = 4, M = 2, fc/fs = 0.125

<details>
<summary>速度と正確さの比較</summary>

###### `cl /std:c++20 /O2 /EHsc /arch:AVX2`, Microsoft (R) C/C++ Optimizing Compiler Version 19.51.36248 for x64

<figure>
<img src="img/quality__4_2_cl.svg" alt="" style="padding-bottom: 12px;"/>
</figure>

###### `clang++ -std=c++20 -O3 -mfma -march=native`, clang 22.1.3 x86_64-pc-windows-msvc

<figure>
<img src="img/quality__4_2_clang++.svg" alt="" style="padding-bottom: 12px;"/>
</figure>

###### `g++ -std=c++20 -O3 -mfma -march=native`, g++ (GCC) 16.1.1 20260515 (Red Hat 16.1.1-2)

<figure>
<img src="img/quality__4_2_g++.svg" alt="" style="padding-bottom: 12px;"/>
</figure>

</details>

<details>
<summary>周波数特性</summary>

###### `cl /std:c++20 /O2 /EHsc /arch:AVX2`, Microsoft (R) C/C++ Optimizing Compiler Version 19.51.36248 for x64

<figure>
<img src="img/response__4_2_cl.svg" alt="" style="padding-bottom: 12px;"/>
</figure>

###### `clang++ -std=c++20 -O3 -mfma -march=native`, clang 22.1.3 x86_64-pc-windows-msvc

<figure>
<img src="img/response__4_2_clang++.svg" alt="" style="padding-bottom: 12px;"/>
</figure>

###### `g++ -std=c++20 -O3 -mfma -march=native`, g++ (GCC) 16.1.1 20260515 (Red Hat 16.1.1-2)

<figure>
<img src="img/response__4_2_g++.svg" alt="" style="padding-bottom: 12px;"/>
</figure>

</details>

##### Butterworth, order = 8, M = 4, fc/fs = 0.03125

<details>
<summary>速度と正確さの比較</summary>

###### `cl /std:c++20 /O2 /EHsc /arch:AVX2`, Microsoft (R) C/C++ Optimizing Compiler Version 19.51.36248 for x64

<figure>
<img src="img/quality__8_4_cl.svg" alt="" style="padding-bottom: 12px;"/>
</figure>

###### `clang++ -std=c++20 -O3 -mfma -march=native`, clang 22.1.3 x86_64-pc-windows-msvc

<figure>
<img src="img/quality__8_4_clang++.svg" alt="" style="padding-bottom: 12px;"/>
</figure>

###### `g++ -std=c++20 -O3 -mfma -march=native`, g++ (GCC) 16.1.1 20260515 (Red Hat 16.1.1-2)

<figure>
<img src="img/quality__8_4_g++.svg" alt="" style="padding-bottom: 12px;"/>
</figure>

</details>

<details>
<summary>周波数特性</summary>

###### `cl /std:c++20 /O2 /EHsc /arch:AVX2`, Microsoft (R) C/C++ Optimizing Compiler Version 19.51.36248 for x64

<figure>
<img src="img/response__8_4_cl.svg" alt="" style="padding-bottom: 12px;"/>
</figure>

###### `clang++ -std=c++20 -O3 -mfma -march=native`, clang 22.1.3 x86_64-pc-windows-msvc

<figure>
<img src="img/response__8_4_clang++.svg" alt="" style="padding-bottom: 12px;"/>
</figure>

###### `g++ -std=c++20 -O3 -mfma -march=native`, g++ (GCC) 16.1.1 20260515 (Red Hat 16.1.1-2)

<figure>
<img src="img/response__8_4_g++.svg" alt="" style="padding-bottom: 12px;"/>
</figure>

</details>

##### Butterworth, order = 12, M = 4, fc/fs = 0.02

<details>
<summary>速度と正確さの比較</summary>

###### `cl /std:c++20 /O2 /EHsc /arch:AVX2`, Microsoft (R) C/C++ Optimizing Compiler Version 19.51.36248 for x64

<figure>
<img src="img/quality__12_4_cl.svg" alt="" style="padding-bottom: 12px;"/>
</figure>

###### `clang++ -std=c++20 -O3 -mfma -march=native`, clang 22.1.3 x86_64-pc-windows-msvc

<figure>
<img src="img/quality__12_4_clang++.svg" alt="" style="padding-bottom: 12px;"/>
</figure>

###### `g++ -std=c++20 -O3 -mfma -march=native`, g++ (GCC) 16.1.1 20260515 (Red Hat 16.1.1-2)

<figure>
<img src="img/quality__12_4_g++.svg" alt="" style="padding-bottom: 12px;"/>
</figure>

</details>

<details>
<summary>周波数特性</summary>

###### `cl /std:c++20 /O2 /EHsc /arch:AVX2`, Microsoft (R) C/C++ Optimizing Compiler Version 19.51.36248 for x64

<figure>
<img src="img/response__12_4_cl.svg" alt="" style="padding-bottom: 12px;"/>
</figure>

###### `clang++ -std=c++20 -O3 -mfma -march=native`, clang 22.1.3 x86_64-pc-windows-msvc

<figure>
<img src="img/response__12_4_clang++.svg" alt="" style="padding-bottom: 12px;"/>
</figure>

###### `g++ -std=c++20 -O3 -mfma -march=native`, g++ (GCC) 16.1.1 20260515 (Red Hat 16.1.1-2)

<figure>
<img src="img/response__12_4_g++.svg" alt="" style="padding-bottom: 12px;"/>
</figure>

</details>

##### Butterworth, order = 16, M = 8, fc/fs = 0.015625

<details>
<summary>速度と正確さの比較</summary>

###### `cl /std:c++20 /O2 /EHsc /arch:AVX2`, Microsoft (R) C/C++ Optimizing Compiler Version 19.51.36248 for x64

<figure>
<img src="img/quality__16_8_cl.svg" alt="" style="padding-bottom: 12px;"/>
</figure>

###### `clang++ -std=c++20 -O3 -mfma -march=native`, clang 22.1.3 x86_64-pc-windows-msvc

<figure>
<img src="img/quality__16_8_clang++.svg" alt="" style="padding-bottom: 12px;"/>
</figure>

###### `g++ -std=c++20 -O3 -mfma -march=native`, g++ (GCC) 16.1.1 20260515 (Red Hat 16.1.1-2)

<figure>
<img src="img/quality__16_8_g++.svg" alt="" style="padding-bottom: 12px;"/>
</figure>

</details>

<details>
<summary>周波数特性</summary>

###### `cl /std:c++20 /O2 /EHsc /arch:AVX2`, Microsoft (R) C/C++ Optimizing Compiler Version 19.51.36248 for x64

<figure>
<img src="img/response__16_8_cl.svg" alt="" style="padding-bottom: 12px;"/>
</figure>

###### `clang++ -std=c++20 -O3 -mfma -march=native`, clang 22.1.3 x86_64-pc-windows-msvc

<figure>
<img src="img/response__16_8_clang++.svg" alt="" style="padding-bottom: 12px;"/>
</figure>

###### `g++ -std=c++20 -O3 -mfma -march=native`, g++ (GCC) 16.1.1 20260515 (Red Hat 16.1.1-2)

<figure>
<img src="img/response__16_8_g++.svg" alt="" style="padding-bottom: 12px;"/>
</figure>

</details>

#### 所見
Butterworth フィルタについては、どの polyphase IIR の実装であっても速度、正確さ共に transposed direct form II (TDF-II) を超える性能と言えそうです。

Hybrid Simple が速いです。 Hybrid Simple は伝達関数の分子を FIR として計算しており、また FIR 内で Kahan summation を使っていない実装です。

SOS Simple が二番目に早いです。 `clang++` では order = 4 かつ f64 のフィルタで Hybrid Simple を上回る性能になっています。

Kahan summation については予想通り遅くなっています。正確さの向上については大きく向上するとは言い難く、一部のケースを除いて SOS Simple のほうが正確に見えます。

f64 と f32 で実行速度がほぼ同じです。 IIR フィルタと Kahan summation は 1 サンプル前の計算結果を利用するので自動ベクトル化が入らないことが原因と考えられます。 Hybrid Simple は FIR 部分で自動ベクトル化が入っていると捉えれば納得できる結果と言えます (`g++` では確認) 。複数の信号を SIMD で並列計算するときは f32 のほうが速くなると予想されます。

f32 では振幅特性がざっくりと 120 dB を下回る範囲で、周波数特性が歪んでいます。周波数特性のプロットで暴れているのはすべて f32 を表す点線です。

### 係数の数
以下は軸 1 のフィルタ係数の表現の違いによる係数の数を比較した表です。ハイブリッド形式のほうがフィルタ係数の数が少なくなる傾向があります。 SOS 形式については $a_0=1$ とすることで、 1 セクションあたりの係数の数を 5 として数えています。

<details>
<summary>ハイブリッド形式の係数の数の表</summary>

↓ Order \ M → | 1 |  2 |  3 |  4 |  5 |  6 |   7 |   8 |  9 |  10 |  11 |  12 |  13 |  14 |  15 |  16
-------------:|--:|---:|---:|---:|---:|---:|----:|----:|---:|----:|----:|----:|----:|----:|----:|---:
1             |   |    |    |    |    |    |     |     |    |     |     |     |     |     |     |    
2             |   |  7 |    |    |    |    |     |     |    |     |     |     |     |     |     |    
3             |   |    | 14 |    |    |    |     |     |    |     |     |     |     |     |     |    
4             |   | 13 |    | 21 |    |    |     |     |    |     |     |     |     |     |     |    
5             |   |    |    |    | 32 |    |     |     |    |     |     |     |     |     |     |    
6             |   | 19 | 25 |    |    | 43 |     |     |    |     |     |     |     |     |     |    
7             |   |    |    |    |    |    |  58 |     |    |     |     |     |     |     |     |    
8             |   | 25 |    | 41 |    |    |     |  73 |    |     |     |     |     |     |     |    
9             |   |    | 38 |    |    |    |     |     | 92 |     |     |     |     |     |     |    
10            |   | 31 |    |    | 61 |    |     |     |    | 111 |     |     |     |     |     |    
11            |   |    |    |    |    |    |     |     |    |     | 134 |     |     |     |     |    
12            |   | 37 | 49 | 61 |    | 85 |     |     |    |     |     | 157 |     |     |     |    
13            |   |    |    |    |    |    |     |     |    |     |     |     | 184 |     |     |    
14            |   | 43 |    |    |    |    | 113 |     |    |     |     |     |     | 211 |     |    
15            |   |    | 62 |    | 92 |    |     |     |    |     |     |     |     |     | 242 |    
16            |   | 49 |    | 81 |    |    |     | 145 |    |     |     |     |     |     |     | 273

</details>

<details>
<summary>SOS 形式の係数の数の表</summary>

Order \ M | 1 |  2 |   3 |   4 |   5 |   6 |   7 |   8 |   9 |  10 |  11 |  12 |  13 |  14 |  15 |  16
---------:|--:|---:|----:|----:|----:|----:|----:|----:|----:|----:|----:|----:|----:|----:|----:|---:
1         |   |    |     |     |     |     |     |     |     |     |     |     |     |     |     |    
2         |   | 10 |     |     |     |     |     |     |     |     |     |     |     |     |     |    
3         |   |    |  30 |     |     |     |     |     |     |     |     |     |     |     |     |    
4         |   | 20 |     |  40 |     |     |     |     |     |     |     |     |     |     |     |    
5         |   |    |     |     |  75 |     |     |     |     |     |     |     |     |     |     |    
6         |   | 30 |  45 |     |     |  90 |     |     |     |     |     |     |     |     |     |    
7         |   |    |     |     |     |     | 140 |     |     |     |     |     |     |     |     |    
8         |   | 40 |     |  80 |     |     |     | 160 |     |     |     |     |     |     |     |    
9         |   |    |  75 |     |     |     |     |     | 225 |     |     |     |     |     |     |    
10        |   | 50 |     |     | 125 |     |     |     |     | 250 |     |     |     |     |     |    
11        |   |    |     |     |     |     |     |     |     |     | 330 |     |     |     |     |    
12        |   | 60 |  90 | 120 |     | 180 |     |     |     |     |     | 360 |     |     |     |    
13        |   |    |     |     |     |     |     |     |     |     |     |     | 455 |     |     |    
14        |   |    |     |     |     |     | 245 |     |     |     |     |     |     | 490 |     |    
15        |   |    | 120 |     | 200 |     |     |     |     |     |     |     |     |     | 600 |    
16        |   |    |     | 160 |     |     |     | 320 |     |     |     |     |     |     |     | 640

</details>

## 実装
以下は Python による設計コードと C++20 による Hybrid Simple のフィルタ実装です。

<details>
<summary>Python によるフィルタ設計の実装例</summary>

次数が高いフィルタは正確に設計できないことがあるので mpmath を使っています。 `tf2sos_mp` を省略しているので、そのままでは動きません。

```python
from mpmath import mp

def design_polyphase_iir(
    z_mp,
    p_mp,
    k_mp,
    M: int,
    output: str = "ba",
    as_float: bool = False,
    workdps: int | None = None,
):
    """
    Decomposes an arbitrary IIR filter (given by its zeros, poles, and gain)
    for polyphase down-sampling.

    Parameters:
    -----------
    z_mp : list or iterable
        Zeros of the filter as mpmath numbers.
    p_mp : list or iterable
        Poles of the filter as mpmath numbers.
    k_mp : mpmath number
        System gain.
    M : int
        Down-sampling factor.
    output : str, optional
        Format of the output filter: "ba", "sos", or "hybrid".
    as_float : bool, optional
        Cast the coefficients to Python float.
    workdps : int, optional
        Working precision in decimal places (dps). If None, the current
        ambient mpmath precision context is used.

    Returns:
    --------
    tuple or list
        The returned structure depends on the value of the `output` parameter:
        * If "ba" (default):
          Returns a tuple `(q_polyphase, a_low_poly)` where:
            - `q_polyphase` (list of lists): Numerator coefficients for each of the `M` polyphase branches.
            - `a_low_poly` (list): Denominator coefficients shared by all branches (expressed in powers of z^-M).
        * If "sos":
          Returns `sos_polyphase` (list of lists of lists): A list of length `M` containing the Second-Order Section (SOS) representations for each polyphase branch.
        * If "hybrid":
          Returns a tuple `(q_polyphase, sos_sections)` where:
            - `q_polyphase` (list of lists): Numerator coefficients for each of the `M` polyphase branches.
            - `sos_sections` (list of lists): Denominator coefficients `[a1, a2]` for each pole pair.

        Coefficients are returned as standard Python floats if `as_float` is True; otherwise, they are `mpmath` types.
    """
    if output not in ("ba", "sos", "hybrid"):
        raise ValueError("output must be one of 'ba', 'sos', or 'hybrid'")

    with mp.workdps(workdps):

        def poly_from_roots(roots):
            """
            Expands (x - r1)(x - r2)... into polynomial coefficients.
            Returns coefficients in descending order of powers (highest power first).
            """
            c = [mp.mpf(1)]
            for r in roots:
                c += [0]
                for i in range(len(c) - 1, 0, -1):
                    c[i] -= r * c[i - 1]
            return c

        def mp_convolve(a, b):
            len_a = len(a)
            len_b = len(b)
            out = [mp.mpf(0.0)] * (len_a + len_b - 1)
            for i in range(len_a):
                for j in range(len_b):
                    out[i + j] += a[i] * b[j]
            return out

        def mp_deconvolve(num, den):
            """
            Polynomial division: num(z) / den(z).
            Returns (quotient, remainder).
            Assumes coefficients are in descending order of powers (highest power first).
            """
            div = list(num)
            sz = len(div) - len(den) + 1
            if sz <= 0:
                return [mp.mpf(0)], div

            q = [mp.mpf(0)] * sz
            for k in range(sz):
                q[k] = div[k] / den[0]
                for j, val in enumerate(den):
                    div[k + j] -= q[k] * val

            return q, div

        b_poly = [mp.re(x * k_mp) for x in poly_from_roots(z_mp)]
        a_poly = [mp.re(x) for x in poly_from_roots(p_mp)]

        p_new = [p**M for p in p_mp]
        a_low_poly = [mp.re(x) for x in poly_from_roots(p_new)]

        deg_low = len(a_low_poly) - 1
        a_high_poly = [mp.mpf(0.0)] * (deg_low * M + 1)
        for i, val in enumerate(a_low_poly):
            a_high_poly[i * M] = val

        s_poly, remainder = mp_deconvolve(a_high_poly, a_poly)
        q_poly = mp_convolve(b_poly, s_poly)

        q_polyphase = []
        for k in range(M):
            q_polyphase.append(q_poly[k::M])

        if output == "ba":
            if as_float:
                q_polyphase = apply(q_polyphase, float)
                a_low_poly = apply(a_low_poly, float)
            return q_polyphase, a_low_poly

        elif output == "sos":
            sos_polyphase = []
            for k in range(M):
                q_k = q_poly[k::M]
                sos_section = tf2sos_mp(q_k, a_low_poly)
                sos_polyphase.append(sos_section)

            if as_float:
                sos_polyphase = apply(sos_polyphase, float)
            return sos_polyphase

        elif output == "hybrid":
            sos_sections = []
            pool = list(p_new)
            while len(pool) > 0:
                p1 = pool.pop(0)

                if abs(mp.im(p1)) < 1e-20:
                    a1 = -mp.re(p1)
                    a2 = mp.mpf(0.0)
                    sos_sections.append([a1, a2])
                    continue

                best_idx = -1
                best_err = mp.mpf("inf")

                for i, p2 in enumerate(pool):
                    err = abs(p1 - mp.conj(p2))
                    if err < best_err:
                        best_err = err
                        best_idx = i

                if best_idx != -1:
                    p2 = pool.pop(best_idx)
                    a1 = -2 * mp.re(p1)
                    a2 = abs(p1) ** 2
                    sos_sections.append([a1, a2])
                else:
                    a1 = -mp.re(p1)
                    a2 = mp.mpf(0.0)
                    sos_sections.append([a1, a2])

            if as_float:
                q_polyphase = apply(q_polyphase, float)
                sos_sections = apply(sos_sections, float)
            return q_polyphase, sos_sections
```

</details>

<details>
<summary>C++20 による Hybrid Simple の実装</summary>

```c++
#pragma once

#include <array>
#include <cmath>
#include <concepts>
#include <cstddef>
#include <ranges>
#include <tuple>
#include <utility>

// Hybrid 形式の polyphase IIR の係数。 (Butterworth, order=8, cutoff=0.3125, M=4)
template<typename T> struct Hybrid {
  static constexpr int nPhase = 4; // 文中の M に対応。

  // 分母の SOS: {a1, a2} for 1 / (1 + a1*z^-1 + a2*z^-2)
  static constexpr std::array<std::array<T, 2>, 4> denom{{
    { T(7.60834734050858330e-01), T(2.70020217408328600e-01) },
    { T(2.45817342985895443e-01), T(1.83487154523652787e-02) },
    { T(6.21709429503566474e-02), T(1.10913984116065879e-03) },
    { T(-8.69548164131942261e-03), T(1.06677160707193438e-04) },
  }};

  // 分子の FIR.
  static constexpr std::array<std::array<T, 9>, nPhase> branches{{
    { T(4.67360371460534638e-04), T(1.78965525913988399e-01), T(2.66265704737198072e-01), T(1.82695165925812308e-01), T(3.59985359464480154e-02), T(3.47770758488648735e-03), T(9.76101633490265749e-05), T(5.13705983926379128e-07), T(5.56796711341444609e-11) },
    { T(5.13399690923943292e-03), T(2.82921645103193231e-01), T(2.25756981966161091e-01), T(1.35987195285923984e-01), T(2.25319762725637897e-02), T(1.58374704330561684e-03), T(3.19637728170234209e-05), T(9.14850632661454864e-08), T(0.00000000000000000e+00) },
    { T(2.62245483078802598e-02), T(3.38960281490826409e-01), T(2.16172793820407716e-01), T(8.99252575894855938e-02), T(1.32404246635418957e-02), T(6.76027869152339453e-04), T(9.30032662397550523e-06), T(1.25270724575848982e-08), T(0.00000000000000000e+00) },
    { T(8.25615432363466933e-02), T(3.22027248938468458e-01), T(2.10495757082402951e-01), T(5.67774689416133613e-02), T(7.09628583197262353e-03), T(2.68501683273136621e-04), T(2.36627742844577689e-06), T(1.17002404623042308e-09), T(0.00000000000000000e+00) },
  }};
};

template<std::floating_point T, auto Coefficients> class FirSection {
private:
  static constexpr std::size_t N = Coefficients.size();

  alignas(64) std::array<T, (N > 0 ? N - 1 : 0)> s{};

public:
  void reset() { s.fill(T(0)); }

  inline T process(T input) {
    if constexpr (N == 0) { return T(0); }

    constexpr auto b0 = Coefficients[0];
    if constexpr (N == 1) { return b0 * input; }

    T y = std::fma(b0, input, s[0]);
    for (std::size_t i = 0; i < N - 2; ++i) {
      s[i] = std::fma(Coefficients[i + 1], input, s[i + 1]);
    }
    s[N - 2] = Coefficients[N - 1] * input;
    return y;
  }
};

template<std::floating_point T, auto Coefficients> class IirSection {
private:
  static constexpr std::size_t nSections = Coefficients.size();
  std::array<std::array<T, 2>, nSections> s{};

  template<std::size_t... I> inline void process_sections(T& val, std::index_sequence<I...>) {
    ((process_one_section<I>(val)), ...);
  }

  template<std::size_t Index> inline void process_one_section(T& val) {
    T y_out = val + s[Index][0];
    s[Index][0] = std::fma(-Coefficients[Index][0], y_out, s[Index][1]);
    s[Index][1] = -Coefficients[Index][1] * y_out;
    val = y_out;
  }

public:
  void reset() {
    for (auto& sec : s) { sec.fill(T(0)); }
  }

  inline T process(T input) {
    if constexpr (nSections == 0) { return input; }
    T y = input;
    process_sections(y, std::make_index_sequence<nSections>{});
    return y;
  }
};

template<std::floating_point T, typename Coefficients> class PolyphaseIir {
public:
  static constexpr int nPhase = Coefficients::nPhase;

private:
  template<std::size_t... I> static auto make_fir_tuple(std::index_sequence<I...>) {
    return std::tuple<FirSection<T, Coefficients::branches[I]>...>{};
  }
  using FirTuple = decltype(make_fir_tuple(std::make_index_sequence<nPhase>{}));

  FirTuple fir_branches;
  IirSection<T, Coefficients::denom> iir_filter;

  template<std::size_t Index> inline T process_branch(T input) {
    return std::get<Index>(fir_branches).process(input);
  }

  template<std::size_t... I>
  inline T sum_branches(const std::array<T, nPhase>& inputs, std::index_sequence<I...>) {
    return (... + process_branch<I>(inputs[I]));
  }

  template<std::size_t... I> inline void reset_branches(std::index_sequence<I...>) {
    ((std::get<I>(fir_branches).reset()), ...);
  }

public:
  void reset() {
    reset_branches(std::make_index_sequence<nPhase>{});
    iir_filter.reset();
  }

  inline T process(const std::array<T, nPhase>& inputs) {
    T fir_sum = sum_branches(inputs, std::make_index_sequence<nPhase>{});
    return iir_filter.process(fir_sum);
  }
};

// 使用例。
template<std::floating_point T> auto down_sample(std::vector<T>& inputs) {
  PolyphaseIir<T, Hybrid<T>> filter;

  constexpr size_t M = Hybrid<T>::nPhase;
  constexpr size_t n_sample = inputs.size() / M;

  std::vector<T> outputs;
  outputs.reserve(n_sample);

  for (size_t n = 0; n < n_sample; n += M) {
    std::array<T, M> frame;
    for (size_t k = 0; k < M; ++k) { frame[k] = inputs[n + k]; }
    outputs.push_back(filter.process(frame));
  }

  return outputs;
}
```

</details>

完全な形のコードを以下のリンク先に掲載しています。

- [filter_notes/polyphase_iir at master · ryukau/filter_notes · GitHub](https://github.com/ryukau/filter_notes/tree/master/polyphase_iir)

Python による完全な設計コードはリンク先の `design.py` と `signal_mp.py` に書かれています。 `design.py` は `design_polyphase_iir` を含むメインの設計コード、 `signal_mp.py` は主に SciPy の実装を mpmath に移植したものです。

C++20 による他のフィルタ形式の実装は `polyphaseiir.hpp` に書かれています。

## 参考文献
- [Msps | Analog Devices](https://www.analog.com/en/resources/glossary/msps.html)
