## 放物線アタックと指数曲線ディケイの組み合わせの別解
位置 $\hat{y}_A(t)$ の式です。

$$
\hat{y}_A(t) = \begin{cases}
e^{\gamma t} a_A \dfrac{t^2}{2}, & \text{if} \enspace 0 \leq t \leq t_A, \\
e^{\gamma t} \left(- b_A \dfrac{t^2}{2} + \hat{v}_A t + \hat{h}_A \right), & \text{if} \enspace t_A \leq t \leq t_p, \\
e^{\gamma t},                    & \text{if} \enspace t_p \leq t. \\
\end{cases}
$$

- $t$ : エンベロープがトリガされてからの経過時間
- $\gamma$ : 指数曲線の減衰の速さ
- $a_A$ : 放物線の加速時の加速度の絶対値
- $b_A$ : 放物線の減速時の加速度の絶対値
- $t_A$ : 放物線の加減速を切り替える時点
- $t_p$ : 放物線がピークに到達する時点
- $t_k$ : $\hat{y}$ がピークに到達する時点
- $\hat{v}_A$ : 時点 $t_A$ での速度
- $\hat{h}_A$ : 時点 $t_A$ での位置

$\gamma$ は 0 以下です。

次の手順でピークに到達する時点 $t_k$ のある範囲を調べます。

- $t_A$ の時点で速度が負の値なら $0 \leq t_k < t_A$ 。
- そうでなければ $t_A \leq t_k$ 。

区間 $0 \leq t \leq t_A$ での速度 $\hat{\dot{y}}_A(t)$ を求めます。 Maxima を使います。

```maxima
y1_d1: diff(exp(γ * t) * a_A * t^2 / 2, t);
```

$$
\begin{aligned}
\dot{\hat{y}}_A(t)
  &= \frac{a_A t^2 \gamma e^{t \gamma}}{2} + a_A t e^{t \gamma} \\
  &= e^{\gamma t} a_A t \left( \frac{\gamma t}{2} + 1 \right), \qquad 0 \leq t \leq t_A.
\end{aligned}
$$

$\hat{v}_A$ と $\hat{h}_A$ の式です。

$$
\begin{aligned}
\hat{v}_A &= \dot{\hat{y}}_A(t_A)
  = e^{\gamma t_A} a_A t_A \left( \frac{\gamma t_A}{2} + 1 \right) \\
\hat{h}_A &= \hat{y}_A(t_A)
  = e^{\gamma t_A} a_A \dfrac{t_A^2}{2} \\
\end{aligned}
$$

$\hat{v}_A$ が負の値なら $\hat{y}$ のピークは $0 \leq t < t_A$ のどこかにあります。負の値になるのは $\gamma$ だけなので、 $\left( \dfrac{\gamma t_A}{2} + 1 \right)$ が 0 以下になるかどうかを調べれば十分です。判定式を変形すると次の形になります。

$$
\gamma t_A \leq -2
$$

$\dot{\hat{y}}_A(t) = 0$ となる時点 $t_{k, A}$ を求めます。

```maxima
solve(0 = y1_d1, t);
```

2 つの解が得られますが、 $t = 0$ は使いません。エンベロープをトリガした瞬間の速度が 0 であることは事前に分かっているからです。よって $t_{k, A}$ は次の式になります。

$$
t_{k, A} = - \frac{2}{\gamma}
$$

区間 $t_A \leq t \leq t_p$ で速度 $\hat{\dot{y}}_A(t)$ が 0 になる時点 $t_{k,p}$ を求めます。

```maxima
y2_d1: diff(exp(γ * t) * (-b_A * t^2 / 2 + v_A * t + h_A), t);
solve(0 = y2_d1, t);
```

2 つの解が得られました。

$$
\begin{aligned}
t &= -\frac{\sqrt{(v_A^2 + 2 b_A h_A) \gamma^2 + {b_A}^2} - v_A \gamma + b_A}{b_A \gamma} \\
t &= \frac{\sqrt{(v_A^2 + 2 b_A h_A) \gamma^2 + {b_A}^2} + v_A \gamma - b_A}{b_A \gamma} \\
\end{aligned}
$$

実装して確認したところ $t_{k,p,2}$ を使うと範囲 $[0, 1]$ の出力が得られました。

実装へのリンクです。ここで紹介した計算方法は `getPeakAlt` に実装されています。

- TODO リンク
