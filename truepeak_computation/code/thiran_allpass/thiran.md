## Thiran オールパスフィルタ
離散オールパスフィルタの伝達関数です。

$$
H(z) = \frac{
  a_N + a_{(N-1)} z^{-1} + \dots + a_1 z^{-(N-1)} + z^{-N}
}{
  1 + a_1 z^{-1} + \dots + a_{N-1} z^{-(N-1)} + a_N z^{-N}
}
$$

$N$ 次の Thiran オールパスフィルタで使われる $a_1, \dots , a_N$ の値です。

$$
a_k = (-1)^k \binom{N}{k} \prod_{n=0}^{N} \frac{D - N + n}{D - N + k + n}
$$

- https://ccrma.stanford.edu/~jos/pasp/Thiran_Allpass_Interpolators.html

Maxima で $a_k$ の式を導出します。

```maxima
displayThiran(N) := for k: 1 while k <= N do display(
  concat(a, k) = (-1)^k * binomial(N, k)
    * product((D - N + n) / (D - N + k + n), n, 0, N)
);

thiranCoefficient(N) := block([result: []],
  for k: 1 while k <= N do
    result: endcons(
      (-1)^k * binomial(N, k) * product((D - N + n) / (D - N + k + n), n, 0, N),
      result
    ),
  result
);

thiranTransferFunction(N) := block([num: 0, den: 1, a],
  a: thiranCoefficient(N),
  num: num + z^-length(a),
  for k: 1 thru length(a) do (
    den: den + a[k] * z^-k,
    num: num + a[k] * z^-(length(a) - k)
  ),
  num / den
);
```

`thiran(4);` の実行結果です。

```
a1=-(4*(D-4))/(D+1)
a2=(6*(D-4)*(D-3))/((D+1)*(D+2))
a3=-(4*(D-4)*(D-3)*(D-2))/((D+1)*(D+2)*(D+3))
a4=((D-4)*(D-3)*(D-2)*(D-1))/((D+1)*(D+2)*(D+3)*(D+4))
```

[Direct form II](https://ccrma.stanford.edu/~jos/fp/Direct_Form_II.html) で実装します。

```javascript
class ThiranAllpass4 {
  constructor() {
    this.reset()
  }

  reset() {
    this.x1 = 0
    this.x2 = 0
    this.x3 = 0
    this.x4 = 0
  }

  push(input) {
    this.x4 = this.x3
    this.x3 = this.x2
    this.x2 = this.x1
    this.x1 = input
  }

  process(input, fraction) {
    var d = fraction + 4

    var q1 = (d - 4) / (d + 1)
    var q2 = (d - 3) / (d + 2) * q1
    var q3 = (d - 2) / (d + 3) * q2
    var a1 = -4 * q1
    var a2 = 6 * q2
    var a3 = -4 * q3
    var a4 = (d - 1) / (d + 4) * q3

    var x0 = input - a1 * this.x1 - a2 * this.x2 - a3 * this.x3 - a4 * this.x4
    var y0 = a4 * x0 + a3 * this.x1 + a2 * this.x2 + a1 * this.x3 + this.x4

    this.push(x0)
    return y0
  }
}
```

---

$$
1
-3 \frac{(D-3)}{(D+1)} z^{-1}
+3 \frac{(D-3)(D-2)}{(D+1)(D+2)} z^{-2}
-\frac{(D-3)(D-2)(D-1)}{(D+1)(D+2)(D+3)} z^{-3}
$$

---

$$
\begin{aligned}
A &= \frac{(D-3)}{(D+1)}\\
B &= \frac{(D-2)}{(D+2)}\\
C &= \frac{(D-1)}{(D+3)}
\end{aligned}
$$

---

$$
1
-3 A z^{-1}
+3 A B z^{-2}
-  A B C z^{-3}
$$

この形は多分無理。

---

$$
1 - 3 x + 3 x^2 - x^3 \equiv (1 - x)^3
$$
