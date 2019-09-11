# PTR の式を求める Maxima のコード
- $f_0$ : 基本周波数。
- $f_s$ : サンプリング周波数。
- $T = f_0/f_s$
- $\phi = n T \bmod 1$
- $n = \phi / T$

数式上は $\phi$ と $n$ が循環参照になっていますが、実装では $\phi$ から $n$ を計算できます。

```python
class PTROscillator:
    def __init__(self, samplerate, frequency):
        self.phi = 0.0
        self.T = frequency / samplerate
        self.h = 1.0

    def process(self):
        self.phi += self.T
        if self.phi >= 1.0:
            self.phi -= 1.0
        return self.PTRSaw3(self.phi, self.T, self.h)

    def PTRSaw3(self, phi, T, h):
        n = phi / T
        if n >= 3.0:
            return 2*T*n-3*T-1
        if 0.0 <= n and n < 1.0:
            return -(h*n**3)/3+2*T*n+2*h-3*T-1
        if 1.0 <= n and n < 2.0:
            return (2*h*n**3)/3-3*h*n**2+3*h*n+2*T*n+h-3*T-1
        if 2.0 <= n and n < 3.0:
            return -(h*n**3)/3+3*h*n**2-9*h*n+2*T*n+9*h-3*T-1
        return 0  # 念のため。
```

以下のコードを wxMaxima にコピペすれば式が得られます。

## Saw
```maxima
DPWSaw(N) := block([result: [], eq: x = 0],
  integration_constant: 'C,
  reset(integration_constant_counter),
  for i: 1 thru N do block([a, c],
    eq: integrate(eq, x),
    a: rhs(solve(eq, concat('C, i))[1]),
    if i > 1 then (
      c: rhs(solve(subst(1, x, a) = subst(-1, x, a), concat('C, i - 1))[1]),
      a: subst(c, concat('C, i - 1), a),
      eq: subst(c, concat('C, i - 1), eq)
    ),
    a: num(a),
    a: a / coeff(a, x, hipow(a, x)),
    result: endcons(expand(a), result)
  ),
  push(x, result)
);

PTRSaw(dpw, N) := block([ptr],
  removeTinverse(eq) := block([r: 0],
    if atom(eq) then r: eq
    else for term in args(eq) do block([den],
      den: denom(term),
      if atom(den) then
        if den # T then r: r + term
    ),
    r
  ),

  phi(n, boundary) :=
    if boundary = 'n or n > boundary then n * T
    else n * T + h,

  s(n, boundary) := 2 * phi(n, boundary) - 1,

  y(dpw, N) := block([M: N - 1, result: []],
    for b: 0 thru M do block([eq: 0],
      for k: 0 thru M do
        eq: eq + (-1)^k * binomial(M, k) * subst(s(n - k, n - b), x, dpw[N]),
      result: endcons(eq, result)
    ),
    result / (2 * T)^(N - 1) / N!
  ),

  ptr: expand(ratsimp(y(dpw, N))),
  for i: 1 thru length(ptr) do ptr[i]: removeTinverse(ptr[i]),
  ptr
);

dispptrsaw(dpwsaw) := block([N: length(dpwsaw)],
  for i: 1 thru N do block([eq],
    eq: PTRSaw(dpwsaw, i),

    disp(concat('PTRSaw, i - 1)),
    disp([(length(eq) - 1) * T <= phi, eq[1]]),
    for j: 2 thru length(eq) do
      disp([(j - 2) * T <= phi and phi < (j - 1) * T, eq[j]])
  )
);

dpwsaw: DPWSaw(10);
dispptrsaw(dpwsaw);

/* format.py で利用。 */
/*
dispptrlist(dpw, ptrfunc) := block([L: []],
  for i: 1 thru length(dpw) do L: endcons(ptrfunc(dpw, i), L),
  L
);
dispptrlist(dpwsaw, PTRSaw);
*/
```

## Triangle
位相が $[0, 1)$ で一周するとき、 $[0, 0.5)$ の範囲の計算式だけを求めています。

導出中の遷移領域の位置 `k` と、波形の不連続点 `b` の関係に応じて絶対値 `abs(x)` を `x` か `-x` に置き換えています。位相が $[0, 1)$ で一周するとき、 $[0, 0.5)$ の範囲の計算式だけを求めています。

Saw とは $\phi$ と $s$ が異なります。スケーリング係数は Saw の 2 倍です。

```maxima
DPWTri(N) := block([result: []],
  tri(N) := block([poly: 0, poly_diff: [], A: [], D],
    for i: 1 thru N do (
      if oddp(i) then
        poly: poly + if i = N then x^i else concat(a, i) * x^i
      elseif i = N then
        poly: poly + x^(i - 1) * abs(x)
    ),

    D: subst(x, abs(x), poly),
    D: diff(D, x),
    poly_diff: endcons(D, poly_diff),
    while hipow(D, x) > 1 do (
      D: diff(D, x, 2),
      poly_diff: endcons(D, poly_diff)
    ),

    for i: length(poly_diff) step -1 thru 1 do (
      D: solve(subst(1 / 2, x, poly_diff[i]) = 0, concat(a, i * 2 - 1)),
      A: append(A, subst(A, D))
    ),

    subst(A, poly)
  ),

  for i: 1 thru N do result: endcons(tri(i), result),
  result
);

PTRTri(dpw, N) := block([ptr: []],
  phi(n, boundary) :=
    if boundary = 'n or n > boundary then n * T
    else n * T + 1,

  s_odd(n, boundary) := 1 / 2 - abs(1 - 2 * phi(n, boundary)),
  s_even(n, boundary) := 1 - 2 * phi(n, boundary),

  replaceabs(eq, sign) := block(
    f(x) := if not atom(x) and op(x) = abs then sign * part(x, 1) else x,
    if atom(eq) then eq else scanmap(f, eq)
  ),

  y(dpw, N, s) := block([M: N - 1, result: []],
    for b: 0 thru M do block([eq: 0, boundary],
      boundary: 'n - b,
      for k: 0 thru M do block([abs_x],
        abs_x: if b = 0 then x elseif k < b then x else -x,
        eq: eq + (-1)^k * binomial(M, k) * replaceabs(
          subst(s(n - k, boundary), x, dpw[N]),
          if b = 0 then -1 elseif k < b then -1 else 1
        )
      ),
      result: endcons(eq, result)
    ),
    result * 2 / (2 * T)^(N - 1) / N!
  ),

  ptr: y(dpw, N, if oddp(N) then s_odd else s_even),
  ptr: expand(ratsimp(ptr)),
  ptr
);

dispptrtri(dpwtri) := block([N: length(dpwtri)],
  for i: 1 thru N do block([eq],
    eq: PTRTri(dpwtri, i),
    disp(concat('PTRTri, i - 1)),
    for j: 1 thru length(eq) do disp([eq[j]])
  )
);

dispptrlist(dpw) := block([L: []],
  for i: 1 thru length(dpw) do L: endcons(PTRTri(dpw, i), L),
  L
);

dpwtri: DPWTri(11);
dispptrtri(dpwtri);
```

## Ramp
- [Ramp function - Wikipedia](https://en.wikipedia.org/wiki/Ramp_function)

Saw を元にして境界の向こう側を 0 に置き直しました。

```maxima
DPWRamp(N) := block([result: [], eq: x + 1 = 0],
  integration_constant: 'C,
  reset(integration_constant_counter),
  for i: 1 thru N do block([a, c],
    eq: integrate(eq, x),
    a: rhs(solve(eq, concat('C, i))[1]),
    if i > 1 then (
      /* solve の式を変更。 */
      c: rhs(solve(subst(-1, x, a) = 0, concat('C, i - 1))[1]),
      a: subst(c, concat('C, i - 1), a),
      eq: subst(c, concat('C, i - 1), eq)
    ),
    a: num(a),
    a: a / coeff(a, x, hipow(a, x)),
    result: endcons(expand(a), result)
  ),
  push(x, result)
);

PTRRamp(dpw, N) := block([ptr],
  removeTinverse(eq) := block([r: 0],
    if atom(eq) then r: eq
    else for term in args(eq) do block([den],
      den: denom(term),
      if atom(den) then
        if den # T then r: r + term
    ),
    r
  ),

  phi(n, boundary) :=
    if boundary = 'n or n > boundary then n * T
    else 0, /* 境界を超えると 0 。 */

  s(n, boundary) := 2 * phi(n, boundary) - 1,

  y(dpw, N) := block([M: N - 1, result: []],
    for b: 0 thru M do block([eq: 0],
      for k: 0 thru M do
        eq: eq + (-1)^k * binomial(M, k) * subst(s(n - k, n - b), x, dpw[N]),
      result: endcons(eq, result)
    ),
    result / (2 * T)^(N - 1) / N!
  ),

  ptr: expand(ratsimp(y(dpw, N))),
  for i: 1 thru length(ptr) do ptr[i]: removeTinverse(ptr[i]),
  ptr
);

dispptrramp(dpwramp) := block([N: length(dpwramp)],
  for i: 1 thru N do block([eq],
    eq: PTRRamp(dpwramp, i),

    disp(concat('PTRRamp, i - 1)),
    disp([(length(eq) - 1) * T <= phi, eq[1]]),
    for j: 2 thru length(eq) do
      disp([(j - 2) * T <= phi and phi < (j - 1) * T, eq[j]])
  )
);

dispptrlist(dpw) := block([L: []],
  for i: 1 thru length(dpw) do L: endcons(PTRRamp(dpw, i), L),
  L
);

dpwramp: DPWRamp(10);
dispptrramp(dpwramp);
```

## Step
Ramp 関数の微分です。

- [Heaviside step function - Wikipedia](https://en.wikipedia.org/wiki/Heaviside_step_function)

```maxima
DPWStep(N) := block([result: [], eq: x + 1 = 0],
  integration_constant: 'C,
  reset(integration_constant_counter),
  for i: 1 thru N do block([a, c],
    eq: integrate(eq, x),
    a: rhs(solve(eq, concat('C, i))[1]),
    if i > 1 then (
      c: rhs(solve(subst(-1, x, a) = 0, concat('C, i - 1))[1]),
      a: subst(c, concat('C, i - 1), a),
      eq: subst(c, concat('C, i - 1), eq)
    ),
    a: num(a),
    a: a / coeff(a, x, hipow(a, x)),
    result: endcons(diff(expand(a), x), result) /* diff を追加。 */
  ),
  push(x, result)
);

PTRStep(dpw, N) := block([ptr],
  removeTinverse(eq) := block([r: 0],
    if atom(eq) then r: eq
    else for term in args(eq) do block([den],
      den: denom(term),
      if atom(den) then
        if den # T then r: r + term
    ),
    r
  ),

  phi(n, boundary) :=
    if boundary = 'n or n > boundary then n * T
    else h, /* 0 でも h でも結果が同じになる。 */

  s(n, boundary) := 2 * phi(n, boundary) - 1,

  y(dpw, N) := block([M: N - 1, result: []],
    for b: 0 thru M do block([eq: 0],
      for k: 0 thru M do
        eq: eq + (-1)^k * binomial(M, k) * subst(s(n - k, n - b), x, dpw[N]),
      result: endcons(eq, result)
    ),
    result / (2 * T)^(N - 1) / N!
  ),

  ptr: expand(ratsimp(y(dpw, N))),
  for i: 1 thru length(ptr) do ptr[i]: removeTinverse(ptr[i]),
  ptr
);

dpwstep: DPWStep(10);
dispptrramp(dpwstep);
```
