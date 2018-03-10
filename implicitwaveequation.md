<style>
body {
  max-width: 704px;
  margin: auto;
  padding: 32px 8px;
}

code {
  overflow-x: scroll;
  overflow-y: hidden;
  white-space: pre;
}

.katex {
  font-size: 1.3em !important;
}
</style>

# Implicit FDM で1次元の波のシミュレーション
<a href="waveequation.html">前回のシミュレーション</a> では [Explicit FDM](https://en.wikipedia.org/wiki/Finite_difference_method#Explicit_method) を用いていたので、パラメータの値によっては発散するという問題がありました。 [Implicit FDM](https://en.wikipedia.org/wiki/Finite_difference_method#Implicit_method) を使えば計算コストと引き換えに発散しなくなります。


## Explicit FDM と ImplicitFDM
Explicit な有限差分と Implicit な有限差分では空間方向の微分の形が変わります。

例として、1次元の波動方程式を有限差分の形にします。

$$
\frac{\partial^2 u(x,\,t)}{\partial t^2} = c^2 \frac{\partial^2 u(x,\,t)}{\partial x^2}
$$

Explicit な有限差分です。

$$
\frac{u{\left (x,t + dt \right )} - 2\,u{\left (x,t \right )} + u{\left (x,t - dt \right )}}{dt^{2}}
=
c^2\,\frac{u{\left (- dx + x,t \right )} - 2\,u{\left (x,t \right )} + u{\left (dx + x,t \right )}}{dx^{2}}
$$

Implicit な有限差分です。

$$
\frac{u{\left (x,t + dt \right )} - 2\,u{\left (x,t \right )} + u{\left (x,t - dt \right )}}{dt^{2}}
=
c^2\,\frac{u{\left (- dx + x,t + dt \right )} - 2\,u{\left (x,t + dt \right )} + u{\left (dx + x,t + dt \right )}}{dx^{2}}
$$

右辺の $u(\,\_\,,t)$ が $u(\,\_\,, t + dt)$ に変わっています。 $u(\,\_\,, t + dt)$ は今から計算する未来のステップの値です。

Explicit な有限差分では $u(\,\_\,, t + dt)$ の形になる項は $u(x,t + dt)$ だけです。したがって出てきた有限差分の式を $u(x,t + dt)$ について解けば計算できる形になります。

Implicit な有限差分では $u(\,\_\,, t + dt)$ の形になる項が空間方向の微分によって複数出てきます。1つの式を整理しても計算できる形にはなりませんが、すべての $x$ について式を立てれば連立方程式として解くことができます。

## 連立方程式を立てる
1次元の波動方程式をImplicit FDMで解く連立方程式の形にします。

1次元の波動方程式を再掲します。

$$
\frac{\partial^2 u(x,\,t)}{\partial t^2} = c^2 \frac{\partial^2 u(x,\,t)}{\partial x^2}
$$

SymPyを使ってImplicitな有限差分に変形します。

```python
from sympy import *

c, t, u, x, dx, dt = symbols('c t u x dx dt')

u_xt1 = u(x, t + dt)
u_dx2 = u_xt1.diff(x, x).as_finite_difference(dx)

u_xt = u(x, t)
u_dt2 = u_xt.diff(t, t).as_finite_difference(dt)

eq = Eq(u_dt2, c**2 * u_dx2)
ans = solve(eq, u(x, t + dt))

latex(expand(ans[0]))
```

出力された式です。

$$
\begin{aligned}
u{\left (x,t \right )}
=& \frac{c^{2} dt^{2}}{dx^{2}} u{\left (x,dt + t \right )}
- \frac{c^{2} dt^{2}}{2 dx^{2}} u{\left (- dx + x,dt + t \right )}
- \frac{c^{2} dt^{2}}{2 dx^{2}} u{\left (dx + x,dt + t \right )}\\
&+ \frac{1}{2} u{\left (x,- dt + t \right )}
+ \frac{1}{2} u{\left (x,dt + t \right )}
\end{aligned}
$$

整理します。

$$
\begin{gathered}
\alpha u{\left (x - dx,t + dt \right )}
+ \beta u{\left (x,t + dt \right )}
+ \alpha u{\left (x + dx,t + dt \right )}
= \space
u{\left (x,t \right )}
- \frac{1}{2} u{\left (x,t - dt \right )}\\

\alpha = \space -\frac{c^{2} dt^{2}}{2 \; dx^{2}}, \quad
\beta = \frac{c^{2} dt^{2}}{dx^{2}} + \space \frac{1}{2}
\end{gathered}
$$

連立方程式をたてます。

$$
\mathbf{A} \mathbf{u}^{t+1} = \mathbf{u}^t - 0.5 \mathbf{u}^{t-1}
$$

$\mathbf{A}$ と $\mathbf{u}^t$ です。

$$
\mathbf{A} =
\begin{bmatrix}
 \beta & l & 0 & 0 & \cdots & 0 & 0\\
 \alpha & \beta & \alpha & 0 & \cdots & 0 & 0\\
 0 & \alpha & \beta & \alpha & \cdots & 0 & 0\\
& \vdots & & & & \vdots &\\
 0 & 0 & 0 & 0 & \cdots & r & \beta
\end{bmatrix}
,\quad
\mathbf{u}^t =
\begin{bmatrix}
u(x_0, t)\\
u(x_1, t)\\
u(x_2, t)\\
\vdots\\
u(x_{n-1}, t)
\end{bmatrix}
$$

$\mathbf{A}$ に含まれる $l$ と $r$ です。

$$
\begin{aligned}
l =& \alpha \left (1 + L \right )\\
r =& \alpha \left (1 + R \right )
\end{aligned}
$$

$L$ と $R$ の値は境界条件によって決まります。固定端のときは $L = R = 0$ 、自由端のときは $L = R = 1$ となります。

立てた連立方程式を $\mathbf{u}^{t+1}$ について解くことでシミュレーションを1ステップ進めることができます。

## 実装
コードを実行すると1次元の波のシミュレーションを行い、結果を `wave1d.png` に画像として書き出します。

連立方程式のソルバに[numpy.linalg.solve](https://docs.scipy.org/doc/numpy/reference/generated/numpy.linalg.solve.html)、画像の書き出しに[Imageio](http://imageio.github.io/)を使っています。

```python
import numpy
import imageio


class Wave1D():
    def __init__(self, length, c, dx, dt, attenuation):
        """
        1次元の波のシミュレータ。

        :param length: 波を表す1次元配列の長さ。
        :param c: 波の速度[m/s]。
        :param dx: 配列の要素間の距離[m]。
        :param dt: シミュレーションの1ステップで進む時間[s]。
        :param attenuation: 厳密でない波の減衰係数。 [0, 1] の範囲。
        """
        # u(x, t) -> self.wave[t][x]
        self.length = length
        self.wave = numpy.zeros((3, self.length))

        self.c = c
        self.dx = dx
        self.dt = dt
        self.attenuation = attenuation

        self.refreshConstants()
        self.initBoundary()
        self.initMatrix()
        self.pick(0, 0)

    def value(self):
        """
        描画用に最新の波を返す。
        """
        return self.wave[0]

    def refreshConstants(self):
        """
        シミュレーションで用いる定数を設定。
        """
        c2_dt2_dx2 = (self.c * self.dt / self.dx)**2
        self.alpha = -c2_dt2_dx2 / 2
        self.beta = 0.5 + c2_dt2_dx2

    def initBoundary(self):
        """
        境界条件を指定。 0 で固定端。 1 で自由端。
        """
        self.L = 0
        self.R = 1

    def initMatrix(self):
        """
        Implicit finite difference method で解く必要のある方程式の設定。
        ax = b の a。
        """
        mat = numpy.zeros((self.length, self.length))

        mat[0][0] = self.beta
        mat[0][1] = self.alpha * (1 + self.L)

        last = self.length - 1
        for i in range(1, last):
            mat[i][i - 1] = self.alpha
            mat[i][i] = self.beta
            mat[i][i + 1] = self.alpha

        mat[last][last - 1] = self.alpha * (1 + self.R)
        mat[last][last] = self.beta

        self.a = mat

    def step(self):
        """
        シミュレーションの1ステップ。
        """
        self.wave = numpy.roll(self.wave, 1, axis=0)
        if self.pickY != 0:
            self.wave[1][self.pickX] = self.pickY
        self.wave[0] = self.attenuation * numpy.linalg.solve(
            self.a,
            self.wave[1] - 0.5 * self.wave[2],
        )

    def reset(self):
        """
        波を 0 で埋めて初期状態に戻す。
        """
        self.wave.fill(0)

    def pick(self, x, y):
        """
        x, y で指定した位置の波をつまむ。

        :param x: 波をつまむ場所。[0, 1] の範囲。
        :param y: 波をつまむ高さ。任意の実数。
        """
        self.pickX = int((self.length - 1) * numpy.clip(x, 0.0, 1.0))
        self.pickY = y


if __name__ == "__main__":
    length = 512
    wave1d = Wave1D(length, 64, 0.1, 0.01, 1)

    result = []
    wave1d.pick(0.5, 1)
    for t in range(0, length):
        wave1d.step()
        result.append(wave1d.value())

    imageio.imwrite("wave1d.png", numpy.array(result))
```

実行結果の画像です。横はあるステップでの波の状態、縦は時間で上から下に向かって進んでいます。

<figure>
<img src="img/waveequationimplicit/wave1d.png" alt="Image of a result of 1d wave simulation." style="width: 512px; padding-bottom: 12px;"/>
</figure>

## デモ
Implicit FDMとExplicit FDMによる波のシミュレーションを比較した動画です。

<video controls>
  <source src="img/waveequationimplicit/wave1d_implicit_vs_explicit.mp4" type="video/mp4">
  Video of 1 dimensional wave simulation that comparing implicit method to explicit method.
</video>

デモのコードは次のリンクの `draw.py` と `wave.py` になります。

<a href="https://github.com/ryukau/filter_notes/tree/master/docs/demo/waveequation_implicit">デモのコードを見る (github.com)</a>

## その他
今回実装したシミュレーションでは減衰係数が1のときでも時間の経過とともに波がなまっていきます。
