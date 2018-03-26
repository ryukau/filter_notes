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

# 分数階微分で1次元の熱伝導-波動方程式
1次元の熱伝導方程式です。

$$
\frac{\partial u}{\partial t} = c^2 \frac{\partial^2 u}{\partial x^2}
$$

1次元の波動方程式です。

$$
\frac{\partial^2 u}{\partial t^2} = c^2 \frac{\partial^2 u}{\partial x^2}
$$

熱伝導方程式と波動方程式を方程式を分数階微分でつなぎます。

$$
\frac{\partial^{1 + \alpha} u}{\partial t^{1 + \alpha}} = c^2 \frac{\partial^2 u}{\partial x^2}
$$

ここでは、この式を熱伝導-波動方程式 (Heat-Wave Equation) と呼びます。拡散-波動方程式 (Diffusion-Wave Equation) と呼ばれることもあります。

## 連立方程式を立てる
Implcit FDM (Finite Difference Method) とGrünwald–Letnikovの分数階微分を組み合わせてシミュレーションします。

Grünwald–Letnikovの分数階微分を再掲します。 $m$ の初期値を明示するために微分演算子 $D$ の表記を少し変えています。

$$
D^{\alpha, k} f(x)
= \lim_{h \to 0} \frac{1}{h^{\alpha}} \sum_{m=k}^{\infty} (-1)^{m}
\frac{\Gamma(\alpha + 1)}{\Gamma(m + 1)\Gamma(\alpha - m + 1)} f(x - mh)
$$

熱伝導-波動方程式を計算できる形にするために、左辺をGrünwald–Letnikovの分数階微分、右辺を有限差分として展開します。今回は $u(\_, t)$ が今から計算する値で $u(\_, t - dt)$ までの値が既に得られているものとしています。

$$
\begin{aligned}
u(x, t) +& \sum_{m=1}^{\infty} (-1)^{m} \binom{1 + \alpha}{m} u(x, t - m\,dt)\\
=&
\frac{c^2 dt^{1 + \alpha}}{dx^2} \left( u(x - dx, t) -2u(x, t) + u(x + dx, t) \right)\\
\end{aligned}
$$

Implicit FDMの形に整理します。

$$
\begin{aligned}
C_1 u(x - dx, t) + C_2 u(x, t) + C_1 u(x + dx, t)
=& \sum_{m=1}^{\infty} (-1)^{m} \binom{1 + \alpha}{m} u(x, t - m\,dt)\\
C_1 = \frac{c^2 dt^{1 + \alpha}}{dx^2},\quad C_2 = -(1 + 2A)
\end{aligned}
$$

連立方程式を立てます。

$$
\mathbf{A} \mathbf{u}^{t} = D^{(1 + \alpha), 1} \mathbf{u}^{t}
$$

$\mathbf{A}$ と $\mathbf{u}^t$ です。

$$
\mathbf{A} =
\begin{bmatrix}
 C_2 & l & 0 & 0 & \cdots & 0 & 0\\
 C_1 & C_2 & C_1 & 0 & \cdots & 0 & 0\\
 0 & C_1 & C_2 & C_1 & \cdots & 0 & 0\\
& \vdots & & & & \vdots &\\
 0 & 0 & 0 & 0 & \cdots & r & C_2
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
l =& C_1 \left (1 + L \right )\\
r =& C_1 \left (1 + R \right )
\end{aligned}
$$

$L$ と $R$ の値は境界条件によって決まります。固定端のときは $L = R = 0$ 、自由端のときは $L = R = 1$ となります。

立てた連立方程式を左辺の $\mathbf{u}^{t}$ について解くことでシミュレーションを1ステップ進めることができます。

## 実装
コードを実行すると1次元の熱-波のシミュレーションを行い、結果を `heat_wave1d.png` に画像として書き出します。

連立方程式のソルバに[numpy.linalg.solve](https://docs.scipy.org/doc/numpy/reference/generated/numpy.linalg.solve.html)、画像の書き出しに[Imageio](http://imageio.github.io/)を使っています。

```python
import numpy
import imageio
from scipy.special import binom


class HeatWave1D():
    def __init__(self, length, c, dx, dt, alpha, attenuation):
        """
        1次元の熱と波のシミュレータ。

        d^(1 + alpha) u / dt^(1 + alpha) = c^2 d^2 u / dx^2

        :param length: 波を表す1次元配列の長さ。
        :param c: 波の速度[m/s]。
        :param dx: 配列の要素間の距離[m]。
        :param dt: シミュレーションの1ステップで進む時間[s]。
        :param alpha: 分数階微分の階数。 [0, 1] の範囲。
        :param attenuation: 厳密でない波の減衰係数。 [0, 1] の範囲。
        """
        # u(x, t) -> self.wave[t][x]
        self.length = length
        self.fracDiffDepth = 256
        self.wave = numpy.zeros((self.fracDiffDepth, self.length))

        self.field = numpy.zeros(self.length)

        self.setParameters(c, dx, dt, alpha, attenuation)
        self.pick(0, 0)

    def value(self):
        """
        描画用に最新の波を返す。
        """
        return self.wave[0]

    def initBoundary(self):
        """
        境界条件を指定。 0 で固定端。 1 で自由端。
        """
        self.L = 0
        self.R = 0

    def setParameters(self, c, dx, dt, alpha, attenuation):
        """
        パラメータを更新する。
        シミュレーションの途中でパラメータが変更された時に呼び出す。
        """
        self.c = c
        self.dx = dx
        self.dt = dt
        self.alpha = alpha
        self.attenuation = attenuation

        self.refreshConstants()
        self.initBoundary()
        self.initMatrix()

    def refreshConstants(self):
        """
        シミュレーションで用いる定数を設定。
        """
        self.C1 = (self.c / self.dx)**2 * self.dt**(1 + self.alpha)
        self.C2 = -(1 + 2 * self.C1)

        self.fracCoefficients = [(-1)**m * binom(1 + self.alpha, m)
                                 for m in range(self.fracDiffDepth)]

    def initMatrix(self):
        """
        Implicit finite difference method で解く必要のある方程式の設定。
        ax = b の a。
        """
        mat = numpy.zeros((self.length, self.length))

        mat[0][0] = self.C2
        mat[0][1] = self.C1 * (1 + self.L)

        last = self.length - 1
        for i in range(1, last):
            mat[i][i - 1] = self.C1
            mat[i][i] = self.C2
            mat[i][i + 1] = self.C1

        mat[last][last - 1] = self.C1 * (1 + self.R)
        mat[last][last] = self.C2

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
            self.fractionalDifference(),
        )

    def fractionalDifference(self):
        """
        分数階微分の計算。方程式 ax = b の bを返す。
        """
        self.field.fill(0)
        for m in range(1, len(self.wave)):
            self.field += self.fracCoefficients[m] * self.wave[m]
        return self.field

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
    wave1d = HeatWave1D(length, 64, 0.1, 0.01, 0.5, 1)

    result = []
    wave1d.pick(0.5, 1)
    for t in range(0, length):
        wave1d.step()
        result.append(wave1d.value())

    imageio.imwrite("heat_wave1d.png", numpy.array(result))
```

実行結果の画像です。横はあるステップでの波の状態、縦は時間で上から下に向かって進んでいます。

<figure>
<img src="img/heat_wave_equation/heat_wave1d.png" alt="Image of a result of 1d heat-wave simulation." style="width: 512px; padding-bottom: 12px;"/>
</figure>

## デモ
$\alpha$ の値を変えて熱伝導-波動をシミュレーションした動画です。

<video width="640px" controls>
  <source src="img/heat_wave_equation/heat_wave1d.mp4" type="video/mp4">
  Video of 1 dimensional heat-wave simulation with variety of alpha.
</video>

デモのコードは次のリンクから読むこともできます。

<a href="https://github.com/ryukau/filter_notes/tree/master/docs/demo/heat_wave_equation">デモのコードを見る (github.com)</a>

## 問題点
自由端にすると発散します。固定端でも $0 < \alpha < 1$ かつ $C1$ が小さいときに発散することがあります。この発散がExplicit FDMの発散と関係がないことを確認するためにImplicit FDMで実装しました。
