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

# 1次元の分数階Zener波動方程式
HolmとNäsholmによる ["A causal and fractional all-frequency wave equation for lossy media"](http://heim.ifi.uio.no/~sverre/papers/2011_HolmNasholm-fractZener-JournAcoustSocAm.pdf) の式(10)で紹介されていた、分数階Zener波動方程式 (fractional Zener model wave equation) を1次元の implicit な形で実装します。

分数階Zener波動方程式です。

$$
\nabla^2 u
- \frac{1}{c_0^2} \frac{\partial^2 u}{\partial t^2}
+ \tau_{\sigma}^{\alpha} \frac{\partial^{\alpha}}{\partial t^{\alpha}} \nabla^2 u
- \frac{\tau_{\varepsilon}^{\beta}}{c_0^2} \frac{\partial^{\beta + 2} u}{\partial t^{\beta + 2}}
= 0
$$

## 連立方程式を立てる
式の簡略化のために演算子 $\hat{D}$ を定義します。

$$
\hat{D}^{\alpha, k} f(x) = \sum_{m=k}^{\infty} (-1)^{m} \binom{\alpha}{m} f(x - mh)
$$

自然数階微分を有限差分に変形します。

$$
\begin{aligned}
\nabla^2 u =&
\frac{1}{dx^{2}} \left(
    u{\left (x - dx,t \right )}
    - 2 u{\left (x,t \right )}
    + u{\left (x + dx,t \right )}
\right)\\

\frac{\partial^2 u}{\partial t^2} =&
\frac{1}{dt^{2}} \left(
    u{\left (x,t \right )}
    - 2 u{\left (x,t - dt \right )}
    + u{\left (x,t - 2dt \right )}
\right)
\end{aligned}
$$

分数階Zener波動方程式に代入します。

$$
\begin{aligned}
0 =
& \frac{1}{dx^{2}} \left(
    u{\left (x - dx,t \right )}
    - 2 u{\left (x,t \right )}
    + u{\left (x + dx,t \right )}
\right)\\
&- \frac{1}{c_0^2 dt^{2}} \left(
    u{\left (x,t \right )}
    - 2 u{\left (x,t - dt \right )}
    + u{\left (x,t - 2dt \right )}
\right)\\
&+ \frac{\tau_{\sigma}^{\alpha}}{dt^{\alpha}} \hat{D}^{\alpha} \left(
    \frac{1}{dx^{2}} \left(
        u{\left (x - dx,t \right )}
        - 2 u{\left (x,t \right )}
        + u{\left (x + dx,t \right )}
    \right)
\right)\\
&- \frac{\tau_{\varepsilon}^{\beta}}{c_0^2 dt^{2 + \beta}} \hat{D}^{2 + \beta} \left(
    u(x, t)
\right)\\
\end{aligned}
$$

移項します。

$$
\begin{aligned}
& \frac{1}{dx^{2}} \left(1 + \frac{\tau_{\sigma}^{\alpha}}{dt^{\alpha}} \right)
\left(
    u{\left (x - dx,t \right )}
    - 2 u{\left (x,t \right )}
    + u{\left (x + dx,t \right )}
\right)\\
&- \frac{1}{c_0^2 dt^{2}} u{\left (x,t \right )}\\
&- \frac{\tau_{\varepsilon}^{\beta}}{c_0^2 dt^{2 + \beta}} u(x, t)\\

=
& \frac{1}{c_0^2 dt^{2}} \left(
    - 2 u{\left (x,t - dt \right )}
    + u{\left (x,t - 2dt \right )}
\right)\\
&- \frac{\tau_{\sigma}^{\alpha}}{dt^{\alpha} dx^{2}} \hat{D}^{\alpha, 1} \left(
    u{\left (x - dx,t \right )}
    - 2 u{\left (x,t \right )}
    + u{\left (x + dx,t \right )}
\right)\\
&+ \frac{\tau_{\varepsilon}^{\beta}}{c_0^2 dt^{2 + \beta}} \hat{D}^{2 + \beta, 1} \left(
    u(x, t)
\right)\\
\end{aligned}
$$

計算しやすい形に整理します。

$$
\begin{aligned}
C_1 u{\left (x - dx,t \right )}
+& C_4 u{\left (x,t \right )}
+ C_1 u{\left (x + dx,t \right )}\\
=&
C_2 \left( u{\left (x,t - 2dt \right )} - 2 u{\left (x,t - dt \right )} \right)\\
&- C_0 \hat{D}^{\alpha, 1} \left(
    u{\left (x - dx,t \right )}
    - 2 u{\left (x,t \right )}
    + u{\left (x + dx,t \right )}
\right)\\
&+C_3 \hat{D}^{2 + \beta, 1} \left( u(x, t) \right)\\
\\

C_0 =& \frac{\tau_{\sigma}^{\alpha}}{dt^{\alpha} dx^2}
,\quad
C_1 = \frac{1}{dx^{2}} + C_0
,\quad
C_2 = \frac{1}{c_0^2 dt^{2}}
\\
C_3 =& C_2 \frac{\tau_{\varepsilon}^{\beta}}{dt^{\beta}}
,\quad
C_4 = -(2 C_1 + C_2 + C_3)
\end{aligned}
$$

連立方程式を立てます。

$$
\mathbf{A} \mathbf{u}^{t} =
C_2 (\mathbf{u}^{t-2} - 2 \mathbf{u}^{t-1})
- C_0 \hat{D}^{\alpha, 1} (\mathbf{u}_{x-1}^{t} - 2 \mathbf{u}_{x}^{t} + \mathbf{u}_{x + 1}^{t})
+ C_3 \hat{D}^{2 + \beta, 1} (\mathbf{u}^{t})
$$

$\mathbf{A}$ と $\mathbf{u}_{x + i}^t$ です。

$$
\mathbf{A} =
\begin{bmatrix}
 C_4 & C_1 & 0 & 0 & \cdots & 0 & 0\\
 C_1 & C_4 & C_1 & 0 & \cdots & 0 & 0\\
 0 & C_1 & C_4 & C_1 & \cdots & 0 & 0\\
& \vdots & & & & \vdots &\\
 0 & 0 & 0 & 0 & \cdots & C_1 & C_4
\end{bmatrix}
,\quad
\mathbf{u}_{x + i}^t =
\begin{bmatrix}
u(x_{0 + i}, t)\\
u(x_{1 + i}, t)\\
u(x_{2 + i}, t)\\
\vdots\\
u(x_{n + i -1}, t)
\end{bmatrix}
$$

$k < 0$ あるいは $n \le k$ のとき $u(x_k, t) = 0$ としています。

左右の端は固定端にしています。自由端にすると発散します。

## 実装
連立方程式のソルバに[numpy.linalg.solve](https://docs.scipy.org/doc/numpy/reference/generated/numpy.linalg.solve.html)を使っています。

```python
import numpy
from scipy.special import binom


class ZenerWave1D():
    def __init__(self, length, c, dx, dt, attenuation, tau_epsilon, tau_sigma,
                 alpha, beta):
        """
        Fractional Zener wave equation のシミュレータ。

        :param length: 波を表す1次元配列の長さ。
        :param c: 波の速度[m/s]。
        :param dx: 配列の要素間の距離[m]。
        :param dt: シミュレーションの1ステップで進む時間[s]。
        :param attenuation: 厳密でない波の減衰係数。 [0, 1] の範囲。
        :param tau_epsilon: Relaxation time.
        :param tau_sigma: Retardation time.
        :param alpha: 分数階微分の階数。 [0, 2] の範囲。
        :param beta: 分数階微分の階数。 [0, 2] の範囲。
        """
        self.length = length
        self.fracDiffDepth = 64
        self.wave = numpy.zeros((self.fracDiffDepth, self.length))

        self.alphaField = numpy.zeros(self.length)
        self.betaField = numpy.zeros(self.length)

        self.setParameters(c, dx, dt, attenuation, tau_epsilon, tau_sigma,
                           alpha, beta)
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
        self.L = 0.0
        self.R = 0.0

    def setParameters(self, c, dx, dt, attenuation, tau_epsilon, tau_sigma,
                      alpha, beta):
        """
        パラメータを更新する。
        シミュレーションの途中でパラメータが変更された時に呼び出す。
        """
        self.c = c
        self.dx = dx
        self.dt = dt
        self.attenuation = attenuation
        self.tau_epsilon = tau_epsilon
        self.tau_sigma = tau_sigma
        self.alpha = alpha
        self.beta = beta

        self.refreshConstants()
        self.initBoundary()
        self.initMatrix()

    def refreshConstants(self):
        """
        シミュレーションで用いる定数を設定。
        """
        dx2 = self.dx**2
        self.C0 = (self.tau_sigma / self.dt)**self.alpha / dx2
        self.C1 = 1 / dx2 + self.C0
        self.C2 = 1 / (self.c * self.dt)**2
        self.C3 = self.C2 * (self.tau_epsilon / self.dt)**self.beta
        self.C4 = -(2 * self.C1 + self.C2 + self.C3)

        self.alphaC = [(-1)**m * binom(self.alpha, m)
                       for m in range(self.fracDiffDepth)]
        self.betaC = [(-1)**m * binom(2 + self.beta, m)
                      for m in range(self.fracDiffDepth)]

    def initMatrix(self):
        """
        Implicit finite difference method で解く必要のある方程式の設定。
        ax = b の a。
        """
        mat = numpy.zeros((self.length, self.length))

        mat[0][0] = self.C4
        mat[0][1] = self.C1 * (1 + self.L)

        last = self.length - 1
        for i in range(1, last):
            mat[i][i - 1] = self.C1
            mat[i][i] = self.C4
            mat[i][i + 1] = self.C1

        mat[last][last - 1] = self.C1 * (1 + self.R)
        mat[last][last] = self.C4

        self.A = mat

    def step(self):
        """
        シミュレーションの1ステップ。
        """
        self.wave = numpy.roll(self.wave, 1, axis=0)
        if self.pickY != 0:
            self.wave[1][self.pickX] = self.pickY
        self.wave[0] = self.attenuation * numpy.linalg.solve(
            self.A, self.rightHandSide())

    def rightHandSide(self):
        """
        方程式 ax = b の bを返す。
        """
        self.alphaField.fill(0)
        self.betaField.fill(0)
        last = self.length - 1
        for m in range(1, len(self.wave)):
            right = numpy.roll(self.wave[m], 1)
            left = numpy.roll(self.wave[m], -1)
            right[0] = 0
            left[last] = 0
            self.alphaField += self.alphaC[m] * (
                left - 2 * self.wave[m] + right)
            self.betaField += self.betaC[m] * self.wave[m]
        return self.C2 * (
            self.wave[2] - 2 * self.wave[1]
        ) - self.C0 * self.alphaField + self.C3 * self.betaField

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
```

## デモ
$\alpha, \beta, \tau_\varepsilon$ の値を変えて分数階Zener波動方程式をシミュレーションした動画です。その他のパラメータは $c_0 = 16,\: dt = 1 / 60,\: dx = 0.0625,\: \tau_\sigma = 1$ で固定しています。

<video width="640px" controls>
  <source src="img/waveequation_fractional_zener/fractional_zener_wave1d.mp4" type="video/mp4">
  Video of simulation of 1 dimensional fractional zener wave equation.
</video>

デモのコードは次のリンクから読むこともできます。

<a href="https://github.com/ryukau/filter_notes/tree/master/docs/demo/waveequation_fractional_zener">デモのコードを見る (github.com)</a>

## その他
パラメータの組み合わせによっては発散します。参考にした論文では $\alpha = \beta$ かつ $\tau_\varepsilon < \tau_\sigma$ の場合のみを扱っています。

分数階微分の計算が重いです。分数階微分の計算はFIRフィルタの畳み込みと同じ形になっているので、フィルタ設計の技術を使って計算コストを下げることができるかもしれません。
