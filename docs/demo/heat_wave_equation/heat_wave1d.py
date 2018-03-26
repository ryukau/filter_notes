import numpy
import imageio
from scipy.special import binom


class HeatWave1D():
    def __init__(self, length, c, dx, dt, alpha, attenuation):
        """
        1次元の熱と波のシミュレータ。

        d^(1 + alpha) u / dt^(1 + alpha) = d^2 u / dx^2

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
        self.R = 1

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
