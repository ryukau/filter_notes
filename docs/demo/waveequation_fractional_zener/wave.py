import math
import numpy
from scipy.special import binom


class Wave1D():
    def __init__(self, length, c, dx, dt, attenuation):
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

    def refreshConstants(self):
        c2_dt2_dx2 = (self.c * self.dt / self.dx)**2
        self.alpha = -c2_dt2_dx2 / 2
        self.beta = 0.5 + c2_dt2_dx2

    def initBoundary(self):
        self.L = 0.0
        self.R = 0.0

    def initMatrix(self):
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

    def value(self):
        return self.wave[0]

    def step(self):
        self.wave = numpy.roll(self.wave, 1, axis=0)
        if self.pickY != 0:
            self.wave[1][self.pickX] = self.pickY
        self.wave[0] = self.attenuation * numpy.linalg.solve(
            self.a,
            self.wave[1] - 0.5 * self.wave[2],
        )

    def reset(self):
        self.wave.fill(0)

    def pick(self, x, y):
        self.pickX = int((self.length - 1) * numpy.clip(x, 0.0, 1.0))
        self.pickY = y

    def setParameters(self, c, dx, dt, attenuation):
        self.c = c
        self.dx = dx
        self.dt = dt
        self.attenuation = attenuation

        self.refreshConstants()
        self.initBoundary()
        self.initMatrix()


class ZenerWave1D():
    def __init__(self, length, c, dx, dt, attenuation, tau_epsilon, tau_sigma,
                 alpha, beta):
        # u(x, t) -> self.wave[t][x]
        self.length = length
        self.fracDiffDepth = 64
        self.wave = numpy.zeros((self.fracDiffDepth, self.length))

        self.alphaField = numpy.zeros(self.length)
        self.betaField = numpy.zeros(self.length)

        self.setParameters(c, dx, dt, attenuation, tau_epsilon, tau_sigma,
                           alpha, beta)
        self.pick(0, 0)

    def refreshConstants(self):
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

    def initBoundary(self):
        self.L = 0.0
        self.R = 0.0

    def initMatrix(self):
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

    def value(self):
        return self.wave[0]

    def step(self):
        self.wave = numpy.roll(self.wave, 1, axis=0)
        if self.pickY != 0:
            self.wave[1][self.pickX] = self.pickY
            # pick_length = 32
            # start = int(self.pickX - pick_length / 2)
            # if start < 1:
            #     start = 1
            # end = int(self.pickX + pick_length / 2)
            # if end > self.length - 2:
            #     end = self.length - 2
            # denom = pick_length - 1
            # for i in range(start, end):
            #     win = (1 - numpy.cos(8 * numpy.pi * (i - start) / denom)) / 2
            #     self.wave[1][i] = self.pickY * win
        self.wave[0] = self.attenuation * numpy.linalg.solve(
            self.A, self.rightHandSide())

    def rightHandSide(self):
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
        self.wave.fill(0)

    def pick(self, x, y):
        self.pickX = int((self.length - 1) * numpy.clip(x, 0.0, 1.0))
        self.pickY = y

    def setParameters(self, c, dx, dt, attenuation, tau_epsilon, tau_sigma,
                      alpha, beta):
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


class ZenerWave1DExplicit(ZenerWave1D):
    def __init__(self, length, c, dx, dt, attenuation, tau_epsilon, tau_sigma,
                 alpha, beta):
        # u(x, t) -> self.wave[t][x]
        self.length = length
        self.fracDiffDepth = 64
        self.wave = numpy.zeros((self.fracDiffDepth, self.length))

        self.alphaField = numpy.zeros(self.length)
        self.betaField = numpy.zeros(self.length)

        self.setParameters(c, dx, dt, attenuation, tau_epsilon, tau_sigma,
                           alpha, beta)
        self.pick(0, 0)

    def refreshConstants(self):
        self.C1 = (self.dt / self.tau_epsilon)**self.beta
        self.C0 = self.C1 * (self.c * self.dt / self.dx)**2
        self.C2 = self.C0 * (self.tau_sigma / self.dt)**self.alpha

        self.alphaC = [(-1)**m * binom(self.alpha, m)
                       for m in range(self.fracDiffDepth)]
        self.betaC = [(-1)**m * binom(2 + self.beta, m)
                      for m in range(self.fracDiffDepth)]

    def step(self):
        self.wave = numpy.roll(self.wave, 1, axis=0)
        if self.pickY != 0:
            self.wave[1][self.pickX] = self.pickY
        self.wave[0] = self.attenuation * self.rightHandSide()

    def rightHandSide(self):
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
        return self.C0 * (left - 2 * self.wave[1] + right) - self.C1 * (
            self.wave[1] - 2 * self.wave[2] + self.wave[3]
        ) + self.C2 * self.alphaField - self.betaField

    def setParameters(self, c, dx, dt, attenuation, tau_epsilon, tau_sigma,
                      alpha, beta):
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


class HeatWave1D(Wave1D):
    def __init__(self, length, c, dx, dt, alpha, attenuation):
        # u(x, t) -> self.wave[t][x]
        self.length = length
        self.fracDiffDepth = 256
        self.wave = numpy.zeros((self.fracDiffDepth, self.length))

        self.field = numpy.zeros(self.length)

        self.setParameters(c, dx, dt, alpha, attenuation)
        self.pick(0, 0)

    def setParameters(self, c, dx, dt, alpha, attenuation):
        self.c = c
        self.dx = dx
        self.dt = dt
        self.alpha = alpha
        self.attenuation = attenuation

        self.refreshConstants()
        self.initBoundary()
        self.initMatrix()

    def refreshConstants(self):
        self.C1 = (self.c / self.dx)**2 * self.dt**(1 + self.alpha)
        self.C2 = -(1 + 2 * self.C1)

        self.fracCoefficients = [(-1)**m * binom(1 + self.alpha, m)
                                 for m in range(self.fracDiffDepth)]

    def initMatrix(self):
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
        self.wave = numpy.roll(self.wave, 1, axis=0)
        if self.pickY != 0:
            self.wave[1][self.pickX] = self.pickY
        self.wave[0] = self.attenuation * numpy.linalg.solve(
            self.a,
            self.fractionalDifference(),
        )

    def step_CN(self):
        self.wave = numpy.roll(self.wave, 1, axis=0)
        if self.pickY != 0:
            self.wave[1][self.pickX] = self.pickY
        right = numpy.roll(self.wave[1], 1)
        right[0] = right[1] if self.L == 1 else 0
        left = numpy.roll(self.wave[1], -1)
        left[len(left) - 1] = left[len(left) - 2] if self.R == 1 else 0
        self.wave[0] = self.attenuation * numpy.linalg.solve(
            self.a,
            self.fractionalDifference() - self.C1 *
            (left - 2 * self.wave[1] + right),
        )

    def fractionalDifference(self):
        self.field.fill(0)
        for m in range(1, len(self.wave)):
            self.field += self.fracCoefficients[m] * self.wave[m]
        return self.field

    def reset(self):
        self.wave.fill(0)

    def pick(self, x, y):
        self.pickX = int((self.length - 1) * numpy.clip(x, 0.0, 1.0))
        self.pickY = y


class Wave1DExplicit(Wave1D):
    def __init__(self, length, c, dx, dt, attenuation):
        self.length = length
        self.wave = numpy.zeros((3, self.length))

        self.c = c
        self.dt = dt
        self.dx = dx
        self.attenuation = attenuation
        self.refreshConstants()

        self.boundaryCondition = 1  # 0: const, 1: free
        self.pick(0, 0)

    def refreshConstants(self):
        self.alpha = (self.c * self.dt / self.dx)**2 / 2
        self.beta = 2 * (1 - self.alpha)

    def step(self):
        self.wave = numpy.roll(self.wave, 1, axis=0)
        self.wave[0][0] = 0
        last = len(self.wave[0]) - 1

        if self.boundaryCondition == 1:
            self.wave[0][last] = self.attenuation * (
                self.alpha * (self.wave[1][last - 1] + self.wave[1][last - 1])
                + self.beta * self.wave[1][last] - self.wave[2][last])
        else:
            self.wave[0][last] = 0

        for x in range(1, last):
            self.wave[0][x] = self.attenuation * (
                self.alpha * (self.wave[1][x + 1] + self.wave[1][x - 1]) +
                self.beta * self.wave[1][x] - self.wave[2][x])

        if self.pickY != 0:
            self.wave[0][self.pickX] = self.pickY


wave1d = Wave1D(512, 40, 0.1, 1 / 60, 1)
