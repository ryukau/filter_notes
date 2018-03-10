import math
import numpy


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
        self.R = 1.0

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
