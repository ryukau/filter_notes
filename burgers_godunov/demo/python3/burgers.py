import numpy

class Burgers1D:
    """
    dx = 1, dt = 1 で固定した inviscid Burgers' equation 。
    """

    def __init__(self, length):
        """length はシミュレーションする波の配列の長さ。"""
        self.wave = numpy.zeros((2, length))
        self.reset()

    def u_star(self, u, v):
        return numpy.where(
            u >= v,
            numpy.where((u + v) / 2 > 0, u, v),
            numpy.where(u > 0, u, numpy.where(v < 0, v, 0)),
        )

    def step(self):
        self.wave = numpy.roll(self.wave, 1, axis=0)

        last = self.wave.shape[1] - 1

        wave_l = numpy.roll(self.wave[1], 1)
        wave_r = numpy.roll(self.wave[1], -1)
        u_star_l = self.u_star(self.wave[1], wave_r)
        u_star_r = self.u_star(wave_l, self.wave[1])
        self.wave[0] = self.wave[1] - (u_star_l * u_star_l - u_star_r * u_star_r) / 2

        self.wave[0][0] = 0
        self.wave[0][last] = 0

        if self.pick_y != 0:
            self.wave[0][self.pick_x] = self.pick_y

    def reset(self):
        self.wave.fill(0)

    def pick(self, x, y):
        """
        x の範囲は [0, 1] 。
        y の範囲は [-1, 1] 。 |y| > 1 で発散する。
        """
        x = max(0, min(x, 1))
        self.pick_x = numpy.int32(x * self.wave.shape[1])
        self.pick_y = y

    def state(self):
        """この関数はできれば使わないほうが速い。"""
        return self.wave[0]
