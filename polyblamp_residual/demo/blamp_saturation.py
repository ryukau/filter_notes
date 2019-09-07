import numpy
import matplotlib.pyplot as pyplot
import soundfile

def naive_clip(sig, limit):
    sig[sig > limit] = limit
    sig[sig < -limit] = -limit
    return sig

class Clipper:
    def __init__(self, limit):
        self.upper = limit
        self.lower = -limit
        self.reset()

    def reset(self):
        self.s0 = 0
        self.s1 = 0
        self.s2 = 0
        self.s3 = 0
        self.state0 = 0  # 0: bypass, 1: positive clip, 2: negative clip
        self.state1 = 0
        self.state2 = 0
        self.state3 = 0

    def push(self, sig):
        self.s0 = self.s1
        self.s1 = self.s2
        self.s2 = self.s3
        self.s3 = sig
        self.state0 = self.state1
        self.state1 = self.state2
        self.state2 = self.state3

    def lagrange4(self, s0, s1, s2, s3):
        return (
            -1 / 6 * s0 + 1 / 2 * s1 - 1 / 2 * s2 + 1 / 6 * s3,
            s0 - 5 / 2 * s1 + 2 * s2 - 1 / 2 * s3,
            -6 / 11 * s0 + 3 * s1 - 3 / 2 * s2 + 1 / 3 * s3,
            s0,
        )

    def estimate(self, rho, s0, s1, s2, s3):
        slope = numpy.abs(s3 - s2)
        d = numpy.abs((s2 - rho) / slope) if slope > 1e-5 else 1
        while d > 1:  # 角の推定に失敗したら適当に発散しない値を代入する。
            d /= 2

        return (d, slope)

    def blamp_residual(self, d, slope):
        return (
            slope * (d**5 / 120),
            slope * (-d**5 / 40 + d**4 / 24 + d**3 / 12 + d**2 / 12 + d / 24 + 1 / 120),
            slope * (d**5 / 40 - d**4 / 12 + d**2 / 3 - d / 2 + 7 / 30),
            slope * (-d**5 / 120 + d**4 / 24 - d**3 / 12 + d**2 / 12 - d / 24 + 1 / 120),
        )

    def process2(self, sig):
        self.push(sig)

        if self.s3 > self.upper:
            self.state3 = 1
            self.s3 = self.upper
        elif self.s3 < self.lower:
            self.state3 = 2
            self.s3 = self.lower
        else:
            self.state3 = 0

        if self.state2 == 1 and self.state3 == 1 and self.state1 != 1:
            d, slope = self.estimate(self.upper, self.s3, self.s2, self.s1, self.s0)
            p0, p1, p2, p3 = self.blamp_residual(1 - d, slope)
            self.s0 -= p0
            self.s1 -= p1
            self.s2 -= p2
            self.s3 -= p3
        elif self.state2 != 1 and self.state1 == 1 and self.state0 == 1:
            d, slope = self.estimate(self.upper, self.s0, self.s1, self.s2, self.s3)
            p3, p2, p1, p0 = self.blamp_residual(d, slope)
            self.s0 -= p0
            self.s1 -= p1
            self.s2 -= p2
            self.s3 -= p3
        elif self.state2 == 2 and self.state3 == 2 and self.state1 != 2:
            d, slope = self.estimate(self.lower, self.s3, self.s2, self.s1, self.s0)
            p0, p1, p2, p3 = self.blamp_residual(1 - d, slope)
            self.s0 += p0
            self.s1 += p1
            self.s2 += p2
            self.s3 += p3
        elif self.state2 != 2 and self.state1 == 2 and self.state0 == 2:
            d, slope = self.estimate(self.lower, self.s0, self.s1, self.s2, self.s3)
            p3, p2, p1, p0 = self.blamp_residual(d, slope)
            self.s0 += p0
            self.s1 += p1
            self.s2 += p2
            self.s3 += p3

        return self.s0

    def processNaive(self, sig):
        self.push(sig)

        if self.s2 > self.upper:
            self.s2 = self.upper
        elif self.s2 < self.lower:
            self.s2 = self.lower

        return self.s0

def to_decibel(data):
    data_abs = numpy.abs(data)
    return 20 * numpy.log10(data_abs / numpy.max(data_abs))

def signal_to_spectrum(sig):
    return to_decibel(numpy.abs(numpy.fft.rfft(sig)))

plot_signal = False
plot_spectrum = False
samplerate = 44100
frequency = 8000  # 周波数が高いとランプ関数の角の推定に失敗して発散しやすい。

phase = numpy.linspace(0, 2 * numpy.pi * frequency, samplerate)
sin = numpy.sin(phase)
sin = numpy.sin(phase + sin)  # FM

clipper = Clipper(0.3)
naive = [clipper.processNaive(s) for s in sin]

clipper.reset()
blamp = [clipper.process2(s) for s in sin]

soundfile.write("snd/hardclip_naive.wav", naive, samplerate, subtype="FLOAT")
soundfile.write("snd/hardclip_blamp.wav", blamp, samplerate, subtype="FLOAT")

if plot_signal:
    pyplot.plot(sin, label="source", color="gray", ls="--", alpha=0.5)
    pyplot.plot(naive, label="naive", color="red", alpha=0.5)
    pyplot.plot(blamp, label="blamp", color="black", alpha=0.75)
    pyplot.ylim((-1.1, 1.1))
    pyplot.grid()
    pyplot.legend()
    pyplot.show()

if plot_spectrum:
    spec_sin = signal_to_spectrum(sin)
    spec_naive = signal_to_spectrum(naive)
    spec_blamp = signal_to_spectrum(blamp)

    pyplot.plot(spec_sin, label="source", color="gray", ls="--", alpha=0.5)
    pyplot.plot(spec_naive, label="naive", alpha=0.5, lw=1, color="red")
    pyplot.plot(spec_blamp, label="blamp", alpha=0.75, lw=1, color="black")
    pyplot.grid()
    pyplot.legend()
    pyplot.ylim((-200, 10))
    pyplot.show()
