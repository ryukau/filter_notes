import numpy as np
import matplotlib.pyplot as plt
import scipy.signal as signal
import soundfile

class SquareOscillator:
    def __init__(self, samplerate, frequnecy):
        self.x0 = 0
        self.x1 = 0
        self.x2 = 0
        self.x3 = 0
        self.x4 = 0
        self.x5 = 0
        self.x6 = 0
        self.x7 = 0

        # g_(n/2) が 0 でないときに PolyBLEP residual の加算をトリガ。
        self.g0 = 0
        self.g1 = 0
        self.g2 = 0
        self.g3 = 0

        # 分数ディレイ時間。
        self.t0 = 0
        self.t1 = 0
        self.t2 = 0
        self.t3 = 0

        # 位相の範囲は [0, 2] 。
        self.phase = 0
        self.tick = 2 * frequnecy / samplerate

    def polyBlep4_0(self, t):
        return -t**4 / 24 + t**3 / 6 - t**2 / 4 + t / 6 - 1 / 24

    def polyBlep4_1(self, t):
        return t**4 / 8 - t**3 / 3 + (2 * t) / 3 - 1 / 2

    def polyBlep4_2(self, t):
        return -t**4 / 8 + t**3 / 6 + t**2 / 4 + t / 6 + 1 / 24

    def polyBlep4_3(self, t):
        return t**4 / 24

    def process4(self):
        self.x3 = self.x2
        self.x2 = self.x1
        self.x1 = self.x0
        self.x0 = 1 if self.phase < 1 else -1

        self.g1 = self.g0

        self.t1 = self.t0

        if self.x0 == self.x1:
            self.g0 = 0
            self.t0 = 0
        else:
            self.g0 = self.x0 - self.x1
            self.t0 = (self.phase - np.floor(self.phase)) / self.tick

        if self.g1 != 0:
            t = max(0, min(self.t1, 1))  # 念のため。
            self.x0 += self.g1 * self.polyBlep4_0(t)
            self.x1 += self.g1 * self.polyBlep4_1(t)
            self.x2 += self.g1 * self.polyBlep4_2(t)
            self.x3 += self.g1 * self.polyBlep4_3(t)

        self.phase += self.tick
        if self.phase >= 2:
            self.phase -= 2

        return self.x3

    def polyBlep6_0(self, t):
        return -t**6 / 720 + t**5 / 120 - t**4 / 48 + t**3 / 36 - t**2 / 48 + t / 120 - 1 / 720

    def polyBlep6_1(self, t):
        return t**6 / 144 - t**5 / 30 + t**4 / 24 + t**3 / 18 - (5 * t**2) / 24 + (
            13 * t) / 60 - 29 / 360

    def polyBlep6_2(self, t):
        return -t**6 / 72 + t**5 / 20 - t**3 / 6 + (11 * t) / 20 - 1 / 2

    def polyBlep6_3(self, t):
        return t**6 / 72 - t**5 / 30 - t**4 / 24 + t**3 / 18 + (5 * t**2) / 24 + (
            13 * t) / 60 + 29 / 360

    def polyBlep6_4(self, t):
        return -t**6 / 144 + t**5 / 120 + t**4 / 48 + t**3 / 36 + t**2 / 48 + t / 120 + 1 / 720

    def polyBlep6_5(self, t):
        return t**6 / 720

    def process6(self):
        self.x5 = self.x4
        self.x4 = self.x3
        self.x3 = self.x2
        self.x2 = self.x1
        self.x1 = self.x0
        self.x0 = 1 if self.phase < 1 else -1

        self.g2 = self.g1
        self.g1 = self.g0

        self.t2 = self.t1
        self.t1 = self.t0

        if self.x0 == self.x1:
            self.g0 = 0
            self.t0 = 0
        else:
            self.g0 = self.x0 - self.x1
            self.t0 = (self.phase - np.floor(self.phase)) / self.tick

        if self.g2 != 0:
            t = max(0, min(self.t2, 1))  # 念のため。
            self.x0 += self.g2 * self.polyBlep6_0(t)
            self.x1 += self.g2 * self.polyBlep6_1(t)
            self.x2 += self.g2 * self.polyBlep6_2(t)
            self.x3 += self.g2 * self.polyBlep6_3(t)
            self.x4 += self.g2 * self.polyBlep6_4(t)
            self.x5 += self.g2 * self.polyBlep6_5(t)

        self.phase += self.tick
        if self.phase >= 2:
            self.phase -= 2

        return self.x5

    def polyBlep8_0(self, t):
        return -t**8 / 40320 + t**7 / 5040 - t**6 / 1440 + t**5 / 720 - t**4 / 576 + t**3 / 720 - t**2 / 1440 + t / 5040 - 1 / 40320

    def polyBlep8_1(self, t):
        return t**8 / 5760 - t**7 / 840 + t**6 / 360 - t**4 / 72 + t**3 / 30 - (
            7 * t**2) / 180 + t / 42 - 31 / 5040

    def polyBlep8_2(self, t):
        return -t**8 / 1920 + t**7 / 336 - t**6 / 288 - t**5 / 80 + (
            19 * t**4) / 576 + t**3 / 48 - (49 * t**2) / 288 + (397 *
                                                                t) / 1680 - 4541 / 40320

    def polyBlep8_3(self, t):
        return t**8 / 1152 - t**7 / 252 + t**5 / 45 - t**3 / 9 + (151 * t) / 315 - 1 / 2

    def polyBlep8_4(self, t):
        return -t**8 / 1152 + t**7 / 336 + t**6 / 288 - t**5 / 80 - (
            19 * t**4) / 576 + t**3 / 48 + (49 * t**2) / 288 + (397 *
                                                                t) / 1680 + 4541 / 40320

    def polyBlep8_5(self, t):
        return t**8 / 1920 - t**7 / 840 - t**6 / 360 + t**4 / 72 + t**3 / 30 + (
            7 * t**2) / 180 + t / 42 + 31 / 5040

    def polyBlep8_6(self, t):
        return -t**8 / 5760 + t**7 / 5040 + t**6 / 1440 + t**5 / 720 + t**4 / 576 + t**3 / 720 + t**2 / 1440 + t / 5040 + 1 / 40320

    def polyBlep8_7(self, t):
        return t**8 / 40320

    def process8(self):
        self.x7 = self.x6
        self.x6 = self.x5
        self.x5 = self.x4
        self.x4 = self.x3
        self.x3 = self.x2
        self.x2 = self.x1
        self.x1 = self.x0
        self.x0 = 1 if self.phase < 1 else -1

        self.g3 = self.g2
        self.g2 = self.g1
        self.g1 = self.g0

        self.t3 = self.t2
        self.t2 = self.t1
        self.t1 = self.t0

        if self.x0 == self.x1:
            self.g0 = 0
            self.t0 = 0
        else:
            self.g0 = self.x0 - self.x1
            self.t0 = (self.phase - np.floor(self.phase)) / self.tick

        if self.g3 != 0:
            t = max(0, min(self.t3, 1))  # 念のため。
            self.x0 += self.g3 * self.polyBlep8_0(t)
            self.x1 += self.g3 * self.polyBlep8_1(t)
            self.x2 += self.g3 * self.polyBlep8_2(t)
            self.x3 += self.g3 * self.polyBlep8_3(t)
            self.x4 += self.g3 * self.polyBlep8_4(t)
            self.x5 += self.g3 * self.polyBlep8_5(t)
            self.x6 += self.g3 * self.polyBlep8_6(t)
            self.x7 += self.g3 * self.polyBlep8_7(t)

        self.phase += self.tick
        if self.phase >= 2:
            self.phase -= 2

        return self.x7

samplerate = 48000
duration = 1
frequency = 1234
nSample = int(duration * samplerate)

time = np.linspace(0, duration, nSample)

phase = 2 * np.pi * frequency * time
naive = signal.square(phase)

square = SquareOscillator(samplerate, frequency)
polyblep4 = np.array([square.process4() for _ in range(nSample)])
polyblep6 = np.array([square.process6() for _ in range(nSample)])
polyblep8 = np.array([square.process8() for _ in range(nSample)])

soundfile.write("naive.wav", 0.2 * naive, samplerate, subtype="FLOAT")
soundfile.write("polyblep4.wav", 0.2 * polyblep4, samplerate, subtype="FLOAT")
soundfile.write("polyblep6.wav", 0.2 * polyblep6, samplerate, subtype="FLOAT")
soundfile.write("polyblep8.wav", 0.2 * polyblep8, samplerate, subtype="FLOAT")

## Time domain plot.
end = int(samplerate / 500)
cmap = plt.get_cmap("viridis")
plt.title(f"Square Wave, {frequency} Hz")
plt.plot(naive[:end], color=cmap(0 / 4), label="Naive")
plt.plot(polyblep4[:end], color=cmap(1 / 4), label="4 pt.")
plt.plot(polyblep6[:end], color=cmap(2 / 4), label="6 pt.")
plt.plot(polyblep8[:end], color=cmap(3 / 4), label="8 pt.")
plt.xlabel("Time [samples]")
plt.ylabel("Amplitude")
plt.grid()
plt.legend()
plt.tight_layout()
plt.show()

## Frequency domain plot.
def powerSpectrum(sig):
    spectrum = np.abs(np.fft.rfft(sig))
    maximum = np.max(spectrum)
    if maximum == 0:
        return spectrum
    return 20 * np.log10(spectrum / maximum)

for name, spectrum in [
    ("Naive", powerSpectrum(naive)),
    ("4 Point", powerSpectrum(polyblep4)),
    ("6 Point", powerSpectrum(polyblep6)),
    ("8 Point", powerSpectrum(polyblep8)),
]:
    plt.figure(figsize=(8, 3))
    plt.plot(spectrum, lw=1, color="black")
    plt.title(name)
    plt.xlabel("Frequency")
    plt.ylabel("Amplitude [dB]")
    plt.ylim((-100, 10))
    plt.grid()
    plt.tight_layout()
    plt.show()
