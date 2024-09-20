import numpy as np
import scipy.signal as signal
import matplotlib.pyplot as plt


def toDecibel(response):
    return 20 * np.log10(np.maximum(np.abs(response), 1e-5))


class LP1:
    def __init__(self, cutoffNormalized) -> None:
        self.b0 = 1
        self.a1 = 1
        self.x1 = 0
        self.y1 = 0
        self.setCutoff(cutoffNormalized)

    def reset(self) -> None:
        self.x1 = 0
        self.y1 = 0

    def setCutoff(self, cutoffNormalized) -> None:
        cut = np.clip(cutoffNormalized, 0.00001, 0.49998)
        k = 1 / np.tan(np.pi * cut)
        a0 = 1 + k
        self.b0 = 1 / a0
        self.a1 = (1 - k) / a0

    def process(self, x0) -> float:
        self.y1 = self.b0 * (x0 + self.x1) - self.a1 * self.y1
        self.x1 = x0
        return self.y1


class HP1:
    def __init__(self, cutoffNormalized) -> None:
        self.b0 = 1
        self.a1 = 1
        self.x1 = 0
        self.y1 = 0
        self.setCutoff(cutoffNormalized)

    def reset(self) -> None:
        self.x1 = 0
        self.y1 = 0

    def setCutoff(self, cutoffNormalized) -> None:
        cut = np.clip(cutoffNormalized, 0.00001, 0.49998)
        k = 1 / np.tan(np.pi * cut)
        a0 = 1 + k
        self.b0 = k / a0
        self.a1 = (1 - k) / a0

    def process(self, x0) -> float:
        self.y1 = self.b0 * (x0 - self.x1) - self.a1 * self.y1
        self.x1 = x0
        return self.y1


class EmaLP:
    def __init__(self, cutoffNormalized) -> None:
        y = 1 - np.cos(2 * np.pi * cutoffNormalized)
        self.k = np.sqrt((y + 2) * y) - y
        self.v = 0

    def process(self, x0) -> float:
        self.v += self.k * (x0 - self.v)
        return self.v


class EmaHP:
    def __init__(self, cutoffNormalized) -> None:
        y = 1 - np.cos(2 * np.pi * cutoffNormalized)
        self.k = np.sqrt((y + 2) * y) - y
        self.v = 0

    def process(self, x0) -> float:
        self.v += self.k * (x0 - self.v)
        return x0 - self.v


sampleRate: int = 48000
cutoff: float = 1000 / sampleRate

impulse = np.zeros(sampleRate)
impulse[0] = 1

lower = np.zeros_like(impulse)
upper = np.zeros_like(impulse)

lowpass = [LP1(cutoff), LP1(cutoff)]
highpass = [HP1(cutoff), HP1(cutoff)]

for i in range(len(impulse)):
    x = impulse[i]
    lower[i] = lowpass[1].process(lowpass[0].process(x))
    upper[i] = highpass[1].process(highpass[0].process(x))

freq, respLower = signal.freqz(lower, worN=2048, fs=sampleRate)
freq, respUpper = signal.freqz(upper, worN=2048, fs=sampleRate)
freq, respSum = signal.freqz(lower - upper, worN=2048, fs=sampleRate)

plt.plot(freq, toDecibel(respLower), alpha=0.75, color="orange", label="Lower")
plt.plot(freq, toDecibel(respUpper), alpha=0.75, color="red", label="Upper")
plt.plot(freq, toDecibel(respSum), alpha=0.75, color="blue", label="Sum")
plt.xscale("log")
plt.xlabel("Frequency [Hz]")
plt.ylabel("Gain [dB]")
plt.grid()
plt.legend()
plt.show()
