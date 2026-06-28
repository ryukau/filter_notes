import numpy as np
import scipy.signal as signal
import matplotlib.pyplot as plt


def toDecibel(response):
    return 20 * np.log10(np.maximum(np.abs(response), 1e-4))


def toPhase(response):
    return (np.angle(response)) / np.pi


def toGroupDelay(response):
    gd = -np.diff(np.unwrap(np.angle(response))) * len(response) / np.pi
    # gd -= gd[0]
    return np.hstack(([gd[0]], gd)) - gd[0]


class ComplexAllpass:
    def __init__(self, cutoffNormalized, R):
        theta = 2 * np.pi * cutoffNormalized
        self.a1 = R * np.exp(1j * theta)
        self.b1 = R / np.exp(1j * theta).conj()
        self.x1 = 0
        self.y1 = 0

    def process(self, x0):
        self.y1 = x0 + self.b1 * self.x1 - self.a1 * self.y1
        self.x1 = x0
        return self.y1


class ComplexLowpass:
    def __init__(self, cutoffNormalized, R):
        theta = 2 * np.pi * cutoffNormalized
        self.a1 = R * np.exp(1j * theta)
        self.x1 = 0
        self.y1 = 0

        gain = (1 - self.a1) / 2
        self.b = gain

    def process(self, x0):
        self.y1 = self.b * (x0 + self.x1) + self.a1 * self.y1
        self.x1 = x0
        return self.y1


class ComplexHighpass:
    def __init__(self, cutoffNormalized, R):
        theta = 2 * np.pi * cutoffNormalized
        self.a1 = R * np.exp(1j * theta)
        self.x1 = 0
        self.y1 = 0

        gain = (1 + self.a1) / 2
        self.b = gain

    def process(self, x0):
        self.y1 = self.b * (x0 - self.x1) + self.a1 * self.y1
        self.x1 = x0
        return self.y1


def complexIR(filterType, cutoffNormalized, phaseNormalized, worN, fs):
    resonator = filterType(cutoffNormalized, phaseNormalized)

    ir = np.zeros(worN, dtype=np.complex128)
    ir[0] = 1
    for idx in range(len(ir)):
        ir[idx] = resonator.process(ir[idx])

    return signal.freqz(ir.real, worN=worN, fs=fs)


def plot(filterType):
    fs = 48000
    # phase = np.linspace(0, 0.9, 10)
    # phase = [np.sqrt(2) / 2]
    phase = [0.9]
    cutoff = 10000 / fs

    cmap = plt.get_cmap("viridis")
    for idx, phs in enumerate(phase):
        freq, resp = complexIR(filterType, cutoff, phs, worN=2**16, fs=fs)
        plt.plot(
            freq, toDecibel(resp), color=cmap(idx / len(phase)), label=f"{phs:.2f}"
        )
    plt.axvline(fs * cutoff, ls="--", lw=1, color="black", alpha=0.1)
    plt.grid(which="both", color="#f0f0f0")
    plt.xlim([10, 20000])
    plt.xscale("log")
    plt.legend()
    plt.tight_layout()
    plt.show()


plot(ComplexLowpass)
