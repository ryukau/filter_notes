import numpy as np
import scipy.signal as signal
import matplotlib.pyplot as plt
import math


class Svf:
    """
    TPT biquad filter described in "The Art of VA Filter Design".

    normalizedFreq = cutoffHz / sampleRate = cutoffRadian / (2 * pi).
    Q is quality factor, or resonance.
    """

    def __init__(self):
        self.reset()

    def reset(self):
        self.ic1eq = 0
        self.ic2eq = 0

    def tick(self, v0, g, k):
        v1 = (self.ic1eq + g * (v0 - self.ic2eq)) / (1 + g * (g + k))
        v2 = self.ic2eq + g * v1
        self.ic1eq = 2 * v1 - self.ic1eq
        self.ic2eq = 2 * v2 - self.ic2eq
        return (v1, v2)

    def processLowpass(self, v0, normalizedFreq, Q):
        g = math.tan(normalizedFreq * math.pi)
        k = 1 / Q
        v1, v2 = self.tick(v0, g, k)
        return v2

    def processBandpass(self, v0, normalizedFreq, Q):
        g = math.tan(normalizedFreq * math.pi)
        k = 1 / Q
        v1, v2 = self.tick(v0, g, k)
        return v1

    def processHighpass(self, v0, normalizedFreq, Q):
        g = math.tan(normalizedFreq * math.pi)
        k = 1 / Q
        v1, v2 = self.tick(v0, g, k)
        return v0 - k * v1 - v2

    def processNotch(self, v0, normalizedFreq, Q):
        g = math.tan(normalizedFreq * math.pi)
        k = 1 / Q
        v1, v2 = self.tick(v0, g, k)
        return v0 - k * v1

    def processPeak(self, v0, normalizedFreq, Q):
        g = math.tan(normalizedFreq * math.pi)
        k = 1 / Q
        v1, v2 = self.tick(v0, g, k)
        return v0 - k * v1 - 2 * v2

    def processAllpass(self, v0, normalizedFreq, Q):
        g = math.tan(normalizedFreq * math.pi)
        k = 1 / Q
        v1, v2 = self.tick(v0, g, k)
        return v0 - 2 * k * v1


class SvfShelving:
    """
    normalizedFreq = cutoffHz / sampleRate
                   = cutoffRadian / (2 * pi).
    Q is quality factor, or resonance.
    """

    def __init__(self, normalizedFreq, Q, shelvingType="low"):
        self.reset()
        self.freq = math.tan(normalizedFreq * math.pi)
        self.k = 1 / Q
        if shelvingType == "bell":
            self.processFunc = self.processBell
        elif shelvingType == "high":
            self.processFunc = self.processHighShelf
        else:  # shelvingType == "low"
            self.processFunc = self.processLowShelf

    def reset(self):
        self.ic1eq = 0
        self.ic2eq = 0

    def tick(self, v0, g, k):
        v1 = (self.ic1eq + g * (v0 - self.ic2eq)) / (1 + g * (g + k))
        v2 = self.ic2eq + g * v1
        self.ic1eq = 2 * v1 - self.ic1eq
        self.ic2eq = 2 * v2 - self.ic2eq
        return (v1, v2)

    def processAmp(self, v0, gainAmp):
        """gainAmp >= 0"""
        return self.processFunc(v0, math.sqrt(gainAmp))

    def processdB(self, v0, gaindB):
        return self.processFunc(v0, math.pow(10, gaindB / 40))

    def processBell(self, v0, A):
        g = self.freq
        kk = self.k / A
        v1, v2 = self.tick(v0, g, kk)
        return v0 + (A * A - 1) * kk * v1

    def processLowShelf(self, v0, A):
        g = self.freq / math.sqrt(A)
        v1, v2 = self.tick(v0, g, self.k)
        return v0 + (A - 1) * self.k * v1 + (A * A - 1) * v2

    def processHighShelf(self, v0, A):
        g = self.freq * math.sqrt(A)
        v1, v2 = self.tick(v0, g, self.k)
        return A * A * (v0 - self.k * v1 - v2) + A * self.k * v1 + v2


def getResponse(sig):
    spec = np.fft.rfft(sig)
    absed = np.abs(spec)
    power = 20 * np.log10(absed)  # Do not normalize amplitude.
    phase = np.angle(spec)
    return (power, phase)


def testSvf():
    """
    Peaking filter gain amplitude is 2 * Q.
    In decibel, `gain = 20 * log10(2 * Q)`.
    Examples:
    - 0.5 when Q = 0.25
    - 1.0 when Q = 0.5
    - 2.0 when Q = 1.0
    """
    sampleRate = 48000
    Q = math.sqrt(2) / 2
    # Q = 0.25
    cutoffHz = 500

    sig = np.zeros(sampleRate)
    sig[0] = 1

    def applySvf(svfMethod, svf, sig, cutoff, Q):
        svf.reset()
        output = np.empty_like(sig)
        for i, value in enumerate(sig):
            output[i] = svfMethod(value, cutoff, Q)
        return output.tolist()

    cutoff = cutoffHz / sampleRate  # Normalized frequency.
    svf = Svf()
    result = {
        "lowpass": applySvf(svf.processLowpass, svf, sig, cutoff, Q),
        "bandpass": applySvf(svf.processBandpass, svf, sig, cutoff, Q),
        "highpass": applySvf(svf.processHighpass, svf, sig, cutoff, Q),
        "notch": applySvf(svf.processNotch, svf, sig, cutoff, Q),
        "peak": applySvf(svf.processPeak, svf, sig, cutoff, Q),
        "allpass": applySvf(svf.processAllpass, svf, sig, cutoff, Q),
    }

    fig, ax = plt.subplots(2, 1)
    cmap = plt.get_cmap("plasma")
    index = 0
    for svfType, sig in result.items():
        power, phase = getResponse(sig)
        color = cmap(index / len(result))
        ax[0].plot(power, lw=1, alpha=0.5, color=color, label=svfType)
        ax[1].plot(phase, lw=1, alpha=0.5, color=color, label=svfType)
        index += 1
    ax[0].set_ylabel("Gain [dB]")
    ax[0].set_ylim((-30, 10))
    ax[1].set_ylabel("Phase [rad/sample]")
    for axis in ax:
        axis.axvline(cutoffHz, lw=1, ls="--", color="black", alpha=0.25, label="cutoff")
        axis.set_xscale("log")
        axis.grid(which="both", color="#f0f0f0")
        axis.legend(ncol=2)
    fig.set_size_inches((9, 6))
    fig.tight_layout()
    plt.show()


def testShelving():
    sampleRate = 48000
    Q = math.sqrt(2) / 2
    cutoffHz = 500
    dB = 6
    amp = math.pow(10, dB / 20)
    # amp = 0.5
    # dB = 20 * math.log10(amp)

    sig = np.zeros(sampleRate)
    sig[0] = 1

    def applySvf(svfMethod, sig, gain):
        output = np.empty_like(sig)
        for i, value in enumerate(sig):
            output[i] = svfMethod(value, gain)
        return output.tolist()

    result = {}
    for shelvingType in ["low", "high", "bell"]:
        svf = SvfShelving(cutoffHz / sampleRate, Q, shelvingType)
        result[shelvingType] = {
            "amp": applySvf(svf.processAmp, sig, amp),
            "dB": applySvf(svf.processdB, sig, dB),
        }

    fig, ax = plt.subplots(2, 1)
    cmap = plt.get_cmap("plasma")
    nElement = len(result["low"]) * len(result)
    index = 0
    for shelvingType, data in result.items():
        for gainType, sig in data.items():
            power, phase = getResponse(sig)
            label = shelvingType + gainType
            c0 = cmap(index / nElement)
            c1 = cmap((index + 1) / nElement)
            ax[0].plot(power, lw=1, alpha=0.5, color=c0, label=label)
            ax[1].plot(phase, lw=1, alpha=0.5, color=c1, label=label)
        index += 2
    ax[0].set_ylabel("Gain [dB]")
    ax[1].set_ylabel("Phase [rad/sample]")
    for axis in ax:
        axis.axvline(cutoffHz, lw=1, ls="--", color="black", alpha=0.25, label="cutoff")
        axis.set_xscale("log")
        axis.grid(which="both", color="#f0f0f0")
        axis.legend()
    fig.set_size_inches((9, 6))
    fig.tight_layout()
    plt.show()


def testPhaseCorrection():
    sampleRate = 48000
    Q = math.sqrt(2) / 2
    # Q = 0.25
    cutoffHz = 500

    sig = np.zeros(sampleRate)
    sig[0] = 1

    def applySvf(svfMethod, svf, sig, cutoff, Q):
        svf.reset()
        output = np.empty_like(sig)
        for i, value in enumerate(sig):
            output[i] = svfMethod(value, cutoff, Q)
        return output.tolist()

    cutoff = cutoffHz / sampleRate  # Normalized frequency.
    svf = Svf()
    allpass = applySvf(svf.processAllpass, svf, sig, cutoff, Q)
    filtered = applySvf(svf.processBandpass, svf, sig, cutoff, Q)
    corrected = applySvf(svf.processAllpass, svf, filtered, cutoff, Q)
    result = {"allpass": allpass, "filtered": filtered, "corrected": corrected}

    fig, ax = plt.subplots(2, 1)
    cmap = plt.get_cmap("plasma")
    index = 0
    for svfType, sig in result.items():
        power, phase = getResponse(sig)
        color = cmap(index / len(result))
        ax[0].plot(power, lw=1, alpha=0.5, color=color, label=svfType)
        ax[1].plot(phase, lw=1, alpha=0.5, color=color, label=svfType)
        index += 1
    ax[0].set_ylabel("Gain [dB]")
    ax[0].set_ylim((-30, 10))
    ax[1].set_ylabel("Phase [rad/sample]")
    for axis in ax:
        axis.axvline(cutoffHz, lw=1, ls="--", color="black", alpha=0.25, label="cutoff")
        axis.set_xscale("log")
        axis.grid(which="both", color="#f0f0f0")
        axis.legend(ncol=2)
    fig.tight_layout()
    fig.set_size_inches((12, 8))
    plt.show()


if __name__ == "__main__":
    testSvf()
    testShelving()
    # testPhaseCorrection()
