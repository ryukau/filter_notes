import numpy as np
import scipy.signal as signal
import matplotlib.pyplot as plt


def toDecibel(response):
    return 20 * np.log10(np.maximum(np.abs(response), 1e-7))


def breakFreqToA(freqNormalized):
    t = np.tan(np.pi * freqNormalized)
    return (t - 1) / (t + 1)


def aToBreakFreq(a):
    return -np.arctan((a + 1) / (a - 1)) / np.pi


def unwrapPhase(phase):
    phase /= 2 * np.pi
    phase -= np.floor(phase)
    phase *= 2 * np.pi
    return phase - np.max(phase)


def generateSine(normalizedFrequency: float, length: int):
    sig = np.empty(length)
    phase = 0
    delta = normalizedFrequency
    for i in range(len(sig)):
        sig[i] = np.sin(2 * np.pi * phase)
        phase += delta
        phase -= np.floor(phase)
    return sig


def generateSweep(normalizedFrequency: float, length: int, freqRatio=2):
    sig = np.empty(length)
    phase = 0
    delta = np.geomspace(freqRatio * normalizedFrequency, normalizedFrequency, length)
    for i in range(len(sig)):
        sig[i] = np.sin(2 * np.pi * phase)
        phase += delta[i]
        phase -= np.floor(phase)
    return sig


def generateStepSine(normalizedFrequency: float, length: int, freqRatio=2):
    nStep = 4
    delta = [
        np.full(length // nStep, freq)
        for freq in np.geomspace(
            freqRatio * normalizedFrequency, normalizedFrequency, nStep
        )
    ]
    remaining = length - nStep * (length // nStep)
    delta.append(np.full(remaining, normalizedFrequency))
    delta = np.hstack(delta)

    sig = np.empty(length)
    phase = 0.0
    for i in range(len(sig)):
        sig[i] = np.sin(2 * np.pi * phase)
        phase += delta[i]
        phase -= np.floor(phase)
    return sig


def testNotch():
    fs = 48000
    ir = np.zeros(48000)
    ir[0] = 1

    rho = 0.99
    a = -2 * np.cos(2 * np.pi * 24000 / fs)
    print(f"a={a}")

    x1 = 0
    x2 = 0
    y1 = 0
    y2 = 0
    for i in range(len(ir)):
        x0 = ir[i]
        y0 = x0 + a * x1 + x2 - rho * a * y1 - rho * rho * y2
        ir[i] = y0

        x2 = x1
        x1 = x0
        y2 = y1
        y1 = y0

    freq, resp = signal.freqz(ir, fs=fs)
    plt.plot(freq, toDecibel(resp), color="black")
    # plt.xscale("log")
    plt.grid()
    plt.tight_layout()
    plt.show()


def testNotch2T():
    """
    Direct-form II, transposed.
    """
    fs = 48000
    ir = np.zeros(48000)
    ir[0] = 1

    rho = 0.9
    a = -2 * np.cos(2 * np.pi * 8000 / fs)

    v1 = 0
    v2 = 0
    for i in range(len(ir)):
        a1 = rho * a
        a2 = rho * rho

        x0 = ir[i]
        y0 = x0 + v1
        v1 = a * x0 - a1 * y0 + v2
        v2 = x0 - a2 * y0
        ir[i] = y0

    freq, resp = signal.freqz(ir, fs=fs)
    plt.plot(freq, toDecibel(resp), color="black")
    # plt.xscale("log")
    plt.grid()
    plt.tight_layout()
    plt.show()


def testNotchGainNormalized():
    """
    `input - resonator`.

    Reference: https://ccrma.stanford.edu/~jos/filters/Constant_Peak_Gain_Resonator.html
    """
    fs = 48000
    ir = np.zeros(48000)
    ir[0] = 1

    rho = 0.9
    a = -2 * np.cos(2 * np.pi * 8000 / fs)
    g = (1 - rho * rho) / 2

    x1 = 0
    x2 = 0
    y1 = 0
    y2 = 0
    for i in range(len(ir)):
        x0 = ir[i]
        y0 = x0 - x2 - rho * a * y1 - rho * rho * y2
        ir[i] = x0 - y0 * g

        x2 = x1
        x1 = x0
        y2 = y1
        y1 = y0

    freq, resp = signal.freqz(ir, fs=fs)
    plt.plot(freq, toDecibel(resp), color="black")
    # plt.xscale("log")
    plt.grid()
    plt.tight_layout()
    plt.show()


def adaptiveNotchCpz1(x, rho=0.99, mu=1, initialGuess=0.5):
    """
    Reference: http://www.apsipa.org/proceedings/2018/pdfs/0001355.pdf

    Direct-form I implementaion.
    """
    out = np.zeros_like(x)
    grad = np.zeros_like(x)

    a = -2 * np.cos(2 * np.pi * initialGuess)

    x1 = 0
    x2 = 0
    y1 = 0
    y2 = 0
    s1 = 0
    s2 = 0
    for i in range(len(x)):
        a1 = rho * a
        a2 = rho * rho

        x0 = x[i]
        y0 = x0 + a * x1 + x2 - a1 * y1 - a2 * y2
        out[i] = y0

        s0 = (1 - rho) * x0 - rho * (1 - rho) * x2 - a1 * s1 - a2 * s2
        a = np.clip(a - 2 * mu * y0 * s0, -2, 2)
        grad[i] = a

        x2 = x1
        x1 = x0
        y2 = y1
        y1 = y0
        s2 = s1
        s1 = s0

    sos = [1, a, 1, 1, rho * a, rho * rho]
    return (out, grad, sos)


def adaptiveNotchCpz2(x, rho=0.99, mu=1, initialGuess=0.5):
    """
    Direct-form II.
    """
    out = np.zeros_like(x)
    grad = np.zeros_like(x)

    a = -2 * np.cos(2 * np.pi * initialGuess)

    v1 = 0
    v2 = 0
    for i in range(len(x)):
        a1 = rho * a
        a2 = rho * rho

        x0 = x[i]
        v0 = x0 - a1 * v1 - a2 * v2
        y0 = v0 + a * v1 + v2
        out[i] = y0

        s0 = (1 - rho) * v0 - rho * (1 - rho) * v2
        a = np.clip(a - 2 * mu * y0 * s0, -2, 2)
        grad[i] = a

        v2 = v1
        v1 = v0

    sos = [1, a, 1, 1, rho * a, rho * rho]
    return (out, grad, sos)


def adaptiveNotchCpz2T(x, rho=0.99, mu=1, initialGuess=0.5):
    """
    Direct-form II, transposed.
    """
    out = np.zeros_like(x)
    grad = np.zeros_like(x)

    a = -2 * np.cos(2 * np.pi * initialGuess)

    v1 = 0
    v2 = 0
    w1 = 0
    w2 = 0
    for i in range(len(x)):
        a1 = rho * a
        a2 = rho * rho

        x0 = x[i]
        y0 = x0 + v1
        v1 = a * x0 - a1 * y0 + v2
        v2 = x0 - a2 * y0
        out[i] = y0

        s0 = (1 - rho) * x0 + w1
        w1 = -a1 * s0 + w2
        w2 = rho * (rho - 1) * x0 - a2 * s0
        a = np.clip(a - 2 * mu * y0 * s0, -2, 2)
        grad[i] = a

    sos = [1, a, 1, 1, rho * a, rho * rho]
    return (out, grad, sos)


def adaptiveNotchAM1(x, rho=0.99, mu=1, initialGuess=0.5):
    out = np.zeros_like(x)
    grad = np.zeros_like(x)

    a = -2 * np.cos(2 * np.pi * initialGuess)

    x1 = 0
    x2 = 0
    y1 = 0
    y2 = 0
    for i in range(len(x)):
        x0 = x[i]
        y0 = x0 + a * x1 + x2 - rho * a * y1 - rho * rho * y2
        out[i] = y0

        a = np.clip(a - 2 * mu * y0 * x1, -2, 2)
        grad[i] = a

        x2 = x1
        x1 = x0
        y2 = y1
        y1 = y0

    sos = [1, a, 1, 1, rho * a, rho * rho]
    return (out, grad, sos)


def adaptiveNotchAM2(x, rho=0.99, mu=1, initialGuess=0.5):
    out = np.zeros_like(x)
    grad = np.zeros_like(x)

    a = -2 * np.cos(2 * np.pi * initialGuess)

    x1 = 0
    x2 = 0
    y1 = 0
    y2 = 0

    q1 = 0
    r1 = 0
    for i in range(len(x)):
        x0 = x[i]
        y0 = x0 + a * x1 + x2 - rho * a * y1 - rho * rho * y2
        out[i] = y0

        omega_a = np.pi - np.arccos(a / 2)
        t = np.tan(omega_a / 2)
        k_ap = (t - 1) / (t + 1)

        r1 = k_ap * (x0 - r1) + q1
        q1 = x0
        a = np.clip(a - mu * np.abs(y0) * x0 * r1, -2, 2)

        grad[i] = a

        x2 = x1
        x1 = x0
        y2 = y1
        y1 = y0

    sos = [1, a, 1, 1, rho * a, rho * rho]
    return (out, grad, sos)


class BiquadNotch:
    def __init__(self):
        self.x1 = 0
        self.x2 = 0
        self.y1 = 0
        self.y2 = 0

    def sos(self, cutoffNormalized, notchWidth):
        theta = 2 * np.pi * cutoffNormalized
        cs = np.cos(theta)
        sn = np.sin(theta)
        alpha = sn * np.sinh((np.log(2) * notchWidth * theta) / (2 * sn))

        a0 = 1 + alpha
        a1 = (-2 * cs) / a0
        a2 = (1 - alpha) / a0
        b0 = (1) / a0
        b1 = (-2 * cs) / a0
        b2 = (1) / a0

        return [b0, b1, b2, 1, a1, a2]

    def process(self, x0, cutoffNormalized, notchWidth):
        sos = self.sos(cutoffNormalized, notchWidth)

        y0 = (
            sos[0] * x0
            + sos[1] * self.x1
            + sos[2] * self.x2
            - sos[4] * self.y1
            - sos[5] * self.y2
        )

        self.x2 = self.x1
        self.x1 = x0
        self.y2 = self.y1
        self.y1 = y0

        return y0


def adaptiveNotchAM3(x, rho=0.99, mu=1, initialGuess=0.5):
    out = np.zeros_like(x)
    grad = np.zeros_like(x)

    minCut = 0.00001
    maxCut = 0.49998

    notch = BiquadNotch()
    cutoffNormalized = np.clip(initialGuess, minCut, maxCut)

    q1 = 0
    r1 = 0
    for i in range(len(x)):
        x0 = x[i]
        out[i] = notch.process(x0, cutoffNormalized, rho)

        k_ap = breakFreqToA(cutoffNormalized)

        r1 = k_ap * (x0 - r1) + q1
        q1 = x0
        cutoffNormalized = np.clip(cutoffNormalized - mu * x0 * r1, minCut, maxCut)

        grad[i] = cutoffNormalized

    return (out, grad, notch.sos(cutoffNormalized, rho))


def adaptiveNotchMix(x, rho=0.99, mu=1, initialGuess=0.5):
    out = np.zeros_like(x)
    grad = np.zeros_like(x)

    a = -2 * np.cos(2 * np.pi * initialGuess)

    v1 = 0
    v2 = 0
    q1 = 0
    r1 = 0
    for i in range(len(x)):
        a1 = rho * a
        a2 = rho * rho

        x0 = x[i]
        v0 = x0 - a1 * v1 - a2 * v2
        y0 = v0 + a * v1 + v2
        out[i] = y0

        notchFreq = (np.pi - np.arccos(a / 2)) / (2 * np.pi)
        k_ap = breakFreqToA(notchFreq)

        r1 = k_ap * (x0 - r1) + q1
        q1 = x0
        a_AM = np.clip(a - mu * np.abs(y0) * x0 * r1, -2, 2)

        s0 = (1 - rho) * v0 - rho * (1 - rho) * v2
        a_CPZ = np.clip(a - 2 * mu * y0 * s0, -2, 2)

        a = 0.5 * (a_AM + a_CPZ)

        grad[i] = a

        v2 = v1
        v1 = v0

    sos = [1, a, 1, 1, rho * a, rho * rho]
    return (out, grad, sos)


def compareAdaptiveNotch():
    seed = 4564
    fs = 48000
    sineFreq = 1000 / fs
    initialGuess = 0.0
    noiseAmp = 0.0
    rho = 0.90
    mu = 1 / 2**9

    length = 20480 * 4

    rng = np.random.default_rng(seed)
    sig = (
        rng.normal(0, noiseAmp / 3, length)
        # + generateStepSine(sineFreq, length, 1 / 10)
        + generateStepSine(sineFreq, length, 10)
    )
    # sig = rng.normal(0, noiseAmp / 3, length) + generateSine(sineFreq, length)
    # outB, gradB, sosB = adaptiveNotchCpz1(sig, rho, mu, initialGuess)
    outA, gradA, sosA = adaptiveNotchCpz2(sig, rho, mu, initialGuess)
    # outB, gradB, sosB = adaptiveNotchCpz2T(sig, rho, mu, initialGuess)
    # outB, gradB, sosB = adaptiveNotchAM1(sig, rho, mu, initialGuess)
    outB, gradB, sosB = adaptiveNotchAM2(sig, rho, mu, initialGuess)
    # outB, gradB, sosB = adaptiveNotchAM3(sig, rho, mu, initialGuess)
    # outB, gradB, sosB = adaptiveNotchMix(sig, rho, mu, initialGuess)

    fig, ax = plt.subplots(2, 2)
    fig.set_size_inches((12, 8))

    ax[0][0].set_title("Output")
    ax[0][0].axhline(noiseAmp, color="black", alpha=0.75, lw=1, ls="--", label="Target")
    ax[0][0].axhline(-noiseAmp, color="black", alpha=0.75, lw=1, ls="--")
    ax[0][0].plot(sig, color="black", alpha=0.5, label="Source")
    ax[0][0].plot(outA, color="red", alpha=0.5, label="A")
    ax[0][0].plot(outB, color="blue", alpha=0.5, label="B")
    ax[0][0].set_xlabel("Time [sample]")
    ax[0][0].set_ylabel("Amplitude")

    ax[1][0].set_title("Gradient")
    ax[1][0].plot(gradA, color="red", alpha=0.5, label="A")
    ax[1][0].plot(gradB, color="blue", alpha=0.5, label="B")
    ax[1][0].set_xlabel("Time [sample]")
    ax[1][0].set_ylabel("Amplitude")

    ax[0][1].set_title("Output Spectrum")
    ax[0][1].axvline(
        sineFreq * fs, color="black", alpha=0.75, lw=1, ls="--", label="Target"
    )
    freq = np.fft.rfftfreq(len(sig), 1 / fs)
    resp = np.fft.rfft(sig) / len(sig)
    ax[0][1].plot(freq, toDecibel(resp), color="black", alpha=0.5, label="Source")
    resp = np.fft.rfft(outA) / len(sig)
    ax[0][1].plot(freq, toDecibel(resp), color="red", alpha=0.5, label="A")
    resp = np.fft.rfft(outB) / len(sig)
    ax[0][1].plot(freq, toDecibel(resp), color="blue", alpha=0.5, label="B")
    ax[0][1].set_xscale("log")
    ax[0][1].set_xlabel("Frequency [Hz]")
    ax[0][1].set_ylabel("Amplitude [dB]")

    ax[1][1].set_title("Last Amplitude Response")
    ax[1][1].axvline(
        sineFreq * fs, color="black", alpha=0.75, lw=1, ls="--", label="Target"
    )
    freq, resp = signal.sosfreqz(sosA, worN=2**16, fs=fs)
    ax[1][1].plot(freq, toDecibel(resp), color="red", alpha=0.5, label="A")
    freq, resp = signal.sosfreqz(sosB, worN=2**16, fs=fs)
    ax[1][1].plot(freq, toDecibel(resp), color="blue", alpha=0.5, label="B")
    ax[1][1].set_xscale("log")
    ax[1][1].set_xlabel("Frequency [Hz]")
    ax[1][1].set_ylabel("Gain [dB]")

    for row in ax:
        for axis in row:
            axis.grid()
            axis.legend()
    plt.tight_layout()
    plt.show()


def testNormalize():
    def normalizeGain(sos):
        zeroGain = (sos[0] + sos[1] + sos[2]) / (1 + sos[4] + sos[5])
        nyquistGain = (sos[0] - sos[1] + sos[2]) / (1 - sos[4] + sos[5])
        gain = zeroGain if zeroGain >= nyquistGain else nyquistGain

        s = sos.copy()
        s[0] /= gain
        s[1] /= gain
        s[2] /= gain
        return s

    def compareGain(sos, a):
        zeroGain = (sos[0] + sos[1] + sos[2]) / (1 + sos[4] + sos[5])
        nyquistGain = (sos[0] - sos[1] + sos[2]) / (1 - sos[4] + sos[5])
        gain0 = zeroGain if zeroGain >= nyquistGain else nyquistGain
        gain1 = zeroGain if a >= 0 else nyquistGain
        print(gain0, gain1, gain0 - gain1)

    fs = 48000
    freqHz = 6000
    rho = 0.1  # Notch narrowness.
    a = -2 * np.cos(2 * np.pi * freqHz / fs)  # Initial guess.

    sos = [1, a, 1, 1, rho * a, rho * rho]
    compareGain(sos, a)
    sos = normalizeGain(sos)
    freq, resp = signal.sosfreqz(sos, fs=fs)
    plt.plot(freq, toDecibel(resp), color="black", alpha=0.5, label="source")
    # plt.xscale("log")
    plt.grid()
    plt.legend()
    plt.tight_layout()
    plt.show()


def multiple_formatter(denominator=2, number=np.pi, latex="\pi"):
    """
    https://stackoverflow.com/questions/40642061/how-to-set-axis-ticks-in-multiples-of-pi-python-matplotlib
    """

    def gcd(a, b):
        while b:
            a, b = b, a % b
        return a

    def _multiple_formatter(x, pos):
        den = denominator
        num = int(np.rint(den * x / number))
        com = gcd(num, den)
        (num, den) = (int(num / com), int(den / com))
        if den == 1:
            if num == 0:
                return r"$0$"
            if num == 1:
                return r"$%s$" % latex
            elif num == -1:
                return r"$-%s$" % latex
            else:
                return r"$%s%s$" % (num, latex)
        else:
            # if num == 1:
            #     return r"$\frac{%s}{%s}$" % (latex, den)
            # elif num == -1:
            #     return r"$\frac{-%s}{%s}$" % (latex, den)
            # else:
            #     return r"$\frac{%s%s}{%s}$" % (num, latex, den)
            if num == 1:
                return r"$%s/%s$" % (latex, den)
            elif num == -1:
                return r"$-%s/%s$" % (latex, den)
            else:
                return r"$%s%s/%s$" % (num, latex, den)

    return _multiple_formatter


class Multiple:
    def __init__(self, denominator=2, number=np.pi, latex="\pi"):
        self.denominator = denominator
        self.number = number
        self.latex = latex

    def locator(self):
        return plt.MultipleLocator(self.number / self.denominator)

    def formatter(self):
        return plt.FuncFormatter(
            multiple_formatter(self.denominator, self.number, self.latex)
        )


def plotNotchResponse():
    cutoffNormalized = 0.1
    rho = 0.8

    a = -2 * np.cos(2 * np.pi * cutoffNormalized)
    sos = [1, a, 1, 1, rho * a, rho * rho]
    freq, resp = signal.sosfreqz(sos, worN=2048, fs=1)

    fig, ax = plt.subplots(3)
    fig.suptitle(f"Notch Filter (cutoff={cutoffNormalized}, ρ={rho})")
    fig.set_size_inches((6, 7))

    ax[0].plot(freq, toDecibel(resp), color="black")
    ax[0].set_ylabel("Gain [dB]")

    ax[1].plot(freq, np.angle(resp), color="black")
    ax[1].set_ylabel("Phase [rad]")
    major = Multiple(2, np.pi, r"\pi")
    minor = Multiple(4, np.pi, r"\pi")
    ax[1].yaxis.set_major_locator(major.locator())
    ax[1].yaxis.set_minor_locator(minor.locator())
    ax[1].yaxis.set_major_formatter(plt.FuncFormatter(major.formatter()))
    ax[1].yaxis.set_minor_formatter(plt.FuncFormatter(minor.formatter()))

    tf = signal.sos2tf([sos])
    freq, gd = signal.group_delay(tf, w=2048, fs=1)
    ax[2].plot(freq, gd, color="black")
    ax[2].set_ylabel("Group Delay [sample]")

    for axis in ax:
        axis.grid(which="both", color="#f0f0f0")
        axis.set_xscale("log")
        axis.axvline(cutoffNormalized, color="black", ls="--", lw=1)
    fig.tight_layout()
    plt.show()


def plotNotchPhaseResponse():
    cutoffNormalized = 0.01
    rho = [0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 0.99]

    fig, axis = plt.subplots(1)
    fig.suptitle(f"Notch Filter Phase Response (Cutoff={cutoffNormalized})")
    fig.set_size_inches((6, 4.5))

    cmap = plt.get_cmap("plasma")
    for index, R in enumerate(rho):
        a = -2 * np.cos(2 * np.pi * cutoffNormalized)
        sos = [1, a, 1, 1, R * a, R * R]
        freq, resp = signal.sosfreqz(sos, worN=2048, fs=1)
        phase = unwrapPhase(np.angle(resp))
        axis.plot(freq, phase, color=cmap(index / len(rho)), alpha=0.75, label=f"ρ={R}")

    axis.set_xlabel("Normalized Frequency [rad/2π]")
    axis.set_ylabel("Phase [rad]")
    major = Multiple(2, np.pi, r"\pi")
    minor = Multiple(4, np.pi, r"\pi")
    axis.yaxis.set_major_locator(major.locator())
    axis.yaxis.set_minor_locator(minor.locator())
    axis.yaxis.set_major_formatter(plt.FuncFormatter(major.formatter()))
    # axis.yaxis.set_minor_formatter(plt.FuncFormatter(minor.formatter()))

    axis.axvline(
        cutoffNormalized, color="black", ls=":", lw=2, alpha=0.25, label="Cutoff"
    )
    axis.axhline(-np.pi / 2, color="black", ls="-", lw=1)
    axis.axhline(-np.pi * 3 / 2, color="black", ls="-", lw=1)

    axis.grid(which="both", color="#f0f0f0")
    axis.set_xscale("log")
    axis.legend()
    fig.tight_layout()
    plt.show()


def plotCpzPhaseResponse():
    cutoffNormalized = 0.01
    rho = 0.8

    fig, axis = plt.subplots(1)
    fig.set_size_inches((6, 4.5))
    axis.set_title("Phase Difference of y and s in CPZ-ANF")

    a = -2 * np.cos(2 * np.pi * cutoffNormalized)
    sosNotch = [1, a, 1, 1, rho * a, rho * rho]
    freq, resp = signal.sosfreqz(sosNotch, worN=2048, fs=1)
    phaseNotch = unwrapPhase(np.angle(resp))
    axis.plot(freq, phaseNotch, color="blue", label="y")

    sosAdapt = [1 - rho, 0, -rho * (1 - rho), 1, rho * a, rho * rho]
    freq, resp = signal.sosfreqz(sosAdapt, worN=2048, fs=1)
    phaseAdapt = np.angle(resp)
    axis.plot(freq, phaseAdapt, color="red", label="s")

    axis.plot(freq, phaseNotch - phaseAdapt, color="orange", label="Diff.")

    axis.set_xlabel("Normalized Frequency [rad/2π]")
    axis.set_ylabel("Phase [rad]")
    major = Multiple(2, np.pi, r"\pi")
    minor = Multiple(4, np.pi, r"\pi")
    axis.yaxis.set_major_locator(major.locator())
    axis.yaxis.set_minor_locator(minor.locator())
    axis.yaxis.set_major_formatter(plt.FuncFormatter(major.formatter()))
    # axis.yaxis.set_minor_formatter(plt.FuncFormatter(minor.formatter()))

    axis.axvline(
        cutoffNormalized, color="black", ls=":", lw=2, alpha=0.25, label="Cutoff"
    )
    axis.axhline(-np.pi / 2, color="black", ls="-", lw=1)
    axis.axhline(-np.pi * 3 / 2, color="black", ls="-", lw=1)

    axis.grid(which="both", color="#f0f0f0")
    axis.set_xscale("log")
    axis.legend()
    fig.tight_layout()
    plt.show()


def plotAllpassPhaseResponse():
    fig, axis = plt.subplots(1)
    fig.set_size_inches((6, 3.5))
    axis.set_title(f"Phase Response of First Order Allpass Filter")

    breakFreq = np.geomspace(5e-3, 0.5, 8)
    aList = breakFreqToA(breakFreq)

    cmap = plt.get_cmap("plasma")
    for index, a in enumerate(aList):
        sos = [a, 1, 0, 1, a, 0]
        freq, resp = signal.sosfreqz(sos, worN=2048, fs=1)
        phase = np.angle(resp)
        axis.plot(
            freq,
            phase,
            color=cmap(index / len(aList)),
            alpha=0.75,
            label=f"k_AP={a:.3f}",
        )
        # axis.axvline(aToBreakFreq(a), color="black", ls=":", lw=2, alpha=0.25)

    axis.set_xlabel("Normalized Frequency [rad/2π]")
    axis.set_ylabel("Phase [rad]")
    major = Multiple(2, np.pi, r"\pi")
    minor = Multiple(4, np.pi, r"\pi")
    axis.yaxis.set_major_locator(major.locator())
    axis.yaxis.set_minor_locator(minor.locator())
    axis.yaxis.set_major_formatter(plt.FuncFormatter(major.formatter()))
    axis.yaxis.set_minor_formatter(plt.FuncFormatter(minor.formatter()))

    axis.axhline(-np.pi / 2, color="black", ls="-", lw=1)
    # axis.axhline(-np.pi * 3 / 2, color="black", ls="-", lw=1)

    axis.grid(which="both", color="#f0f0f0")
    axis.set_xscale("log")
    axis.legend(ncol=1, loc=3)
    fig.tight_layout()
    plt.show()


def plotPhaseDiff():
    def fixGain(sig):
        v = np.max(np.abs(sig))
        if v != 0:
            return sig / v
        return sig

    fs = 48000
    length = 1024
    notchFreqHz = 100
    oscFreqHz = 1000
    rho = 0.8  # Notch narrowness.

    a = -2 * np.cos(2 * np.pi * notchFreqHz / fs)
    sos = [1, a, 1, 1, rho * a, rho * rho]

    source = generateSine(oscFreqHz / fs, length)
    output = signal.sosfilt(sos, source)

    z = np.exp(2j * np.pi * oscFreqHz / fs)
    response = (sos[0] * z * z + sos[1] * z + sos[2]) / (
        sos[3] * z * z + sos[4] * z + sos[5]
    )
    print(f"Gain : {np.abs(response)}\nPhase: {np.angle(response)}")

    trim = 128
    source = source[trim:]
    output = output[trim:]

    plt.plot(source, color="blue", alpha=0.5, label="source")
    plt.plot(fixGain(output), color="red", alpha=0.5, label="output")
    plt.plot(fixGain(source * output), color="orange", alpha=0.5, label="AM")
    plt.grid()
    plt.legend()
    plt.tight_layout()
    plt.show()


def plotBias():
    fs = 48000
    length = 2**16
    notchFreqHz = 100
    rho = 0.8  # Notch narrowness.

    a = -2 * np.cos(2 * np.pi * notchFreqHz / fs)
    sos = [1, a, 1, 1, rho * a, rho * rho]

    oscFreqNormalized = np.geomspace(100, fs / 4, 64) / fs
    bias = []
    for freq in oscFreqNormalized:
        source = generateSine(freq, length)
        output = signal.sosfilt(sos, source)
        am = source * output
        bias.append(np.average(am))

    plt.axvline(notchFreqHz, color="black", ls="--", alpha=0.5, label="Cutoff")
    plt.plot(oscFreqNormalized * fs, bias, color="blue", alpha=0.5, label="Average")
    plt.xscale("log")
    plt.grid(which="both", color="#f0f0f0")
    plt.legend()
    plt.tight_layout()
    plt.show()


def plotWidthNormalizedNotchResponse():
    cutoffNormalized = np.geomspace(1e-3, 0.1, 10)
    rhoBase = 0.9
    worN = np.geomspace(1e-4, 0.5, 2**16)

    fig, ax = plt.subplots(3)
    fig.suptitle(f"Width Normalized Notch Filter (ρ={rhoBase})")
    fig.set_size_inches((6, 7))

    cmap = plt.get_cmap("plasma")
    for index, cutoff in enumerate(cutoffNormalized):
        theta = 2 * np.pi * cutoff
        a = -2 * np.cos(theta)
        sn = np.sin(theta)
        rho = sn * rhoBase**sn
        sos = [1, a, 1, 1, rho * a, rho * rho]

        color = cmap(index / len(cutoffNormalized))

        freq, resp = signal.sosfreqz(sos, worN=worN, fs=1)
        ax[0].plot(freq, toDecibel(resp), color=color, label=f"{cutoff:.2e}")
        ax[1].plot(freq, np.angle(resp), color=color)

        tf = signal.sos2tf([sos])
        freq, gd = signal.group_delay(tf, w=worN, fs=1)
        ax[2].plot(freq, gd, color=color)

    ax[0].set_ylabel("Gain [dB]")
    # ax[0].legend()

    ax[1].set_ylabel("Phase [rad]")
    major = Multiple(2, np.pi, r"\pi")
    minor = Multiple(4, np.pi, r"\pi")
    ax[1].yaxis.set_major_locator(major.locator())
    ax[1].yaxis.set_minor_locator(minor.locator())
    ax[1].yaxis.set_major_formatter(plt.FuncFormatter(major.formatter()))
    ax[1].yaxis.set_minor_formatter(plt.FuncFormatter(minor.formatter()))

    ax[2].set_ylabel("Group Delay [sample]")

    for axis in ax:
        axis.grid(which="both", color="#f0f0f0")
        axis.set_xscale("log")
    fig.tight_layout()
    plt.show()


def plotWidthNormalizedNotchResponse2():
    cutoffNormalized = np.geomspace(1e-3, 0.1, 10)
    rhoBase = 0.9
    worN = np.geomspace(1e-4, 0.5, 2**16)

    fig, ax = plt.subplots(3)
    fig.suptitle(f"Width Normalized Notch Filter (ρ={rhoBase})")
    fig.set_size_inches((6, 7))

    cmap = plt.get_cmap("plasma")
    for index, cutoff in enumerate(cutoffNormalized):

        theta = 2 * np.pi * cutoff
        cs = np.cos(theta)
        sn = np.sin(theta)
        alpha = sn * np.sinh((np.log(2) * rhoBase * theta) / (2 * sn))

        b0 = 1
        b1 = -2 * cs
        b2 = 1
        a0 = 1 + alpha
        a1 = -2 * cs
        a2 = 1 - alpha
        sos = np.array([b0 / a0, b1 / a0, b2 / a0, 1, a1 / a0, a2 / a0])

        color = cmap(index / len(cutoffNormalized))

        freq, resp = signal.sosfreqz(sos, worN=worN, fs=1)
        ax[0].plot(freq, toDecibel(resp), color=color, label=f"{cutoff:.2e}")
        ax[1].plot(freq, np.angle(resp), color=color)

        tf = signal.sos2tf([sos])
        freq, gd = signal.group_delay(tf, w=worN, fs=1)
        ax[2].plot(freq, gd, color=color)

    ax[0].set_ylabel("Gain [dB]")
    # ax[0].legend()

    ax[1].set_ylabel("Phase [rad]")
    major = Multiple(2, np.pi, r"\pi")
    minor = Multiple(4, np.pi, r"\pi")
    ax[1].yaxis.set_major_locator(major.locator())
    ax[1].yaxis.set_minor_locator(minor.locator())
    ax[1].yaxis.set_major_formatter(plt.FuncFormatter(major.formatter()))
    ax[1].yaxis.set_minor_formatter(plt.FuncFormatter(minor.formatter()))

    ax[2].set_ylabel("Group Delay [sample]")

    for axis in ax:
        axis.grid(which="both", color="#f0f0f0")
        axis.set_xscale("log")
    fig.tight_layout()
    plt.show()


if __name__ == "__main__":
    # testNotch()
    # testNotchGainNormalized()
    # testNotch2T()
    # testNormalize()
    # plotNotchResponse()
    # plotNotchPhaseResponse()
    # plotCpzPhaseResponse()
    # plotAllpassPhaseResponse()
    # plotPhaseDiff()
    # plotBias()
    # plotWidthNormalizedNotchResponse()
    # plotWidthNormalizedNotchResponse2()
    compareAdaptiveNotch()
