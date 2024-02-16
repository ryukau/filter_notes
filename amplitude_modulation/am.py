import scipy.signal as signal
import numpy as np
import matplotlib.pyplot as plt


def checkFreq(length, sampleRate):
    freq = np.fft.rfftfreq(length, 1 / sampleRate)
    print(freq[length // 8])


def applySpectrumBandpass(sig, lowerCutoff, upperCutoff):
    """
    `lowerCutoff` and `upperCutoff` are normalized frequencies in [0, 0.5).
    """
    spectrum = np.fft.rfft(sig)

    # Highpass.
    lowerIndex = int(len(spectrum) * 2 * lowerCutoff)
    spectrum[: min(lowerIndex, len(spectrum))] = 0

    # Lowpass.
    upperIndex = int(len(spectrum) * 2 * upperCutoff)
    spectrum[max(0, upperIndex) :] = 0

    return np.fft.irfft(spectrum)


def generateSin(nFrame: int, normalizedFreq):
    return generateLinearSweep(nFrame, normalizedFreq, normalizedFreq)


def generateLinearSweep(nFrame: int, startFreq=0, endFreq=0.5):
    freq = np.linspace(startFreq, endFreq, nFrame)
    phase = np.empty(nFrame)
    phi = 0
    for i in range(nFrame):
        phi += freq[i]
        phi -= np.floor(phi)
        phase[i] = phi
    return np.sin(2 * np.pi * phase)


def upperAM(carrior, modulator):
    c0 = signal.hilbert(carrior)
    m0 = signal.hilbert(modulator)
    return c0.real * m0.real - c0.imag * m0.imag


def lowerAM(carrior, modulator):
    c0 = signal.hilbert(carrior)
    m0 = signal.hilbert(modulator)
    return c0.real * m0.real + c0.imag * m0.imag


def frequencyShift(sig, normalizedFreq):
    analytic = signal.hilbert(sig)
    norm = np.abs(analytic)
    theta = np.angle(analytic)
    time = np.linspace(0, len(analytic), len(analytic))
    return norm * np.cos(theta + 2 * np.pi * normalizedFreq * time)


def antiAliasedUpperAM(carrior, modulator):
    """Remove aliasing using 2-fold up-sampling."""
    upCar = signal.resample(carrior, 2 * len(carrior))
    upMod = signal.resample(modulator, 2 * len(modulator))
    return signal.resample(upperAM(upCar, upMod), len(carrior))


def antiAliasedLowerAM(carrior, modulator):
    """
    Steps are below. `[lower, upper)` is a range of normalized frequencies:

    1. 2-fold up-sampling.
    2. Double the frequencies of carrior using frequency shifter.
    3. Apply lower side-band AM.
    4. Apply highpass with cutoff at 0.25.
    5. Shift frequencies of AM signal to [0, 0.25).
    6. 2-fold down-sampling.
    """
    upCar = signal.resample(carrior, 2 * len(carrior))
    upMod = signal.resample(modulator, 2 * len(modulator))

    shiftedCar = frequencyShift(upCar, 0.25)
    am = lowerAM(shiftedCar, upMod)
    filtered = applySpectrumBandpass(am, 0.25, 0.5)
    result = frequencyShift(filtered, -0.25)

    return signal.resample(result, len(carrior))


def antiAliasedAM(carrior, modulator):
    """Same as `antiAliasedUpperAM`."""
    upCar = signal.resample(carrior, 2 * len(carrior))
    upMod = signal.resample(modulator, 2 * len(modulator))
    return signal.resample(upCar * (1 + upMod), len(carrior))


def antiAliasedAMFull(carrior, modulator):
    """
    Similar to `antiAliasedLowerAM`. Main difference is that 3-fold up or down-sampling is used because of upper side-band.
    """
    upCar = signal.resample(carrior, 3 * len(carrior))
    upMod = signal.resample(modulator, 3 * len(modulator))

    shiftedCar = frequencyShift(upCar, 1 / 6)
    am = shiftedCar * (1 + upMod)
    filtered = applySpectrumBandpass(am, 1 / 6, 2 / 6)
    result = frequencyShift(filtered, -1 / 6)

    return signal.resample(result, len(carrior))


class FirstOrderAllpass:
    def __init__(self, a):
        self.x1 = 0
        self.y1 = 0
        self.a = a

    def process(self, x0):
        self.y1 = self.a * (x0 - self.y1) + self.x1
        self.x1 = x0
        return self.y1


def getHalfbandAllpass():
    h0_a = [
        0.0765690656031399,
        0.264282270318935,
        0.47939467893641907,
        0.661681722389424,
        0.7924031566294969,
        0.8776927911111817,
        0.9308500986629166,
        0.9640156636878193,
        0.9862978287283355,
    ]
    h1_a = [
        0.019911761024506557,
        0.16170648261075027,
        0.37320978687920564,
        0.5766558985008232,
        0.7334355636406803,
        0.8399227128761151,
        0.9074601780285125,
        0.9492937701934973,
        0.9760539731706528,
        0.9955323321150525,
    ]
    return (
        [FirstOrderAllpass(a) for a in h0_a],
        [FirstOrderAllpass(a) for a in h1_a],
    )


def applyAllpass(sig, allpasses):
    for ap in allpasses:
        sig = ap.process(sig)
    return sig


def applyHalfBandUpsampler(sig):
    ap0, ap1 = getHalfbandAllpass()
    output = np.zeros(2 * len(sig))
    for idx, x0 in enumerate(sig):
        output[2 * idx + 0] = applyAllpass(x0, ap1)
        output[2 * idx + 1] = applyAllpass(x0, ap0)
    return output


def applyHalfBandDownsampler(sig):
    ap0, ap1 = getHalfbandAllpass()
    output = np.zeros(len(sig) // 2)
    for idx in range(len(output)):
        y0 = applyAllpass(sig[2 * idx + 0], ap0)
        y1 = applyAllpass(sig[2 * idx + 1], ap1)
        output[idx] = 0.5 * (y0 + y1)
    return output


def toDecibel(x):
    """-144.5 dB = 2^-24, that is the machine epsilon of float32."""
    absed = np.abs(x)
    dB = 20 * np.log10(absed / np.max(absed))
    return np.where(dB < -144.5, -144.5, dB)


def plotSig(axis, sig):
    axis[0].set_ylabel("Amplitude")
    axis[0].set_xlabel("Time [sample]")
    axis[0].plot(sig)

    window = signal.get_window("hann", 511)
    sft = signal.ShortTimeFFT(window, len(sig) // 1000, 1)
    stspc = sft.stft(sig)
    axis[1].set_ylabel("Frequency [rad/2π]")
    axis[1].set_xlabel("Time [sample]")
    axis[1].imshow(
        toDecibel(stspc),
        cmap="magma",
        origin="lower",
        aspect="auto",
        extent=sft.extent(len(sig)),
    )

    spectrum = np.fft.rfft(sig, norm="ortho")
    absed = np.abs(spectrum)
    freq = np.fft.rfftfreq(len(sig), 1)
    axis[2].set_ylabel("Gain [dB]")
    axis[2].set_xlabel("Frequency [rad/2π]")
    # axis[2].set_xscale("log")
    axis[2].plot(freq, 20 * np.log10(absed))


def testHalfBandFilters():
    nFrame = 20000
    sig = generateLinearSweep(nFrame)

    up = applyHalfBandUpsampler(sig)
    down = applyHalfBandDownsampler(sig)

    fig, ax = plt.subplots(3, 3)
    ax = list(map(list, zip(*ax)))
    fig.set_size_inches((15, 9))
    fig.suptitle(
        "Half-band Filter Test (Top: Signal, Mid: Spectrogram, Bottom: Abs. Spectrum)"
    )

    plotSig(ax[0], sig)
    plotSig(ax[1], up)
    plotSig(ax[2], down)

    ax[0][0].set_title("Source")
    ax[1][0].set_title("Up-sampled")
    ax[2][0].set_title("Down-sampled")

    fig.tight_layout()
    plt.show()


def antiAliasedAMLowLatency(carrior, modulator):
    """Only upper side-band is antialiased."""
    upCar = applyHalfBandUpsampler(carrior)
    upMod = applyHalfBandUpsampler(modulator)
    return applyHalfBandDownsampler(upCar * (1 + upMod))


def printSos(sos):
    print([v[:3] + v[4:] for v in sos.tolist()])


def antiAliasedAMLowLatencyFull(carrior, modulator):
    def expand3(sig):
        y = np.zeros(3 * len(sig))
        y[0::3] = sig
        return y

    # Highpass stopband attenuation is -60 dB to achieve sharper fall off.
    lowCut = (1 + 60 / 48000) / 6
    highCut = 1.925 / 6
    sosBp = np.vstack(
        [
            signal.ellip(16, 0.01, 60, lowCut, "highpass", output="sos", fs=1),
            signal.ellip(16, 0.01, 140, highCut, "lowpass", output="sos", fs=1),
        ]
    )
    sosLp = signal.ellip(16, 0.01, 140, 0.925 / 6, "lowpass", output="sos", fs=1)

    printSos(sosBp)
    # printSos(sosLp)

    upCar = signal.sosfilt(sosLp, expand3(carrior))
    upMod = signal.sosfilt(sosLp, expand3(modulator))

    shiftedCar = frequencyShift(upCar, 1 / 6)
    am = shiftedCar * (1 + upMod)
    # filtered = applySpectrumBandpass(am, 1 / 6, 2 / 6)
    filtered = signal.sosfilt(sosBp, am)
    result = frequencyShift(filtered, -1 / 6)

    return signal.sosfilt(sosLp, result)[::3]


def antiAliasedAMLowQuality(carrior, modulator):
    """
    Light weight antialiasing of AM. Only upper side-band is antialiased.

    - 2-fold up-sampling using linear interpolation.
    - Down-sampling lowpass is Butterworth.
    - Cutoff frequency of lowpasses are far below the Nyquist frequency.
    """
    x = np.linspace(0, len(carrior), 2 * len(carrior), endpoint=False)
    xp = np.arange(len(carrior))
    upCar = np.interp(x, xp, carrior)
    upMod = np.interp(x, xp, modulator)

    preLowpass = signal.butter(2, 3 / 16, output="sos", fs=2)
    upCarLp = signal.sosfilt(preLowpass, upCar)
    upModLp = signal.sosfilt(preLowpass, upMod)

    am = upCarLp * (1 + upModLp)
    sos = signal.butter(2, 1 / 3, output="sos", fs=2)
    result = signal.sosfilt(sos, am)
    return result[::2]


def testAm():
    sampleRate = 48000
    durationSeconds = 1
    nFrame = durationSeconds * sampleRate

    # rng = np.random.default_rng(6814384)
    # noise = rng.uniform(-0.1, 0.1, nFrame)
    # filtered = applySpectrumBandpass(noise, 0.0, 0.05)

    # car = filtered
    car = generateLinearSweep(nFrame)
    mod = generateSin(nFrame, 0.25)

    # result = car
    # result = car * (1 + mod)
    # result = car + upperAM(car, mod)
    # result = car + lowerAM(car, mod)
    # result = frequencyShift(car, 0.1)
    # result = antiAliasedUpperAM(car, mod)
    # result = antiAliasedLowerAM(car, mod)
    # result = antiAliasedAM(car, mod)
    # result = antiAliasedAMFull(car, mod)
    # result = antiAliasedAMLowLatency(car, mod)
    result = antiAliasedAMLowLatencyFull(car, mod)
    # result = antiAliasedAMLowQuality(car, mod)

    fig, ax = plt.subplots(3, 1)
    fig.set_size_inches((8, 8))

    plotSig(ax, result)
    ax[0].set_title("AM Signal")
    ax[1].set_title("Spectrogram")
    ax[2].set_title("Absolute Spectrum")

    fig.tight_layout()
    plt.show()


if __name__ == "__main__":
    # testHalfBandFilters()
    testAm()
