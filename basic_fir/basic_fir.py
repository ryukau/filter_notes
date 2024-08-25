import numpy as np
import matplotlib.pyplot as plt
import scipy.signal as signal


def plotResponse(fir, cutoffLow=None, cutoffHigh=None, name=""):
    sampleRate = 1

    fig, ax = plt.subplots(2, 2)

    ω, response = signal.freqz(fir, fs=sampleRate, worN=2048)
    gain = 20 * np.log10(np.abs(response))
    phase = np.unwrap(np.angle(response))
    _, delay = signal.group_delay((fir, 1), fs=sampleRate, w=2048)

    ax[0][0].plot(ω, gain, color="black")
    ax[0][1].plot(ω, delay, color="black")
    ax[1][0].plot(ω, phase, color="black")
    ax[1][1].plot(fir, color="black")

    ax[0][0].set_ylabel("Gain [dB]")
    ax[0][0].set_xlabel(f"Normalized Frequency [rad/2π]")
    ax[0][0].set_ylim((-60, 10))
    # ax[0][0].set_xscale("log")

    ax[0][1].set_ylabel("Group Delay [sample]")
    ax[0][1].set_xlabel(f"Normalized Frequency [rad/2π]")
    ax[0][1].set_ylim((-0.5, 64.5))
    # ax[0][1].set_xscale("log")

    ax[1][0].set_ylabel("Phase [rad]")
    ax[1][0].set_xlabel(f"Normalized Frequency [rad/2π]")
    # ax[1][0].set_xscale("log")

    ax[1][1].set_ylabel("Amplitude of FIR Coefficient")
    ax[1][1].set_xlabel(f"Time [sample]")
    # ax[1][1].set_xticks(np.arange(len(fir)))

    def setCutoffLine(axis):
        if cutoffLow is not None:
            axis.axvline(cutoffLow, color="red", alpha=0.5, label="Cutoff Low")
        if cutoffHigh is not None:
            axis.axvline(cutoffHigh, color="red", alpha=0.5, label="Cutoff High")
        axis.legend(ncol=1)

    setCutoffLine(ax[0][0])
    setCutoffLine(ax[0][1])
    setCutoffLine(ax[1][0])

    for row in ax:
        for axis in row:
            # axis.set_xlim((sampleRate / 8, sampleRate / 2))
            axis.grid(which="both", color="#f0f0f0")

    fig.set_size_inches((10, 6))
    plt.tight_layout()
    plt.show()
    # plt.savefig(f"img/fir_filter_{name}_response.svg")
    # plt.close()


def sinc(x):
    if x == 0:
        return 1
    return np.sin(np.pi * x) / (np.pi * x)


def modifiedSinc(x, cutoff):
    if x == 0:
        return 2 * cutoff
    return np.sin(np.pi * 2 * cutoff * x) / (np.pi * x)


def lowpassFir(length: int, cutoff: float):
    length -= (length + 1) % 2  # 係数の数を奇数にする。
    mid = length // 2

    fir = np.zeros(length)
    for i in range(length):
        x = i - mid
        fir[i] = modifiedSinc(x, cutoff)
    return fir


def highpassFir(length: int, cutoff: float):
    fir = -lowpassFir(length, cutoff)
    mid = length // 2
    fir[mid] -= np.sum(fir)
    return fir


def bandpassFir(length: int, cutoffLow: float, cutoffHigh: float):
    length -= (length + 1) % 2  # 係数の数を奇数にする。
    mid = length // 2

    fir = np.zeros(length)
    for i in range(length):
        x = i - mid
        fir[i] = modifiedSinc(x, cutoffLow) - modifiedSinc(x, cutoffHigh)
    return fir


def bandrejectFir(length: int, cutoffLow: float, cutoffHigh: float):
    fir = lowpassFir(length, cutoffLow) + highpassFir(length, cutoffHigh)
    return fir


bandLow = 0.125
bandHigh = 0.25
plotResponse(lowpassFir(63, bandHigh), bandHigh, None, "lowpass")
plotResponse(highpassFir(63, bandLow), None, bandLow, "highpass")
plotResponse(bandpassFir(63, bandLow, bandHigh), bandLow, bandHigh, "bandpass")
plotResponse(bandrejectFir(63, bandLow, bandHigh), bandLow, bandHigh, "bandreject")
