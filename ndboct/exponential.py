import numpy as np
import matplotlib.pyplot as plt
import scipy.signal as signal
import sympy


def toDecibel(response):
    return 20 * np.log10(np.maximum(np.abs(response), 1e-7))


def plotReferenceSlope(axis, startFreq, slopeDB, step, color="red", label="reference"):
    refX = []
    refY = []
    for idx in range(0, step):
        refX.append(startFreq * 2**idx)
        refY.append(slopeDB * idx)
    axis.plot(refX, refY, label=label, color=color, ls="--")


def plotResponse(fir, slopeDecibel):
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
    # ax[0][0].set_xlim((0.5 / 100, 0.5))
    ax[0][0].set_ylim((-60, 10))
    ax[0][0].set_xscale("log")

    plotReferenceSlope(ax[0][0], sampleRate / 2**8, slopeDecibel, 8)

    ax[0][1].set_ylabel("Group Delay [sample]")
    ax[0][1].set_xlabel(f"Normalized Frequency [rad/2π]")
    # ax[0][1].set_ylim((-0.5, 512.5))
    ax[0][1].set_xscale("log")

    ax[1][0].set_ylabel("Phase [rad]")
    ax[1][0].set_xlabel(f"Normalized Frequency [rad/2π]")
    ax[1][0].set_xscale("log")

    ax[1][1].set_ylabel("Amplitude of FIR Coefficient")
    ax[1][1].set_xlabel(f"Time [sample]")
    # ax[1][1].set_xticks(np.arange(len(fir)))

    for row in ax:
        for axis in row:
            # axis.set_xlim((sampleRate / 8, sampleRate / 2))
            axis.grid(which="both", color="#f0f0f0")
            # axis.legend(ncol=1)

    fig.set_size_inches((10, 6))
    plt.tight_layout()
    plt.show()
    # plt.savefig(f"img/fir_filter_{name}_response.svg")
    # plt.close()


def plotVaryingSlopeResponse(firFunc, length, targetSlopes):
    sampleRate = 48000
    cmapMain = plt.get_cmap("viridis")
    cmapRef = plt.get_cmap("gray")
    fig, ax = plt.subplots(len(targetSlopes), 2)
    for index, slopeDecibel in enumerate(targetSlopes):
        fir = firFunc(length, slopeDecibel)
        ω, response = signal.freqz(fir, fs=sampleRate, worN=2048)
        gain = 20 * np.log10(np.abs(response))
        ax[index][0].plot(
            ω,
            gain,
            color=cmapMain(index / len(targetSlopes)),
            label=f"{slopeDecibel} dB/oct",
        )
        refFreq = 11 - index
        plotReferenceSlope(
            ax[index][0],
            sampleRate / 2**refFreq,
            slopeDecibel,
            refFreq,
            "red",
            f"Ref. {slopeDecibel} dB/oct",
        )
        ax[index][1].plot(
            fir,
            color=cmapMain(index / len(targetSlopes)),
            label=f"{slopeDecibel} dB/oct",
        )
    for row in ax:
        row[0].set_ylabel("Gain [dB]")
        row[0].set_xscale("log")
        row[0].set_xlim([10, sampleRate // 2])
        row[0].set_ylim([-60, 10])

        row[1].set_ylabel("Amplitude")

        for axis in row:
            axis.grid(which="both", color="#f0f0f0")
            axis.legend()
    ax[0][0].set_title("Amplitude Reponse")
    ax[0][1].set_title("FIR Coefficient")
    ax[len(targetSlopes) - 1][0].set_xlabel(f"Frequency [Hz]")
    ax[len(targetSlopes) - 1][1].set_xlabel(f"Time [sample]")
    fig.set_size_inches((10, 8))
    plt.tight_layout()
    plt.show()


def iftExpCurve():
    ω = sympy.Symbol("ω", real=True)
    x = sympy.Symbol("x", real=True, positive=True)
    k = sympy.Symbol("k", real=True, positive=True)
    result = sympy.integrate(
        sympy.exp(-k * sympy.Abs(ω)) * sympy.exp(sympy.I * ω * x),
        (ω, -sympy.oo, sympy.oo),
    )
    result = sympy.simplify(result)
    print(sympy.latex(result))


def slopeFir(length: int, k: float):
    length -= (length + 1) % 2  # 係数の数を奇数にする。
    mid = length // 2

    fir = np.zeros(length)
    for i in range(length):
        x = i - mid
        fir[i] = 2 * k / (k * k + x * x)
    return fir / np.sum(fir)


def exponentialWindow(length: int, slopeDecibel: float):
    mid = length // 2

    fir = np.zeros(length)
    tau = (length / 2) * (np.e * 10 ** (slopeDecibel / 10))
    for i in range(length):
        x = i - mid
        fir[i] = np.exp(-np.abs(x) / tau)
    return fir / np.sum(fir)


def exponentialWindow2(length: int, slopeDecibel: float):
    tau = (length / 2) * (np.e / 10 ** (slopeDecibel / 10))
    fir = signal.get_window(("exponential", None, tau), length)
    return fir / np.sum(fir)


def testExponentialWindow():
    slopeDecibel = -6
    plotResponse(exponentialWindow(511, slopeDecibel), slopeDecibel)


def pinkNoiseFilter(length: int, slopeDecibel: float):
    N = length + length % 2  # Make it even number.
    mid = N // 2
    x = np.arange(1 - mid, 1 + mid)

    alpha = slopeDecibel / (40 * np.log10(1 / np.sqrt(2)))

    fir = np.ones(N)
    fir += np.cos(np.pi * x) / mid**alpha
    for k in range(1, mid):
        fir += 2 * np.cos(2 * np.pi * k * x / N) / k**alpha
    fir /= N
    if length % 2 == 1:
        return fir[0:-1]
    return fir


# plotResponse(slopeFir(63, 2), -9)
# plotResponse(exponentialWindow(64, -12), -12)
# plotResponse(pinkNoiseFilter(64, -12), -12)

length = 2048
targetSlopes = np.arange(-3, -13, -3)
# plotVaryingSlopeResponse(exponentialWindow, length, targetSlopes)
plotVaryingSlopeResponse(pinkNoiseFilter, length, targetSlopes)
