import numpy as np
import scipy.signal as signal
import matplotlib.pyplot as plt


def toDecibel(r):
    db = 20 * np.log10(np.abs(r))
    return np.where(db < -150, -150, db)


def toPhase(r):
    return np.angle(r) / np.pi


def getSosLP(sampleRate, cutoffHz):
    k = 1 / np.tan(np.pi * cutoffHz / sampleRate)

    a0 = 1 + k
    a1 = (1 - k) / a0
    g = 1 / a0

    return [[g, g, 0, 1, a1, 0]]  # `sos` (second order series) format in scipy.signal.


def getSosHP(sampleRate, cutoffHz):
    k = 1 / np.tan(np.pi * cutoffHz / sampleRate)

    a0 = 1 + k
    b0 = k / a0
    b1 = -b0
    a1 = (1 - k) / a0

    return [[b0, b1, 0, 1, a1, 0]]


def getSosEmaLP(sampleRate, cutoffHz):
    """
                    y0 = y1 + kp * (x0 - y1)
                    y0 = kp * x0 + (1 - kp) * y1
    y0 - (1 - kp) * y1 = kp * x0
    """
    y = 1 - np.cos(2 * np.pi * cutoffHz / sampleRate)
    kp = np.sqrt(y * y + 2 * y) - y
    return [[kp, 0, 0, 1, kp - 1, 0]]


def getSosEmaHP(sampleRate, cutoffHz):
    """
    H(z) = 1 - kp / (1 + (kp - 1)z^-1)
         = (1 - kp + (kp - 1)z^-1) / (1 + (kp - 1)z^-1)
    """

    y = 1 - np.cos(2 * np.pi * cutoffHz / sampleRate)
    kp = np.sqrt(y * y + 2 * y) - y
    return [[1 - kp, kp - 1, 0, 1, kp - 1, 0]]


def plotSosResponse(sosFunc, sampleRate, cutoffHz):
    sos = sosFunc(sampleRate, cutoffHz)
    freq, response = signal.sosfreqz(sos, worN=2048, fs=sampleRate)

    plt.rcParams["figure.figsize"] = (6, 3)
    plt.title(f"Bilinear LP1 Gain, Accurate Cutoff, fs={sampleRate}, cutoff={cutoffHz}")
    plt.plot(freq, toDecibel(response), label="gain")
    # plt.plot(freq, toPhase(response), label="phase")
    plt.axvline(cutoffHz, color="black", ls="--", label="cutoff")
    plt.xscale("log")
    plt.xlabel("Frequency [Hz]")
    plt.ylim([-6, 3])
    plt.yticks(np.arange(-6, 3 + 1))
    plt.ylabel("Gain [dB]")
    plt.grid(which="both", color="#f8f8f8")
    plt.legend()
    plt.tight_layout()
    plt.show()


def plotSosResponses(sosFunc, sampleRate):
    cmap = plt.get_cmap("plasma")
    cutoffs = [10, 30, 100, 300, 1000, 3000, 10000]

    fig, ax = plt.subplots(2, 1)
    fig.set_size_inches((8, 6))
    fig.suptitle(f"Bilinear Transformed 1-pole Lowpass, fs={sampleRate}")
    for index, cutoffHz in enumerate(cutoffs):
        sos = sosFunc(sampleRate, cutoffHz)
        freq, response = signal.sosfreqz(sos, worN=2**16, fs=sampleRate)
        ax[0].plot(
            freq,
            toDecibel(response),
            alpha=0.75,
            color=cmap(index / len(cutoffs)),
            label=f"$f_c$={cutoffHz:.0f}",
        )
        ax[1].plot(
            freq,
            toPhase(response),
            alpha=0.75,
            color=cmap(index / len(cutoffs)),
            label=f"$f_c$={cutoffHz:.0f}",
        )
    ax[0].scatter(
        cutoffs,
        [-3 for _ in range(len(cutoffs))],
        color="black",
        marker="+",
        zorder=2,
        label="cutoff",
    )
    ax[1].scatter(
        cutoffs,
        [-0.25 for _ in range(len(cutoffs))],
        # [0.25 for _ in range(len(cutoffs))],
        color="black",
        marker="+",
        zorder=2,
        label="cutoff",
    )
    for axis in ax:
        axis.set_xscale("log")
        axis.grid(which="both", color="#f4f4f4")
        axis.legend()
    ax[0].set_ylim([-6, 2])
    ax[0].set_yticks(np.arange(-6, 2 + 1))
    ax[0].set_ylabel("Gain [dB]")
    ax[1].set_yticks(np.linspace(-0.5, 0, 3))
    # ax[1].set_yticks(np.linspace(+0.5, 0, 3))
    ax[1].set_ylabel("Phase [dB]")
    ax[1].set_xlabel("Frequency [Hz]")
    fig.tight_layout()
    plt.show()
    # plt.savefig("../img/blt_one_pole.svg")


def plotSosIR(sosFunc, sampleRate, cutoffHz):
    """IR == Impulse Response."""

    impulse = np.zeros(2**8)
    impulse[0] = 1

    sos = sosFunc(sampleRate, cutoffHz)
    ir = signal.sosfilt(sos, impulse)

    plt.plot(impulse)
    plt.plot(ir)
    plt.grid()
    plt.show()


def testNyquist(sampleRate, cutoffHz):
    x = np.ones(2**8)
    x[1::2] = -1
    x *= signal.get_window("hann", len(x))

    y_blt = signal.sosfilt(getSosLP(sampleRate, cutoffHz), x)
    y_ema = signal.sosfilt(getSosEmaLP(sampleRate, cutoffHz), x)

    plt.title("Nyquist Frequency Attenuation (Hann Window)")
    plt.plot(y_ema, color="black", alpha=0.25, label="EMA")
    plt.plot(y_blt, color="red", alpha=0.5, label="bilinear")
    plt.xlabel("Time [sample]")
    plt.ylabel("Amplitude")
    plt.legend()
    plt.grid()
    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    sampleRate = 48000
    cutoffHz = 200

    plotSosResponse(getSosLP, sampleRate, cutoffHz)
    plotSosResponses(getSosLP, sampleRate)
    plotSosIR(getSosLP, sampleRate, cutoffHz)
    testNyquist(sampleRate, cutoffHz)
