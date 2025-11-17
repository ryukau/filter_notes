import numpy as np
import scipy.signal as signal
import matplotlib.pyplot as plt
import json


def plotReset(data):
    latency = data["latency"]
    low0 = data["testIrLow"]
    low1 = data["testResetLow"]
    high0 = data["testIrHigh"]
    high1 = data["testResetHigh"]

    firLength = 0
    for i in range(len(low0) - 1, -1, -1):
        if low0[i] != 0:
            firLength = i
            break

    fig, ax = plt.subplots(2, 1)
    fig.set_size_inches((6, 6))

    ax[0].plot(low0, color="red", alpha=0.5, label="first")
    ax[0].plot(low1, color="blue", alpha=0.5, label="reset")

    ax[1].plot(high0, color="red", alpha=0.5, label="first")
    ax[1].plot(high1, color="blue", alpha=0.5, label="reset")

    for axis in ax:
        axis.axvline(latency, color="black", ls="--", lw=1, label="latency")
        axis.grid()
        axis.legend()
        axis.set_xlim([0, firLength + 1])
    fig.tight_layout()
    plt.show()


def plotResponse(data):
    cutoff = data["cutoff"]
    low0 = data["testIrLow"]
    low1 = data["testResetLow"]
    high0 = data["testIrHigh"]
    high1 = data["testResetHigh"]

    def getGain(x):
        freq, resp = signal.freqz(x, 1, worN=2048, fs=1)
        absed = np.abs(resp)
        return [freq, 20 * np.log10(absed / np.max(absed))]

    fig, ax = plt.subplots(2, 1)
    fig.set_size_inches((6, 6))
    ax[0].plot(*getGain(low0), color="red", alpha=0.5, label="first")
    ax[0].plot(*getGain(low1), color="blue", alpha=0.5, label="reset")

    ax[1].plot(*getGain(high0), color="red", alpha=0.5, label="first")
    ax[1].plot(*getGain(high1), color="blue", alpha=0.5, label="reset")

    for axis in ax:
        axis.axvline(cutoff, color="black", ls="--", lw=1, label="cutoff")
        axis.grid()
        axis.legend()

    fig.tight_layout()
    plt.show()


if __name__ == "__main__":
    with open("cpp_test_output.json", "r", encoding="utf-8") as fp:
        data = json.load(fp)
    plotReset(data)
    plotResponse(data)
