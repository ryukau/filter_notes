import numpy as np
import matplotlib.pyplot as plt
import scipy.signal as signal
import soundfile
from pathlib import Path

def ampTodB(amp):
    return 20 * np.log10(amp)

def loadSound(pathStr, upRate):
    data, fs = soundfile.read(pathStr, always_2d=True)
    ch0 = data.T[0]
    up = signal.resample(ch0, upRate * len(ch0))
    return (ch0, up)

def plot(source, result):
    upRate = 16

    sigSource, upSource = loadSound(source, upRate)
    sigResult, upResult = loadSound(result, upRate)

    print(f"Source true-peak: {ampTodB(np.max(np.abs(upSource)))} dB")
    print(f"Result true-peak: {ampTodB(np.max(np.abs(upResult)))} dB")

    length = 512  #len(sigSource)
    offset = 0

    sigSource = sigSource[offset:offset + length]
    upSource = upSource[offset:offset + length * upRate]
    sigResult = sigResult[offset:offset + length]
    upResult = upResult[offset:offset + length * upRate]

    timeUp = np.linspace(0, len(sigSource), len(upSource), endpoint=False)

    fig, ax = plt.subplots(1, 2)
    ax[0].set_title(f"Source")
    ax[0].plot(timeUp, upSource, lw=1, color="red", alpha=0.5, label=f"{upRate}x Source")
    ax[0].plot(sigSource, lw=1, color="black", alpha=1, label="Source")

    ax[1].set_title(f"True-peak Limited")
    ax[1].plot(
        timeUp,
        upResult,
        lw=1,
        color="red",
        alpha=0.5,
        zorder=11,
        label=f"{upRate}x TP Limited",
    )
    ax[1].plot(sigResult, lw=1, color="black", alpha=1, zorder=11, label="TP Limited")

    peakIndex = np.where(np.abs(upResult) > 1)[0]
    if len(peakIndex) > 0:
        print(peakIndex)
        print(upResult[peakIndex])
        ax[1].scatter(
            peakIndex.astype(np.float64) / upRate,
            upResult[peakIndex],
            marker="s",
            zorder=10,
            label="Overshoot",
        )

    for axis in ax:
        # axis.axhline(1, color="yellow", alpha=0.5, lw=1, ls="--", label="Threshold")
        # axis.axhline(-1, color="yellow", alpha=0.5, lw=1, ls="--")
        axis.set_ylim((-2.5, 2.5))
        axis.set_ylabel("Amplitude")
        axis.set_xlabel("Time [sample]")
        axis.legend(loc=1)
        axis.grid(color="#f0f0f0")
    fig.set_size_inches((10, 4))
    fig.tight_layout()
    plt.show()

if __name__ == "__main__":
    plot("data/flip204568972046735.wav", "data/flip_current.wav")
