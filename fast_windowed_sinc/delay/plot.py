import numpy as np
import scipy.signal as signal
import matplotlib.pyplot as plt
import json
from pathlib import Path
from mpl_toolkits.axes_grid1 import ImageGrid


imgPath = Path("img")


def plotDelayTime(name, data):
    """
    Max delay time is the half of signal length. See `delay.cpp`.
    """
    cmap = plt.get_cmap("plasma")
    for key, value in data.items():
        delayTime = int(key)
        plt.plot(value, lw=1, alpha=0.75, color=cmap(delayTime / len(data)), label=key)
    plt.title(f"Delay Time Test: {name}")
    plt.grid()
    # plt.show()
    plt.savefig(imgPath / Path(f"testDelayTime_{name}.svg"))
    plt.close()


def plotSpectrogram(name, data):
    sampleRateHz = float(data["sampleRate"])
    sig = np.array(data["sinMod"])
    logMagOffset = 1e-5

    frameSize = 1024
    win = signal.get_window("hann", frameSize)
    sft = signal.ShortTimeFFT(
        win, len(win) // 128, fs=sampleRateHz, scale_to="magnitude"
    )
    mag = sft.stft(sig)

    plt.figure(figsize=(10, 5))
    plt.title(f"Modulation Test: {name}")
    im1 = plt.imshow(
        20 * np.log10(abs(mag) + logMagOffset),
        origin="lower",
        aspect="auto",
        extent=sft.extent(len(sig)),
        cmap="magma",
    )
    plt.xlabel("Time [s]")
    plt.ylabel("Frequency [Hz]")
    # plt.ylim([10, 20000])
    # plt.yscale("log")
    offsetStr = np.format_float_scientific(logMagOffset, trim="-", exp_digits=1)
    plt.colorbar(im1, label=f"20*log10(abs(Magnitude) + {offsetStr}) ≈ Magnitude in dB")
    plt.tight_layout()
    # plt.show()
    plt.savefig(imgPath / Path(f"testAntialiasing_{name}.svg"))
    plt.close()


def collageSpectrogram(fullData):
    """
    Reference on plotting:
    - [python - Matplotlib - Tight layout of multiple subplots with colorbar - Stack Overflow](https://stackoverflow.com/questions/74267861/matplotlib-tight-layout-of-multiple-subplots-with-colorbar)
    """
    sources = [
        "DelayInt",
        "DelayLinear",
        "DelayLagrange3",
        "DelayAntialiasedBiquad2",
        "DelayTriangleBiquadSine",
        "DelayBlackmanHarrisBiquadSine",
    ]
    logMagOffset = 1e-5

    nSources = len(sources)
    nCol = int(np.ceil(np.sqrt(nSources)))
    nRow = nCol - (nCol * nCol - nSources) // nCol
    fig, ax = plt.subplots(nRow, nCol, layout="compressed")
    fig.set_size_inches((16, 8))
    im = None
    for idx, key in enumerate(sources):
        data = fullData[key]["testAntialiasing"]
        sampleRateHz = float(data["sampleRate"])
        sig = np.array(data["sinMod"])

        frameSize = 1024
        win = signal.get_window("hann", frameSize)
        sft = signal.ShortTimeFFT(
            win, len(win) // 128, fs=sampleRateHz, scale_to="magnitude"
        )
        mag = sft.stft(sig)

        axis = ax[idx // nCol][idx % nCol]
        axis.set_title(key)
        im = axis.imshow(
            20 * np.log10(abs(mag) + logMagOffset),
            origin="lower",
            aspect="auto",
            extent=sft.extent(len(sig)),
            cmap="magma",
        )
        axis.set_xlabel("Time [s]")
        axis.set_ylabel("Frequency [Hz]")
    if im is not None:
        offsetStr = np.format_float_scientific(logMagOffset, trim="-", exp_digits=1)
        plt.colorbar(
            im,
            ax=ax.ravel().tolist(),
            aspect=40,
            pad=0.025,
            label=f"20*log10(abs(Magnitude) + {offsetStr}) ≈ Magnitude in dB",
        )
    plt.savefig(imgPath / Path(f"aa_comparison_on_delaytime_mod.svg"))
    # plt.show()
    plt.close()


if __name__ == "__main__":
    imgPath.mkdir(parents=True, exist_ok=True)

    with open("cpp_test_output.json", "r", encoding="utf-8") as fp:
        fullData = json.load(fp)

    for name, data in fullData.items():
        plotDelayTime(name, data["testDelayTime"])
        plotSpectrogram(name, data["testAntialiasing"])
    collageSpectrogram(fullData)
