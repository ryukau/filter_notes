import matplotlib.pyplot as plt
import numpy as np
import scipy.signal as signal
import soundfile
from pathlib import Path


def plot(sig, samplerate, name):
    freq, time, spec = signal.spectrogram(
        sig,
        samplerate,
        "blackmanharris",
        nperseg=512,
        noverlap=512 * 7 // 8,
    )

    absed = np.abs(spec)
    spec = 20 * np.log10(absed / np.max(absed))
    spec = np.where(spec < -200, -200, spec)

    plt.figure(figsize=(8, 4), tight_layout=True)
    plt.pcolormesh(time, freq, spec, cmap="magma", shading="gouraud")
    # plt.yscale("log")
    # plt.ylim((20, samplerate / 2))
    plt.ylabel("Frequency [Hz]")
    plt.xlabel("Time [s]")
    plt.savefig("img/" + name + ".png")


def plotAll():
    Path("img").mkdir(parents=True, exist_ok=True)
    for path in Path("snd").glob("*.wav"):
        print(f"{path}")
        data, samplerate = soundfile.read(str(path))
        plot(data, samplerate, path.stem)


def getSpectrum(path):
    sig, samplerate = soundfile.read(str(path))

    freq, time, spec = signal.spectrogram(
        sig,
        samplerate,
        "blackmanharris",
        nperseg=512,
        noverlap=512 * 7 // 8,
    )

    absed = np.abs(spec)
    spec = 20 * np.log10(absed / np.max(absed))
    spec = np.where(spec < -200, -200, spec)

    return (time, freq, spec)


def plotDownSampling():
    data = [
        getSpectrum("snd/Source.wav"),
        getSpectrum("snd/Array_Down.wav"),
    ]

    fig, axes = plt.subplots(2)
    fig.set_size_inches(8, 6)

    for idx, ax in enumerate(axes):
        ax.pcolormesh(*data[idx], cmap="magma", shading="gouraud")
        ax.set_ylabel("Frequency [Hz]")
        ax.set_xticks([i / 2 for i in range(3)])
    axes[-1].set_xlabel("Time [s]")

    axes[0].set_title("Source (Sine Sweep from $0$ to $f_s$ Hz)")
    axes[1].set_title("2-fold Down-sampled")

    plt.tight_layout()
    plt.savefig("img/Result_Down.png")


def plotUpSampling():
    data = [
        getSpectrum("snd/Source.wav"),
        getSpectrum("snd/Array_Up.wav"),
    ]

    fig, axes = plt.subplots(2)
    fig.set_size_inches(8, 6)

    for idx, ax in enumerate(axes):
        ax.pcolormesh(*data[idx], cmap="magma", shading="gouraud")
        ax.set_ylabel("Frequency [Hz]")
        ax.set_xticks([i / 2 for i in range(3)])
    axes[-1].set_xlabel("Time [s]")

    axes[0].set_title("Source (Sine Sweep from $0$ to $f_s/2$ Hz)")
    axes[1].set_title("2-fold Up-sampled")

    plt.tight_layout()
    plt.savefig("img/Result_Up.png")


if __name__ == "__main__":
    print("--- Plotting")

    plotAll()
    plotDownSampling()
    plotUpSampling()
